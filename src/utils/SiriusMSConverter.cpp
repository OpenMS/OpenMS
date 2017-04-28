// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

//not sure if more #include directives are needed

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

#include <fstream>

using namespace OpenMS;
using namespace std;



//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
@page UTILS_SiriusMSConverter

@brief Tool for the conversion of mzML files to (Sirius).ms files

       Needed for the interal data structure of the Sirius command line tool

<B>The command line parameters of this tool are:</B>
@verbinclude UTILS_SiriusMSConverter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude UTILS_SiriusMSConverter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSiriusMSConverter :
    public TOPPBase
{
public:
  TOPPSiriusMSConverter() :
    TOPPBase("SiriusMSConverter", "Tool for the conversion of mzML files to (Sirius).ms files", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerOutputFile_("out", "<file>", "", "Output of (Sirius).ms files");
    setValidFormats_("out", ListUtils::create<String>("csv"));

    registerIntOption_("batch_size", "<value>", 20,"The number of compounds used in one output file", false);
  }

  //Precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz, const MSSpectrum<Peak1D>& spectrum1, double left_tolerance, double right_tolerance)
  {
    MSSpectrum<Peak1D>::ConstIterator left = spectrum1.MZBegin(test_mz - left_tolerance);
    MSSpectrum<Peak1D>::ConstIterator right = spectrum1.MZEnd(test_mz + right_tolerance);

    // no MS1 precursor peak in +- tolerance window found
    if (left == right || left->getMZ() > test_mz + right_tolerance)
    {
      return -1; //not sure if that is allright
    }

    MSSpectrum<Peak1D>::ConstIterator max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

    if (max_intensity_it == right) return -1;

    return max_intensity_it - spectrum1.begin();
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    Int batch_size = getIntOption_("batch_size");

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------

    // laod spectra
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, spectra);
    int count = 0;

    //check for all spectra at the beginning if spectra are centroided
    // determine type of spectral data (profile or centroided) - only checking first spectrum (could be ms2 spectrum)
    SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

    if (spectrum_type == SpectrumSettings::RAWDATA)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided spectra are needed. Please use PeakPicker to convert the spectra.");
    }

    // loop over all spectra
    ofstream os;

    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      // process only MS2 spectra
      if (s_it->getMSLevel() != 2) continue;

      const MSSpectrum<Peak1D>& spectrum = *s_it;

      int scan_index = s_it - spectra.begin();

      const vector<Precursor>& precursor = spectrum.getPrecursors();
      double collision = precursor[0].getActivationEnergy(); //not working?

      IonSource::Polarity p = spectrum.getInstrumentSettings().getPolarity(); //charge

      // needed later for writing in ms file
      int int_charge(1);

      if (p == IonSource::Polarity::POSITIVE)
      {
        int_charge = +1;
      }
      else
      {
        LOG_WARN << "SiriusMSConverter (due to Sirius) only support positive ion mode and mono charged analytes." << endl;
        continue;
      }

      //there should be only one precursor and MS2 should contain peaks to be considered
      if (precursor.size() == 1 && !spectrum.empty())
      {
        //read charge annotated to MS2
        int precursor_charge = precursor[0].getCharge();

        //sirius only supports +1 charge so far
        if (precursor_charge > 1 || precursor_charge <= -1)
        {
          LOG_WARN << "SiriusMSConverter (due to Sirius) only support positively mono charged analytes." << endl;
        }

        //get m/z and intensity of precursor != MS1 spectrum
        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();

        //find corresponding ms1 spectra (precursor)
        PeakMap::ConstIterator s_it2 = spectra.getPrecursorSpectrum(s_it);

        double test_mz = precursor_mz;

        vector<Peak1D> isotopes;
        isotopes.clear();

        if (s_it2->getMSLevel() != 1)
        {
          LOG_WARN << "Error: no MS1 spectrum for this precursor. No isotopes considered in sirius." << endl;
        }
        //get the precursor in the ms1 spectrum (highest intensity in the range of the precursor mz +- 0.1 Da
        else
        {
          const MSSpectrum<Peak1D>& spectrum1 = *s_it2;

          Int mono_index = getHighestIntensityPeakInMZRange(test_mz, spectrum1, 0.2, 0.2);

          if (mono_index != -1)
          {
            const Peak1D& max_mono_peak = spectrum1[mono_index];
            isotopes.push_back(max_mono_peak);

            // make sure the 13C isotopic peak is picked up by doubling the (fractional) mass difference (approx. 1.0066)
            const double C13_dd = 2.0 * (Constants::C13C12_MASSDIFF_U - 1.0);
            Int iso1_index = getHighestIntensityPeakInMZRange(max_mono_peak.getMZ() + 1.0, spectrum1, 0, C13_dd);

            if (iso1_index != -1)
            {
              const Peak1D& iso1_peak = spectrum1[iso1_index];
              isotopes.push_back(iso1_peak);
              Int iso2_index = getHighestIntensityPeakInMZRange(iso1_peak.getMZ() + 1.0, spectrum1, 0, C13_dd);

              if (iso2_index != -1)
              {
                const Peak1D& iso2_peak = spectrum1[iso2_index];
                isotopes.push_back(iso2_peak);
              }
            }
          }
        }

        String query_id = String("unknown") + String(scan_index);

        streamsize prec(0);

        if (count == 0)
        {
            // store data
            String unique_name =  String(File::getUniqueName()).toQString(); //if not done this way - always new "unique name"
            String dir = getStringOption_("out");
            String filename = dir  + "/" + unique_name.toQString() + "_" + query_id.toQString() + ".ms";

            // close previous (.ms) file
            if (os.is_open()) os.close();

            // create temporary input file (.ms)
            os.open(filename.c_str());

            if (!os)
            {
              throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
            }
            prec = os.precision();
            os.precision(12);
        }


        //TODO: Collision energy optional for MS2 - something wrong with .getActivationEnergy() (Precursor)
        //TODO: MS2 instensity cutoff? -> to reduce the interference of low intensity peaks (?) - can do that for specific spectra (hard coded/soft coded?)

        //write internal unique .ms data as sirius input
        os << fixed;
        os << ">compound " << query_id << "\n";
        if (isotopes.empty() == false)
        {
          os << ">parentmass " << isotopes[0].getMZ() << fixed << "\n";
        }
        else
        {
          os << ">parentmass " << precursor_mz << fixed << "\n";
        }
        os << ">charge " << int_charge << "\n\n";

        // Use precursor m/z & int and no ms1 spectra is available else use values from ms1 spectrum

        if (isotopes.empty() == false) //if ms1 spectrum was present
        {
          os << ">ms1" << "\n";
          //m/z and intensity have to be higher than 1e-10
          //the intensity of the peaks of the isotope pattern have to be smaller than the one before
          if (isotopes[0].getMZ() > 1e-10 && isotopes[0].getIntensity() > 1e-10){ os << isotopes[0].getMZ() << " " << isotopes[0].getIntensity() << "\n";}
          if (isotopes[1].getMZ() > 1e-10 && isotopes[1].getIntensity() > 1e-10 && isotopes[1].getIntensity() < isotopes[0].getIntensity() > 1e-10){ os << isotopes[1].getMZ() << " " << isotopes[1].getIntensity() << "\n";}
          if (isotopes[2].getMZ() > 1e-10 && isotopes[2].getIntensity() > 1e-10 && isotopes[2].getIntensity() < isotopes[1].getIntensity() > 1e-10){ os << isotopes[2].getMZ() << " " << isotopes[2].getIntensity() << "\n";}
          os << "\n";
        }
        else
        {
          if(precursor_int != 0) // if no ms1 spectrum was present but precursor intensity is known
          {
            os << ">ms1" << "\n"
               << precursor_mz << " " << precursor_int << "\n\n";
          }
        }

        //if collision energy was given - write it into .ms file if not use ms2 instead
        if (collision == 0.0)
        {
          os << ">ms2" << "\n";
        }
        else
        {
          os << ">collision" << " " << collision << "\n";
        }

        //single spectrum peaks
        for (Size i = 0; i != spectrum.size(); ++i)
        {
          const Peak1D& peak = spectrum[i];
          double mz = peak.getMZ();
          float intensity = peak.getIntensity();

          //intensity has to be higher than zero
          if (intensity != 0)
          {
            os << mz << " " << intensity << "\n";
          }
        }
        os << "\n";
      }

      // increase count and reset to zero if batch size reached
      count = (count + 1) % batch_size;

    }

    // close previous (.ms) file
    if (os.is_open()) os.close();

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPSiriusMSConverter tool;
  return tool.main(argc, argv);
}

/// @endcond
