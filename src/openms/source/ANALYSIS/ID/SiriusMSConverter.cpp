// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>

#include <fstream>

#include <QDir>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  //Precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz, const MSSpectrum& spectrum1, double left_tolerance, double right_tolerance)
  {
    MSSpectrum::ConstIterator left = spectrum1.MZBegin(test_mz - left_tolerance);
    MSSpectrum::ConstIterator right = spectrum1.MZEnd(test_mz + right_tolerance);

    // no MS1 precursor peak in +- tolerance window found
    if (left == right || left->getMZ() > test_mz + right_tolerance)
    {
      return -1;
    }

    MSSpectrum::ConstIterator max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

    if (max_intensity_it == right)
    {
      return -1;
    }

    return max_intensity_it - spectrum1.begin();
  }

  void SiriusMSFile::store(const PeakMap &spectra, const OpenMS::String & msfile)
  {

    int count_skipped_spectra = 0; // spectra skipped due to precursor charge
    int count_to_pos = 0; // count if charge 0 -> +1
    int count_to_neg = 0; // count if charge 0 -> -1
    int count_no_ms1 = 0; // count if no precursor was found

    //check for all spectra at the beginning if spectra are centroided
    //determine type of spectral data (profile or centroided) - only checking first spectrum (could be ms2 spectrum)
    SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

    if (spectrum_type == SpectrumSettings::RAWDATA)
    {
      throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided spectra are needed. Please use PeakPicker to convert the spectra.");
    }

    // loop over all spectra in file and write data to ofstream
    ofstream os;

    // create temporary input file (.ms)
    os.open(msfile.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msfile);
    }
    os.precision(12);

    for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      // process only MS2 spectra
      if (s_it->getMSLevel() != 2)
      {
        continue;
      }

      const MSSpectrum& spectrum = *s_it;

      int scan_index = s_it - spectra.begin();

      const vector<Precursor>& precursor = spectrum.getPrecursors();
      double collision = precursor[0].getActivationEnergy(); //extract collision energy - this function does not work

      IonSource::Polarity p = spectrum.getInstrumentSettings().getPolarity(); //charge

      // there should be only one precursor and MS2 should contain peaks to be considered
      if (precursor.size() == 1 && !spectrum.empty())
      {
        // needed later for writing in ms file
        int int_charge = 0;

        // read precursor charge
        int precursor_charge = precursor[0].getCharge();

        // sirius supports only single charged ions (+1; -1)
        // if charge = 0, it will be allocted to +1; -1 depending on Polarity
        if (precursor_charge > 1 || precursor_charge <= -1)
        {
          count_skipped_spectra = count_skipped_spectra + 1;
          continue;
        }
        // set charge value for msfile
        if (p == IonSource::Polarity::POSITIVE && precursor_charge == +1)
        {
          int_charge = +1;
        }
        if (p == IonSource::Polarity::NEGATIVE && precursor_charge == -1)
        {
          int_charge = -1;
        }
        if (p == IonSource::Polarity::POSITIVE && precursor_charge == 0)
        {
          int_charge = +1;
          count_to_pos = count_to_pos + 1;
        }
        if (p == IonSource::Polarity::NEGATIVE && precursor_charge == 0)
        {
          int_charge = -1;
          count_to_neg = count_to_neg +1;
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
          count_no_ms1 = count_no_ms1 +1;
        }
        //get the precursor in the ms1 spectrum (highest intensity in the range of the precursor mz +- 0.1 Da)
        else
        {
          const MSSpectrum& spectrum1 = *s_it2;

          // 0.2 Da left and right of the precursor m/z - we do not expect metabolites with charge 5 or higher.
          Int mono_index = getHighestIntensityPeakInMZRange(test_mz, spectrum1, 0.2, 0.2);

          if (mono_index != -1)
          {
            const Peak1D& max_mono_peak = spectrum1[mono_index];
            isotopes.push_back(max_mono_peak);

            // make sure the 13C isotopic peak is picked up by doubling the (fractional) mass difference (approx. 1.0066)
            const double C13_dd = 2.0 * (Constants::C13C12_MASSDIFF_U - 1.0);
            Int iso1_index = getHighestIntensityPeakInMZRange(max_mono_peak.getMZ() + Constants::C13C12_MASSDIFF_U, spectrum1, 0, C13_dd);

            if (iso1_index != -1)
            {
              const Peak1D& iso1_peak = spectrum1[iso1_index];
              isotopes.push_back(iso1_peak);
              Int iso2_index = getHighestIntensityPeakInMZRange(iso1_peak.getMZ() + Constants::C13C12_MASSDIFF_U, spectrum1, 0, C13_dd);

              if (iso2_index != -1)
              {
                const Peak1D& iso2_peak = spectrum1[iso2_index];
                isotopes.push_back(iso2_peak);
              }
            }
          }
        }

        String query_id = String("unknown") + String(scan_index);

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
        Size no_isotopes = isotopes.size();

        if ( no_isotopes > 0) //if ms1 spectrum was present
        {
          os << ">ms1" << endl;
          //m/z and intensity have to be higher than 1e-10
          //the intensity of the peaks of the isotope pattern have to be smaller than the one before

          double threshold = 1e-10;

          double first_mz = isotopes[0].getMZ();
          double first_intensity = isotopes[0].getIntensity();

          if (first_mz > threshold && first_intensity > threshold)
          {
            os << first_mz << " " << first_intensity << endl;
          }

          if (no_isotopes > 1)
          {
            double second_mz = isotopes[1].getMZ();
            double second_intensity = isotopes[1].getIntensity();

            if (second_mz > threshold && second_intensity > threshold && second_intensity < first_intensity && first_intensity > threshold)
            {
              os << second_mz << " " << second_intensity << endl;
            }

            if (no_isotopes > 2)
            {
              double third_mz = isotopes[2].getMZ();
              double third_intensity = isotopes[2].getIntensity();

              if (third_mz > threshold && third_intensity > threshold && third_intensity < second_intensity && second_intensity > threshold)
              {
                os << third_mz << " " << third_intensity << endl;
              }
            }
          }
          os << endl;
        }
        else
        {
          if (precursor_int != 0) // if no ms1 spectrum was present but precursor intensity is known
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
        for (Size i = 0; i < spectrum.size(); ++i)
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
        os << endl;
      }
    }
    os.close();

    LOG_WARN << "No MS1 spectrum for this precursor. Occurred " << count_no_ms1 << " times." << endl;
    LOG_WARN << count_skipped_spectra << " spectra were skipped due to precursor charge below -1 and above +1." << endl;
    LOG_WARN << "Charge of 0 was set to +1 due to positive polarity " << count_to_pos << " times."<< endl;
    LOG_WARN << "Charge of 0 was set to -1 due to negative polarity " << count_to_neg << " times." << endl;

  }
}

/// @endcond
