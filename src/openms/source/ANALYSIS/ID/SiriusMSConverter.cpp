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

#include <OpenMS/CONCEPT/LogStream.h>
#include <fstream>

#include <QDir>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  // precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz,
                                       const MSSpectrum& spectrum1,
                                       double left_tolerance,
                                       double right_tolerance)
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

  // extract precursor isotope pattern if no feature information is available
  vector<Peak1D> extractPrecursorC13IsotopePattern(const double& precursor_mz,
                                                   const MSSpectrum& precursor_spectrum,
                                                   int& iterations,
                                                   int& charge)
  {
    vector<Peak1D> isotopes;
    int peak_index;
    Peak1D peak;

    // monoisotopic_trace
    peak_index = getHighestIntensityPeakInMZRange(precursor_mz, precursor_spectrum, 0.2, 0.2);
    if (peak_index != -1)
    {
      peak = precursor_spectrum[peak_index];
      isotopes.push_back(peak);
    }

    // further isotope_traces: make sure the 13C isotopic peak is picked up by doubling the (fractional) mass difference (approx. 1.0066)
    double C13_dd = 2.0 * (Constants::C13C12_MASSDIFF_U - 1.0);
    double massdiff = Constants::C13C12_MASSDIFF_U;

    // depending on the charge different MASSDIFF
    if (charge != 0)
    {
      massdiff = massdiff/charge;
    }

    while (peak_index != -1 && iterations > 0)
    {
      peak_index = getHighestIntensityPeakInMZRange(peak.getMZ() + massdiff, precursor_spectrum, 0, C13_dd);
      if (peak_index != -1)
      {
        peak = precursor_spectrum[peak_index];
        isotopes.push_back(peak);
      }
      iterations = iterations - 1;
    }
    return isotopes;
  }

  void SiriusMSFile::store(const PeakMap &spectra,
                           const OpenMS::String & msfile,
                           const std::map< const size_t, const BaseFeature* > & feature_ms2_spectra_map,
                           const bool& feature_only,
                           const int& isotope_pattern_iterations,
                           const bool no_masstrace_info_isotope_pattern)
  {
    int count_skipped_spectra = 0; // spectra skipped due to precursor charge
    int count_to_pos = 0; // count if charge 0 -> +1
    int count_to_neg = 0; // count if charge 0 -> -1
    int count_no_ms1 = 0; // count if no precursor was found

    // check for all spectra at the beginning if spectra are centroided
    // determine type of spectral data (profile or centroided) - only checking first spectrum (could be ms2 spectrum)
    SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

    // extract native id type accession (e.g. MS:1000768) corresponding to native id type (Thermo nativeID format)
    const String native_id_type_accession = spectra.getExperimentalSettings().getSourceFiles()[0].getNativeIDTypeAccession();
    LOG_DEBUG << "native_id_type_accession: " << native_id_type_accession << std::endl;
  
    if (spectrum_type == SpectrumSettings::PROFILE)
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

      // if no feature information is available for the ms2 spectrum and only ms2 with features should be used - continue
      if (!(feature_ms2_spectra_map.find(scan_index) != feature_ms2_spectra_map.end()) && feature_only)
      {
        continue;
      }

      const String native_id = spectrum.getNativeID();
      int scan_number = SpectrumLookup::extractScanNumber(native_id, native_id_type_accession);
      StringList adducts;
      String feature_id;
      int feature_charge = 0;
      vector<std::pair<double, double> > f_isotopes;
      f_isotopes.clear();

      // extract metadata for specific ms2 spectrum
      // if the index exists extract the meta information (if available)
      if (feature_ms2_spectra_map.find(scan_index) != feature_ms2_spectra_map.end())
      {
        const BaseFeature* feature = feature_ms2_spectra_map.at(scan_index);

        feature_id = feature->getUniqueId();
        feature_charge = feature->getCharge();

        if (feature->metaValueExists("adducts"))
        {
          adducts = feature->getMetaValue("adducts");
        }
        if (feature->metaValueExists("masstrace_centroid_mz") && feature->metaValueExists("masstrace_intensity"))
        {
          vector<double> masstrace_centroid_mz = feature->getMetaValue("masstrace_centroid_mz");
          vector<double> masstrace_intensity = feature->getMetaValue("masstrace_intensity");
          if (masstrace_centroid_mz.size() == masstrace_intensity.size())
          {
            for (Size i = 0; i < masstrace_centroid_mz.size(); ++i)
            {
              std::pair<double, double> masstrace_mz_int(masstrace_centroid_mz[i],masstrace_intensity[i]);
              f_isotopes.push_back(masstrace_mz_int);
            }
          }
        }
      }

      const vector<Precursor>& precursor = spectrum.getPrecursors();

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

        // get m/z and intensity of precursor != MS1 spectrum
        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();

        // extract collision energy
        double collision = precursor[0].getActivationEnergy(); 
        
        // find corresponding ms1 spectra (precursor)
        PeakMap::ConstIterator s_it2 = spectra.getPrecursorSpectrum(s_it);

        double test_mz = precursor_mz;

        vector<Peak1D> isotopes;
        isotopes.clear();
        vector<Peak1D> precursor_spec;

        if (s_it2->getMSLevel() != 1)
        {
          count_no_ms1 = count_no_ms1 +1;
        }
        // get the precursor in the ms1 spectrum (highest intensity in the range of the precursor mz +- 0.1 Da)
        else
        {
          const MSSpectrum& precursor_spectrum = *s_it2;
          int interations = isotope_pattern_iterations;
          // extract precursor isotope pattern via C13 isotope distance
          isotopes = extractPrecursorC13IsotopePattern(test_mz, precursor_spectrum, interations, feature_charge);

          for (Size i = 0; i < precursor_spectrum.size(); ++i)
          {
            const Peak1D& peak = precursor_spectrum[i];
            precursor_spec.push_back(peak);
          }

        }

        String query_id = String("_" + feature_id) + String("-" + String(scan_number) + "-") + String("unknown") + String(scan_index);

        // write internal unique .ms data as sirius input
        os << fixed;
        os << ">compound " << query_id << "\n";

        if (!adducts.empty())
        {
          os << ">ionization " << ListUtils::concatenate(adducts, ',') << "\n";
        }

        if (!f_isotopes.empty() && !no_masstrace_info_isotope_pattern)
        {
          os << ">parentmass " << f_isotopes[0].first << fixed << "\n";
        }
        else if (!isotopes.empty())
        {
          os << ">parentmass " << isotopes[0].getMZ() << fixed << "\n";
        }
        else
        {
          os << ">parentmass " << precursor_mz << fixed << "\n";
        }

        if (feature_charge != 0)
        {
          os << ">charge " << feature_charge << "\n\n";
        }
        else
        {
          os << ">charge " << int_charge << "\n\n";
        }

        // Use precursor m/z & int and no ms1 spectra is available else use values from ms1 spectrum
        Size no_isotopes = isotopes.size();
        Size no_f_isotopes = f_isotopes.size();

        if (no_f_isotopes > 0 && !no_masstrace_info_isotope_pattern)
        {
          os << ">ms1merged" << endl;

          // m/z and intensity have to be higher than 1e-10
          for (auto it = f_isotopes.begin(); it != f_isotopes.end(); ++it)
          {
            os << it->first << " " << it->second <<"\n";
          }
        }
        else if (no_isotopes > 0) // if ms1 spectrum was present
        {
          os << ">ms1merged" << endl;

          // m/z and intensity have to be higher than 1e-10
          for (auto it = isotopes.begin(); it != isotopes.end(); ++it)
          {
              os << it->getMZ() << " " << it->getIntensity() << "\n";
          }
        }
        else
        {
          if (precursor_int != 0) // if no ms1 spectrum was present but precursor intensity is known
          {
            os << ">ms1merged" << "\n" << precursor_mz << " " << precursor_int << "\n\n";
          }
        }

        if (!precursor_spec.empty())
        {
          os << ">ms1peaks" << endl;
          for (auto iter = precursor_spec.begin(); iter != precursor_spec.end(); ++iter)
          {
            os << iter->getMZ() << " " << iter->getIntensity() << "\n";
          }
        }

        // if collision energy was given - write it into .ms file if not use ms2 instead
        if (collision == 0.0)
        {
          os << ">ms2peaks" << "\n";
        }
        else
        {
          os << ">collision" << " " << collision << "\n";
        }

        // single spectrum peaks
        for (Size i = 0; i < spectrum.size(); ++i)
        {
          const Peak1D& peak = spectrum[i];
          double mz = peak.getMZ();
          float intensity = peak.getIntensity();

          // intensity has to be higher than zero
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
