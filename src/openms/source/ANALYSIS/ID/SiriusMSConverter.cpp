// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
#include <cstdint>
#include <QDir>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  // precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz,
                                       const MSSpectrum& spectrum1,
                                       double tolerance,
                                       bool ppm)
  {

    // get tolerance window and left/right iterator
    pair<double, double> tolerance_window = Math::getTolWindow(test_mz, tolerance, ppm);

    // Here left has to be smaller than right
    OPENMS_PRECONDITION(tolerance_window.first < tolerance_window.second, "Left has to be smaller than right");

    MSSpectrum::ConstIterator left = spectrum1.MZBegin(tolerance_window.first);
    MSSpectrum::ConstIterator right = spectrum1.MZBegin(tolerance_window.second);

    // no MS1 precursor peak in +- tolerance window found
    if (left == right)
    {
      return -1;
    }

    MSSpectrum::ConstIterator max_intensity_it = max_element(left, right, Peak1D::IntensityLess());

    return max_intensity_it - spectrum1.begin();
  }

  // extract precursor isotope pattern if no feature information is available
  vector<Peak1D> extractPrecursorIsotopePattern(const double& precursor_mz,
                                                const MSSpectrum& precursor_spectrum,
                                                int& iterations,
                                                const int& charge)
  {
    vector<Peak1D> isotopes;
    int peak_index;
    Peak1D peak;

    // monoisotopic_trace
    peak_index = getHighestIntensityPeakInMZRange(precursor_mz, precursor_spectrum, 10, true);
    if (peak_index != -1)
    {
      peak = precursor_spectrum[peak_index];
      isotopes.push_back(peak);
    }

    // further isotope_traces with the mass error of 1 ppm
    double massdiff = Constants::C13C12_MASSDIFF_U;

    // depending on the charge different MASSDIFF
    if (charge != 0)
    {
      massdiff = massdiff/std::abs(charge);
    }

    while (peak_index != -1 && iterations > 0)
    {
      // check for isotope trace with one ppm error
      peak_index = getHighestIntensityPeakInMZRange(peak.getMZ() + massdiff, precursor_spectrum, 1, true);
      if (peak_index != -1)
      {
        peak = precursor_spectrum[peak_index];
        isotopes.push_back(peak);
      }
      --iterations;
    }
    return isotopes;
  }

  void SiriusMSFile::writeMsFile_(ofstream& os,
                                  const MSExperiment& spectra,
                                  const vector<size_t>& ms2_spectra_index,
                                  const SiriusMSFile::AccessionInfo& ainfo,
                                  const StringList& adducts,
                                  const vector<String>& v_description,
                                  const vector<String>& v_sumformula,
                                  const vector<pair<double,double>>& f_isotopes,
                                  int& feature_charge,
                                  uint64_t& feature_id,
                                  const double& feature_rt,
                                  const double& feature_mz,
                                  bool& writecompound,
                                  const bool& no_masstrace_info_isotope_pattern,
                                  const int& isotope_pattern_iterations,
                                  int& count_skipped_spectra,
                                  int& count_assume_mono,
                                  int& count_no_ms1,
                                  std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo)
  {
    // if multiple identifications present for one MS1 and MS2 use all of them and
    // let SIRIUS sort it out using fragment annotation
    for (unsigned int k = 0; k != v_description.size(); ++k)
    {
      if (v_description.size() > 1) { writecompound = true; }
      SiriusMSFile::CompoundInfo cmpinfo;

      for (const size_t& ind : ms2_spectra_index)
      {
        // construct compound info structure
        const MSSpectrum &current_ms2 = spectra[ind];
        const double current_rt = current_ms2.getRT();

        const String native_id = current_ms2.getNativeID();
        const int scan_number = SpectrumLookup::extractScanNumber(native_id, ainfo.native_id_accession);

        const vector<Precursor> &precursor = current_ms2.getPrecursors();

        // get m/z and intensity of precursor
        if (precursor.empty())
        {
          throw Exception::MissingInformation(__FILE__,
                                              __LINE__,
                                              OPENMS_PRETTY_FUNCTION,
                                              "Precursor for MS/MS spectrum was not found.");
        }

        IonSource::Polarity p = current_ms2.getInstrumentSettings().getPolarity(); //charge

        // there should be only one precursor and MS2 should contain peaks to be considered
        if (precursor.size() == 1 && !current_ms2.empty())
        {

          // read precursor charge
          int precursor_charge = precursor[0].getCharge();

          // sirius supports only single charged ions (+1; -1)
          // if charge = 0, it will be changed to +1; -1 depending on Polarity
          if (precursor_charge > 1 || precursor_charge < -1)
          {
            ++count_skipped_spectra;
            continue;
          }

          // set precursor charge for msfile
          // no charge annotated - assume mono-charged
          if (precursor_charge == 0)
          {
            precursor_charge = 1;
            ++count_assume_mono;
          }
          // negative mode - make sure charges are < 0
          if (p == IonSource::Polarity::NEGATIVE) { precursor_charge = -(std::abs(precursor_charge)); }

          // set feature_charge for msfile if feature information is available
          // no charge annotated - assume mono-charged
          if (feature_id != 0 && feature_charge == 0)
          {
            feature_charge = 1;
            ++count_assume_mono;
          }
          // negative mode - make sure charges are < 0
          if (p == IonSource::Polarity::NEGATIVE) { feature_charge = -(std::abs(feature_charge)); }

          // get m/z and intensity of precursor != MS1 spectrum
          double precursor_mz = precursor[0].getMZ();
          float precursor_int = precursor[0].getIntensity();

          // extract collision energy
          double collision = precursor[0].getActivationEnergy();

          // find corresponding ms1 spectra (precursor)
          PeakMap::ConstIterator s_it2 = spectra.getPrecursorSpectrum((spectra.begin() + ind));

          double test_mz = precursor_mz;
          double precursor_rt = 0.0;

          vector<Peak1D> isotopes;
          isotopes.clear();
          vector<Peak1D> precursor_spec;

          // getPrecursorSpectrum returns past-the-end iterator if spectrum is not found.
          if (s_it2 == spectra.end() || s_it2->getMSLevel() != 1)
          {
            ++count_no_ms1;
          }
          // get the precursor in the ms1 spectrum (highest intensity in the range of the precursor mz +- 0.1 Da)
          else
          {
            const MSSpectrum &precursor_spectrum = *s_it2;
            precursor_rt = precursor_spectrum.getRT();
            int interations = isotope_pattern_iterations;
            // extract precursor isotope pattern via C13 isotope distance
            if (feature_id != 0 && feature_charge != 0)
            {
              isotopes = extractPrecursorIsotopePattern(test_mz, precursor_spectrum, interations, feature_charge);
            }
            else
            {
              isotopes = extractPrecursorIsotopePattern(test_mz, precursor_spectrum, interations, precursor_charge);
            }
            for (Size i = 0; i < precursor_spectrum.size(); ++i)
            {
              const Peak1D &peak = precursor_spectrum[i];
              precursor_spec.push_back(peak);
            }
          }

          // construct query_id; remove spaces from string
          // use first
          std::string des_wo_space = v_description[k];
          des_wo_space.erase(std::remove_if(des_wo_space.begin(), des_wo_space.end(), ::isspace), des_wo_space.end());

          String query_id = "_" + String(feature_id) +
                            String("-" + String(scan_number) + "-") +
                            String(des_wo_space) +
                            String(ind);

          if (writecompound)
          {
            // write internal unique .ms data as sirius input
            os << fixed;
            os << ">compound " << query_id << "\n";
            cmpinfo.cmp = query_id;

            if (!f_isotopes.empty() && !no_masstrace_info_isotope_pattern)
            {
              os << ">parentmass " << f_isotopes[0].first << fixed << "\n";
              cmpinfo.pmass = f_isotopes[0].first;
            }
            else if (!isotopes.empty())
            {
              os << ">parentmass " << isotopes[0].getMZ() << fixed << "\n";
              cmpinfo.pmass = isotopes[0].getMZ();
            }
            else
            {
              os << ">parentmass " << precursor_mz << fixed << "\n";
              cmpinfo.pmass = precursor_mz;
            }

            if (!adducts.empty())
            {
              os << ">ionization " << adducts[k] << "\n";
              cmpinfo.ionization = adducts[k];
            }

            // first
            if (v_sumformula[k] != "UNKNOWN")
            {
              os << ">formula " << v_sumformula[k] << "\n";
              cmpinfo.formula = v_sumformula[k];
            }

            if (feature_charge != 0)
            {
              os << ">charge " << feature_charge << "\n";
              cmpinfo.charge = feature_charge;
            }
            else
            {
              os << ">charge " << precursor_charge << "\n";
              cmpinfo.charge = precursor_charge;
            }

            if (feature_rt != 0)
            {
              os << ">rt " << feature_rt << "\n";
              cmpinfo.rt = feature_rt;
            }
            else if (precursor_rt != 0.0)
            {
              os << ">rt " << precursor_rt << "\n";
              cmpinfo.rt = precursor_rt;
            }
            else
            {
              os << ">rt " << current_rt << "\n";
              cmpinfo.rt = current_rt;
            }

            if (feature_mz != 0 && feature_id != 0)
            {
              os << "##fmz " << String(feature_mz) << "\n";
              os << "##fid " << String(feature_id) << "\n";
              cmpinfo.fmz = feature_mz;
              cmpinfo.fid = feature_id;
            }
            os << "##des " << String(v_description[k]) << "\n";
            os << "##specref_format " << "[MS, " << ainfo.native_id_accession <<", "<< ainfo.native_id_type << "]" << endl;
            os << "##source file " << ainfo.sf_path << endl;
            os << "##source format " << "[MS, " << ainfo.sf_accession << ", "<< ainfo.sf_type << ",]" << endl;
            cmpinfo.des = String(v_description[k]);
            cmpinfo.specref_format = String("[MS, " + ainfo.native_id_accession + ", " + ainfo.native_id_type + "]");
            cmpinfo.source_file = ainfo.sf_path;
            cmpinfo.source_format = String("[MS, " + ainfo.sf_accession + ", "+ ainfo.sf_type + ",]" );

            // use precursor m/z & int and no ms1 spectra is available else use values from ms1 spectrum
            Size num_isotopes = isotopes.size();
            Size num_f_isotopes = f_isotopes.size();

            if (num_f_isotopes > 0 && !no_masstrace_info_isotope_pattern)
            {
              os << ">ms1merged" << endl;
              // m/z and intensity have to be higher than 1e-10
              for (auto it = f_isotopes.begin(); it != f_isotopes.end(); ++it)
              {
                os << it->first << " " << it->second << "\n";
              }
              cmpinfo.pint_mono = f_isotopes[0].second;
            }
            else if (num_isotopes > 0) // if ms1 spectrum was present
            {
              os << ">ms1merged" << endl;
              for (auto it = isotopes.begin(); it != isotopes.end(); ++it)
              {
                os << it->getMZ() << " " << it->getIntensity() << "\n";
              }
              cmpinfo.pint_mono = isotopes[0].getIntensity();
            }
            else
            {
              if (precursor_int != 0) // if no ms1 spectrum was present but precursor intensity is known
              {
                os << ">ms1merged" << "\n" << precursor_mz << " " << precursor_int << "\n\n";
                cmpinfo.pint_mono = precursor_int;
              }
            }
          }

          // if a feature_id is present compound should only be written once
          // since every ms2 belongs to the same feature with a different description
          if (feature_id != 0)
          {
            writecompound = false;
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
          os << "##nid " << native_id<< endl;
          // "mid" annotation for multiple possible identifications (native_id_k)
          // fragment mapping will be done using the mid
          String mid = native_id + "_" + k;
          os << "##mid " << mid << endl;
          os << "##scan " << ind << endl;
          os << "##specref " << "ms_run[1]:" << native_id << endl;

          cmpinfo.native_ids.push_back(native_id);
          cmpinfo.mids.push_back(mid);
          cmpinfo.scan_indices.push_back(ind);
          cmpinfo.specrefs.push_back(String("ms_run[1]:" + native_id));

          // single spectrum peaks
          for (Size i = 0; i < current_ms2.size(); ++i)
          {
            const Peak1D &peak = current_ms2[i];
            double mz = peak.getMZ();
            float intensity = peak.getIntensity();

            // intensity has to be higher than zero
            if (intensity != 0)
            {
              os << mz << " " << intensity << "\n";
            }
          }
        }
      }
      v_cmpinfo.push_back(std::move(cmpinfo));
    }
  }

  void SiriusMSFile::store(const MSExperiment& spectra,
                           const OpenMS::String& msfile,
                           const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                           const bool& feature_only,
                           const int& isotope_pattern_iterations,
                           const bool no_masstrace_info_isotope_pattern,
                           std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo)
  {
    const Map<const BaseFeature*, vector<size_t>>& assigned_ms2 = feature_mapping.assignedMS2;
    const vector<size_t> & unassigned_ms2 = feature_mapping.unassignedMS2;

    bool use_feature_information = false;
    bool use_unassigend_ms2 = false;
    bool no_feautre_information = false;

    // Three different possible .ms formats
    // feature information is used (adduct, masstrace_information (FFM+MAD || FFM+AMS || FMM+MAD+AMS [AMS preferred])
    if (!assigned_ms2.empty()) use_feature_information = true;
    // feature information was provided and unassigend ms2 should be used (feature only parameter)
    if (!unassigned_ms2.empty() && !feature_only) use_unassigend_ms2 = true;
    // no feature information was provided (mzml input only)
    if (assigned_ms2.empty() && unassigned_ms2.empty()) no_feautre_information = true;

    int count_skipped_spectra = 0; // spectra skipped due to precursor charge
    int count_assume_mono = 0; // count if mono charge was assumend and set to current ion mode
    int count_no_ms1 = 0; // count if no precursor was found
    int count_skipped_features = 0; // features skipped due to charge

    // check for all spectra at the beginning if spectra are centroided
    // determine type of spectral data (profile or centroided) - only checking first spectrum (could be ms2 spectrum)
    SpectrumSettings::SpectrumType spectrum_type = spectra[0].getType();

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

    AccessionInfo ainfo;

    // sourcefile 
    ainfo.sf_path = spectra.getSourceFiles()[0].getPathToFile();
    ainfo.sf_type = spectra.getSourceFiles()[0].getFileType();
 
    // extract accession by name
    std::set<String> terms;
    ControlledVocabulary cv;
    cv.loadFromOBO("MS", File::find("/CV/psi-ms.obo"));
    cv.getAllChildTerms(terms, "MS:1000560");
    for (std::set<String>::const_iterator it = terms.begin(); it != terms.end(); ++it)
    {
      if (cv.getTerm(*it).name == ainfo.sf_type)
      {
          cv.getTerm(*it);
          ainfo.sf_accession = cv.getTerm(*it).id;
      }
    }  
    // native_id
    ainfo.native_id_accession = spectra.getSourceFiles()[0].getNativeIDTypeAccession();
    ainfo.native_id_type = spectra.getSourceFiles()[0].getNativeIDType();

    vector<String> adducts;
    String description;
    String sumformula;
    vector<String> v_description;
    vector<String> v_sumformula;

    uint64_t feature_id;
    int feature_charge;
    double feature_rt;
    double feature_mz;
    vector<pair<double, double>> f_isotopes;
    f_isotopes.clear();

    // if feature information is available to this first (write features in one compound)
    if (use_feature_information)
    { 
      for (auto it = assigned_ms2.begin();
                it != assigned_ms2.end();
                ++it)
      {
        const BaseFeature* feature = it->first;
        const vector<size_t> feature_associated_ms2 = it->second;
        
        // reset feature information with each iteration
        f_isotopes.clear();

        feature_id = feature->getUniqueId();
        feature_charge = feature->getCharge();
        feature_rt = feature->getRT();
        feature_mz = feature->getMZ();

        // multiple charged compounds are not allowed in sirius
        if (feature_charge > 1 || feature_charge < -1)
        {
          ++count_skipped_features;
          continue;
        }

        // ffm featureXML
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
                pair<double, double> masstrace_mz_int(masstrace_centroid_mz[i],masstrace_intensity[i]);
                f_isotopes.push_back(masstrace_mz_int);
              }
            }
        }

        // prefer adducts from AccurateMassSearch if MetaboliteAdductDecharger and AccurateMassSearch were performed
        // if multiple PeptideHits / identifications occur - use all for SIRIUS
        v_description.clear();
        v_sumformula.clear();
        if (!feature->getPeptideIdentifications().empty() && !feature->getPeptideIdentifications()[0].getHits().empty())
        {
          adducts.clear();

          for (unsigned int j = 0; j != feature->getPeptideIdentifications()[0].getHits().size(); ++j)
          {
           String adduct;
           description = feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("description");
           sumformula = feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("chemical_formula");
           adduct = feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("modifications");

           // change format of description [name] to name
           description.erase(remove_if(begin(description),
                                       end(description),
                                       [](char c) { return c == '[' || c == ']'; }), end(description));

           // change format of adduct information M+H;1+ -> [M+H]1+
           String adduct_prefix = adduct.prefix(';').trim();
           String adduct_suffix = adduct.suffix(';').trim();
           adduct = "[" + adduct_prefix + "]" + adduct_suffix;
           adducts.insert(adducts.begin(), adduct);
           v_description.push_back(description);
           v_sumformula.push_back(sumformula);
          }
        }
        else
        {
          // initialization with UNKNOWN in case no feature information is available.
          v_description.push_back("UNKNOWN");
          v_sumformula.push_back("UNKNOWN");
        }

        bool writecompound = true;
        // call function to writeMsFile to os
        writeMsFile_(os,
                    spectra,
                    feature_associated_ms2,
                    ainfo,
                    adducts,
                    v_description,
                    v_sumformula,
                    f_isotopes,
                    feature_charge,
                    feature_id,
                    feature_rt,
                    feature_mz,
                    writecompound,
                    no_masstrace_info_isotope_pattern,
                    isotope_pattern_iterations,
                    count_skipped_spectra,
                    count_assume_mono,
                    count_no_ms1,
                    v_cmpinfo);
        }
    }

    // if not mappend information available (e.g. empty featurexml or only a few features)
    if (use_unassigend_ms2)
    {
      // no feature information was provided
      bool writecompound = true;
      f_isotopes.clear();
      adducts.clear();
      feature_charge = 0;
      feature_id = 0;
      feature_mz = 0;
      feature_rt = 0;

      writeMsFile_(os,
                   spectra,
                   unassigned_ms2,
                   ainfo,
                   adducts,
                   v_description,
                   v_sumformula,
                   f_isotopes,
                   feature_charge,
                   feature_id,
                   feature_rt,
                   feature_mz,
                   writecompound,
                   no_masstrace_info_isotope_pattern,
                   isotope_pattern_iterations,
                   count_skipped_spectra,
                   count_assume_mono,
                   count_no_ms1,
                   v_cmpinfo);
    }

    if (no_feautre_information)
    {
      // no feature information was provided
      bool writecompound = true;
      v_description.push_back("UNKNOWN");
      v_sumformula.push_back("UNKNOWN");
      f_isotopes.clear();
      adducts.clear();
      feature_charge = 0;
      feature_id = 0;
      feature_mz = 0;
      feature_rt = 0;

      // fill vector with index of all ms2 of the mzml
      vector<size_t> all_ms2;

      for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
      {
        // process only MS2 spectra
        if (s_it->getMSLevel() != 2)
        {
          continue;
        }

        int scan_index = s_it - spectra.begin();

        all_ms2.push_back(scan_index);
      }

      writeMsFile_(os,
                   spectra,
                   all_ms2,
                   ainfo,
                   adducts,
                   v_description,
                   v_sumformula,
                   f_isotopes,
                   feature_charge,
                   feature_id,
                   feature_rt,
                   feature_mz,
                   writecompound,
                   no_masstrace_info_isotope_pattern,
                   isotope_pattern_iterations,
                   count_skipped_spectra,
                   count_assume_mono,
                   count_no_ms1,
                   v_cmpinfo);
    }

    os.close();

    LOG_WARN << "No MS1 spectrum for this precursor. Occurred " << count_no_ms1 << " times." << endl;
    LOG_WARN << count_skipped_spectra << " spectra were skipped due to precursor charge below -1 and above +1." << endl;
    LOG_WARN << "Mono charge assumed and set to charge 1 with respect to current polarity " << count_assume_mono << " times."<< endl;
    LOG_WARN << count_skipped_features << " features were skipped due to feature charge below -1 and above +1." << endl;

  }
} // namespace OpenMS

/// @endcond
