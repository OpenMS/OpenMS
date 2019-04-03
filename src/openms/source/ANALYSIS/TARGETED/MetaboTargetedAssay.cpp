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


#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  bool MetaboTargetedAssay::intensityLess_(Peak1D a, Peak1D b)
  {
    return (a.getIntensity() < b.getIntensity());
  }

  // method to extract a potential transistions based on the ms/ms based of the highest intensity precursor or a consensus spectrum
  std::vector <MetaboTargetedAssay> MetaboTargetedAssay::extractMetaboTargetedAssay(const PeakMap& spectra,
                                                                                    const map<BaseFeature const *, vector < size_t>>& feature_ms2_spectra_map,
                                                                                    const double& precursor_rt_tol,
                                                                                    const double& precursor_mz_distance,
                                                                                    const double& cosine_sim_threshold,
                                                                                    const double& transition_threshold,
                                                                                    const bool& method_consensus_spectrum,
                                                                                    const bool& exclude_ms2_precursor,
                                                                                    const unsigned int& file_counter)
  {
    int transition_group_counter = 0;
    vector <MetaboTargetedAssay> v_mta;

    for (auto it = feature_ms2_spectra_map.begin();
              it != feature_ms2_spectra_map.end();
              ++it)
      {
        TargetedExperiment::Compound cmp;
        cmp.clearMetaInfo();

        vector <ReactionMonitoringTransition> v_rmt;

        String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");
        StringList v_description, v_sumformula, v_adduct;

        double feature_rt;
        const BaseFeature *min_distance_feature = it->first;
        feature_rt = min_distance_feature->getRT();

        // extract metadata from featureXML
        if (!(min_distance_feature-> getPeptideIdentifications().empty()) && !(min_distance_feature->getPeptideIdentifications()[0].getHits().empty()))
        {
          // accurate mass search may provide multiple possible Hits
          // for heuristics use the identification with the smallest mz error (ppm)
          double min_id_mz_error = std::numeric_limits<double>::max();
          for (unsigned int j = 0; j != min_distance_feature->getPeptideIdentifications()[0].getHits().size(); ++j)
          {
            double current_id_mz_error = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("mz_error_ppm");
            // compare the absolute error absolute error
            if (abs(current_id_mz_error) < min_id_mz_error)
            {
              description = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("description");
              sumformula = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("chemical_formula");
              adduct = min_distance_feature->getPeptideIdentifications()[0].getHits()[j].getMetaValue("modifications");

              // change format of description [name] to name
              description.erase(remove_if(begin(description), end(description), [](char c) { return c == '[' || c == ']'; }), end(description));

              // change format of adduct information M+H;1+ -> [M+H]1+
              String adduct_prefix = adduct.prefix(';').trim();
              String adduct_suffix = adduct.suffix(';').trim();
              adduct = "[" + adduct_prefix + "]" + adduct_suffix;

              min_id_mz_error = abs(current_id_mz_error);
             }
          }
          // use the identification with the lowest mass deviation
          v_description.push_back(description);
          v_sumformula.push_back(sumformula);
          v_adduct.push_back(adduct);
        }
        else
        {
          // count UNKNOWN via transition group counter
          v_description.push_back(String(description + "_" + transition_group_counter));
          v_sumformula.push_back(sumformula);
          v_adduct.push_back(adduct);
        }

      double highest_precursor_mz = 0.0;
      float highest_precursor_int = 0.0;
      int highest_precursor_charge = 0;
      MSSpectrum highest_precursor_int_spectrum;
      MSSpectrum transition_spectrum;
      String native_id;

      // find precursor/spectrum with highest intensity precursor
      vector <size_t> index = it->second;

      for (auto index_it = index.begin(); index_it != index.end(); ++index_it)
      {
        const MSSpectrum &spectrum = spectra[*index_it];

        // check if MS2 spectrum is empty
        if (spectrum.empty())
        {
          LOG_WARN << "Empty MSMS spectrum was provided. Please hava a look. Index: " << *index_it << std::endl;
          continue;
        }

        const vector <Precursor> &precursor = spectrum.getPrecursors();

        // get m/z and intensity of precursor
        // only spectra with precursors are in the map, therefore no need to check for their presence
        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();
        int precursor_charge = precursor[0].getCharge();

        native_id = spectrum.getNativeID();

        // spectrum with highest intensity precursor
        if (precursor_int > highest_precursor_int)
        {
          highest_precursor_int = precursor_int;
          highest_precursor_mz = precursor_mz;
          highest_precursor_int_spectrum = spectrum;
          highest_precursor_charge = precursor_charge;
        }
          transition_spectrum = highest_precursor_int_spectrum;
      }

      // if only one MS2 is available and the consensus method is used - jump right to the transition list calculation
      // fallback: highest intensity precursor
      if (method_consensus_spectrum &&index.size()>= 2)
      {
        // transform to binned spectra
        vector <BinnedSpectrum> binned;
        vector <MSSpectrum> similar_spectra;
        MSExperiment exp;
        const BinnedSpectrum binned_highest_int(highest_precursor_int_spectrum,
                                                BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,
                                                false,
                                                1,
                                                BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

        // calculation of contrast angle (cosine simiarity)
        for (auto index_it = index.begin(); index_it != index.end(); ++index_it)
        {
          const MSSpectrum &spectrum = spectra[*index_it];

          const BinnedSpectrum binned_spectrum(spectrum,
                                               BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES,
                                               false,
                                               1,
                                               BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES);

          BinnedSpectralContrastAngle bspa;
          double cosine_sim = bspa(binned_highest_int, binned_spectrum);
          if (cosine_sim > cosine_sim_threshold)
          {
            similar_spectra.push_back(spectrum);
            exp.addSpectrum(spectrum);
          }
        }

        // at least 2 spectra with high consine similarity necessary
        // fallback to highest precursor intensity spectrum (see above)
        if (similar_spectra.size()> 1)
        {
          // calculate consensus spectrum
          exp.sortSpectra();

          SpectraMerger merger;
          Param p;
          p.insert("",SpectraMerger().getDefaults());
          p.setValue("precursor_method:mz_tolerance", precursor_mz_distance);
          p.setValue("precursor_method:rt_tolerance", precursor_rt_tol * 2);
          merger.setParameters(p);

          // all MS spectra should have the same precursor
          merger.mergeSpectraPrecursors(exp);

          // check if all precursors have been merged if not use highest intensity precursor
          if (exp.getSpectra().size()< 2)
          {
            transition_spectrum = exp.getSpectra()[0];
          }
        }
      }

      // transition calculations
      // calculate max intensity peak and threshold
      float max_int = 0.0;
      float min_int = std::numeric_limits<float>::max();

      // sort intensity in MS2 spectrum to extract transitions
      transition_spectrum.sortByIntensity(true);

      // filter out the precursors if they are in the ms2 spectrum;
      if (exclude_ms2_precursor)
      {
        for (auto spec_it = transition_spectrum.begin(); spec_it != transition_spectrum.end(); ++spec_it)
        {
          if (transition_spectrum.getPrecursors()[0].getMZ() == spec_it->getMZ())
          {
            transition_spectrum.erase(spec_it);
            break;
          }
        }
      }

      // find max and min intensity peak
      max_int = max_element(transition_spectrum.begin(), transition_spectrum.end(), intensityLess_)->getIntensity();
      min_int = min_element(transition_spectrum.begin(), transition_spectrum.end(), intensityLess_)->getIntensity();

      // no peaks or all peaks have same intensity (single peak / noise)
      if (min_int >= max_int)
      {
        continue;
      }

      vector <TargetedExperimentHelper::RetentionTime> v_cmp_rt;
      TargetedExperimentHelper::RetentionTime cmp_rt;
      cmp_rt.setRT(feature_rt);
      v_cmp_rt.push_back(cmp_rt);
      cmp.rts = v_cmp_rt;

      description = ListUtils::concatenate(v_description, ",");
      cmp.id = String(transition_group_counter) + "_" + description + "_" + file_counter;
      cmp.setMetaValue("CompoundName", description);
      cmp.smiles_string = "NA";

      sumformula = ListUtils::concatenate(v_sumformula, ",");
      cmp.molecular_formula = sumformula;

      adduct = ListUtils::concatenate(v_adduct, ",");
      cmp.setMetaValue("Adducts", adduct);

      // threshold should be at x % of the maximum intensity
      // hard minimal threshold of min_int * 1.1
      float threshold_transition = max_int * (transition_threshold / 100);
      float threshold_noise = min_int * 1.1;

      int transition_counter = 0;
      // here ms2 spectra information is used
      for (auto spec_it = transition_spectrum.begin(); spec_it != transition_spectrum.end(); ++spec_it)
      {
        ReactionMonitoringTransition rmt;
        rmt.clearMetaInfo();

        float current_int = spec_it->getIntensity();
        double current_mz = spec_it->getMZ();

        // write row for each transition
        // current int has to be higher than transition thresold and should not be smaller than threshold noise
        if (current_int > threshold_transition && current_int > threshold_noise)
        {
          float rel_int = current_int / max_int;

          rmt.setPrecursorMZ(highest_precursor_mz);
          rmt.setProductMZ(current_mz);
          rmt.setLibraryIntensity(rel_int);

          description = ListUtils::concatenate(v_description, ",");
          rmt.setCompoundRef (String(transition_group_counter) + "_" + description + "_" + file_counter);
          rmt.setNativeID (String(transition_group_counter)+ "_" + String(transition_counter)+ "_" + description + "_" + file_counter);

          v_rmt.push_back(rmt);
          transition_counter += 1;
        }
      }
      transition_group_counter += 1;
      MetaboTargetedAssay mta;
      mta.precursor_mz = highest_precursor_mz;
      mta.precursor_rt = feature_rt;
      mta.precursor_int = highest_precursor_int;
      mta.precursor_charge = highest_precursor_charge;
      //mta.transition_quality =
      mta.compound_name = description;
      mta.compound_adduct = adduct;
      mta.potential_cmp = cmp;
      mta.potential_rmts = v_rmt;
      v_mta.push_back(mta);
    }
    return v_mta;
  }

  // method to extract a potential transitions based on the ms/ms based of the highest intensity precursor with fragment annotation using SIRIUS
  std::vector <MetaboTargetedAssay> MetaboTargetedAssay::extractMetaboTargetedAssayFragmentAnnotation(const vector <pair<SiriusMSFile::CompoundInfo, MSSpectrum>>& v_cmp_spec,
                                                                                                      const double& transition_threshold,
                                                                                                      const bool& use_exact_mass,
                                                                                                      const bool& exclude_ms2_precursor,
                                                                                                      const unsigned int& file_counter)
  {
    int transition_group_counter = 0;
    vector <MetaboTargetedAssay> v_mta;

    for (auto it = v_cmp_spec.begin();
              it != v_cmp_spec.end();
              ++it)
    {
      // check if annotated spectrum exists
      MSSpectrum transition_spectrum;
      transition_spectrum = it->second;
      if (transition_spectrum.empty())
      {
        continue;
      }

      TargetedExperiment::Compound cmp;
      cmp.clearMetaInfo();
      vector <ReactionMonitoringTransition> v_rmt;

      String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");

      double feature_rt;
      feature_rt = it->first.rt;
      description = it->first.des;
      int charge = it->first.charge;

      // use annotated metadata
      sumformula = it->second.getMetaValue("annotated_sumformula");
      adduct = it->second.getMetaValue("annotated_adduct");

      // transition calculations
      // calculate max intensity peak and threshold
      float max_int = 0.0;
      float min_int = std::numeric_limits<float>::max();

      // sort intensity in MS2 spectrum to extract transitions
      transition_spectrum.sortByIntensity(true);

      // have to remove ms2 precursor peak before min/max
      double exact_mass_precursor = 0.0;
      for (auto spec_it = transition_spectrum.begin();
           spec_it != transition_spectrum.end();
           ++spec_it)
      {
        int spec_index = spec_it - transition_spectrum.begin();
        OpenMS::DataArrays::StringDataArray explanation_array = transition_spectrum.getStringDataArrays()[0];
        if (explanation_array.getName() != "explanation")
        {
          LOG_WARN << "Fragment explanation was not found. Please check if you annotation works properly." << std::endl;
        }
        else
        {
          // precursor in fragment annotation has the same sumformula as MS1 Precursor
          if (explanation_array[spec_index] == sumformula)
          {
            // save exact mass
            if (use_exact_mass)
            {
              exact_mass_precursor = spec_it->getMZ();
            }
            // remove precursor ms2 entry
            if (exclude_ms2_precursor)
            {
              transition_spectrum.erase(transition_spectrum.begin() + spec_index);
              transition_spectrum.getStringDataArrays()[0]
                  .erase(transition_spectrum.getStringDataArrays()[0].begin() + spec_index);
              transition_spectrum.getFloatDataArrays()[0]
                  .erase(transition_spectrum.getFloatDataArrays()[0].begin() + spec_index);
              break; // if last element have to break if not iterator will go out of range
            }
          }
        }
      }

      // find max and min intensity peak
      max_int = max_element(transition_spectrum.begin(), transition_spectrum.end(), intensityLess_)->getIntensity();
      min_int = min_element(transition_spectrum.begin(), transition_spectrum.end(), intensityLess_)->getIntensity();

      // no peaks or all peaks have same intensity (single peak / noise)
      if (min_int >= max_int)
      {
        continue;
      }

      vector <TargetedExperimentHelper::RetentionTime> v_cmp_rt;
      TargetedExperimentHelper::RetentionTime cmp_rt;
      cmp_rt.setRT(feature_rt);
      v_cmp_rt.push_back(cmp_rt);
      cmp.rts = v_cmp_rt;
      if (description == "UNKNOWN")
      {
        description = String(description + "_" + transition_group_counter);
      }
      cmp.id = String(transition_group_counter) + "_" + description + "_" + file_counter;
      cmp.setMetaValue("CompoundName", description);
      cmp.smiles_string = "NA";

      cmp.molecular_formula = sumformula;
      cmp.setMetaValue("Adducts", adduct);

      // threshold should be at x % of the maximum intensity
      // hard minimal threshold of min_int * 1.1
      float threshold_transition = max_int * (transition_threshold / 100);
      float threshold_noise = min_int * 1.1;
      int transition_counter = 0;

      // extract current StringDataArry with annotations/explanations;
      OpenMS::DataArrays::StringDataArray explanation_array = transition_spectrum.getStringDataArrays()[0];

      // check which entry is saved in the FloatDataArry (control for "use_exact_mass")
      LOG_DEBUG << transition_spectrum.getFloatDataArrays()[0].getName() << " is not used to build the assay library."
                << std::endl;

      // here ms2 spectra information is used
      for (auto spec_it = transition_spectrum.begin();
           spec_it != transition_spectrum.end();
           ++spec_it)
      {
        ReactionMonitoringTransition rmt;
        rmt.clearMetaInfo();
        int peak_index = spec_it - transition_spectrum.begin();

        float current_int = spec_it->getIntensity();
        double current_mz = spec_it->getMZ();
        String
        current_explanation = explanation_array[peak_index];

        // write row for each transition
        // current int has to be higher than transition thresold and should not be smaller than threshold noise
        if (current_int > threshold_transition && current_int > threshold_noise)
        {
          float rel_int = current_int / max_int;

          rmt.setPrecursorMZ((use_exact_mass && exact_mass_precursor != 0.0) ? exact_mass_precursor : it->first.pmass);
          rmt.setProductMZ(current_mz);
          rmt.setLibraryIntensity(rel_int);

          rmt.setCompoundRef(String(transition_group_counter) + "_" + description + "_" + file_counter);
          rmt.setNativeID(String(transition_group_counter) + "_" + String(transition_counter) + "_" + description + "_" +
                          file_counter);

          rmt.setMetaValue("annotation", DataValue(current_explanation));

          v_rmt.push_back(rmt);
          transition_counter += 1;
        }
      }

      transition_group_counter += 1;
      MetaboTargetedAssay mta;
      mta.precursor_mz = (use_exact_mass && exact_mass_precursor != 0.0) ? exact_mass_precursor : it->first.pmass;
      mta.precursor_rt = it->first.rt;
      mta.precursor_int = 0;
      mta.precursor_charge = charge;
      //mta.transition_quality = ;
      mta.compound_name = description;
      mta.compound_adduct = adduct;
      mta.potential_cmp = cmp;
      mta.potential_rmts = v_rmt;
      v_mta.push_back(mta);
    }
    return v_mta;
  }

} // namespace OpenMS

/// @endcond



