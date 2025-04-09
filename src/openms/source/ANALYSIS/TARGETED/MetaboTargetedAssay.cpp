// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <regex>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

  // get charge from adduct in standard format [M+H]+ or [M+H]1+
  // only for singly charged species
  int MetaboTargetedAssay::getChargeFromAdduct_(const String& adduct)
  {
    int adduct_charge;
    String adduct_suffix = adduct.suffix(']').trim();
    // charge one
    if (adduct_suffix == "+")
    {
      adduct_suffix = "1" + adduct_suffix;
    }
    else if (adduct_suffix == "-")
    {
      adduct_suffix = "1" + adduct_suffix;
    }
    else if (adduct_suffix != "1-" && adduct_suffix != "1+")
    {
      OpenMS_Log_warn << "The adduct had the suffix '" << adduct_suffix << "', but only singly positive or singly negative charged adducts are supported." << std::endl;
    }
    String sign = adduct.back();
    adduct_suffix.resize(adduct_suffix.size()-1);
    if (sign == "+")
    {
      adduct_charge = String(adduct_suffix).toInt();
    }
    else
    {
      adduct_charge = String(sign + adduct_suffix).toInt();
    }
    return adduct_charge;
  }

  bool MetaboTargetedAssay::intensityLess_(const Peak1D& a, const Peak1D& b)
  {
    return (a.getIntensity() < b.getIntensity());
  }

  void MetaboTargetedAssay::filterBasedOnTotalOccurrence_(std::vector<MetaboTargetedAssay>& mta, double total_occurrence_filter, size_t in_files_size)
  {
    if (in_files_size > 1 && !mta.empty())
    {
      double total_occurrence = double(mta.size())/double(in_files_size);
      if (!(total_occurrence >= total_occurrence_filter))
      {
        mta.clear(); // return empty vector
      }
    }
  }

  void MetaboTargetedAssay::sortByPrecursorInt(std::vector<MetaboTargetedAssay>& vec_mta)
  {
    sort(vec_mta.begin(),
         vec_mta.end(),
         [](const MetaboTargetedAssay &a, const MetaboTargetedAssay &b) -> bool
         {
           return a.precursor_int > b.precursor_int;
         });
  }

  void MetaboTargetedAssay::filterBasedOnMolFormAdductOccurrence_(std::vector<MetaboTargetedAssay>& mta)
  {
    std::map<std::pair<String, String>, int> occ_map;
    if (!mta.empty())
    {
      for (const auto &t_it : mta)
      {
        auto [it, success] = occ_map.emplace(std::make_pair(t_it.molecular_formula, t_it.compound_adduct), 1);
        if (!success)
        {
          it->second++;
        }
      }

      // find max element in map
      using pair_type = decltype(occ_map)::value_type;
      auto pr = std::max_element(std::begin(occ_map),
                                 std::end(occ_map),
                                 [](const pair_type &p1, const pair_type &p2) { return p1.second < p2.second; });

      // filter vector down to the compound with mol. formula and adduct based on the highest occurrence
      mta.erase(remove_if(mta.begin(),
                          mta.end(),
                          [&pr](const MetaboTargetedAssay& assay)
                          {
                            return assay.molecular_formula != pr->first.first ||
                                   assay.compound_adduct != pr->first.second;
                          }), mta.end());
    }
  }

  // method to extract potential transitions based on the ms/ms of the highest intensity precursor or a consensus spectrum
  std::vector <MetaboTargetedAssay> MetaboTargetedAssay::extractMetaboTargetedAssay(const MSExperiment& spectra,
                                                                                    const FeatureMapping::FeatureToMs2Indices& feature_ms2_index,
                                                                                    const double& precursor_rt_tol,
                                                                                    const double& precursor_mz_distance,
                                                                                    const double& cosine_sim_threshold,
                                                                                    const double& transition_threshold,
                                                                                    const double& min_fragment_mz,
                                                                                    const double& max_fragment_mz,
                                                                                    const bool& method_consensus_spectrum,
                                                                                    const bool& exclude_ms2_precursor,
                                                                                    const unsigned int& file_counter)
  {
    int transition_group_counter = 0;
    vector <MetaboTargetedAssay> v_mta;
    const std::map<BaseFeature const *, vector < size_t>>& feature_ms2_spectra_map = feature_ms2_index.assignedMS2;

    for (const auto& it : feature_ms2_spectra_map)
      {
        TargetedExperiment::Compound cmp;
        cmp.clearMetaInfo();

        vector <ReactionMonitoringTransition> v_rmt;

        String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");
        StringList v_description, v_sumformula, v_adduct;

        double feature_rt;
        int feature_charge;
        const BaseFeature *min_distance_feature = it.first;
        feature_rt = min_distance_feature->getRT();
        feature_charge = min_distance_feature->getCharge();

        // extract metadata from featureXML
        auto metaboliteIdentifications = min_distance_feature->getPeptideIdentifications();
        if (!(metaboliteIdentifications.empty()) && !(metaboliteIdentifications[0].getHits().empty()))
        {
          // accurate mass search may provide multiple possible Hits
          // for heuristics use the identification with the smallest mz error (ppm)
          double min_id_mz_error = std::numeric_limits<double>::max();
          for (const auto& mhit : min_distance_feature->getPeptideIdentifications()[0].getHits())
          {
            double current_id_mz_error = mhit.getMetaValue("mz_error_ppm");
            // compare the absolute error absolute error
            if (abs(current_id_mz_error) < min_id_mz_error)
            {
              description = mhit.getMetaValue("description");
              sumformula = mhit.getMetaValue("chemical_formula");
              adduct = mhit.getMetaValue("modifications");

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
      vector <size_t> index = it.second;

      for (auto index_it = index.begin(); index_it != index.end(); ++index_it)
      {
        const MSSpectrum &spectrum = spectra[*index_it];

        // check if MS2 spectrum is empty
        if (spectrum.empty())
        {
          OPENMS_LOG_WARN << "Empty MS/MS spectrum was provided. Please manually investigate at index: " << *index_it << std::endl;
          continue;
        }

        const vector <Precursor> &precursor = spectrum.getPrecursors();

        // get m/z and intensity of precursor
        if (precursor.empty())
        {
          throw Exception::MissingInformation(__FILE__,
                                              __LINE__,
                                              OPENMS_PRETTY_FUNCTION,
                                              "Precursor for MS/MS spectrum was not found.");
        }

        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();
        int precursor_charge = precursor[0].getCharge();

        // if precursor charge is not annotated - use feature charge
        if (precursor_charge == 0)
        {
          precursor_charge = feature_charge;
        }

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

      // check if more than one MS/MS spectrum is available to use the consensus method
      // the MS/MS spectrum of the highest intensity precursor is used as reference and compared
      // to the other MS/MS spectrum of a specific feature.
      // if the cosine similarity is over the manually set threshold these are merged via SpectraMerger
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

        // calculation of contrast angle (cosine similarity)
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
          if (abs(transition_spectrum.getPrecursors()[0].getMZ() - spec_it->getMZ()) < 1e-3)
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
      v_cmp_rt.push_back(std::move(cmp_rt));
      cmp.rts = std::move(v_cmp_rt);

      cmp.setChargeState(highest_precursor_charge);

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
      float threshold_noise = min_int * noise_threshold_constant_;

      int transition_counter = 0;
      // here ms2 spectra information is used
      for (auto spec_it = transition_spectrum.begin(); spec_it != transition_spectrum.end(); ++spec_it)
      {
        ReactionMonitoringTransition rmt;
        rmt.clearMetaInfo();

        float current_int = spec_it->getIntensity();
        double current_mz = spec_it->getMZ();

        // write row for each transition
        // current int has to be higher than transition threshold and should not be smaller than threshold noise
        // current_mz has to be higher than min_fragment_mz and lower than max_fragment_mz
        if (current_int > threshold_transition && current_int > threshold_noise && current_mz > min_fragment_mz && current_mz < max_fragment_mz)
        {
          float rel_int = current_int / max_int;

          rmt.setPrecursorMZ(highest_precursor_mz);
          rmt.setProductMZ(current_mz);
          TargetedExperimentHelper::TraMLProduct product;
          product.setMZ(current_mz);
          // charge state from adduct
          if (!adduct.empty() && adduct != "UNKNOWN")
          {
            product.setChargeState(getChargeFromAdduct_(adduct));
          }
          rmt.setProduct(product);
          rmt.setLibraryIntensity(rel_int);
          description = ListUtils::concatenate(v_description, ",");
          rmt.setCompoundRef (String(transition_group_counter) + "_" + description + "_" + file_counter);
          rmt.setNativeID (String(transition_group_counter)+ "_" + String(transition_counter)+ "_" + description + "_" + file_counter);
          rmt.setDecoyTransitionType(ReactionMonitoringTransition::DecoyTransitionType::TARGET); // no decoys are generated without SIRIUS

          v_rmt.push_back(std::move(rmt));
          transition_counter += 1;
        }
      }
      transition_group_counter += 1;
      MetaboTargetedAssay mta;
      mta.precursor_int = highest_precursor_int;
      mta.compound_name = description;
      mta.compound_adduct = adduct;
      mta.precursor_mz = highest_precursor_mz;
      mta.molecular_formula = sumformula;
      mta.compound_rt = feature_rt;
      mta.compound_file = file_counter;
      mta.potential_cmp = cmp;
      mta.potential_rmts = v_rmt;

      // do not report if no valid transitions are found after filtering
      if (!mta.potential_rmts.empty())
      {
        v_mta.push_back(std::move(mta));
      }
    }
    return v_mta;
  }

  // method to extract potential transitions based on the ms/ms based of the highest intensity precursor with fragment annotation using SIRIUS
  std::vector <MetaboTargetedAssay> MetaboTargetedAssay::extractMetaboTargetedAssayFragmentAnnotation(const vector < CompoundTargetDecoyPair >& v_cmp_spec,
                                                                                                      const double& transition_threshold,
                                                                                                      const double& min_fragment_mz,
                                                                                                      const double& max_fragment_mz,
                                                                                                      const bool& use_exact_mass,
                                                                                                      const bool& exclude_ms2_precursor)
  {
    int entry_counter = 0; // counts each entry - to ensure the same count for targets, decoys from the same sirius workspace
    vector <MetaboTargetedAssay> v_mta;

    for (const auto& it : v_cmp_spec)
    {
      // check if annotated object exists
      const MetaboTargetedAssay::CompoundTargetDecoyPair &csp = it;

      vector<MSSpectrum> non_empty_spectra;
      if (!csp.target_decoy_spectra.target.empty())
      {
        MSSpectrum target = csp.target_decoy_spectra.target;
        non_empty_spectra.push_back(target);
      }
      if (!csp.target_decoy_spectra.decoy.empty())
      {
        MSSpectrum decoy = csp.target_decoy_spectra.decoy;
        non_empty_spectra.push_back(decoy);
      }

      // iterate over both entries - targets and decoys (SiriusTargetDecoySpectra)
      // count target and decoy as one entry - to ensure same numbering of targets and decoys
      for (auto& transition_spectrum : non_empty_spectra)
      {
        TargetedExperiment::Compound cmp;
        cmp.clearMetaInfo();
        vector <ReactionMonitoringTransition> v_rmt;

        String description("UNKNOWN"), sumformula("UNKNOWN"), adduct("UNKNOWN");

        double feature_rt = csp.compound_info.rt;
        description = csp.compound_info.des;
        int charge = csp.compound_info.charge;
        double precursor_int = csp.compound_info.pint_mono;

        // use annotated metadata
        sumformula = transition_spectrum.getMetaValue("annotated_sumformula");
        adduct = transition_spectrum.getMetaValue("annotated_adduct");
        int decoy = transition_spectrum.getMetaValue("decoy");

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

          OpenMS::DataArrays::StringDataArray explanation_array;
          if (!transition_spectrum.getStringDataArrays().empty())
          {
            explanation_array = transition_spectrum.getStringDataArrays()[0];
            if (explanation_array.getName() != "explanation")
            {
              OPENMS_LOG_WARN << "Fragment explanation was not found. Please check if your annotation works properly." << std::endl;
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
                  if (decoy == 0) // second mass FloatDataArray only available for targets
                  {
                    transition_spectrum.getFloatDataArrays()[0]
                      .erase(transition_spectrum.getFloatDataArrays()[0].begin() + spec_index);
                  }
                  break; // break to not increment when erase es called.
                }
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
          OPENMS_LOG_DEBUG << "The annotated spectrum does not have any peaks after the intensity filter step, or all peaks have the same intensity: " << csp.compound_info.cmp << std::endl;
          continue;
        }

        vector <TargetedExperimentHelper::RetentionTime> v_cmp_rt;
        TargetedExperimentHelper::RetentionTime cmp_rt;
        cmp_rt.setRT(feature_rt);
        v_cmp_rt = {cmp_rt};
        cmp.rts = {v_cmp_rt};
        cmp.setChargeState(charge);
        String identifier_suffix = adduct + "_" + int(feature_rt) + "_" + csp.compound_info.file_index;

        if (description == "UNKNOWN")
        {
          description = String(description + "_" + entry_counter);
        }
        // compoundID has to be unique over all the files
        // feature_rt if the same ID was detected twice at different retention times in the same file
        if (decoy == 0)
        {
          cmp.id = String(entry_counter) + "_" + description + "_" + identifier_suffix;
          cmp.setMetaValue("CompoundName", description);
        }
        else if (decoy == 1)
        {
          description = String(description + "_decoy");
          cmp.id = String(entry_counter) + "_" + description + "_" + identifier_suffix;
          cmp.setMetaValue("CompoundName", description);
        }

        OPENMS_LOG_DEBUG << "Processed annotated Spectra - mapping of the description and the SIRIUS identifier." << std::endl;
        OPENMS_LOG_DEBUG << "Description: " << description << std::endl;
        OPENMS_LOG_DEBUG << "SIRIUS_workspace_identifier: " << csp.compound_info.cmp << std::endl;
        OPENMS_LOG_DEBUG << "Compound identifier: " << cmp.id << std::endl;

        cmp.smiles_string = "NA";
        cmp.molecular_formula = sumformula;
        cmp.setMetaValue("Adducts", adduct);
        cmp.setMetaValue("decoy", decoy);
        if (!csp.compound_info.native_ids_id.empty())
        {
          cmp.setMetaValue("native_ids_id", csp.compound_info.native_ids_id);
        }
        if (!csp.compound_info.m_ids_id.empty())
        {
          cmp.setMetaValue("m_ids_id", csp.compound_info.m_ids_id);
        }
        if (!csp.compound_info.cmp.empty())
        {
          cmp.setMetaValue("sirius_workspace_identifier", csp.compound_info.cmp);
        }

        // threshold should be at x % of the maximum intensity
        // hard minimal threshold of min_int * 1.1
        float threshold_transition = max_int * (transition_threshold / 100);
        float threshold_noise = min_int * noise_threshold_constant_;
        int transition_counter = 0;

        // extract current StringDataArray with annotations/explanations;
        OpenMS::DataArrays::StringDataArray explanation_array = transition_spectrum.getStringDataArrays()[0];

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
          String current_explanation = explanation_array[peak_index];

          // write row for each transition
          // current int has to be higher than transition threshold and should not be smaller than threshold noise
          // current_mz has to be higher than min_fragment_mz and lower than max_fragment_mz
          if (current_int > threshold_transition && current_int > threshold_noise && current_mz > min_fragment_mz && current_mz < max_fragment_mz)
          {
            float rel_int = current_int / max_int;
            rmt.setPrecursorMZ((use_exact_mass && exact_mass_precursor != 0.0) ? exact_mass_precursor : csp.compound_info.pmass);
            rmt.setProductMZ(current_mz);
            TargetedExperimentHelper::TraMLProduct product;
            product.setMZ(current_mz);
            // charge state from adduct
            if (!adduct.empty() && adduct != "UNKNOWN")
            {
              product.setChargeState(getChargeFromAdduct_(adduct));
            }
            rmt.setProduct(product);
            rmt.setLibraryIntensity(rel_int);
            rmt.setCompoundRef(cmp.id);
            rmt.setNativeID(String(entry_counter) + "_" + String(transition_counter) + "_" + description + "_" + identifier_suffix);
            rmt.setMetaValue("annotation", DataValue(current_explanation));
            if (!csp.compound_info.native_ids_id.empty())
            {
              rmt.setMetaValue("native_ids_id", csp.compound_info.native_ids_id);
            }
            if (!csp.compound_info.m_ids_id.empty())
            {
              rmt.setMetaValue("m_ids_id", csp.compound_info.m_ids_id);
            }
            if (decoy == 1)
            {
              rmt.setDecoyTransitionType(ReactionMonitoringTransition::DecoyTransitionType::DECOY);
            }
            else
            {
              rmt.setDecoyTransitionType(ReactionMonitoringTransition::DecoyTransitionType::TARGET);
            }
            v_rmt.push_back(std::move(rmt));
            transition_counter += 1;
          }
        }
        MetaboTargetedAssay mta;
        mta.precursor_int = precursor_int;
        mta.compound_name = description;
        mta.compound_adduct = adduct;

        if (use_exact_mass)
        {
          mta.precursor_mz = exact_mass_precursor;
        }
        else
        {
          mta.precursor_mz = csp.compound_info.pmass;
        }

        mta.molecular_formula = sumformula;
        mta.compound_rt = feature_rt;
        mta.compound_file = csp.compound_info.file_index;

        mta.potential_cmp = cmp;
        mta.potential_rmts = v_rmt;

        if (!mta.potential_rmts.empty())
        {
          v_mta.push_back(std::move(mta));
        }
      }
      entry_counter += 1;
    }
    return v_mta;
  }

  // method to pair compound information (SiriusMSFile) with the annotated target spectrum from Sirius based on the m_id (unique identifier)
  std::vector< MetaboTargetedAssay::CompoundTargetDecoyPair > MetaboTargetedAssay::pairCompoundWithAnnotatedTDSpectraPairs(const std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                                                                                                                    const std::vector<SiriusFragmentAnnotation::SiriusTargetDecoySpectra>& annotated_spectra)
  {
    vector< MetaboTargetedAssay::CompoundTargetDecoyPair > v_cmp_spec;
    for (const auto& cmp : v_cmpinfo)
    {
      for (const auto& spectra : annotated_spectra)
      {
        if (cmp.m_ids_id == spectra.target.getName()) // the m_id is saved at MSSpectrum level as its name
        {
          v_cmp_spec.emplace_back(cmp, spectra);
        }
      }
    }
    return v_cmp_spec;
  }

  std::unordered_map< UInt64 , vector<MetaboTargetedAssay> > MetaboTargetedAssay::buildAmbiguityGroup(const vector<MetaboTargetedAssay>& v_mta,const double& ar_mz_tol, const double& ar_rt_tol, const String& ar_mz_tol_unit_res, size_t in_files_size)
  {
    String decoy_suffix = "_decoy";

    // group target and decoy position in vector based on the unique CompoundID/TransitionGroupID
    std::map<String, MetaboTargetedAssay::TargetDecoyGroup> target_decoy_groups;
    for (Size i = 0; i < v_mta.size(); ++i)
    {
      MetaboTargetedAssay current_entry = v_mta[i];
      if (!current_entry.potential_rmts.empty()) // should never be empty
      {
        // remove "decoy" tag from compound id for correct mapping
        if (current_entry.potential_rmts[0].getDecoyTransitionType() ==
            ReactionMonitoringTransition::DecoyTransitionType::DECOY)
        {
          String compoundId = current_entry.potential_cmp.id;
          compoundId.erase(compoundId.find(decoy_suffix), decoy_suffix.size());
          auto [it , success] = target_decoy_groups.emplace(compoundId, MetaboTargetedAssay::TargetDecoyGroup());
          if (success)
          {
            it->second.decoy_index = i;
            it->second.decoy_mz = current_entry.precursor_mz;
            it->second.decoy_rt = current_entry.compound_rt;
            it->second.decoy_file_number = current_entry.compound_file;
          }
          else
          {
            it->second.decoy_index = i;
            it->second.decoy_mz = current_entry.precursor_mz;
            it->second.decoy_rt = current_entry.compound_rt;
            it->second.decoy_file_number = current_entry.compound_file;
          }
        }
        if (current_entry.potential_rmts[0].getDecoyTransitionType() ==
            ReactionMonitoringTransition::DecoyTransitionType::TARGET)
        {
          auto [it , success] = target_decoy_groups.emplace(current_entry.potential_cmp.id, MetaboTargetedAssay::TargetDecoyGroup());
          if (success)
          {
            it->second.target_index = i;
            it->second.target_mz = current_entry.precursor_mz;
            it->second.target_rt = current_entry.compound_rt;
            it->second.target_file_number = current_entry.compound_file;
          }
          else
          {
            it->second.target_index = i;
            it->second.target_mz = current_entry.precursor_mz;
            it->second.target_rt = current_entry.compound_rt;
            it->second.target_file_number = current_entry.compound_file;
          }
        }
      }
    }

    std::unordered_map< UInt64 , vector<MetaboTargetedAssay> > ambiguity_groups;
    vector <FeatureMap> feature_maps;

    size_t loop_size;
    if ( in_files_size > 1)
    {
      loop_size = in_files_size;
    }
    else
    {
      loop_size = 2; // needs at least two FeatureMaps to compare - if not available add an empty one.
    }

    for (size_t i = 0; i < loop_size; i++)
    {
      FeatureMap fmap;
      fmap.setUniqueId();
      fmap.ensureUniqueId();
      String internal_file_path = "File" + std::to_string(i) + ".mzML";
      fmap.setPrimaryMSRunPath({internal_file_path});
      feature_maps.emplace_back(fmap);
    }

    for (const auto& it : target_decoy_groups)
    {
      // create minimal feature (mz,rt) with reference back to vector
      Feature f;
      f.setUniqueId();
      f.ensureUniqueId();
      PeptideIdentification pep;
      vector<PeptideIdentification> v_pep;

      // check - no target and decoy available
      if (it.second.target_mz == 0.0 && it.second.decoy_mz == 0.0)
      {
        continue;
      }
      // target and decoy available - check correspondence
      else if (it.second.target_mz != 0.0 && it.second.decoy_mz != 0.0)
      {
        if (!(it.second.target_mz == it.second.decoy_mz &&
              it.second.target_rt == it.second.decoy_rt &&
              it.second.target_file_number == it.second.decoy_file_number))
        {
          OPENMS_LOG_DEBUG << "The decoy and target do not correspond: " <<
                           " target_mz: " << it.second.target_mz <<
                           " decoy_mz: " << it.second.decoy_mz <<
                           " target_rt: " << it.second.target_rt  <<
                           " decoy_rt: " << it.second.decoy_rt  <<
                           " target_file_number: " << it.second.target_file_number <<
                           " decoy_file_number: " << it.second.decoy_file_number << std::endl;
          continue;
        }
      }

      // feature is based on target
      // check if target is valid before generating feature
      if (it.second.target_mz == 0.0 && it.second.target_rt == 0.0)
      {
        continue;
      }

      DPosition<2> pt(it.second.target_rt, it.second.target_mz);
      f.setPosition(pt);

      if (it.second.target_index != -1)
      {
        pep.setMetaValue("v_mta_target_index", DataValue(it.second.target_index));
      }

      if (it.second.decoy_index != -1)
      {
        pep.setMetaValue("v_mta_decoy_index", DataValue(it.second.decoy_index));
      }
      v_pep.push_back(pep);
      f.setPeptideIdentifications(v_pep);

      size_t cfile = it.second.target_file_number;
      feature_maps[cfile].push_back(f);
    }

    ConsensusMap c_map;
    FeatureGroupingAlgorithmQT fgaqt;
    Param param = fgaqt.getDefaults();
    param.setValue("ignore_charge", "true");
    param.setValue("distance_RT:max_difference", ar_rt_tol);
    param.setValue("distance_MZ:max_difference", ar_mz_tol);
    param.setValue("distance_MZ:unit", ar_mz_tol_unit_res);
    fgaqt.setParameters(param);

    // build ambiguity groups based on FeatureGroupingAlgorithmQt
    fgaqt.group(feature_maps, c_map);

    // assign unique ids
    c_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

    // build ambiguity groups based on consensus entries
    for (const auto& c_it : c_map)
    {
      vector <PeptideIdentification> v_pep;
      v_pep = c_it.getPeptideIdentifications();
      vector <MetaboTargetedAssay> ambi_group;
      for (const auto& p_it : v_pep)
      {
        if (p_it.metaValueExists("v_mta_target_index"))
        {
          int index = p_it.getMetaValue("v_mta_target_index");
          ambi_group.push_back(v_mta[index]);
        }
        if (p_it.metaValueExists("v_mta_decoy_index"))
        {
          int index = p_it.getMetaValue("v_mta_decoy_index");
          ambi_group.push_back(v_mta[index]);
        }
      }

      UInt64 entry = c_it.getUniqueId();
      auto mit = ambiguity_groups.find(entry);

      // allow to add targets and decoys, since they are grouped independently
      // targets and decoys have the exact mz and rt
      if (!(mit == ambiguity_groups.end()))
      {
        mit->second.insert(mit->second.end(), ambi_group.begin(), ambi_group.end());
      }
      else
      {
        ambiguity_groups.emplace(entry, ambi_group);
      }
    }
    return ambiguity_groups;
  }

  // resolve IDs based on the consensusXML
  // use the one with the highest intensity
  void MetaboTargetedAssay::resolveAmbiguityGroup(std::unordered_map< UInt64, vector<MetaboTargetedAssay> >& map_mta_filter, const double& total_occurrence_filter, size_t in_files_size)
  {
    vector<UInt64> empty_keys;
    for (auto& map_it : map_mta_filter)
    {
      // split the vector in targets and decoys
      vector<MetaboTargetedAssay> targets;
      vector<MetaboTargetedAssay> decoys;
      vector<MetaboTargetedAssay> targetdecoy;
      for (const auto& it : map_it.second)
      {
        if (it.potential_rmts[0].getDecoyTransitionType() == OpenMS::ReactionMonitoringTransition::TARGET)
        {
          targets.push_back(it);
        }
        else
        {
          decoys.push_back(it);
        }
      }
      // filter based on occurrence in samples (e.g. at least in 20% of the samples)
      MetaboTargetedAssay::filterBasedOnTotalOccurrence_(targets, total_occurrence_filter, in_files_size);
      MetaboTargetedAssay::filterBasedOnTotalOccurrence_(decoys, total_occurrence_filter, in_files_size);

      // if multiple possible identifications are reported within one ambiguity group
      // use the one with the highest occurrence
      MetaboTargetedAssay::filterBasedOnMolFormAdductOccurrence_(targets);
      MetaboTargetedAssay::filterBasedOnMolFormAdductOccurrence_(decoys);

      // From the resolved groups, use the target and decoy with the highest precursor intensity
      if (targets.empty() && decoys.empty())
      {
        empty_keys.push_back(map_it.first);
      }
      if (!targets.empty())
      {
        sortByPrecursorInt(targets);
        targetdecoy.push_back(targets[0]);
      }
      if (!decoys.empty())
      {
        sortByPrecursorInt(decoys);
        targetdecoy.push_back(decoys[0]);
      }
      map_it.second = targetdecoy;
    }

    for (const auto& key : empty_keys)
    {
      map_mta_filter.erase(key);
    }
  }
} // namespace OpenMS

/// @endcond
