// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  MRMFeatureFilter::MRMFeatureFilter() :
    DefaultParamHandler("MRMFeatureFilter")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  MRMFeatureFilter::~MRMFeatureFilter() = default;

  void MRMFeatureFilter::getDefaultParameters(Param& params) const
  {
    params.clear();

    params.setValue("flag_or_filter", "flag", "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC.", {"advanced"});
    params.setValidStrings("flag_or_filter", {"flag","filter"});
  }

  void MRMFeatureFilter::updateMembers_()
  {
    flag_or_filter_ = param_.getValue("flag_or_filter").toString();
  }

  void MRMFeatureFilter::FilterFeatureMap(FeatureMap& features,
    const MRMFeatureQC& filter_criteria,
    const TargetedExperiment& transitions
  )
  {
    // initialize QC variables
    FeatureMap features_filtered;

    // iterate through each component_group/feature
    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      String component_group_name = (String)features.at(feature_it).getMetaValue("PeptideRef");

      std::map<String, int> labels_and_transition_types = countLabelsAndTransitionTypes(features.at(feature_it), transitions);

      // initialize the new feature and subordinates
      std::vector<Feature> subordinates_filtered;
      bool cg_qc_pass = true;
      StringList cg_qc_fail_message_vec;
      UInt cg_tests_count{ 0 };

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features.at(feature_it).getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)features.at(feature_it).getSubordinates().at(sub_it).getMetaValue("native_id");
        bool c_qc_pass = true;
        StringList c_qc_fail_message_vec;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs.size(); ++cg_qc_it)
        {
          if (filter_criteria.component_group_qcs.at(cg_qc_it).component_group_name == component_group_name)
          {
            const double rt = features.at(feature_it).getRT();
            if (!checkRange(rt,
              filter_criteria.component_group_qcs.at(cg_qc_it).retention_time_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).retention_time_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("retention_time");
            }

            const double intensity = features.at(feature_it).getIntensity();
            if (!checkRange(intensity,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("intensity");
            }

            const double quality = features.at(feature_it).getOverallQuality();
            if (!checkRange(quality,
              filter_criteria.component_group_qcs.at(cg_qc_it).overall_quality_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).overall_quality_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("overall_quality");
            }
            // labels and transition counts QC
            if (!checkRange(labels_and_transition_types["n_heavy"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_heavy_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_heavy_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_heavy");
            }
            if (!checkRange(labels_and_transition_types["n_light"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_light_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_light_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_light");
            }
            if (!checkRange(labels_and_transition_types["n_detecting"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_detecting_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_detecting_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_detecting");
            }
            if (!checkRange(labels_and_transition_types["n_quantifying"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_quantifying_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_quantifying_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_quantifying");
            }
            if (!checkRange(labels_and_transition_types["n_identifying"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_identifying_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_identifying_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_identifying");
            }
            if (!checkRange(labels_and_transition_types["n_transitions"],
              filter_criteria.component_group_qcs.at(cg_qc_it).n_transitions_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).n_transitions_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_transitions");
            }

            cg_tests_count += 9;

            // ion ratio QC
            for (size_t sub_it2 = 0; sub_it2 < features.at(feature_it).getSubordinates().size(); ++sub_it2)
            {
              String component_name2 = (String)features.at(feature_it).getSubordinates().at(sub_it2).getMetaValue("native_id");
              // find the ion ratio pair
              if (!filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1.empty()
               && !filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2.empty()
               && filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1 == component_name
               && filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2 == component_name2)
              {
                double ion_ratio = calculateIonRatio(features.at(feature_it).getSubordinates().at(sub_it), features.at(feature_it).getSubordinates().at(sub_it2), filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_feature_name);

                if (!checkRange(ion_ratio,
                  filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_l,
                  filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_u))
                {
                  cg_qc_pass = false;
                  cg_qc_fail_message_vec.push_back("ion_ratio_pair[" + component_name + "/" + component_name2 + "]");
                }
                ++cg_tests_count;
              }
            }

            //std::pair<const String, std::pair<double, double>>
            for (const auto& kv : filter_criteria.component_group_qcs.at(cg_qc_it).meta_value_qc)
            {
              bool metavalue_exists{ false };
              if (!checkMetaValue(features.at(feature_it), kv.first, kv.second.first, kv.second.second, metavalue_exists))
              {
                cg_qc_pass = false;
                cg_qc_fail_message_vec.push_back(kv.first);
              }
              if (metavalue_exists) ++cg_tests_count;
            }
          }
        }

        UInt c_tests_count{ 0 };
        // iterate through feature/sub-feature QCs/filters
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs.at(c_qc_it).component_name == component_name)
          {
            // RT check
            const double rt = features.at(feature_it).getSubordinates().at(sub_it).getRT();
            if (!checkRange(rt,
              filter_criteria.component_qcs.at(c_qc_it).retention_time_l,
              filter_criteria.component_qcs.at(c_qc_it).retention_time_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("retention_time");
            }

            // intensity check
            double intensity = features.at(feature_it).getSubordinates().at(sub_it).getIntensity();
            if (!checkRange(intensity,
              filter_criteria.component_qcs.at(c_qc_it).intensity_l,
              filter_criteria.component_qcs.at(c_qc_it).intensity_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("intensity");
            }

            // overall quality check getQuality
            double quality = features.at(feature_it).getSubordinates().at(sub_it).getOverallQuality();
            if (!checkRange(quality,
              filter_criteria.component_qcs.at(c_qc_it).overall_quality_l,
              filter_criteria.component_qcs.at(c_qc_it).overall_quality_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("overall_quality");
            }

            c_tests_count += 3;

            // metaValue checks
            for (auto const& kv : filter_criteria.component_qcs.at(c_qc_it).meta_value_qc)
            {
              bool metavalue_exists{ false };
              if (!checkMetaValue(features.at(feature_it).getSubordinates().at(sub_it), kv.first, kv.second.first, kv.second.second, metavalue_exists))
              {
                c_qc_pass = false;
                c_qc_fail_message_vec.push_back(kv.first);
              }
              if (metavalue_exists) ++c_tests_count;
            }
          }
        }

        const double c_score = c_tests_count ? 1.0 - c_qc_fail_message_vec.size() / (double)c_tests_count : 1.0;
        features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_score", c_score);

        // Copy or Flag passing/failing subordinates
        if (c_qc_pass && flag_or_filter_ == "filter")
        {
          subordinates_filtered.push_back(features.at(feature_it).getSubordinates().at(sub_it));
        }
        else if (c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_pass", true);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_message", StringList());
        }
        else if (!c_qc_pass && flag_or_filter_ == "filter")
        {
          // do nothing
        }
        else if (!c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_pass", false);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_message", getUniqueSorted(c_qc_fail_message_vec));
        }
      }

      const double cg_score = cg_tests_count ? 1.0 - cg_qc_fail_message_vec.size() / (double)cg_tests_count : 1.0;
      features.at(feature_it).setMetaValue("QC_transition_group_score", cg_score);

      // Copy or Flag passing/failing Features
      if (cg_qc_pass && flag_or_filter_ == "filter" && !subordinates_filtered.empty())
      {
        Feature feature_filtered(features.at(feature_it));
        feature_filtered.setSubordinates(subordinates_filtered);
        features_filtered.push_back(feature_filtered);
      }
      else if (cg_qc_pass && flag_or_filter_ == "filter" && subordinates_filtered.empty())
      {
        // do nothing
      }
      else if (cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_pass", true);
        features.at(feature_it).setMetaValue("QC_transition_group_message", StringList());
      }
      else if (!cg_qc_pass && flag_or_filter_ == "filter")
      {
        // do nothing
      }
      else if (!cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_pass", false);
        features.at(feature_it).setMetaValue("QC_transition_group_message", getUniqueSorted(cg_qc_fail_message_vec));
      }
    }

    // replace with the filtered featureMap
    if (flag_or_filter_ == "filter")
    {
      features = features_filtered;
    }
  }

  void MRMFeatureFilter::FilterFeatureMapPercRSD(FeatureMap& features, const MRMFeatureQC& filter_criteria, const MRMFeatureQC& filter_values)
  {
    // initialize QC variables
    FeatureMap features_filtered;

    // iterate through each component_group/feature
    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      String component_group_name = (String)features.at(feature_it).getMetaValue("PeptideRef");

      // initialize the new feature and subordinates
      std::vector<Feature> subordinates_filtered;
      bool cg_qc_pass = true;
      StringList cg_qc_fail_message_vec;
      UInt cg_tests_count{ 0 };

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features.at(feature_it).getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)features.at(feature_it).getSubordinates().at(sub_it).getMetaValue("native_id");
        bool c_qc_pass = true;
        StringList c_qc_fail_message_vec;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs.size(); ++cg_qc_it)
        {
          if (filter_criteria.component_group_qcs.at(cg_qc_it).component_group_name == component_group_name)
          {
            if (!checkRange(filter_values.component_group_qcs.at(cg_qc_it).retention_time_u,
              filter_criteria.component_group_qcs.at(cg_qc_it).retention_time_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).retention_time_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("retention_time");
            }

            if (!checkRange(filter_values.component_group_qcs.at(cg_qc_it).intensity_u,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("intensity");
            }

            if (!checkRange(filter_values.component_group_qcs.at(cg_qc_it).overall_quality_u,
              filter_criteria.component_group_qcs.at(cg_qc_it).overall_quality_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).overall_quality_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("overall_quality");
            }

            cg_tests_count += 3;

            // ion ratio QC
            if (!filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1.empty()
             && !filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2.empty()) {
              if (!checkRange(filter_values.component_group_qcs.at(cg_qc_it).ion_ratio_u,
                filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_l,
                filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_u))
              {
                cg_qc_pass = false;
                cg_qc_fail_message_vec.push_back("ion_ratio_pair[" + filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1 + "/" + filter_criteria.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2 + "]");
              }
              ++cg_tests_count;
            }

            for (const auto& kv : filter_criteria.component_group_qcs.at(cg_qc_it).meta_value_qc)
            {
              if (!checkRange(filter_values.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second,
                kv.second.first,
                kv.second.second))
              {
                cg_qc_pass = false;
                cg_qc_fail_message_vec.push_back(kv.first);
              }
              ++cg_tests_count;
            }
          }
        }

        UInt c_tests_count{ 0 };
        // iterate through feature/sub-feature QCs/filters
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs.at(c_qc_it).component_name == component_name)
          {
            // RT check
            if (!checkRange(filter_values.component_qcs.at(c_qc_it).retention_time_u,
              filter_criteria.component_qcs.at(c_qc_it).retention_time_l,
              filter_criteria.component_qcs.at(c_qc_it).retention_time_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("retention_time");
            }

            // intensity check
            if (!checkRange(filter_values.component_qcs.at(c_qc_it).intensity_u,
              filter_criteria.component_qcs.at(c_qc_it).intensity_l,
              filter_criteria.component_qcs.at(c_qc_it).intensity_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("intensity");
            }

            // overall quality check getQuality
            if (!checkRange(filter_values.component_qcs.at(c_qc_it).overall_quality_u,
              filter_criteria.component_qcs.at(c_qc_it).overall_quality_l,
              filter_criteria.component_qcs.at(c_qc_it).overall_quality_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("overall_quality");
            }

            c_tests_count += 3;

            // metaValue checks
            for (auto const& kv : filter_criteria.component_qcs.at(c_qc_it).meta_value_qc)
            {
              if (!checkRange(filter_values.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second,
                kv.second.first,
                kv.second.second))
              {
                c_qc_pass = false;
                c_qc_fail_message_vec.push_back(kv.first);
              }
              ++c_tests_count;
            }
          }
        }

        const double c_score = c_tests_count ? 1.0 - c_qc_fail_message_vec.size() / (double)c_tests_count : 1.0;
        features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%RSD_score", c_score);

        // Copy or Flag passing/failing subordinates
        if (c_qc_pass && flag_or_filter_ == "filter")
        {
          subordinates_filtered.push_back(features.at(feature_it).getSubordinates().at(sub_it));
        }
        else if (c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%RSD_pass", true);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%RSD_message", StringList());
        }
        else if (!c_qc_pass && flag_or_filter_ == "filter")
        {
          // do nothing
        }
        else if (!c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%RSD_pass", false);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%RSD_message", getUniqueSorted(c_qc_fail_message_vec));
        }
      }

      const double cg_score = cg_tests_count ? 1.0 - cg_qc_fail_message_vec.size() / (double)cg_tests_count : 1.0;
      features.at(feature_it).setMetaValue("QC_transition_group_%RSD_score", cg_score);

      // Copy or Flag passing/failing Features
      if (cg_qc_pass && flag_or_filter_ == "filter" && !subordinates_filtered.empty())
      {
        Feature feature_filtered(features.at(feature_it));
        feature_filtered.setSubordinates(subordinates_filtered);
        features_filtered.push_back(feature_filtered);
      }
      else if (cg_qc_pass && flag_or_filter_ == "filter" && subordinates_filtered.empty())
      {
        // do nothing
      }
      else if (cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_%RSD_pass", true);
        features.at(feature_it).setMetaValue("QC_transition_group_%RSD_message", StringList());
      }
      else if (!cg_qc_pass && flag_or_filter_ == "filter")
      {
        // do nothing
      }
      else if (!cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_%RSD_pass", false);
        features.at(feature_it).setMetaValue("QC_transition_group_%RSD_message", getUniqueSorted(cg_qc_fail_message_vec));
      }
    }

    // replace with the filtered featureMap
    if (flag_or_filter_ == "filter")
    {
      features = features_filtered;
    }
  }

  void MRMFeatureFilter::FilterFeatureMapBackgroundInterference(FeatureMap& features, const MRMFeatureQC& filter_criteria, const MRMFeatureQC& filter_values)
  {
    // initialize QC variables
    FeatureMap features_filtered;

    // iterate through each component_group/feature
    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      String component_group_name = (String)features.at(feature_it).getMetaValue("PeptideRef");

      // initialize the new feature and subordinates
      std::vector<Feature> subordinates_filtered;
      bool cg_qc_pass = true;
      StringList cg_qc_fail_message_vec;
      UInt cg_tests_count{ 0 };

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features.at(feature_it).getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)features.at(feature_it).getSubordinates().at(sub_it).getMetaValue("native_id");
        bool c_qc_pass = true;
        StringList c_qc_fail_message_vec;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs.size(); ++cg_qc_it)
        {
          if (filter_criteria.component_group_qcs.at(cg_qc_it).component_group_name == component_group_name)
          {
            // intensity check
            const double perc_background_interference = filter_values.component_group_qcs.at(cg_qc_it).intensity_u / features.at(feature_it).getIntensity() * 100;
            if (!checkRange(perc_background_interference,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_l,
              filter_criteria.component_group_qcs.at(cg_qc_it).intensity_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("intensity");
            }
            ++cg_tests_count;
          }
        }

        UInt c_tests_count{ 0 };
        // iterate through feature/sub-feature QCs/filters
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs.at(c_qc_it).component_name == component_name)
          {
            // intensity check
            const double perc_background_interference = filter_values.component_qcs.at(c_qc_it).intensity_u / features.at(feature_it).getSubordinates().at(sub_it).getIntensity() * 100;
            if (!checkRange(perc_background_interference,
              filter_criteria.component_qcs.at(c_qc_it).intensity_l,
              filter_criteria.component_qcs.at(c_qc_it).intensity_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("intensity");
            }
            ++c_tests_count;
          }
        }

        const double c_score = c_tests_count ? 1.0 - c_qc_fail_message_vec.size() / (double)c_tests_count : 1.0;
        features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%BackgroundInterference_score", c_score);

        // Copy or Flag passing/failing subordinates
        if (c_qc_pass && flag_or_filter_ == "filter")
        {
          subordinates_filtered.push_back(features.at(feature_it).getSubordinates().at(sub_it));
        }
        else if (c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%BackgroundInterference_pass", true);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%BackgroundInterference_message", StringList());
        }
        else if (!c_qc_pass && flag_or_filter_ == "filter")
        {
          // do nothing
        }
        else if (!c_qc_pass && flag_or_filter_ == "flag")
        {
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%BackgroundInterference_pass", false);
          features.at(feature_it).getSubordinates().at(sub_it).setMetaValue("QC_transition_%BackgroundInterference_message", getUniqueSorted(c_qc_fail_message_vec));
        }
      }

      const double cg_score = cg_tests_count ? 1.0 - cg_qc_fail_message_vec.size() / (double)cg_tests_count : 1.0;
      features.at(feature_it).setMetaValue("QC_transition_group_%BackgroundInterference_score", cg_score);

      // Copy or Flag passing/failing Features
      if (cg_qc_pass && flag_or_filter_ == "filter" && !subordinates_filtered.empty())
      {
        Feature feature_filtered(features.at(feature_it));
        feature_filtered.setSubordinates(subordinates_filtered);
        features_filtered.push_back(feature_filtered);
      }
      else if (cg_qc_pass && flag_or_filter_ == "filter" && subordinates_filtered.empty())
      {
        // do nothing
      }
      else if (cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_%BackgroundInterference_pass", true);
        features.at(feature_it).setMetaValue("QC_transition_group_%BackgroundInterference_message", StringList());
      }
      else if (!cg_qc_pass && flag_or_filter_ == "filter")
      {
        // do nothing
      }
      else if (!cg_qc_pass && flag_or_filter_ == "flag")
      {
        features.at(feature_it).setMetaValue("QC_transition_group_%BackgroundInterference_pass", false);
        features.at(feature_it).setMetaValue("QC_transition_group_%BackgroundInterference_message", getUniqueSorted(cg_qc_fail_message_vec));
      }
    }

    // replace with the filtered featureMap
    if (flag_or_filter_ == "filter")
    {
      features = features_filtered;
    }
  }

  void MRMFeatureFilter::EstimateDefaultMRMFeatureQCValues(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions, const bool& init_template_values) const
  {
    // iterate through each sample and accumulate the min/max values in the samples in the filter_template
    for (size_t sample_it = 0; sample_it < samples.size(); sample_it++) {

      // iterate through each component_group/feature
      for (size_t feature_it = 0; feature_it < samples.at(sample_it).size(); ++feature_it)
      {
        String component_group_name = (String)samples.at(sample_it).at(feature_it).getMetaValue("PeptideRef");
        std::map<String, int> labels_and_transition_types = countLabelsAndTransitionTypes(samples.at(sample_it).at(feature_it), transitions);

        // iterate through each component/sub-feature
        for (size_t sub_it = 0; sub_it < samples.at(sample_it).at(feature_it).getSubordinates().size(); ++sub_it)
        {
          String component_name = (String)samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getMetaValue("native_id");

          // iterate through multi-feature/multi-sub-feature QCs/filters
          // iterate through component_groups
          for (size_t cg_qc_it = 0; cg_qc_it < filter_template.component_group_qcs.size(); ++cg_qc_it)
          {
            if (filter_template.component_group_qcs.at(cg_qc_it).component_group_name == component_group_name)
            {
              const double rt = samples.at(sample_it).at(feature_it).getRT();
              if (sample_it == 0 && init_template_values) {
                initRange(rt,
                  filter_template.component_group_qcs.at(cg_qc_it).retention_time_l,
                  filter_template.component_group_qcs.at(cg_qc_it).retention_time_u);
              } else {
                updateRange(rt,
                  filter_template.component_group_qcs.at(cg_qc_it).retention_time_l,
                  filter_template.component_group_qcs.at(cg_qc_it).retention_time_u);
              }

              const double intensity = samples.at(sample_it).at(feature_it).getIntensity();
              if (sample_it == 0 && init_template_values) {
                initRange(intensity,
                  filter_template.component_group_qcs.at(cg_qc_it).intensity_l,
                  filter_template.component_group_qcs.at(cg_qc_it).intensity_u);
              } else {
                updateRange(intensity,
                  filter_template.component_group_qcs.at(cg_qc_it).intensity_l,
                  filter_template.component_group_qcs.at(cg_qc_it).intensity_u);
              }

              const double quality = samples.at(sample_it).at(feature_it).getOverallQuality();
              if (sample_it == 0 && init_template_values) {
                initRange(quality,
                  filter_template.component_group_qcs.at(cg_qc_it).overall_quality_l,
                  filter_template.component_group_qcs.at(cg_qc_it).overall_quality_u);
              } else {
                updateRange(quality,
                  filter_template.component_group_qcs.at(cg_qc_it).overall_quality_l,
                  filter_template.component_group_qcs.at(cg_qc_it).overall_quality_u);
              }

              // labels and transition counts QC
              if (sample_it == 0 && init_template_values) {
                initRange(labels_and_transition_types["n_heavy"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_heavy_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_heavy_u);
                initRange(labels_and_transition_types["n_light"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_light_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_light_u);
                initRange(labels_and_transition_types["n_detecting"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_detecting_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_detecting_u);
                initRange(labels_and_transition_types["n_quantifying"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_quantifying_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_quantifying_u);
                initRange(labels_and_transition_types["n_identifying"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_identifying_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_identifying_u);
                initRange(labels_and_transition_types["n_transitions"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_transitions_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_transitions_u);
              } else {
                updateRange(labels_and_transition_types["n_heavy"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_heavy_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_heavy_u);
                updateRange(labels_and_transition_types["n_light"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_light_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_light_u);
                updateRange(labels_and_transition_types["n_detecting"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_detecting_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_detecting_u);
                updateRange(labels_and_transition_types["n_quantifying"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_quantifying_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_quantifying_u);
                updateRange(labels_and_transition_types["n_identifying"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_identifying_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_identifying_u);
                updateRange(labels_and_transition_types["n_transitions"],
                  filter_template.component_group_qcs.at(cg_qc_it).n_transitions_l,
                  filter_template.component_group_qcs.at(cg_qc_it).n_transitions_u);
              }

              // ion ratio QC
              for (size_t sub_it2 = 0; sub_it2 < samples.at(sample_it).at(feature_it).getSubordinates().size(); ++sub_it2)
              {
                String component_name2 = (String)samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it2).getMetaValue("native_id");
                // find the ion ratio pair
                if (!filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1.empty()
                 && !filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2.empty()
                 && filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1 == component_name
                 && filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2 == component_name2)
                {
                  double ion_ratio = calculateIonRatio(samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it), samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it2), filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_feature_name);
                  if (sample_it == 0 && init_template_values) {
                    initRange(ion_ratio,
                      filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_l,
                      filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_u);
                  } else {
                    updateRange(ion_ratio,
                      filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_l,
                      filter_template.component_group_qcs.at(cg_qc_it).ion_ratio_u);
                  }
                }
              }

              for (auto& kv : filter_template.component_group_qcs.at(cg_qc_it).meta_value_qc)
              {
                bool metavalue_exists{ false };
                if (sample_it == 0 && init_template_values) {
                  updateMetaValue(samples.at(sample_it).at(feature_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
                } else {
                  initMetaValue(samples.at(sample_it).at(feature_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
                }
              }
            }
          }

          // iterate through feature/sub-feature QCs/filters
          for (size_t c_qc_it = 0; c_qc_it < filter_template.component_qcs.size(); ++c_qc_it)
          {
            if (filter_template.component_qcs.at(c_qc_it).component_name == component_name)
            {
              // RT check
              const double rt = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getRT();
              if (sample_it == 0 && init_template_values) {
                initRange(rt,
                  filter_template.component_qcs.at(c_qc_it).retention_time_l,
                  filter_template.component_qcs.at(c_qc_it).retention_time_u);
              } else {
                updateRange(rt,
                  filter_template.component_qcs.at(c_qc_it).retention_time_l,
                  filter_template.component_qcs.at(c_qc_it).retention_time_u);
              }

              // intensity check
              double intensity = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getIntensity();
              if (sample_it == 0 && init_template_values) {
                initRange(intensity,
                  filter_template.component_qcs.at(c_qc_it).intensity_l,
                  filter_template.component_qcs.at(c_qc_it).intensity_u);
              } else {
                updateRange(intensity,
                  filter_template.component_qcs.at(c_qc_it).intensity_l,
                  filter_template.component_qcs.at(c_qc_it).intensity_u);
              }

              // overall quality check getQuality
              double quality = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getOverallQuality();
              if (sample_it == 0 && init_template_values) {
                initRange(quality,
                  filter_template.component_qcs.at(c_qc_it).overall_quality_l,
                  filter_template.component_qcs.at(c_qc_it).overall_quality_u);
              } else {
                updateRange(quality,
                  filter_template.component_qcs.at(c_qc_it).overall_quality_l,
                  filter_template.component_qcs.at(c_qc_it).overall_quality_u);
              }

              // metaValue checks
              for (auto& kv : filter_template.component_qcs.at(c_qc_it).meta_value_qc)
              {
                bool metavalue_exists{ false };
                if (sample_it == 0 && init_template_values) {
                  initMetaValue(samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
                } else {
                  updateMetaValue(samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
                }
              }
            }
          }
        }
      }
    }
  }

  void MRMFeatureFilter::TransferLLOQAndULOQToCalculatedConcentrationBounds(const std::vector<AbsoluteQuantitationMethod>& quantitation_method, MRMFeatureQC& filter_template)
  {
    // Iterate through the quantitation method and update the MetaValue for `calculated_concentration` in the filter_template
    for (const AbsoluteQuantitationMethod& quant_method : quantitation_method) {
      if (quant_method.getLLOQ() == 0 && quant_method.getULOQ() == 0) continue;

      // iterate through feature/sub-feature QCs/filters
      for (size_t c_qc_it = 0; c_qc_it < filter_template.component_qcs.size(); ++c_qc_it) {
        if (filter_template.component_qcs.at(c_qc_it).component_name == quant_method.getComponentName()) {

          // update the lower/upper bound for the `calculated_concentration` metaValue
          filter_template.component_qcs.at(c_qc_it).meta_value_qc.at("calculated_concentration").first = quant_method.getLLOQ();
          filter_template.component_qcs.at(c_qc_it).meta_value_qc.at("calculated_concentration").second = quant_method.getULOQ();
        }
      }
    }
  }

  void MRMFeatureFilter::EstimatePercRSD(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const
  {
    // iterate through each sample and accumulate the values in the filter_values
    std::vector<MRMFeatureQC> filter_values;
    accumulateFilterValues(filter_values, samples, filter_template, transitions);

    // Determine the AVE for each filter_template value
    MRMFeatureQC filter_mean;
    calculateFilterValuesMean(filter_mean, filter_values, filter_template);

    // Determine the STD for each filter_template value
    MRMFeatureQC filter_var;
    calculateFilterValuesVar(filter_var, filter_values, filter_mean, filter_template);

    // Determine the %RSD for each filter_template value
    calculateFilterValuesPercRSD(filter_template, filter_mean, filter_var);
  }

  void MRMFeatureFilter::EstimateBackgroundInterferences(const std::vector<FeatureMap>& samples, MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const
  {
    // iterate through each sample and accumulate the values in the filter_values
    std::vector<MRMFeatureQC> filter_values;
    accumulateFilterValues(filter_values, samples, filter_template, transitions);

    // Determine the AVE for each filter_template value
    calculateFilterValuesMean(filter_template, filter_values, filter_template);
  }

  std::map<String, int> MRMFeatureFilter::countLabelsAndTransitionTypes(
    const Feature& component_group,
    const TargetedExperiment& transitions) const
  {
    int n_heavy(0), n_light(0), n_quant(0), n_detect(0), n_ident(0), n_trans(0);
    std::map<String, int> output;

    for (size_t cg_it = 0; cg_it < component_group.getSubordinates().size(); ++cg_it)
    {

      // extract out the matching transition
      ReactionMonitoringTransition transition;
      for (size_t trans_it = 0; trans_it < transitions.getTransitions().size(); ++trans_it)
      {
        if (transitions.getTransitions()[trans_it].getNativeID() == component_group.getSubordinates()[cg_it].getMetaValue("native_id"))
        {
          transition = transitions.getTransitions()[trans_it];
          break;
        }
      }

      // count labels and transition types
      String label_type = (String)component_group.getSubordinates()[cg_it].getMetaValue("LabelType");
      if (label_type == "Heavy")
      {
        ++n_heavy;
      } 
      else if (label_type == "Light")
      {
        ++n_light;
      }
      if (transition.isQuantifyingTransition())
      {
        ++n_quant;
      }
      if (transition.isIdentifyingTransition())
      {
        ++n_ident;
      }
      if (transition.isDetectingTransition())
      {
        ++n_detect;
      }
      ++n_trans;
    }

    // record
    output["n_heavy"] = n_heavy;
    output["n_light"] = n_light;
    output["n_quantifying"] = n_quant;
    output["n_identifying"] = n_ident;
    output["n_detecting"] = n_detect;
    output["n_transitions"] = n_trans;

    return output;
  }

  double MRMFeatureFilter::calculateIonRatio(const Feature& component_1, const Feature& component_2, const String& feature_name) const
  {
    double ratio = 0.0;
    // member feature_name access
    if (feature_name == "intensity")
    {
      if (component_1.metaValueExists("native_id")&& component_2.metaValueExists("native_id"))
      {
        const double feature_1 = component_1.getIntensity();
        const double feature_2 = component_2.getIntensity();
        ratio = feature_1 / feature_2;
      }
      else if (component_1.metaValueExists("native_id"))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = component_1.getIntensity();
        ratio = feature_1;
      }
    }
    // metaValue feature_name access
    else
    {
      if (component_1.metaValueExists(feature_name)&& component_2.metaValueExists(feature_name))
      {
        const double feature_1 = component_1.getMetaValue(feature_name);
        const double feature_2 = component_2.getMetaValue(feature_name);
        ratio = feature_1 / feature_2;
      }
      else if (component_1.metaValueExists(feature_name))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
        const double feature_1 = component_1.getMetaValue(feature_name);
        ratio = feature_1;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
      }
    }

    return ratio;
  }

  bool MRMFeatureFilter::checkMetaValue(
    const Feature& component,
    const String& meta_value_key,
    const double& meta_value_l,
    const double& meta_value_u,
    bool& key_exists
  ) const
  {
    bool check = true;
    if (component.metaValueExists(meta_value_key)) {
      key_exists = true;
      const double meta_value = (double)component.getMetaValue(meta_value_key);
      check = checkRange(meta_value, meta_value_l, meta_value_u);
    } else {
      key_exists = false;
      OPENMS_LOG_DEBUG << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }
    return check;
  }

  void MRMFeatureFilter::updateMetaValue(const Feature& component, const String& meta_value_key, double& meta_value_l, double& meta_value_u, bool& key_exists) const
  {
    if (component.metaValueExists(meta_value_key))
    {
      key_exists = true;
      const double meta_value = (double)component.getMetaValue(meta_value_key);
      updateRange(meta_value, meta_value_l, meta_value_u);
    }
    else
    {
      key_exists = false;
      OPENMS_LOG_DEBUG << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }
  }

  void MRMFeatureFilter::setMetaValue(const Feature& component, const String& meta_value_key, double& meta_value_l, double& meta_value_u, bool& key_exists) const
  {
    if (component.metaValueExists(meta_value_key))
    {
      key_exists = true;
      const double meta_value = (double)component.getMetaValue(meta_value_key);
      setRange(meta_value, meta_value_l, meta_value_u);
    }
    else
    {
      key_exists = false;
      OPENMS_LOG_DEBUG << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }
  }

  void MRMFeatureFilter::initMetaValue(const Feature& component, const String& meta_value_key, double& meta_value_l, double& meta_value_u, bool& key_exists) const
  {
    if (component.metaValueExists(meta_value_key))
    {
      key_exists = true;
      const double meta_value = (double)component.getMetaValue(meta_value_key);
      initRange(meta_value, meta_value_l, meta_value_u);
    }
    else
    {
      key_exists = false;
      OPENMS_LOG_DEBUG << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }
  }

  StringList MRMFeatureFilter::getUniqueSorted(const StringList& messages) const
  {
    StringList unique{ messages };
    std::sort(unique.begin(), unique.end());
    unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    return unique;
  }

  void MRMFeatureFilter::accumulateFilterValues(std::vector<MRMFeatureQC>& filter_values, const std::vector<FeatureMap>& samples, const MRMFeatureQC& filter_template, const TargetedExperiment& transitions) const
  {
    // iterate through each sample and accumulate the values in the filter_values
    for (size_t sample_it = 0; sample_it < samples.size(); sample_it++) {
      MRMFeatureQC filter_value = filter_template;

      // iterate through each component_group/feature
      for (size_t feature_it = 0; feature_it < samples.at(sample_it).size(); ++feature_it)
      {
        String component_group_name = (String)samples.at(sample_it).at(feature_it).getMetaValue("PeptideRef");
        std::map<String, int> labels_and_transition_types = countLabelsAndTransitionTypes(samples.at(sample_it).at(feature_it), transitions);

        // iterate through each component/sub-feature
        for (size_t sub_it = 0; sub_it < samples.at(sample_it).at(feature_it).getSubordinates().size(); ++sub_it)
        {
          String component_name = (String)samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getMetaValue("native_id");

          // iterate through multi-feature/multi-sub-feature QCs/filters
          // iterate through component_groups
          for (size_t cg_qc_it = 0; cg_qc_it < filter_value.component_group_qcs.size(); ++cg_qc_it)
          {
            if (filter_value.component_group_qcs.at(cg_qc_it).component_group_name == component_group_name)
            {
              const double rt = samples.at(sample_it).at(feature_it).getRT();
              setRange(rt,
                filter_value.component_group_qcs.at(cg_qc_it).retention_time_l,
                filter_value.component_group_qcs.at(cg_qc_it).retention_time_u);

              const double intensity = samples.at(sample_it).at(feature_it).getIntensity();
              setRange(intensity,
                filter_value.component_group_qcs.at(cg_qc_it).intensity_l,
                filter_value.component_group_qcs.at(cg_qc_it).intensity_u);

              const double quality = samples.at(sample_it).at(feature_it).getOverallQuality();
              setRange(quality,
                filter_value.component_group_qcs.at(cg_qc_it).overall_quality_l,
                filter_value.component_group_qcs.at(cg_qc_it).overall_quality_u);

              // labels and transition counts QC
              setRange(labels_and_transition_types["n_heavy"],
                filter_value.component_group_qcs.at(cg_qc_it).n_heavy_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_heavy_u);
              setRange(labels_and_transition_types["n_light"],
                filter_value.component_group_qcs.at(cg_qc_it).n_light_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_light_u);
              setRange(labels_and_transition_types["n_detecting"],
                filter_value.component_group_qcs.at(cg_qc_it).n_detecting_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_detecting_u);
              setRange(labels_and_transition_types["n_quantifying"],
                filter_value.component_group_qcs.at(cg_qc_it).n_quantifying_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_quantifying_u);
              setRange(labels_and_transition_types["n_identifying"],
                filter_value.component_group_qcs.at(cg_qc_it).n_identifying_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_identifying_u);
              setRange(labels_and_transition_types["n_transitions"],
                filter_value.component_group_qcs.at(cg_qc_it).n_transitions_l,
                filter_value.component_group_qcs.at(cg_qc_it).n_transitions_u);

              // ion ratio QC
              for (size_t sub_it2 = 0; sub_it2 < samples.at(sample_it).at(feature_it).getSubordinates().size(); ++sub_it2)
              {
                String component_name2 = (String)samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it2).getMetaValue("native_id");
                // find the ion ratio pair
                if (!filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1.empty()
                 && !filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2.empty()
                 && filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_1 == component_name
                 && filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_pair_name_2 == component_name2)
                {
                  double ion_ratio = calculateIonRatio(samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it), samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it2), filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_feature_name);

                  setRange(ion_ratio,
                    filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_l,
                    filter_value.component_group_qcs.at(cg_qc_it).ion_ratio_u);
                }
              }

              for (auto& kv : filter_value.component_group_qcs.at(cg_qc_it).meta_value_qc)
              {
                bool metavalue_exists{ false };
                setMetaValue(samples.at(sample_it).at(feature_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
              }
            }
          }

          // iterate through feature/sub-feature QCs/filters
          for (size_t c_qc_it = 0; c_qc_it < filter_value.component_qcs.size(); ++c_qc_it)
          {
            if (filter_value.component_qcs.at(c_qc_it).component_name == component_name)
            {
              // RT check
              const double rt = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getRT();
              setRange(rt,
                filter_value.component_qcs.at(c_qc_it).retention_time_l,
                filter_value.component_qcs.at(c_qc_it).retention_time_u);

              // intensity check
              double intensity = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getIntensity();
              setRange(intensity,
                filter_value.component_qcs.at(c_qc_it).intensity_l,
                filter_value.component_qcs.at(c_qc_it).intensity_u);

              // overall quality check getQuality
              double quality = samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it).getOverallQuality();
              setRange(quality,
                filter_value.component_qcs.at(c_qc_it).overall_quality_l,
                filter_value.component_qcs.at(c_qc_it).overall_quality_u);

              // metaValue checks
              for (auto& kv : filter_value.component_qcs.at(c_qc_it).meta_value_qc)
              {
                bool metavalue_exists{ false };
                setMetaValue(samples.at(sample_it).at(feature_it).getSubordinates().at(sub_it), kv.first, kv.second.first, kv.second.second, metavalue_exists);
              }
            }
          }
        }
      }
      filter_values.push_back(filter_value);
    }
  }

  void MRMFeatureFilter::zeroFilterValues(MRMFeatureQC& filter_zeros, const MRMFeatureQC& filter_template) const
  {
    // Create a zero filter template for subsequent AVE and STD calculations
    filter_zeros = filter_template;
    for (size_t cg_qc_it = 0; cg_qc_it < filter_zeros.component_group_qcs.size(); ++cg_qc_it) {
      filter_zeros.component_group_qcs.at(cg_qc_it).retention_time_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).retention_time_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).intensity_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).intensity_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).overall_quality_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).overall_quality_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_heavy_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_heavy_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_light_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_light_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_detecting_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_detecting_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_quantifying_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_quantifying_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_identifying_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_identifying_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_transitions_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).n_transitions_u = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).ion_ratio_l = 0;
      filter_zeros.component_group_qcs.at(cg_qc_it).ion_ratio_u = 0;
      for (auto& kv : filter_zeros.component_group_qcs.at(cg_qc_it).meta_value_qc) {
        kv.second.first = 0;
        kv.second.second = 0;
      }
    }
    for (size_t c_qc_it = 0; c_qc_it < filter_zeros.component_qcs.size(); ++c_qc_it) {
      filter_zeros.component_qcs.at(c_qc_it).retention_time_l = 0;
      filter_zeros.component_qcs.at(c_qc_it).retention_time_u = 0;
      filter_zeros.component_qcs.at(c_qc_it).intensity_l = 0;
      filter_zeros.component_qcs.at(c_qc_it).intensity_u = 0;
      filter_zeros.component_qcs.at(c_qc_it).overall_quality_l = 0;
      filter_zeros.component_qcs.at(c_qc_it).overall_quality_u = 0;
      for (auto& kv : filter_zeros.component_qcs.at(c_qc_it).meta_value_qc) {
        kv.second.first = 0;
        kv.second.second = 0;
      }
    }
  }

  void MRMFeatureFilter::calculateFilterValuesMean(MRMFeatureQC& filter_mean, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_template) const
  {
    // Determine the AVE for each filter_template value
    zeroFilterValues(filter_mean, filter_template);
    for (const MRMFeatureQC& filter : filter_values) { // Accumulate the sum
      for (size_t cg_qc_it = 0; cg_qc_it < filter_mean.component_group_qcs.size(); ++cg_qc_it) {
        filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l += filter.component_group_qcs.at(cg_qc_it).retention_time_l;
        filter_mean.component_group_qcs.at(cg_qc_it).retention_time_u += filter.component_group_qcs.at(cg_qc_it).retention_time_u;
        filter_mean.component_group_qcs.at(cg_qc_it).intensity_l += filter.component_group_qcs.at(cg_qc_it).intensity_l;
        filter_mean.component_group_qcs.at(cg_qc_it).intensity_u += filter.component_group_qcs.at(cg_qc_it).intensity_u;
        filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l += filter.component_group_qcs.at(cg_qc_it).overall_quality_l;
        filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u += filter.component_group_qcs.at(cg_qc_it).overall_quality_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l += filter.component_group_qcs.at(cg_qc_it).n_heavy_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u += filter.component_group_qcs.at(cg_qc_it).n_heavy_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_light_l += filter.component_group_qcs.at(cg_qc_it).n_light_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_light_u += filter.component_group_qcs.at(cg_qc_it).n_light_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l += filter.component_group_qcs.at(cg_qc_it).n_detecting_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u += filter.component_group_qcs.at(cg_qc_it).n_detecting_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l += filter.component_group_qcs.at(cg_qc_it).n_quantifying_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u += filter.component_group_qcs.at(cg_qc_it).n_quantifying_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l += filter.component_group_qcs.at(cg_qc_it).n_identifying_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u += filter.component_group_qcs.at(cg_qc_it).n_identifying_u;
        filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l += filter.component_group_qcs.at(cg_qc_it).n_transitions_l;
        filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u += filter.component_group_qcs.at(cg_qc_it).n_transitions_u;
        filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l += filter.component_group_qcs.at(cg_qc_it).ion_ratio_l;
        filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u += filter.component_group_qcs.at(cg_qc_it).ion_ratio_u;
        for (auto& kv : filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc) {
          kv.second.first += filter.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first;
          kv.second.second += filter.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second;
        }
      }
      for (size_t c_qc_it = 0; c_qc_it < filter_mean.component_qcs.size(); ++c_qc_it) {
        filter_mean.component_qcs.at(c_qc_it).retention_time_l += filter.component_qcs.at(c_qc_it).retention_time_l;
        filter_mean.component_qcs.at(c_qc_it).retention_time_u += filter.component_qcs.at(c_qc_it).retention_time_u;
        filter_mean.component_qcs.at(c_qc_it).intensity_l += filter.component_qcs.at(c_qc_it).intensity_l;
        filter_mean.component_qcs.at(c_qc_it).intensity_u += filter.component_qcs.at(c_qc_it).intensity_u;
        filter_mean.component_qcs.at(c_qc_it).overall_quality_l += filter.component_qcs.at(c_qc_it).overall_quality_l;
        filter_mean.component_qcs.at(c_qc_it).overall_quality_u += filter.component_qcs.at(c_qc_it).overall_quality_u;
        for (auto& kv : filter_mean.component_qcs.at(c_qc_it).meta_value_qc) {
          kv.second.first += filter.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first;
          kv.second.second += filter.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second;
        }
      }
    }
    for (size_t cg_qc_it = 0; cg_qc_it < filter_mean.component_group_qcs.size(); ++cg_qc_it) {// Divide by the size (performed separately due to int types...)
      filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l = filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).retention_time_u = filter_mean.component_group_qcs.at(cg_qc_it).retention_time_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).intensity_l = filter_mean.component_group_qcs.at(cg_qc_it).intensity_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).intensity_u = filter_mean.component_group_qcs.at(cg_qc_it).intensity_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l = filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u = filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l = filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u = filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_light_l = filter_mean.component_group_qcs.at(cg_qc_it).n_light_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_light_u = filter_mean.component_group_qcs.at(cg_qc_it).n_light_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l = filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u = filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l = filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u = filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l = filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u = filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l = filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u = filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l = filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l / filter_values.size();
      filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u = filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u / filter_values.size();
      for (auto& kv : filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc) {
        kv.second.first = filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first / filter_values.size();
        kv.second.second = filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second / filter_values.size();
      }
    }
    for (size_t c_qc_it = 0; c_qc_it < filter_mean.component_qcs.size(); ++c_qc_it) {
      filter_mean.component_qcs.at(c_qc_it).retention_time_l = filter_mean.component_qcs.at(c_qc_it).retention_time_l / filter_values.size();
      filter_mean.component_qcs.at(c_qc_it).retention_time_u = filter_mean.component_qcs.at(c_qc_it).retention_time_u / filter_values.size();
      filter_mean.component_qcs.at(c_qc_it).intensity_l = filter_mean.component_qcs.at(c_qc_it).intensity_l / filter_values.size();
      filter_mean.component_qcs.at(c_qc_it).intensity_u = filter_mean.component_qcs.at(c_qc_it).intensity_u / filter_values.size();
      filter_mean.component_qcs.at(c_qc_it).overall_quality_l = filter_mean.component_qcs.at(c_qc_it).overall_quality_l / filter_values.size();
      filter_mean.component_qcs.at(c_qc_it).overall_quality_u = filter_mean.component_qcs.at(c_qc_it).overall_quality_u / filter_values.size();
      for (auto& kv : filter_mean.component_qcs.at(c_qc_it).meta_value_qc) {
        kv.second.first = filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first / filter_values.size();
        kv.second.second = filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second / filter_values.size();
      }
    }
  }

  void MRMFeatureFilter::calculateFilterValuesVar(MRMFeatureQC& filter_var, const std::vector<MRMFeatureQC>& filter_values, const MRMFeatureQC& filter_mean, const MRMFeatureQC& filter_template) const
  {
    // Determine the STD for each filter_template value
    zeroFilterValues(filter_var, filter_template);
    for (const MRMFeatureQC& filter : filter_values) { // Accumulate the squared sum of the difference from the mean
      for (size_t cg_qc_it = 0; cg_qc_it < filter_var.component_group_qcs.size(); ++cg_qc_it) {
        filter_var.component_group_qcs.at(cg_qc_it).retention_time_l += std::pow(filter.component_group_qcs.at(cg_qc_it).retention_time_l - filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).retention_time_u += std::pow(filter.component_group_qcs.at(cg_qc_it).retention_time_u - filter_mean.component_group_qcs.at(cg_qc_it).retention_time_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).intensity_l += std::pow(filter.component_group_qcs.at(cg_qc_it).intensity_l - filter_mean.component_group_qcs.at(cg_qc_it).intensity_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).intensity_u += std::pow(filter.component_group_qcs.at(cg_qc_it).intensity_u - filter_mean.component_group_qcs.at(cg_qc_it).intensity_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).overall_quality_l += std::pow(filter.component_group_qcs.at(cg_qc_it).overall_quality_l - filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).overall_quality_u += std::pow(filter.component_group_qcs.at(cg_qc_it).overall_quality_u - filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_heavy_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_heavy_l - filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_heavy_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_heavy_u - filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_light_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_light_l - filter_mean.component_group_qcs.at(cg_qc_it).n_light_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_light_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_light_u - filter_mean.component_group_qcs.at(cg_qc_it).n_light_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_detecting_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_detecting_l - filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_detecting_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_detecting_u - filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_quantifying_l - filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_quantifying_u - filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_identifying_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_identifying_l - filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_identifying_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_identifying_u - filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_transitions_l += std::pow(filter.component_group_qcs.at(cg_qc_it).n_transitions_l - filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).n_transitions_u += std::pow(filter.component_group_qcs.at(cg_qc_it).n_transitions_u - filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u, 2);
        filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_l += std::pow(filter.component_group_qcs.at(cg_qc_it).ion_ratio_l - filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l, 2);
        filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_u += std::pow(filter.component_group_qcs.at(cg_qc_it).ion_ratio_u - filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u, 2);
        for (auto& kv : filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc) {
          kv.second.first += std::pow(filter.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first - filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first, 2);
          kv.second.second += std::pow(filter.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second - filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second, 2);
        }
      }
      for (size_t c_qc_it = 0; c_qc_it < filter_var.component_qcs.size(); ++c_qc_it) {
        filter_var.component_qcs.at(c_qc_it).retention_time_l += std::pow(filter.component_qcs.at(c_qc_it).retention_time_l - filter_mean.component_qcs.at(c_qc_it).retention_time_l, 2);
        filter_var.component_qcs.at(c_qc_it).retention_time_u += std::pow(filter.component_qcs.at(c_qc_it).retention_time_u - filter_mean.component_qcs.at(c_qc_it).retention_time_u, 2);
        filter_var.component_qcs.at(c_qc_it).intensity_l += std::pow(filter.component_qcs.at(c_qc_it).intensity_l - filter_mean.component_qcs.at(c_qc_it).intensity_l, 2);
        filter_var.component_qcs.at(c_qc_it).intensity_u += std::pow(filter.component_qcs.at(c_qc_it).intensity_u - filter_mean.component_qcs.at(c_qc_it).intensity_u, 2);
        filter_var.component_qcs.at(c_qc_it).overall_quality_l += std::pow(filter.component_qcs.at(c_qc_it).overall_quality_l - filter_mean.component_qcs.at(c_qc_it).overall_quality_l, 2);
        filter_var.component_qcs.at(c_qc_it).overall_quality_u += std::pow(filter.component_qcs.at(c_qc_it).overall_quality_u - filter_mean.component_qcs.at(c_qc_it).overall_quality_u, 2);
        for (auto& kv : filter_var.component_qcs.at(c_qc_it).meta_value_qc) {
          kv.second.first += std::pow(filter.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first - filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first, 2);
          kv.second.second += std::pow(filter.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second - filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second, 2);
        }
      }
    }
    for (size_t cg_qc_it = 0; cg_qc_it < filter_var.component_group_qcs.size(); ++cg_qc_it) {// Divide by the size (performed separately due to int types...)
      filter_var.component_group_qcs.at(cg_qc_it).retention_time_l = filter_var.component_group_qcs.at(cg_qc_it).retention_time_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).retention_time_u = filter_var.component_group_qcs.at(cg_qc_it).retention_time_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).intensity_l = filter_var.component_group_qcs.at(cg_qc_it).intensity_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).intensity_u = filter_var.component_group_qcs.at(cg_qc_it).intensity_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).overall_quality_l = filter_var.component_group_qcs.at(cg_qc_it).overall_quality_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).overall_quality_u = filter_var.component_group_qcs.at(cg_qc_it).overall_quality_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_heavy_l = filter_var.component_group_qcs.at(cg_qc_it).n_heavy_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_heavy_u = filter_var.component_group_qcs.at(cg_qc_it).n_heavy_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_light_l = filter_var.component_group_qcs.at(cg_qc_it).n_light_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_light_u = filter_var.component_group_qcs.at(cg_qc_it).n_light_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_detecting_l = filter_var.component_group_qcs.at(cg_qc_it).n_detecting_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_detecting_u = filter_var.component_group_qcs.at(cg_qc_it).n_detecting_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_l = filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_u = filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_identifying_l = filter_var.component_group_qcs.at(cg_qc_it).n_identifying_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_identifying_u = filter_var.component_group_qcs.at(cg_qc_it).n_identifying_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_transitions_l = filter_var.component_group_qcs.at(cg_qc_it).n_transitions_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).n_transitions_u = filter_var.component_group_qcs.at(cg_qc_it).n_transitions_u / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_l = filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_l / (filter_values.size() - 1);
      filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_u = filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_u / (filter_values.size() - 1);
      for (auto& kv : filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc) {
        kv.second.first = filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first / (filter_values.size() - 1);
        kv.second.second = filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second / (filter_values.size() - 1);
      }
    }
    for (size_t c_qc_it = 0; c_qc_it < filter_var.component_qcs.size(); ++c_qc_it) {
      filter_var.component_qcs.at(c_qc_it).retention_time_l = filter_var.component_qcs.at(c_qc_it).retention_time_l / (filter_values.size() - 1);
      filter_var.component_qcs.at(c_qc_it).retention_time_u = filter_var.component_qcs.at(c_qc_it).retention_time_u / (filter_values.size() - 1);
      filter_var.component_qcs.at(c_qc_it).intensity_l = filter_var.component_qcs.at(c_qc_it).intensity_l / (filter_values.size() - 1);
      filter_var.component_qcs.at(c_qc_it).intensity_u = filter_var.component_qcs.at(c_qc_it).intensity_u / (filter_values.size() - 1);
      filter_var.component_qcs.at(c_qc_it).overall_quality_l = filter_var.component_qcs.at(c_qc_it).overall_quality_l / (filter_values.size() - 1);
      filter_var.component_qcs.at(c_qc_it).overall_quality_u = filter_var.component_qcs.at(c_qc_it).overall_quality_u / (filter_values.size() - 1);
      for (auto& kv : filter_var.component_qcs.at(c_qc_it).meta_value_qc) {
        kv.second.first = filter_var.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first / (filter_values.size() - 1);
        kv.second.second = filter_var.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second / (filter_values.size() - 1);
      }
    }
  }

  void MRMFeatureFilter::calculateFilterValuesPercRSD(MRMFeatureQC& filter_rsd, const MRMFeatureQC& filter_mean, const MRMFeatureQC& filter_var) const
  {
    // Determine the %RSD for each filter_rsd value
    for (size_t cg_qc_it = 0; cg_qc_it < filter_rsd.component_group_qcs.size(); ++cg_qc_it) {
      filter_rsd.component_group_qcs.at(cg_qc_it).retention_time_l = (filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).retention_time_l) / filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).retention_time_u = (filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).retention_time_l) / filter_mean.component_group_qcs.at(cg_qc_it).retention_time_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).intensity_l = (filter_mean.component_group_qcs.at(cg_qc_it).intensity_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).intensity_l) / filter_mean.component_group_qcs.at(cg_qc_it).intensity_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).intensity_u = (filter_mean.component_group_qcs.at(cg_qc_it).intensity_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).intensity_u) / filter_mean.component_group_qcs.at(cg_qc_it).intensity_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).overall_quality_l = (filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).overall_quality_l) / filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).overall_quality_u = (filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).overall_quality_u) / filter_mean.component_group_qcs.at(cg_qc_it).overall_quality_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_heavy_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_heavy_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_heavy_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_heavy_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_heavy_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_light_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_light_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_light_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_light_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_light_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_light_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_light_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_light_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_detecting_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_detecting_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_detecting_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_detecting_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_detecting_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_quantifying_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_quantifying_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_quantifying_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_quantifying_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_identifying_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_identifying_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_identifying_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_identifying_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_identifying_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_transitions_l = (filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_transitions_l) / filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).n_transitions_u = (filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).n_transitions_u) / filter_mean.component_group_qcs.at(cg_qc_it).n_transitions_u * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).ion_ratio_l = (filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_l) / filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_l * 100 : 0;
      filter_rsd.component_group_qcs.at(cg_qc_it).ion_ratio_u = (filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).ion_ratio_u) / filter_mean.component_group_qcs.at(cg_qc_it).ion_ratio_u * 100 : 0;
      for (auto& kv : filter_rsd.component_group_qcs.at(cg_qc_it).meta_value_qc) {
        kv.second.first = (filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first) / filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).first * 100 : 0;
        kv.second.second = (filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second != 0) ? std::sqrt(filter_var.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second) / filter_mean.component_group_qcs.at(cg_qc_it).meta_value_qc.at(kv.first).second * 100 : 0;
      }
    }
    for (size_t c_qc_it = 0; c_qc_it < filter_rsd.component_qcs.size(); ++c_qc_it) {
      filter_rsd.component_qcs.at(c_qc_it).retention_time_l = (filter_mean.component_qcs.at(c_qc_it).retention_time_l != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).retention_time_l) / filter_mean.component_qcs.at(c_qc_it).retention_time_l * 100 : 0;
      filter_rsd.component_qcs.at(c_qc_it).retention_time_u = (filter_mean.component_qcs.at(c_qc_it).retention_time_u != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).retention_time_u) / filter_mean.component_qcs.at(c_qc_it).retention_time_u * 100 : 0;
      filter_rsd.component_qcs.at(c_qc_it).intensity_l = (filter_mean.component_qcs.at(c_qc_it).intensity_l != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).intensity_l) / filter_mean.component_qcs.at(c_qc_it).intensity_l * 100 : 0;
      filter_rsd.component_qcs.at(c_qc_it).intensity_u = (filter_mean.component_qcs.at(c_qc_it).intensity_u != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).intensity_u) / filter_mean.component_qcs.at(c_qc_it).intensity_u * 100 : 0;
      filter_rsd.component_qcs.at(c_qc_it).overall_quality_l = (filter_mean.component_qcs.at(c_qc_it).overall_quality_l != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).overall_quality_l) / filter_mean.component_qcs.at(c_qc_it).overall_quality_l * 100 : 0;
      filter_rsd.component_qcs.at(c_qc_it).overall_quality_u = (filter_mean.component_qcs.at(c_qc_it).overall_quality_u != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).overall_quality_u) / filter_mean.component_qcs.at(c_qc_it).overall_quality_u * 100 : 0;
      for (auto& kv : filter_rsd.component_qcs.at(c_qc_it).meta_value_qc) {
        kv.second.first = (filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first) / filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).first * 100 : 0;
        kv.second.second = (filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second != 0) ? std::sqrt(filter_var.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second) / filter_mean.component_qcs.at(c_qc_it).meta_value_qc.at(kv.first).second * 100 : 0;
      }
    }
  }

  template <typename T>
  bool MRMFeatureFilter::checkRange(const T& value, const T& value_l, const T& value_u) const
  {
    return value >= value_l&& value <= value_u;
  }
  template<typename T>
  void MRMFeatureFilter::updateRange(const T& value, T& value_l, T& value_u) const
  {
    if (value < value_l) value_l = value;
    if (value > value_u) value_u = value;
  }
  template<typename T>
  void MRMFeatureFilter::setRange(const T& value, T& value_l, T& value_u) const
  {
    if (value >= T(0)) {
      value_l = T(0);
      value_u = value;
    } else {
      value_l = value;
      value_u = T(0);
    }
  }
  template<typename T>
  void MRMFeatureFilter::initRange(const T& value, T& value_l, T& value_u) const
  {
    value_l = value;
    value_u = value;
  }
}
