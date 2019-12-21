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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFilter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>


#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/QcMLFile.h>

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

  MRMFeatureFilter::~MRMFeatureFilter()
  {
  }

  void MRMFeatureFilter::getDefaultParameters(Param& params) const
  {
    params.clear();

    params.setValue("flag_or_filter", "flag", "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC.", ListUtils::create<String>("advanced"));
    params.setValidStrings("flag_or_filter", ListUtils::create<String>("flag,filter"));

    //TODO: future implementation for QcML reporting
    params.setValue("report_xic", "false", "Embed an image of the XIC in the QC report.", ListUtils::create<String>("advanced"));
    params.setValidStrings("report_xic", ListUtils::create<String>("true,false"));

    //TODO: future implementation for QcML reporting
    params.setValue("report_tic", "false", "Embed an image of the TIC in the QC report.", ListUtils::create<String>("advanced"));
    params.setValidStrings("report_tic", ListUtils::create<String>("true,false"));
  }

  void MRMFeatureFilter::updateMembers_()
  {
    flag_or_filter_ = param_.getValue("flag_or_filter").toString();
    report_xic_ = param_.getValue("report_xic").toBool();
    report_tic_ = param_.getValue("report_tic").toBool();
  }

  void MRMFeatureFilter::FilterFeatureMap(FeatureMap& features,
    const MRMFeatureQC& filter_criteria,
    const TargetedExperiment & transitions
  )
  {
    // initialize QC variables
    FeatureMap features_filtered;

    // bool qc_pass;
    String concentration_units;// iterate through each component_group/feature

    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      String component_group_name = (String)features[feature_it].getMetaValue("PeptideRef");
      // std::cout << "component_group_name" << component_group_name << std::endl; //debugging

      std::map<String,int> labels_and_transition_types = countLabelsAndTransitionTypes(features[feature_it], transitions);

      // initialize the new feature and subordinates
      std::vector<Feature> subordinates_filtered;
      bool cg_qc_pass = true;
      StringList cg_qc_fail_message_vec;
      UInt cg_tests_count{0};

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features[feature_it].getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)features[feature_it].getSubordinates()[sub_it].getMetaValue("native_id");
        // std::cout << "component_name" << component_name << std::endl; //debugging
        bool c_qc_pass = true;
        StringList c_qc_fail_message_vec;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs.size(); ++cg_qc_it)
        {
          if (filter_criteria.component_group_qcs[cg_qc_it].component_group_name == component_group_name)
          {
            const double rt = features[feature_it].getRT();
            if (!checkRange(rt,
              filter_criteria.component_group_qcs[cg_qc_it].retention_time_l,
              filter_criteria.component_group_qcs[cg_qc_it].retention_time_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("retention_time");
            }

            const double intensity = features[feature_it].getIntensity();
            if (!checkRange(intensity,
              filter_criteria.component_group_qcs[cg_qc_it].intensity_l,
              filter_criteria.component_group_qcs[cg_qc_it].intensity_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("intensity");
            }

            const double quality = features[feature_it].getOverallQuality();
            if (!checkRange(quality,
              filter_criteria.component_group_qcs[cg_qc_it].overall_quality_l,
              filter_criteria.component_group_qcs[cg_qc_it].overall_quality_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("overall_quality");
            }
            // labels and transition counts QC
            // std::cout << "n_heavy" << std::endl; //debugging
            if (!checkRange(labels_and_transition_types["n_heavy"],
              filter_criteria.component_group_qcs[cg_qc_it].n_heavy_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_heavy_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_heavy");
            }
            // std::cout << "n_light" << std::endl; //debugging
            if (! checkRange(labels_and_transition_types["n_light"],
              filter_criteria.component_group_qcs[cg_qc_it].n_light_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_light_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_light");
            }
            // std::cout << "n_detecting" << std::endl; //debugging
            if (! checkRange(labels_and_transition_types["n_detecting"],
              filter_criteria.component_group_qcs[cg_qc_it].n_detecting_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_detecting_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_detecting");
            }
            // std::cout << "n_quantifying" << std::endl; //debugging
            if (! checkRange(labels_and_transition_types["n_quantifying"],
              filter_criteria.component_group_qcs[cg_qc_it].n_quantifying_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_quantifying_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_quantifying");
            }
            // std::cout << "n_identifying" << std::endl; //debugging
            if (! checkRange(labels_and_transition_types["n_identifying"],
              filter_criteria.component_group_qcs[cg_qc_it].n_identifying_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_identifying_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_identifying");
            }
            // std::cout << "n_transitions" << std::endl; //debugging
            if (! checkRange(labels_and_transition_types["n_transitions"],
              filter_criteria.component_group_qcs[cg_qc_it].n_transitions_l,
              filter_criteria.component_group_qcs[cg_qc_it].n_transitions_u))
            {
              cg_qc_pass = false;
              cg_qc_fail_message_vec.push_back("n_transitions");
            }

            cg_tests_count += 9;

            // ion ratio QC
            for (size_t sub_it2 = 0; sub_it2 < features[feature_it].getSubordinates().size(); ++sub_it2)
            {
              String component_name2 = (String)features[feature_it].getSubordinates()[sub_it2].getMetaValue("native_id");
              // find the ion ratio pair
              if (filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_pair_name_1 != ""
                && filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_pair_name_2 != ""
                && filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_pair_name_1 == component_name
                && filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_pair_name_2 == component_name2)
              {
                double ion_ratio = calculateIonRatio(features[feature_it].getSubordinates()[sub_it], features[feature_it].getSubordinates()[sub_it2], filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_feature_name);

                // std::cout << "ion_ratio" << std::endl; //debugging
                if (! checkRange(ion_ratio,
                  filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_l,
                  filter_criteria.component_group_qcs[cg_qc_it].ion_ratio_u))
                {
                  cg_qc_pass = false;
                  cg_qc_fail_message_vec.push_back("ion_ratio_pair[" + component_name + "/" + component_name2 + "]");
                }
                ++cg_tests_count;
              }
            }

            for (const std::pair<String,std::pair<double,double>>& kv : filter_criteria.component_group_qcs[cg_qc_it].meta_value_qc)
            {
              bool metavalue_exists {false};
              if (!checkMetaValue(features[feature_it], kv.first, kv.second.first, kv.second.second, metavalue_exists))
              {
                cg_qc_pass = false;
                cg_qc_fail_message_vec.push_back(kv.first);
              }
              if (metavalue_exists) ++cg_tests_count;
            }
          }
        }

        UInt c_tests_count{0};
        // iterate through feature/sub-feature QCs/filters
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs[c_qc_it].component_name == component_name)
          {
            // RT check
            const double rt = features[feature_it].getSubordinates()[sub_it].getRT();
            // std::cout << "RT" << std::endl; //debugging
            if (!checkRange(rt,
              filter_criteria.component_qcs[c_qc_it].retention_time_l,
              filter_criteria.component_qcs[c_qc_it].retention_time_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("retention_time");
            }

            // intensity check
            // std::cout << "Intensity" << std::endl; //debugging
            double intensity = features[feature_it].getSubordinates()[sub_it].getIntensity();
            if (!checkRange(intensity,
              filter_criteria.component_qcs[c_qc_it].intensity_l,
              filter_criteria.component_qcs[c_qc_it].intensity_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("intensity");
            }

            // overall quality check getQuality
            double quality = features[feature_it].getSubordinates()[sub_it].getOverallQuality();
            // std::cout << "Quality" << std::endl; //debugging
            if (!checkRange(quality,
              filter_criteria.component_qcs[c_qc_it].overall_quality_l,
              filter_criteria.component_qcs[c_qc_it].overall_quality_u))
            {
              c_qc_pass = false;
              c_qc_fail_message_vec.push_back("overall_quality");
            }

            c_tests_count += 3;

            // metaValue checks
            for (auto const& kv : filter_criteria.component_qcs[c_qc_it].meta_value_qc)
            {
              // std::cout << "MetaData" << std::endl; //debugging
              bool metavalue_exists{false};
              if (!checkMetaValue(features[feature_it].getSubordinates()[sub_it], kv.first, kv.second.first, kv.second.second, metavalue_exists))
              {
                c_qc_pass = false;
                c_qc_fail_message_vec.push_back(kv.first);
              }
              if (metavalue_exists) ++c_tests_count;
            }
          }
        }

        const double c_score = c_tests_count ? 1.0 - c_qc_fail_message_vec.size() / (double)c_tests_count : 1.0;
        features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_score", c_score);

        // Copy or Flag passing/failing subordinates
        if (c_qc_pass && flag_or_filter_ == "filter")
        {
          // std::cout << "copied passing subordinate" << std::endl; //debugging
          subordinates_filtered.push_back(features[feature_it].getSubordinates()[sub_it]);
        }
        else if (c_qc_pass && flag_or_filter_ == "flag")
        {
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_pass", true);
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_message", StringList());
        }
        else if (!c_qc_pass && flag_or_filter_ == "filter")
        {
          // do nothing
          // std::cout << "omitted failing subordinate" << std::endl; //debugging
        }
        else if (!c_qc_pass && flag_or_filter_ == "flag")
        {
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_pass", false);
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_message", getUniqueSorted(c_qc_fail_message_vec));
        }
      }

      const double cg_score = cg_tests_count ? 1.0 - cg_qc_fail_message_vec.size() / (double)cg_tests_count : 1.0;
      features[feature_it].setMetaValue("QC_transition_group_score", cg_score);

      // Copy or Flag passing/failing Features
      if (cg_qc_pass && flag_or_filter_ == "filter" && subordinates_filtered.size() > 0)
      {
        // std::cout << "copied passing feature" << std::endl; //debugging
        Feature feature_filtered(features[feature_it]);
        feature_filtered.setSubordinates(subordinates_filtered);
        features_filtered.push_back(feature_filtered);
      }
      else if (cg_qc_pass && flag_or_filter_ == "filter" && subordinates_filtered.size() == 0)
      {
        // do nothing
        // std::cout << "omitted failing feature" << std::endl; //debugging
      }
      else if (cg_qc_pass && flag_or_filter_ == "flag")
      {
        features[feature_it].setMetaValue("QC_transition_group_pass", true);
        features[feature_it].setMetaValue("QC_transition_group_message", StringList());
      }
      else if (!cg_qc_pass && flag_or_filter_ == "filter")
      {
        // do nothing
        // std::cout << "omitted failing feature" << std::endl; //debugging
      }
      else if (!cg_qc_pass && flag_or_filter_ == "flag")
      {
        features[feature_it].setMetaValue("QC_transition_group_pass", false);
        features[feature_it].setMetaValue("QC_transition_group_message", getUniqueSorted(cg_qc_fail_message_vec));
      }
    }

    // replace with the filtered featureMap
    if (flag_or_filter_ == "filter")
    {
      features = features_filtered;
    }
  }

  std::map<String,int> MRMFeatureFilter::countLabelsAndTransitionTypes(
    const Feature & component_group,
    const TargetedExperiment & transitions)
  {
    int n_heavy(0), n_light(0), n_quant(0), n_detect(0), n_ident(0), n_trans(0);
    std::map<String,int> output;

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

  double MRMFeatureFilter::calculateIonRatio(const Feature & component_1, const Feature & component_2, const String & feature_name)
  {
    double ratio = 0.0;
    // member feature_name access
    if (feature_name == "intensity")
    {
      if (component_1.metaValueExists("native_id") && component_2.metaValueExists("native_id"))
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
      if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
      {
        const double feature_1 = component_1.getMetaValue(feature_name);
        const double feature_2 = component_2.getMetaValue(feature_name);
        ratio = feature_1/feature_2;
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
    const Feature & component,
    const String & meta_value_key,
    const double & meta_value_l,
    const double & meta_value_u,
    bool & key_exists
  ) const
  {
    bool check = true;
    if (component.metaValueExists(meta_value_key))
    {
      key_exists = true;
      const double meta_value = (double)component.getMetaValue(meta_value_key);
      check = checkRange(meta_value, meta_value_l, meta_value_u);
    }
    else
    {
      key_exists = false;
      OPENMS_LOG_DEBUG << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }
    return check;
  }

  StringList MRMFeatureFilter::getUniqueSorted(const StringList& messages) const
  {
    StringList unique {messages};
    std::sort(unique.begin(), unique.end());
    unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
    return unique;
  }

  template <typename T>
  bool MRMFeatureFilter::checkRange(const T& value, const T& value_l, const T& value_u) const
  {
    // std::cout << "value: " << (String)value << " lb: " << (String)value_l << " ub: " << (String)value_u << std::endl; //debugging
    return value >= value_l && value <= value_u;
  }

  //TODO: Future addition to allow for generating a QcML attachment for QC reporting
  // void MRMFeatureFilter::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  // {
  //   //TODO
  // }

}

