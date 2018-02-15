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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskeyt $
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
    defaults_.setValue("flag_or_filter", "flag", "Flag or Filter (i.e., remove) Components or transitions that do not pass the QC.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("flag_or_filter", ListUtils::create<String>("flag,filter"));
    
    //TODO:  future implementation for QcML reporting
    defaults_.setValue("report_xic", "false", "Embed an image of the XIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_xic", ListUtils::create<String>("true,false"));

    //TODO:  future implementation for QcML reporting
    defaults_.setValue("report_tic", "false", "Embed an image of the TIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_tic", ListUtils::create<String>("true,false"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMFeatureFilter::~MRMFeatureFilter()
  {
  }

  void MRMFeatureFilter::updateMembers_()
  {
    flag_or_filter_ = (String)param_.getValue("flag_or_filter");
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

    String delim = ";";

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
      std::vector<String> cg_qc_fail_message_vec;

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features[feature_it].getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)features[feature_it].getSubordinates()[sub_it].getMetaValue("native_id"); 
        // std::cout << "component_name" << component_name << std::endl; //debugging
        bool c_qc_pass = true;
         std::vector<String> c_qc_fail_message_vec;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs.size(); ++cg_qc_it)
        {
          if (filter_criteria.component_group_qcs[cg_qc_it].component_group_name == component_group_name)
          {
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
              }
            }
          }
        }
        
        // iterate through feature/sub-feature QCs/filters        
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs[c_qc_it].component_name == component_name)
          {
            // RT check
            double rt = features[feature_it].getSubordinates()[sub_it].getRT(); //check!
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

            // metaValue checks
            for (auto const& kv : filter_criteria.component_qcs[c_qc_it].meta_value_qc)
            {
              // std::cout << "MetaData" << std::endl; //debugging
              if (!checkMetaValue(features[feature_it].getSubordinates()[sub_it], kv.first, kv.second.first, kv.second.second))
              {
                c_qc_pass = false;
                c_qc_fail_message_vec.push_back("metaValue[" + kv.first + "]");
              }
            }
          }
        }

        // Copy or Flag passing/failing subordinates
        if (c_qc_pass && flag_or_filter_ == "filter")
        {
          // std::cout << "copied passing subordinate" << std::endl; //debugging
          subordinates_filtered.push_back(features[feature_it].getSubordinates()[sub_it]);
        }
        else if (c_qc_pass && flag_or_filter_ == "flag")
        {
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_pass", true);
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_message", "");
        }
        else if (!c_qc_pass && flag_or_filter_ == "filter")
        {
          // do nothing
          // std::cout << "omitted failing subordinate" << std::endl; //debugging
        }
        else if (!c_qc_pass && flag_or_filter_ == "flag")
        {
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_pass", false);
          String c_qc_fail_message = uniqueJoin(c_qc_fail_message_vec, delim);
          features[feature_it].getSubordinates()[sub_it].setMetaValue("QC_transition_message", c_qc_fail_message);
        }
      }

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
        features[feature_it].setMetaValue("QC_transition_group_message", "");
      }
      else if (!cg_qc_pass && flag_or_filter_ == "filter")
      {
        // do nothing
        // std::cout << "omitted failing feature" << std::endl; //debugging
      }   
      else if (!cg_qc_pass && flag_or_filter_ == "flag")
      {
        features[feature_it].setMetaValue("QC_transition_group_pass", false);
        String cg_qc_fail_message = uniqueJoin(cg_qc_fail_message_vec, delim);
        features[feature_it].setMetaValue("QC_transition_group_message", cg_qc_fail_message);
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
    if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
    {
      double feature_1 = component_1.getMetaValue(feature_name);
      double feature_2 = component_2.getMetaValue(feature_name);
      ratio = feature_1/feature_2;
      
    } 
    else if (component_1.metaValueExists(feature_name))
    {
      LOG_INFO << "Warning: no ion pair found for transition_id " << component_1.getMetaValue("native_id") << ".";
      double feature_1 = component_1.getMetaValue(feature_name);
      ratio = feature_1;
    } 
    else
    {
      LOG_INFO << "Feature metaValue " << feature_name << " not found for transition_ids " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
    }

    return ratio;
  }
  
  bool MRMFeatureFilter::checkMetaValue(const Feature & component, const String & meta_value_key, const double & meta_value_l, const double & meta_value_u)
  {
    bool check = true;
    if (component.metaValueExists(meta_value_key))
    {
      double meta_value = (double)component.getMetaValue(meta_value_key);
      check = checkRange(meta_value,
        meta_value_l,
        meta_value_u);
    }
    else 
    {
      LOG_INFO << "Warning: no metaValue found for transition_id " << component.getMetaValue("native_id") << " for metaValue key " << meta_value_key << ".";
    }

    return check;
  }

  String MRMFeatureFilter::uniqueJoin(std::vector<String>& str_vec, String& delim)
  {
    //remove duplicates
    std::sort(str_vec.begin(), str_vec.end());
    str_vec.erase(std::unique(str_vec.begin(), str_vec.end()), str_vec.end());
    //concatenate
    String str_cat = "";
    for (auto str : str_vec)
    {
      str_cat = str_cat + str + delim;
    }
    // std::cout << str_cat << std::endl; //debugging
    //remove trailing delimm
    if (str_cat != "")
    {
      str_cat = str_cat.substr(0, str_cat.length() - delim.length());
    }
    // std::cout << str_cat << std::endl; //debugging
    return str_cat;
  }

  template <typename T>
  bool MRMFeatureFilter::checkRange(T const& value, T const& value_l, T const& value_u)
  {
    bool range_check = true;
    if (value < value_l
      || value > value_u)
    {
      range_check = false;
    }
    // std::cout << "value: " << (String)value << " lb: " << (String)value_l << " ub: " << (String)value_u << std::endl; //debugging
    return range_check;
  }
  
  //TODO: Future addition to allow for generating a QcML attachment for QC reporting
  // void MRMFeatureFilter::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  // {
  //   //TODO
  // }

}

