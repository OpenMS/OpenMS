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
    
    defaults_.setValue("report_xic", "false", "Embed an image of the XIC in the QC report.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("report_xic", ListUtils::create<String>("true,false"));

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

  void MRMFeatureFilter::FilterFeatureMap(FeatureMap& features, MRMFeatureQC& filter_criteria)
  { 
    // initialize the new feature map
    if (flag_or_filter_ == "filter")
    {
      FeatureMap features_filtered;
    }
    
    // initialize QC variables
    std::map<String,MRMFeatureQC>::iterator feature_qc_it;

    // initialize variables
    String component_name; //i.e., transition_id
    String IS_component_name; //i.e., internal standard transition_id
    String component_group_name; //i.e., peptideRef
    double calculated_concentration;
    bool qc_pass;
    String concentration_units;// iterate through each component_group/feature     

    for (size_t feature_it = 0; feature_it < features.size(); ++feature_it)
    {
      component_group_name = (String)features[feature_it].getMetaValue("PeptideRef");

      std::map<String,int> labels_and_transition_types = countLabelsAndTransitionTypes(features[feature_it]);

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < features[feature_it].getSubordinates().size(); ++sub_it)
      {
        component_name = (String)features[feature_it].getSubordinates()[sub_it].getMetaValue("native_id"); 
        qc_pass = true;

        // iterate through multi-feature/multi-sub-feature QCs/filters
        // iterate through component_groups
        if (sub_it == 0)
        {
          for (size_t cg_qc_it = 0; cg_qc_it < filter_criteria.component_group_qcs_.size(); ++cg_qc_it)
          {
            if (filter_criteria.component_group_qcs_[cg_qc_it].component_group_name_ == component_group_name)
            {
              // labels and transition counts QC
              if (labels_and_transition_types["n_heavy"] < filter_criteria.component_group_qcs_[cg_qc_it].n_heavy_l_
                || labels_and_transition_types["n_heavy"] > filter_criteria.component_group_qcs_[cg_qc_it].n_heavy_u_)
              {
                qc_pass = false;
              }
              if (labels_and_transition_types["n_light"] < filter_criteria.component_group_qcs_[cg_qc_it].n_light_l_
                || labels_and_transition_types["n_light"] > filter_criteria.component_group_qcs_[cg_qc_it].n_light_u_)
              {
                qc_pass = false;
              }
              if (labels_and_transition_types["n_detecting"] < filter_criteria.component_group_qcs_[cg_qc_it].n_detecting_l_
                || labels_and_transition_types["n_detecting"] > filter_criteria.component_group_qcs_[cg_qc_it].n_detecting_u_)
              {
                qc_pass = false;
              }
              if (labels_and_transition_types["n_quantifying"] < filter_criteria.component_group_qcs_[cg_qc_it].n_quantifying_l_
                || labels_and_transition_types["n_quantifying"] > filter_criteria.component_group_qcs_[cg_qc_it].n_quantifying_u_)
              {
                qc_pass = false;
              }
              if (labels_and_transition_types["n_identifying"] < filter_criteria.component_group_qcs_[cg_qc_it].n_identifying_l_
                || labels_and_transition_types["n_identifying"] > filter_criteria.component_group_qcs_[cg_qc_it].n_identifying_u_)
              {
                qc_pass = false;
              }
              if (labels_and_transition_types["n_transitions"] < filter_criteria.component_group_qcs_[cg_qc_it].n_transitions_l_
                || labels_and_transition_types["n_transitions"] > filter_criteria.component_group_qcs_[cg_qc_it].n_transitions_u_)
              {
                qc_pass = false;
              }

              // ion ratio QC
              for (size_t sub_it2 = 0; sub_it2 < features[feature_it].getSubordinates().size(); ++sub_it2)
              {
                component_name2 = (String)features[feature_it].getSubordinates()[sub_it2].getMetaValue("native_id"); 

                // find the ion ratio pair
                if (filter_criteria.component_group_qcs_[cg_qc_it].ion_ratio_pair_name_1_ == component_name
                  && filter_criteria.component_group_qcs_[cg_qc_it].ion_ratio_pair_name_2_ == component_name2)
                {
                  double ion_ratio = calculateIonRatio(features[feature_it].getSubordinates()[sub_it], features[feature_it].getSubordinates()[sub_it2], filter_criteria.component_group_qcs_[cg_qc_it].ion_ratio_feature_name_);
                  
                  if (ion_ratio < filter_criteria.component_group_qcs_[cg_qc_it].ion_ratio_l_
                  || ion_ratio > filter_criteria.component_group_qcs_[cg_qc_it].ion_ratio_u_)
                  {
                    qc_pass = false;
                  }
                }
              }
            }
          }
        }
        // iterate through feature/sub-feature QCs/filters        
        for (size_t c_qc_it = 0; c_qc_it < filter_criteria.component_qcs_.size(); ++c_qc_it)
        {
          if (filter_criteria.component_qcs_[c_qc_it].component_name_ == component_name)
          {
            // RT check
            double rt = features[feature_it].getSubordinates()[sub_it].getRT(); //check!
            if (rt < filter_criteria.component_qcs_[c_qc_it].retention_time_l_
              && rt > filter_criteria.component_qcs_[c_qc_it].retention_time_u_)
              {
                qc_pass = false;
              }

            // intensity check
            double intensity = features[feature_it].getSubordinates()[sub_it].getIntensity();
            if (intensity < filter_criteria.component_qcs_[c_qc_it].intensity_l_
              && intensity > filter_criteria.component_qcs_[c_qc_it].intensity_u_)
              {
                qc_pass = false;
              }

            // overall quality check getQuality
            double quality = features[feature_it].getSubordinates()[sub_it].getQuality();
            if (quality < filter_criteria.component_qcs_[c_qc_it].overall_quality_l_
              && quality > filter_criteria.component_qcs_[c_qc_it].overall_quality_u_)
              {
                qc_pass = false;
              }

            // metaValue checks
            for (auto const& kv : filter_criteria.component_qcs_[c_qc_it].meta_value_qc_)
            {
              if (!checkMetaValue(features[feature_it].getSubordinates()[sub_it], kv.first, kv.second.first, kv.second.second))
              {
                qc_pass = false;
              }
            }
          }
        }
      }


    }
  }
  
  std::map<String,int> MRMFeatureFilter::countLabelsAndTransitionTypes(Feature & component_group)
  {
    //TODO
  }
  
  double MRMFeatureFilter::calculateIonRatio(Feature & component_1, Feature & component_2, String & feature_name)
  {
    //TODO
  }
  
  bool MRMFeatureFilter::checkMetaValue(Feature & component, String & meta_value_key, String & meta_value_l, String & meta_value_u)
  {
    //TODO
  }
  
  void MRMFeatureFilter::FeatureMapToAttachment(FeatureMap& features, QcMLFile::Attachment& attachment)
  {
    //TODO
  }

}

