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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <algorithm>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>

namespace OpenMS
{
  // bool comp(const std::pair<double, String>& a, const std::pair<double, String>& b) {
  //   return a.first < b.first;
  // }

        

  MRMFeatureSelector::MRMFeatureSelector() :
    DefaultParamHandler("MRMFeatureSelector")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  MRMFeatureSelector::~MRMFeatureSelector() {}

  void MRMFeatureSelector::optimize_Tr() {
  }

  void MRMFeatureSelector::optimize_score() {}

  void MRMFeatureSelector::select_MRMFeature_qmip(
    FeatureMap& features,
    TargetedExperiment& targeted_exp
  )
  {
    std::cout << "=======START========" << std::endl;
    std::vector<std::pair<double, String>> time_to_name;
    std::map< String, std::vector<std::map<String, DataValue>> > feature_name_map;
    size_t feature_count = 0;
    for (FeatureMap::const_iterator it = features.begin(); it != features.end(); ++it) {
      String component_group_name = it->getMetaValue("PeptideRef").toString();
      double retention_time = it->getRT();
      double assay_retention_time = it->getMetaValue("assay_rt");
      UInt64 transition_id = it->getUniqueId();
      std::vector<String> keys;
      it->getKeys(keys);
      std::map<String, DataValue> feature_properties {
        {"retention_time", retention_time},
        {"transition_id", transition_id},
        {"component_group_name", component_group_name},
        {"component_name", component_group_name},
      };
      for (Size i = 0; i < keys.size(); i++) {
        feature_properties[keys[i]] = it->getMetaValue(keys[i]);
      }
      time_to_name.push_back(std::make_pair(assay_retention_time, component_group_name));
      if (feature_name_map.find(component_group_name) == feature_name_map.end()) {
        feature_name_map[component_group_name] = std::vector<std::map<String, DataValue>>();
      }
      feature_name_map[component_group_name].push_back(feature_properties);
      ++feature_count;
      if (!getSelectTransitionGroup()) {
        for (std::vector<Feature>::const_iterator sub_it = it->getSubordinates().begin();
            sub_it != it->getSubordinates().end(); ++sub_it) {
          String component_name = sub_it->getMetaValue("native_id").toString();
          time_to_name.push_back(std::make_pair(assay_retention_time, component_name));
          std::map<String, DataValue> subfeature_properties;
          subfeature_properties = feature_properties;
          subfeature_properties["component_name"] = component_name;
          std::vector<String> subkeys;
          sub_it->getKeys(subkeys);
          for (Size i = 0; i < subkeys.size(); i++) {
            subfeature_properties[subkeys[i]] = sub_it->getMetaValue(subkeys[i]);
          }
          if (feature_name_map.find(component_name) == feature_name_map.end()) {
            feature_name_map[component_name] = std::vector<std::map<String, DataValue>>();
          }
          feature_name_map[component_name].push_back(subfeature_properties);
          ++feature_count;
        }
      }
    }
    std::cout << feature_count << " features detected" << std::endl;
    sort(time_to_name.begin(), time_to_name.end());
    double window_length = getSegmentWindowLength();
    double step_length = getSegmentWindowLength();
    if (window_length == -1 && step_length == -1) {
      window_length = time_to_name.size();
      step_length = time_to_name.size();
    }
    size_t n_segments = std::ceil(time_to_name.size() / step_length) ;
    std::cout << n_segments << " SEGMENTS" << std::endl;
    for (size_t i=0; i < n_segments; ++i) {
      size_t start = step_length*i;
      size_t end = std::min(start + window_length, (double)time_to_name.size());
      std::vector<std::pair<double, String>> time_slice(time_to_name.begin() + start, time_to_name.begin() + end);
      for (const std::pair<double, String>& j : time_slice) {
        std::cout << j.first << " " << j.second << std::endl;
      }
      std::cout << i << " SEGMENT!!!" << std::endl;
    }
  }

  void MRMFeatureSelector::select_MRMFeature_score() {}

  double MRMFeatureSelector::make_score(
    Feature& feature,
    String& metaValue,
    double& weight // maybe others
  )
  {}

  void MRMFeatureSelector::setNNThreshold(const double& nn_threshold)
  {
    nn_threshold_ = nn_threshold;
  }

  double MRMFeatureSelector::getNNThreshold() const
  {
    return nn_threshold_;
  }

  void MRMFeatureSelector::setLocalityWeight(const bool& locality_weight)
  {
    locality_weight_ = locality_weight;
  }

  bool MRMFeatureSelector::getLocalityWeight() const
  {
    return locality_weight_;
  }

  void MRMFeatureSelector::setSelectTransitionGroup(const bool& select_transition_group)
  {
    select_transition_group_ = select_transition_group;
  }

  bool MRMFeatureSelector::getSelectTransitionGroup() const
  {
    return select_transition_group_;
  }

  void MRMFeatureSelector::setSegmentWindowLength(const double& segment_window_length)
  {
    segment_window_length_ = segment_window_length;
  }

  double MRMFeatureSelector::getSegmentWindowLength() const
  {
    return segment_window_length_;
  }

  void MRMFeatureSelector::setSegmentStepLength(const double& segment_step_length)
  {
    segment_step_length_ = segment_step_length;
  }

  double MRMFeatureSelector::getSegmentStepLength() const
  {
    return segment_step_length_;
  }

  void MRMFeatureSelector::setSelectHighestCount(const bool& select_highest_count)
  {
    select_highest_count_ = select_highest_count;
  }

  bool MRMFeatureSelector::getSelectHighestCount() const
  {
    return select_highest_count_;
  }

  void MRMFeatureSelector::setVariableType(const String& variable_type)
  {
    variable_type_ = variable_type;
  }

  String MRMFeatureSelector::getVariableType() const
  {
    return variable_type_;
  }

  void MRMFeatureSelector::setOptimalThreshold(const double& optimal_threshold)
  {
    optimal_threshold_ = optimal_threshold;
  }

  double MRMFeatureSelector::getOptimalThreshold() const
  {
    return optimal_threshold_;
  }

  void MRMFeatureSelector::getDefaultParameters(Param& params)
  {
    params.clear();
    // TODO Adjust defaults
    // TODO set limits on parameters
    params.setValue("nn_threshold", 4.0);
    params.setValue("locality_weight", "false");
    params.setValue("select_transition_group", "true");
    params.setValue("segment_window_length", 8.0);
    params.setValue("segment_step_length", 4.0);
    params.setValue("select_highest_count", "false");
    params.setValue("variable_type", "continuous");
    params.setValue("optimal_threshold", 0.5);
  }

  void MRMFeatureSelector::updateMembers_()
  {
    nn_threshold_ = (double)param_.getValue("nn_threshold");
    locality_weight_ = param_.getValue("locality_weight").toBool();
    select_transition_group_ = param_.getValue("select_transition_group").toBool();
    segment_window_length_ = (double)param_.getValue("segment_window_length");
    segment_step_length_ = (double)param_.getValue("segment_step_length");
    select_highest_count_ = param_.getValue("select_highest_count").toBool();
    variable_type_ = (String)param_.getValue("variable_type");
    optimal_threshold_ = (double)param_.getValue("optimal_threshold");
  }
}
