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
#include <unordered_set>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>


namespace OpenMS
{
  MRMFeatureSelector::MRMFeatureSelector() :
    DefaultParamHandler("MRMFeatureSelector")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  MRMFeatureSelector::~MRMFeatureSelector() {}

  Int MRMFeatureSelector::_addVariable (LPWrapper& problem, String& name) {
    Int index = problem.addColumn();
    problem.setColumnBounds(index, 0., 1., LPWrapper::DOUBLE_BOUNDED);
    problem.setColumnName(index, name);
    problem.setColumnType(index, LPWrapper::CONTINUOUS);
    problem.setObjective(index, 1.);
    return index;
  }

  void MRMFeatureSelector::_addConstraint(LPWrapper& problem, size_t size, Int *indices_array, double *values_array, String name, double lb, double ub, LPWrapper::Type param) {
    std::vector<Int> indices(indices_array, indices_array + size / sizeof(Int) );
    std::vector<double> values(values_array, values_array + size / sizeof(double) );
    problem.addRow(indices, values, name, lb, ub, param);
  }

  void MRMFeatureSelector::optimize_Tr(
      std::vector<std::pair<double, String>> time_to_name, 
      std::map< String, std::vector<Feature> > feature_name_map
  ) {
    std::cout << "=======START OPTIMIZE TR========" << std::endl;
    std::unordered_set<std::string> variables;
    LPWrapper problem;
    problem.setObjectiveSense(LPWrapper::MIN);
    for (size_t cnt1=0; cnt1 < time_to_name.size(); ++cnt1) {
      std::vector<Feature> feature_row1 = feature_name_map[time_to_name[cnt1].second];
      for (size_t i=0; i < feature_row1.size(); ++i) {
        String name1 = time_to_name[cnt1].second + "_" + (String)feature_row1[i].getUniqueId();
        if (variables.find(name1) == variables.end()) {
            _addVariable(problem, name1);
            variables.insert(name1);
        }
        // TODO: nearest neighbours, not all the components
        for (size_t cnt2=0; cnt2 < time_to_name.size(); ++cnt2) {
          if (cnt1 == cnt2) continue;
          std::vector<Feature> feature_row2 = feature_name_map[time_to_name[cnt2].second];
          for (size_t j=0; j < feature_row2.size(); ++j) {
            String name2 = time_to_name[cnt2].second + "_" + (String)feature_row2[j].getUniqueId();
            if (variables.find(name2) == variables.end()) {
                _addVariable(problem, name2);
                variables.insert(name2);
            }
            String var_qp_name = time_to_name[cnt1].second + "_" + (String)i + "-" + time_to_name[cnt2].second + "_" + (String)j;
            Int index_var_qp = _addVariable(problem, var_qp_name);
            Int index1 = problem.getColumnIndex(name1);
            Int index2 = problem.getColumnIndex(name2);
            Int indices1[] = {index1, index_var_qp};
            double values[] = {1., -1.};
            _addConstraint(problem, 2, indices1, values, var_qp_name + "-QP1", 0., 1., LPWrapper::LOWER_BOUND_ONLY);
            Int indices2[] = {index2, index_var_qp};
            _addConstraint(problem, 2, indices2, values, var_qp_name + "-QP2", 0., 1., LPWrapper::LOWER_BOUND_ONLY);
            Int indices3[] = {index1, index2, index_var_qp};
            double values3[] = {1., 1., -1};
            _addConstraint(problem, 3, indices3, values3, var_qp_name + "-QP3", 0., 1., LPWrapper::UPPER_BOUND_ONLY);
          }
        }
      }
    }
    std::cout << "=======END OPTIMIZE TR========" << std::endl;
  }

  void MRMFeatureSelector::optimize_score() {}

  void MRMFeatureSelector::select_MRMFeature_qmip(
    FeatureMap& features,
    TargetedExperiment& targeted_exp
  )
  {
    std::cout << "=======START========" << std::endl;
    std::unordered_set<std::string> names;
    std::vector<std::pair<double, String>> time_to_name;
    // std::map< String, std::vector<std::map<String, DataValue>> > feature_name_map;
    std::map< String, std::vector<Feature> > feature_name_map;
    size_t feature_count = 0;
    for (FeatureMap::iterator it = features.begin(); it != features.end(); ++it) {
      String component_group_name = it->getMetaValue("PeptideRef").toString();
      double assay_retention_time = it->getMetaValue("assay_rt");
      if (names.find(component_group_name) != names.end()) {
        continue;
      }
      names.insert(component_group_name);
      time_to_name.push_back(std::make_pair(assay_retention_time, component_group_name));
      if (feature_name_map.find(component_group_name) == feature_name_map.end()) {
        feature_name_map[component_group_name] = std::vector<Feature>();
      }
      feature_name_map[component_group_name].push_back(*it);
      ++feature_count;
      if (!getSelectTransitionGroup()) {
        for (std::vector<Feature>::const_iterator sub_it = it->getSubordinates().begin();
            sub_it != it->getSubordinates().end(); ++sub_it) {
          String component_name = sub_it->getMetaValue("native_id").toString();
          if (names.find(component_name) != names.end()) {
            continue;
          }
          names.insert(component_name);
          time_to_name.push_back(std::make_pair(assay_retention_time, component_name));
          if (feature_name_map.find(component_name) == feature_name_map.end()) {
            feature_name_map[component_name] = std::vector<Feature>();
          }
          feature_name_map[component_name].push_back(*sub_it);
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
      optimize_Tr(time_slice, feature_name_map);
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
