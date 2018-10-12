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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureScheduler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  MRMFeatureScheduler::MRMFeatureScheduler() :
      DefaultParamHandler("MRMFeatureScheduler")
  {}

  MRMFeatureScheduler::~MRMFeatureScheduler() {}

  void MRMFeatureScheduler::schedule_MRMFeatures(MRMFeatureSelector& feature_selector, const FeatureMap& features, FeatureMap& output_features) {
    FeatureMap features_mutable = features;
    for (size_t i = 0; i < segment_window_lengths_.size(); ++i) {
      Param param;
      param.setValue("nn_threshold", nn_thresholds_[i]);
      param.setValue("locality_weight", locality_weights_[i]);
      param.setValue("select_transition_group", select_transition_groups_[i]);
      param.setValue("segment_window_length", segment_window_lengths_[i]);
      param.setValue("segment_step_length", segment_step_lengths_[i]);
      param.setValue("select_highest_count", select_highest_counts_[i]);
      param.setValue("variable_type", variable_types_[i]);
      param.setValue("optimal_threshold", optimal_thresholds_[i]);
      feature_selector.setParameters(param);
      feature_selector.select_MRMFeature(features_mutable, output_features);
      features_mutable = output_features;
    }
  }

  void MRMFeatureScheduler::schedule_MRMFeaturesQMIP(const FeatureMap& features, FeatureMap& output_features) {
    MRMFeatureSelectorQMIP feature_selector;
    schedule_MRMFeatures(feature_selector, features, output_features);
  }

  void MRMFeatureScheduler::schedule_MRMFeatures_score(const FeatureMap& features, FeatureMap& output_features) {
    MRMFeatureSelectorScore feature_selector;
    schedule_MRMFeatures(feature_selector, features, output_features);
  }

  void MRMFeatureScheduler::setNNThresholds(const std::vector<double>& nn_thresholds)
  {
    nn_thresholds_ = nn_thresholds;
  }

  void MRMFeatureScheduler::setLocalityWeights(const std::vector<String>& locality_weights)
  {
    locality_weights_ = locality_weights;
  }

  void MRMFeatureScheduler::setSelectTransitionGroups(const std::vector<String>& select_transition_groups)
  {
    select_transition_groups_ = select_transition_groups;
  }

  void MRMFeatureScheduler::setSegmentWindowLengths(const std::vector<double>& segment_window_lengths)
  {
    segment_window_lengths_ = segment_window_lengths;
  }

  void MRMFeatureScheduler::setSegmentStepLengths(const std::vector<double>& segment_step_lengths)
  {
    segment_step_lengths_ = segment_step_lengths;
  }

  void MRMFeatureScheduler::setSelectHighestCounts(const std::vector<String>& select_highest_counts)
  {
    select_highest_counts_ = select_highest_counts;
  }

  void MRMFeatureScheduler::setVariableTypes(const std::vector<String>& variable_types)
  {
    variable_types_ = variable_types;
  }

  void MRMFeatureScheduler::setOptimalThresholds(const std::vector<double>& optimal_thresholds)
  {
    optimal_thresholds_ = optimal_thresholds;
  }
}
