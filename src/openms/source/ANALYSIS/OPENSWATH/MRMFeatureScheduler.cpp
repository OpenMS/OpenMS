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
  MRMFeatureScheduler::MRMFeatureScheduler() {}

  MRMFeatureScheduler::~MRMFeatureScheduler() {}

  FeatureMap MRMFeatureScheduler::schedule_MRMFeatures(MRMFeatureSelector& feature_selector, FeatureMap& features) {
    int nn_thresholds[] = {2, 4, 6, 4};
    bool locality_weights[] = {false, true, true, true};
    bool select_transition_groups[] = {true, true, true, true};
    int segment_window_lengths[] = {12, 24, 48, -1};
    int segment_step_lengths[] = {2, 6, 12, -1};
    bool select_highest_counts[] = {false, false, false, false};
    String variable_types[] = {"continous", "continous", "continous", "continous"};
    double optimal_thresholds[] = {0.5, 0.5, 0.5, 0.5};
    FeatureMap optimal_features = features;
    for (size_t i = 0; i < 2; ++i) {
      feature_selector.setNNThreshold(nn_thresholds[i]);
      feature_selector.setLocalityWeight(locality_weights[i]);
      feature_selector.setSelectTransitionGroup(select_transition_groups[i]);
      feature_selector.setSegmentWindowLength(segment_window_lengths[i]);
      feature_selector.setSegmentStepLength(segment_step_lengths[i]);
      feature_selector.setSelectHighestCount(select_highest_counts[i]);
      feature_selector.setVariableType(variable_types[i]);
      feature_selector.setOptimalThreshold(optimal_thresholds[i]);
      optimal_features = feature_selector.select_MRMFeature(optimal_features);
      std::cout << "optimal_features " << optimal_features.size() << std::endl;
    }
    return optimal_features;
  }

  FeatureMap MRMFeatureScheduler::schedule_MRMFeaturesQMIP(FeatureMap& features) {
    MRMFeatureSelectorQMIP feature_selector;
    return schedule_MRMFeatures(feature_selector, features);
  }

  FeatureMap MRMFeatureScheduler::schedule_MRMFeatures_score(FeatureMap& features) {
    MRMFeatureSelectorScore feature_selector;
    return schedule_MRMFeatures(feature_selector, features);
  }
}
