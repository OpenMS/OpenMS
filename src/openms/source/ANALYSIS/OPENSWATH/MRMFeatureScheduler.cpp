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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureScheduler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  void MRMFeatureScheduler::schedule_MRMFeatures(
    MRMFeatureSelector& feature_selector,
    const FeatureMap& features,
    FeatureMap& output_features
  ) const
  {
    FeatureMap features_mutable = features;
    for (const SelectorParameters& params : parameters_) {
      feature_selector.setNNThreshold(params.nn_threshold);
      feature_selector.setLocalityWeight(params.locality_weight);
      feature_selector.setSelectTransitionGroup(params.select_transition_group);
      feature_selector.setSegmentWindowLength(params.segment_window_length);
      feature_selector.setSegmentStepLength(params.segment_step_length);
      feature_selector.setVariableType(params.variable_type);
      feature_selector.setOptimalThreshold(params.optimal_threshold);
      feature_selector.setScoreWeights(params.score_weights);

      feature_selector.select_MRMFeature(features_mutable, output_features);
      features_mutable = output_features;
    }
  }

  void MRMFeatureScheduler::schedule_MRMFeaturesQMIP(const FeatureMap& features, FeatureMap& output_features) const
  {
    MRMFeatureSelectorQMIP feature_selector;
    schedule_MRMFeatures(feature_selector, features, output_features);
  }

  void MRMFeatureScheduler::schedule_MRMFeaturesScore(const FeatureMap& features, FeatureMap& output_features) const
  {
    MRMFeatureSelectorScore feature_selector;
    schedule_MRMFeatures(feature_selector, features, output_features);
  }

  void MRMFeatureScheduler::setSchedulerParameters(const std::vector<SelectorParameters>& parameters)
  {
    parameters_ = parameters;
  }

  std::vector<MRMFeatureScheduler::SelectorParameters>& MRMFeatureScheduler::getSchedulerParameters()
  {
    return parameters_;
  }
}
