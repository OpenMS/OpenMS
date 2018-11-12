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

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <map>
#include <vector>

namespace OpenMS
{
  /**
    Class used to schedule multiple calls to `MRMFeatureSelector`
  */
  class OPENMS_DLLAPI MRMFeatureScheduler
  {
public:
    MRMFeatureScheduler() = default;
    ~MRMFeatureScheduler() = default;

    /**
      Structure to easily feed the parameters to the `MRMFeatureSelector` derived classes
    */
    struct SelectorParameters
    {
      SelectorParameters() = default;

      SelectorParameters(
        Int nn,
        bool lw,
        bool stg,
        Int swl,
        Int ssl,
        MRMFeatureSelector::VariableType vt,
        double ot,
        std::map<String, MRMFeatureSelector::LambdaScore>& sw
      ) :
        nn_threshold(nn),
        locality_weight(lw),
        select_transition_group(stg),
        segment_window_length(swl),
        segment_step_length(ssl),
        variable_type(vt),
        optimal_threshold(ot),
        score_weights(sw) {}

      Int    nn_threshold            = 4; ///< Nearest neighbor threshold: the number of components or component groups to the left and right to include in the optimization problem (i.e. number of nearest compounds by Tr to include in network)
      bool   locality_weight         = false; ///< Weight compounds with a nearer Tr greater than compounds with a further Tr
      bool   select_transition_group = true; ///< Use components groups instead of components for retention time optimization
      Int    segment_window_length   = 8; ///< Number of components or component groups to include in the network
      Int    segment_step_length     = 4; ///< Number of of components or component groups to shift the `segment_window_length` at each loop
      MRMFeatureSelector::VariableType variable_type = MRMFeatureSelector::VariableType::CONTINUOUS; ///< INTEGER or CONTINUOUS
      double optimal_threshold       = 0.5; ///< Value above which the transition group or transition is considered optimal (0 < x < 1)
      std::map<String, MRMFeatureSelector::LambdaScore> score_weights; ///< Weights for the scores
    };

    /**
      Calls `feature_selector.select_MRMFeature()` feeding it the parameters found in `parameters_`.
      It calls said method `parameters_.size()` times, using the result of each cycle as input
      for the next cycle.

      @param[in] feature_selector Base class for the feature selector to use
      @param[in] features Input features
      @param[out] selected_features Selected features
    */
    void scheduleMRMFeatures(MRMFeatureSelector& feature_selector, const FeatureMap& features, FeatureMap& selected_features) const;

    /// Calls `scheduleMRMFeatures()` using a `MRMFeatureSelectorScore` selector
    void scheduleMRMFeaturesScore(const FeatureMap& features, FeatureMap& selected_features) const;

    /// Calls `scheduleMRMFeatures()` using a `MRMFeatureSelectorQMIP` selector
    void scheduleMRMFeaturesQMIP(const FeatureMap& features, FeatureMap& selected_features) const;

    /// Setter for the scheduler's parameters
    void setSchedulerParameters(const std::vector<SelectorParameters>& parameters);

    /// Getter for the scheduler's parameters
    std::vector<SelectorParameters>& getSchedulerParameters(void);

private:
    /// Parameters for a single call to the scheduler. All elements will be consumed.
    std::vector<SelectorParameters> parameters_;
  };
}
