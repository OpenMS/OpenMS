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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>

namespace OpenMS
{
  class OPENMS_DLLAPI MRMFeatureSelector
  {
    friend class MRMFeatureSelectorQMIP;
    friend class MRMFeatureSelectorScore;
public:
    MRMFeatureSelector() = default;
    virtual ~MRMFeatureSelector() = default;

    /// To test private methods
    friend class MRMFeatureSelector_test;

    virtual void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map< String, std::vector<Feature> >& feature_name_map,
      std::vector<String>& result
    ) = 0;

    void select_MRMFeature(const FeatureMap& features, FeatureMap& features_filtered);

    void setNNThreshold(const Int nn_threshold);
    Int getNNThreshold() const;

    void setLocalityWeight(const bool locality_weight);
    bool getLocalityWeight() const;

    void setSelectTransitionGroup(const bool select_transition_group);
    bool getSelectTransitionGroup() const;

    void setSegmentWindowLength(const Int segment_window_length);
    Int getSegmentWindowLength() const;

    void setSegmentStepLength(const Int segment_step_length);
    Int getSegmentStepLength() const;

    void setVariableType(const String& variable_type);
    String getVariableType() const;

    void setOptimalThreshold(const double optimal_threshold);
    double getOptimalThreshold() const;

    void setScoreWeights(const std::map<String, String>& score_weights);
    std::map<String, String> getScoreWeights() const;

private:
    Int    nn_threshold_            = 4;
    bool   locality_weight_         = false;
    bool   select_transition_group_ = true;
    Int    segment_window_length_   = 8;
    Int    segment_step_length_     = 4;
    String variable_type_           = "continuous";
    double optimal_threshold_       = 0.5;
    std::map<String, String> score_weights_;

    Int addVariable(LPWrapper& problem, const String& name, const bool bounded = true, const double obj = 1.0) const;

    void addConstraint(
      LPWrapper& problem,
      std::vector<Int> indices,
      std::vector<double> values,
      const String& name,
      const double lb,
      const double ub,
      const LPWrapper::Type param
    ) const;

    void constructToList(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map
    ) const;

    double weight_func(const double score, const String& lambda_score) const;

    double make_score(const Feature& feature) const;

    String remove_spaces(String str) const;
  };

  class OPENMS_DLLAPI MRMFeatureSelectorQMIP : public MRMFeatureSelector
  {
public:
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map< String, std::vector<Feature> >& feature_name_map,
      std::vector<String>& result
    );
  };

  class OPENMS_DLLAPI MRMFeatureSelectorScore : public MRMFeatureSelector
  {
public:
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map< String, std::vector<Feature> >& feature_name_map,
      std::vector<String>& result
    );
  };

  class MRMFeatureSelector_test : public MRMFeatureSelectorQMIP
  {
public:
    MRMFeatureSelector_test() = default;
    ~MRMFeatureSelector_test() = default;

    void constructToList(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map
    ) const
    {
      selector_.constructToList(features, time_to_name, feature_name_map);
    }

    double weight_func(const double score, const String& lambda_score) const
    {
      return selector_.weight_func(score, lambda_score);
    }

    double make_score(const Feature& feature) const
    {
      return selector_.make_score(feature);
    }

    String remove_spaces(String str) const
    {
      return selector_.remove_spaces(str);
    }

    void setScoreWeights(const std::map<String, String>& score_weights)
    {
      selector_.setScoreWeights(score_weights);
    }

    MRMFeatureSelectorQMIP selector_;
  };
}
