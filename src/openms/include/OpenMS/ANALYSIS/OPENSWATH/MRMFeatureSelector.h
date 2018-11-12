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
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  class OPENMS_DLLAPI MRMFeatureSelector
  {
public:
    MRMFeatureSelector() = default;
    virtual ~MRMFeatureSelector() = default;

    /// To test private and protected methods
    friend class MRMFeatureSelector_test;

    /**
      Derived classes implement this pure virtual method.

      It setups the linear programming problem and solves it.

      @param[in] time_to_name Pairs representing a mapping of retention times to transition names
      @param[in] feature_name_map Transitions' names to their features objects
      @param[out] result Transitions' names filtered out of the LP problem
    */
    virtual void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map<String, std::vector<Feature>>& feature_name_map,
      std::vector<String>& result
    ) = 0;

    /**
      Splits the features into time segments. `optimize()` method is run on each of these segments.

      @param[in] features Input features
      @param[out] selected_filtered Output features
    */
    void selectMRMFeature(const FeatureMap& features, FeatureMap& selected_filtered);

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

protected:
    /**
      Add variable to the LP problem instantiated in `optimize()`

      @param[in,out] problem LPWrapper object
      @param[in] name Column name
      @param[in] bounded Double bounded if true, otherwise Unbounded.
      @param[in] obj Objective value

      @return The variable's column index
    */
    Int addVariable_(
      LPWrapper& problem,
      const String& name,
      const bool bounded = true,
      const double obj = 1.0
    ) const;

    /**
      Scoring method used by the optimizer. Metavalues to use are decided by
      the `score_weights_` class' member.
      The returned value is used in the LP problems' variables and contraints.

      @param[in] feature Input feature

      @return Computed score
    */
    double computeScore_(const Feature& feature) const;

    /**
      Add constraint to the LP problem instantiated in `optimize()`

      @param[in,out] problem LPWrapper object
      @param[in] indices LP matrix indices
      @param[in] values LP matrix values
      @param[in] name Row name
      @param[in] lb Lower bound
      @param[in] ub Upper bound
      @param[in] param Row type
    */
    void addConstraint_(
      LPWrapper& problem,
      std::vector<Int> indices,
      std::vector<double> values,
      const String& name,
      const double lb,
      const double ub,
      const LPWrapper::Type param
    ) const;

private:
    Int    nn_threshold_            = 4;
    bool   locality_weight_         = false;
    bool   select_transition_group_ = true;
    Int    segment_window_length_   = 8;
    Int    segment_step_length_     = 4;
    String variable_type_           = "continuous";
    double optimal_threshold_       = 0.5;
    std::map<String, String> score_weights_;

    /**
      Construct the target transition's or transition group's retention times that
      will be used to score candidate features based on their deviation from the
      relative distance between the target transition's or transition group's times

      @param[in] features Input features
      @param[out] time_to_name Pairs representing a mapping of retention times to transition names
      @param[out] feature_name_map Transitions' names to their features objects
    */
    void constructTargTransList_(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map
    ) const;

    /**
      Transform the given score through the chosen lambda function

      Possible values for `lambda_score` are:
      - "lambda score: score*1.0"
      - "lambda score: 1/score"
      - "lambda score: log(score)"
      - "lambda score: 1/log(score)"
      - "lambda score: 1/log10(score)"

      @throw Exception::IllegalArgument When an invalid `lambda_score` is passed

      @param[in] score Value to transform
      @param[in] lambda_score A string representing the desired transformation

      @return The weighted value
    */
    double weightScore_(const double score, const String& lambda_score) const;

    /// Removes spaces from the given string, not-in-place.
    String removeSpaces_(String str) const;
  };

  /**
    Class used to select MRMFeatures based on relative retention time using a
    quadratic mixed integer programming (QMIP) formulation.
    The method is described in [TODO: update when published]
  */
  class OPENMS_DLLAPI MRMFeatureSelectorQMIP : public MRMFeatureSelector
  {
public:
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map< String, std::vector<Feature> >& feature_name_map,
      std::vector<String>& result
    );
  };

  /**
    Class used to select MRMFeatures based on a linear programming where each
    possible transition is weighted by a user defined score (most often retention
    time and peak intensity). The method is described in [TODO: update when published].
  */
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

    void constructTargTransList_(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map
    ) const
    {
      selector_.constructTargTransList_(features, time_to_name, feature_name_map);
    }

    double weightScore_(const double score, const String& lambda_score) const
    {
      return selector_.weightScore_(score, lambda_score);
    }

    double computeScore_(const Feature& feature) const
    {
      return selector_.computeScore_(feature);
    }

    String removeSpaces_(String str) const
    {
      return selector_.removeSpaces_(str);
    }

    void setScoreWeights(const std::map<String, String>& score_weights)
    {
      selector_.setScoreWeights(score_weights);
    }

    MRMFeatureSelectorQMIP selector_;
  };
}
