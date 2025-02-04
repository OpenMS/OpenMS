// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  /** @class MRMFeatureSelector

    A Base class (it contains a pure virtual function named `optimize()`) for
    selection of MRM Features through Linear Programming.
  */
  class OPENMS_DLLAPI MRMFeatureSelector
  {
public:
    MRMFeatureSelector() = default;
    virtual ~MRMFeatureSelector() = default;

    enum class VariableType
    {
      INTEGER = 1,
      CONTINUOUS
    };

    enum class LambdaScore
    {
      LINEAR = 1,
      INVERSE,
      LOG,
      INVERSE_LOG,
      INVERSE_LOG10
    };

    /// To test private and protected methods
    friend class MRMFeatureSelector_test;

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
      Derived classes implement this pure virtual method.

      It sets up the linear programming problem and solves it.

      @param[in] time_to_name Pairs representing a mapping of retention times to transition names
      @param[in] feature_name_map Transitions' names to their features objects
      @param[out] result Transitions' names filtered out of the LP problem
      @param[in] parameters Parameters
    */
    virtual void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map<String, std::vector<Feature>>& feature_name_map,
      std::vector<String>& result,
      const SelectorParameters& parameters
    ) const = 0;

    /**
      The features are sorted by retention time and splitted into segments with
      the given step and window length. The features are then selected based on
      the results of `optimize()` method applied to each segment. The segments
      may overlap.

      @param[in] features Input features
      @param[out] selected_filtered Output features
      @param[in] parameters Parameters
    */
    void selectMRMFeature(
      const FeatureMap& features,
      FeatureMap& selected_filtered,
      const SelectorParameters& parameters
    ) const;

protected:
    /**
      Add variable to the LP problem instantiated in `optimize()`

      @param[in,out] problem LPWrapper object
      @param[in] name Column name
      @param[in] bounded Double bounded if true, otherwise Unbounded.
      @param[in] obj Objective value
      @param[in] variableType Either integer or continuous

      @return The variable's column index
    */
    Int addVariable_(
      LPWrapper& problem,
      const String& name,
      const bool bounded,
      const double obj,
      const VariableType variableType
    ) const;

    /**
      Scoring method used by the optimizer. Metavalues to use are decided by
      the `score_weights` argument.
      The returned value is used in the LP problems' variables and constraints.

      @param[in] feature Input feature
      @param[in] score_weights Score weights

      @return Computed score
    */
    double computeScore_(const Feature& feature, const std::map<String, LambdaScore>& score_weights) const;

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
      const std::vector<Int>& indices,
      const std::vector<double>& values,
      const String& name,
      const double lb,
      const double ub,
      const LPWrapper::Type param
    ) const;

private:
    /**
      Construct the target transition's or transition group's retention times that
      will be used to score candidate features based on their deviation from the
      relative distance between the target transition's or transition group's times

      @param[in] features Input features
      @param[out] time_to_name Pairs representing a mapping of retention times to transition names
      @param[out] feature_name_map Transitions' names to their features objects
      @param[in] select_transition_group Transition group selection
    */
    void constructTargTransList_(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map,
      const bool select_transition_group
    ) const;

    /**
      Transform the given score through the chosen lambda function

      Possible values for `lambda_score` are:
      - LambdaScore::LINEAR
      - LambdaScore::INVERSE
      - LambdaScore::LOG
      - LambdaScore::INVERSE_LOG
      - LambdaScore::INVERSE_LOG10

      @throw Exception::IllegalArgument When an invalid `lambda_score` is passed

      @param[in] score Value to transform
      @param[in] lambda_score A string representing the desired transformation

      @return The weighted value
    */
    double weightScore_(const double score, const LambdaScore lambda_score) const;

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
    /**
      Set up the linear programming problem and solve it.

      @param[in] time_to_name Pairs representing a mapping of retention times to transition names
      @param[in] feature_name_map Transitions' names to their features objects
      @param[out] result Transitions' names filtered out of the LP problem
      @param[in] parameters Parameters
    */
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map<String, std::vector<Feature>>& feature_name_map,
      std::vector<String>& result,
      const SelectorParameters& parameters
    ) const override;
  };

  /**
    Class used to select MRMFeatures based on a linear programming where each
    possible transition is weighted by a user defined score (most often retention
    time and peak intensity). The method is described in [TODO: update when published].
  */
  class OPENMS_DLLAPI MRMFeatureSelectorScore : public MRMFeatureSelector
  {
public:
    /**
      Set up the linear programming problem and solve it.

      @param[in] time_to_name Pairs representing a mapping of retention times to transition names
      @param[in] feature_name_map Transitions' names to their features objects
      @param[out] result Transitions' names filtered out of the LP problem
      @param[in] parameters Parameters
    */
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map<String, std::vector<Feature>>& feature_name_map,
      std::vector<String>& result,
      const SelectorParameters& parameters
    ) const override;
  };

  class MRMFeatureSelector_test : public MRMFeatureSelectorQMIP
  {
public:
    MRMFeatureSelector_test() = default;
    ~MRMFeatureSelector_test() override = default;

    void constructTargTransList_(
      const FeatureMap& features,
      std::vector<std::pair<double, String>>& time_to_name,
      std::map<String, std::vector<Feature>>& feature_name_map,
      const bool select_transition_group
    ) const
    {
      selector_.constructTargTransList_(features, time_to_name, feature_name_map, select_transition_group);
    }

    double weightScore_(const double score, const LambdaScore lambda_score) const
    {
      return selector_.weightScore_(score, lambda_score);
    }

    double computeScore_(const Feature& feature, const std::map<String, LambdaScore>& score_weights) const
    {
      return selector_.computeScore_(feature, score_weights);
    }

    String removeSpaces_(String str) const
    {
      return selector_.removeSpaces_(str);
    }

    MRMFeatureSelectorQMIP selector_;
  };
}
