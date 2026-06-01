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

  /**
    @brief A base class for selection of MRM Features through Linear Programming optimization.

    This class provides the framework for optimal feature selection in MRM/SRM experiments
    using Linear Programming (LP). The key idea is to select the best peak for each transition
    while maintaining consistency with neighboring transitions based on retention time relationships.

    Two derived implementations are provided:
    - @ref MRMFeatureSelectorQMIP - Uses Quadratic Mixed Integer Programming based on relative RT
    - @ref MRMFeatureSelectorScore - Uses score-weighted linear programming

    @section MRMFeatureSelector_algorithm Algorithm Overview

    1. Features are sorted by retention time
    2. The RT range is divided into overlapping segments (sliding window)
    3. For each segment, an LP problem is formulated and solved
    4. Solutions from segments are merged to produce final selection

    @section MRMFeatureSelector_params Key Parameters

    - `nn_threshold`: Number of nearest neighbors to include in optimization
    - `segment_window_length`: Size of sliding window
    - `segment_step_length`: Step size between windows
    - `variable_type`: INTEGER (exact) or CONTINUOUS (relaxed) LP
    - `optimal_threshold`: Cutoff for considering a feature as selected (0-1)
    - `score_weights`: Weights for different scoring functions (LINEAR, LOG, INVERSE, etc.)

    @see MRMBatchFeatureSelector for iterative batch processing
    @see LPWrapper for the underlying LP solver interface

    @ingroup TargetedQuantitation
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
      Int    segment_step_length     = 4; ///< Number of components or component groups to shift the `segment_window_length` at each loop
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
    @brief @ref MRMFeatureSelector implementation that selects MRM features via a quadratic-programming formulation over relative retention time.

    For each transition, considers neighbouring transitions inside
    @ref MRMFeatureSelector::SelectorParameters::nn_threshold and adds locality-weighted
    pairwise terms keyed on the expected retention-time delta. Per-feature score is
    normalised by the @c nth root of the number of @c score_weights entries (with @c n
    the number of weights) so multi-weight configurations don't dominate the QP.

    Companion to @ref MRMFeatureSelectorScore (linear, no nearest-neighbour coupling).

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI MRMFeatureSelectorQMIP : public MRMFeatureSelector
  {
public:
    /**
      @brief Build the quadratic LP problem for @p time_to_name and write the names of selected features to @p result.

      Uses @c LPWrapper with @c LPWrapper::MIN objective sense. For each transition,
      registers one binary (or continuous, depending on
      @ref MRMFeatureSelector::SelectorParameters::variable_type) variable per feature,
      then for every neighbour within @c nn_threshold adds a locality-weighted pairwise
      term keyed on the expected RT delta. After solving, every column whose value is
      @c >= @c parameters.optimal_threshold contributes its name to @p result.
      @p result is cleared on entry.

      @param[in]  time_to_name      Pairs of (retention time, transition name); the order of this vector defines the neighbour sliding window.
      @param[in]  feature_name_map  Transition name → its candidate features.
      @param[out] result            Names of the selected features (one per column whose solver value clears @c optimal_threshold). Cleared on entry.
      @param[in]  parameters        Algorithm parameters (@c nn_threshold, @c locality_weight, @c variable_type, @c score_weights, @c optimal_threshold, ...).
    */
    void optimize(
      const std::vector<std::pair<double, String>>& time_to_name,
      const std::map<String, std::vector<Feature>>& feature_name_map,
      std::vector<String>& result,
      const SelectorParameters& parameters
    ) const override;
  };

  /**
    @brief @ref MRMFeatureSelector implementation that selects MRM features via a linear program with score-weighted per-feature variables.

    Simpler than @ref MRMFeatureSelectorQMIP — no nearest-neighbour coupling, no
    pairwise quadratic terms. Each transition contributes one binary (or continuous)
    variable per candidate feature plus an equality constraint forcing exactly one
    feature to be picked per transition (sum @c == @c 1). The objective is the sum of
    per-feature scores produced by @c computeScore_ from
    @ref MRMFeatureSelector::SelectorParameters::score_weights.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI MRMFeatureSelectorScore : public MRMFeatureSelector
  {
public:
    /**
      @brief Build the linear program for @p time_to_name and write the names of selected features to @p result.

      Uses @c LPWrapper with @c LPWrapper::MIN objective sense. For each transition,
      registers one variable per candidate feature (binary or continuous, per
      @c parameters.variable_type) scored by @c computeScore_, and adds a
      @c DOUBLE_BOUNDED constraint with both bounds equal to @c 1.0 (i.e. exactly one
      feature must be picked per transition). After solving, every column whose value
      is @c >= @c parameters.optimal_threshold contributes its name to @p result.
      @p result is cleared on entry.

      @param[in]  time_to_name      Pairs of (retention time, transition name). Order is irrelevant — unlike the QMIP variant, no neighbour window is used.
      @param[in]  feature_name_map  Transition name → its candidate features.
      @param[out] result            Names of the selected features (one per column whose solver value clears @c optimal_threshold). Cleared on entry.
      @param[in]  parameters        Algorithm parameters (@c variable_type, @c score_weights, @c optimal_threshold; @c nn_threshold and @c locality_weight are not consulted).
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
