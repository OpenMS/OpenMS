// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Batch processing wrapper for MRMFeatureSelector.

    This class enables iterative refinement of feature selection by applying
    MRMFeatureSelector multiple times with different parameter sets. Each iteration
    uses the output of the previous iteration as input, allowing for progressive
    filtering and optimization.

    @section MRMBatchFeatureSelector_usage Usage

    @code
    FeatureMap features, selected;
    std::vector<MRMFeatureSelector::SelectorParameters> params_list;

    // First pass: relaxed parameters
    params_list.push_back(MRMFeatureSelector::SelectorParameters());

    // Second pass: stricter parameters
    MRMFeatureSelector::SelectorParameters strict;
    strict.optimal_threshold = 0.8;
    params_list.push_back(strict);

    MRMBatchFeatureSelector::batchMRMFeaturesScore(features, selected, params_list);
    @endcode

    @see MRMFeatureSelector for single-pass feature selection
    @see MRMFeatureSelectorQMIP for QMIP-based selection
    @see MRMFeatureSelectorScore for score-based selection

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI MRMBatchFeatureSelector
  {
public:
    MRMBatchFeatureSelector() = delete;
    ~MRMBatchFeatureSelector() = delete;

    /**
      Calls `feature_selector.selectMRMFeature()` feeding it the parameters found in `parameters`.
      It calls said method `parameters.size()` times, using the result of each cycle as input
      for the next cycle.

      @param[in] feature_selector Base class for the feature selector to use
      @param[in] features Input features
      @param[out] selected_features Selected features
      @param[in] parameters Vector of parameters for the multiple calls to the selector
    */
    static void batchMRMFeatures(
      const MRMFeatureSelector& feature_selector,
      const FeatureMap& features,
      FeatureMap& selected_features,
      const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
    );

    /// Calls `batchMRMFeatures()` using a `MRMFeatureSelectorScore` selector
    static void batchMRMFeaturesScore(
      const FeatureMap& features,
      FeatureMap& selected_features,
      const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
    );

    /// Calls `batchMRMFeatures()` using a `MRMFeatureSelectorQMIP` selector
    static void batchMRMFeaturesQMIP(
      const FeatureMap& features,
      FeatureMap& selected_features,
      const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
    );
  };
}
