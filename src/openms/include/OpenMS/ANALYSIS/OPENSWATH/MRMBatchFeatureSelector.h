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
    Class used to batch multiple calls to `MRMFeatureSelector`'s methods.
    The class offers a generic batch method (where the user is supposed to pass
    a `MRMFeatureSelector` derived object) and two specialized versions (Score and QMIP).
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
