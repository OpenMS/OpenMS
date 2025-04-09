// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  void MRMBatchFeatureSelector::batchMRMFeatures(
    const MRMFeatureSelector& feature_selector,
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    FeatureMap input_features = features;
    selected_features.clear();
    for (const MRMFeatureSelector::SelectorParameters& params : parameters)
    {
      feature_selector.selectMRMFeature(input_features, selected_features, params);
      input_features = selected_features;
    }
  }

  void MRMBatchFeatureSelector::batchMRMFeaturesQMIP(
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    MRMFeatureSelectorQMIP feature_selector;
    batchMRMFeatures(feature_selector, features, selected_features, parameters);
  }

  void MRMBatchFeatureSelector::batchMRMFeaturesScore(
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    MRMFeatureSelectorScore feature_selector;
    batchMRMFeatures(feature_selector, features, selected_features, parameters);
  }
}
