// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  class OPENMS_DLLAPI FeatureMapping
  {
    public:

    /// Stores information required for preprocessing
    class FeatureMappingInfo
    {
    public:
      std::vector<FeatureMap> feature_maps; ///< feature data
      KDTreeFeatureMaps kd_tree; ///< KDTree references into feature_maps to provides fast spatial queries
    };

    /// Stores preprocessed feature mapping information
    class FeatureToMs2Indices
    {
    public:
       std::map<const BaseFeature*, std::vector<size_t>> assignedMS2;
       std::vector<size_t> unassignedMS2;
    };

    /**
      @brief Allocate ms2 spectra to feature within the minimal distance

      @return FeatureToMs2Indices

      @param spectra: Input of PeakMap/MSExperiment with spectra information
      @param fm_info: KDTree used for query and match spectra with features
      @param precursor_mz_tolerance: mz_tolerance used for query
      @param precursor_rt_tolerance: rt tolerance used for query
      @param ppm: mz tolerance window calculation in ppm or Da

    */
    static FeatureToMs2Indices assignMS2IndexToFeature(const MSExperiment& spectra,
                                                       const FeatureMappingInfo& fm_info,
                                                       const double& precursor_mz_tolerance,
                                                       const double& precursor_rt_tolerance,
                                                       bool ppm);

  };
} // namespace OpenMS
