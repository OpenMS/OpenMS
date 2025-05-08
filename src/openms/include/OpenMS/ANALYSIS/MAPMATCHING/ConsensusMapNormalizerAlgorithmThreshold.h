// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
   * @brief Algorithms of ConsensusMapNormalizer
   *
   */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmThreshold
  {
private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold(const ConsensusMapNormalizerAlgorithmThreshold & copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold & operator=(const ConsensusMapNormalizerAlgorithmThreshold & rhs);

public:
    /// default constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold();

    /// destructor is not implemented -> private
    virtual ~ConsensusMapNormalizerAlgorithmThreshold();

    /**
     * @brief determines the ratio of all maps to the map with the most features
     * @param map ConsensusMap
     * @param ratio_threshold threshold for the ratio
     * @param acc_filter string describing the regular expression for filtering accessions
     * @param desc_filter string describing the regular expression for filtering descriptions
     */
    static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold, const String & acc_filter, const String & desc_filter);

    /**
     * @brief applies the given ratio to the maps of the consensusMap
     * @param map ConsensusMap
     * @param ratios ratios for the normalization
     */
    static void normalizeMaps(ConsensusMap & map, const std::vector<double> & ratios);
  };

} // namespace OpenMS

