// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
  /**
    @brief A feature grouping algorithm using Wasserstein optimal transport.

    Groups corresponding features across multiple maps by solving a minimum-cost
    network flow problem. Each pair of maps is aligned independently using
    WNetAligner<2> (2D: m/z and RT), and pairwise matchings are merged into
    consensus groups via union-find.

    Features are matched by minimizing the total transport cost between their
    (m/z, RT) positions, with unmatched features assigned to a "trash" node
    at a configurable cost. This produces a 1-to-1 consensus matching per pair.

    @htmlinclude OpenMS_FeatureGroupingAlgorithmWNet.parameters

    @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmWNet :
    public FeatureGroupingAlgorithm
  {
  public:
    /// Default constructor
    FeatureGroupingAlgorithmWNet();

    /// Destructor
    ~FeatureGroupingAlgorithmWNet() override;

    /**
        @brief Applies the algorithm to feature maps

        @pre The data ranges of the input maps have to be up-to-date (use FeatureMap::updateRanges).

        @param[in]  maps  Input feature maps (at least two required).
        @param[out] out   Resulting consensus map with grouped features.

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<FeatureMap>& maps, ConsensusMap& out) override;

    /**
        @brief Applies the algorithm to consensus maps

        @pre The data ranges of the input maps have to be up-to-date (use ConsensusMap::updateRanges).

        @param[in]  maps  Input consensus maps (at least two required).
        @param[out] out   Resulting consensus map with grouped features.

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out) override;

  private:
    FeatureGroupingAlgorithmWNet(const FeatureGroupingAlgorithmWNet&) = delete;
    FeatureGroupingAlgorithmWNet& operator=(const FeatureGroupingAlgorithmWNet&) = delete;
    FeatureGroupingAlgorithmWNet(FeatureGroupingAlgorithmWNet&&) = delete;
    FeatureGroupingAlgorithmWNet& operator=(FeatureGroupingAlgorithmWNet&&) = delete;

    template <typename MapType>
    void group_(const std::vector<MapType>& maps, ConsensusMap& out);
  };

} // namespace OpenMS
