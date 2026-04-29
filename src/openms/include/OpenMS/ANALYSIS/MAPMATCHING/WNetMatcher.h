// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Pairwise point-set matching using Wasserstein optimal transport.

    Matches two sets of 2D points (with associated intensities) by solving
    a minimum-cost network flow problem. Returns 1-to-1 matched index pairs
    selected greedily by descending transport flow.

    This provides a minimal, FeatureMap-independent interface to the WNetAlign
    algorithm. For feature-level grouping across multiple maps, use
    FeatureGroupingAlgorithmWNet instead.

    @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI WNetMatcher
  {
  public:
    /// Distance metric for comparing point positions
    enum class DistanceMetric { L1, L2, LINF };

    /// Result of a pairwise matching
    struct OPENMS_DLLAPI MatchResult
    {
      std::vector<std::pair<Size, Size>> matched_pairs; ///< (index_a, index_b) of matched points
      double cost = 0.0; ///< Total transport cost
    };

    /**
      @brief Match two sets of 2D points using optimal transport.

      Points are matched by minimizing the total transport cost between their
      positions, weighted by intensity. Unmatched points are sent to a "trash"
      node at the specified cost.

      @param[in] positions_a    Positions of first point set (each {dim0, dim1})
      @param[in] intensities_a  Intensities/weights of first point set
      @param[in] positions_b    Positions of second point set (each {dim0, dim1})
      @param[in] intensities_b  Intensities/weights of second point set
      @param[in] metric         Distance metric (L1, L2, or LINF)
      @param[in] max_distance   Maximum distance to consider a match
      @param[in] trash_cost     Cost of leaving a point unmatched

      @return MatchResult with paired indices and total cost
    */
    static MatchResult match(
      const std::vector<std::array<double, 2>>& positions_a,
      const std::vector<double>& intensities_a,
      const std::vector<std::array<double, 2>>& positions_b,
      const std::vector<double>& intensities_b,
      DistanceMetric metric = DistanceMetric::LINF,
      double max_distance = 100.0,
      double trash_cost = 100.0
    );

    /**
      @brief Convert a string to a DistanceMetric enum value.

      Recognizes "L1", "L2", and "LINF" (maps to DistanceMetric::L1/L2/LINF
      respectively). Any other input emits a warning via OPENMS_LOG_WARN and
      returns DistanceMetric::LINF as the fallback.
    */
    static DistanceMetric metricFromString(const std::string& s);
  };

} // namespace OpenMS
