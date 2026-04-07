// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/OpenMSConfig.h>

#include <map>
#include <utility>

namespace OpenMS
{
  class ConsensusMap;

  /**
    @brief Post-link filter that removes unreliable match-between-runs (MBR) transferred
           feature handles from a ConsensusMap based on pairwise RT residual statistics.

    For each ConsensusFeature, handles whose source feature carries a peptide
    identification act as "anchors". Pairs of handles in the same ConsensusFeature
    that share the same modified peptide sequence and charge are used to populate
    a per-map-pair distribution of RT residuals (median + robust sigma derived from MAD).

    In a second pass, every handle that does NOT carry a peptide identification
    (typical MBR transfer) is checked against all ID-bearing handles in its
    ConsensusFeature: if its RT is statistically inconsistent with the per-pair
    anchor distribution for the corresponding map pair, the handle is removed
    and the consensus position is recomputed.

    A pooled global distribution is used as fallback for map pairs with too few
    anchors.

    The filter operates strictly post-link and does not require access to the
    source FeatureMaps: peptide identifications carry a "map_index" MetaValue
    that ConsensusFeature::insert(map_index, BaseFeature) attaches at link time.

    @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI ConsensusMapMBRFilter :
    public DefaultParamHandler
  {
public:
    /// Unordered pair of map indices, smaller index first
    using PairKey = std::pair<UInt64, UInt64>;

    /// Per-pair RT residual statistics
    struct PairStats
    {
      Size   n         = 0;
      double median    = 0.0;        ///< signed: rt(map_a) - rt(map_b), with a < b
      double mad_sigma = 0.0;        ///< 1.4826 * MAD, floored to min_pair_sigma
      Size   removals  = 0;          ///< number of unidentified handles removed in this pair
    };

    /// Per-map-pair report
    struct Report
    {
      std::map<PairKey, PairStats> pair_stats;
      PairStats                    pooled;
      Size                         handles_removed = 0;
    };

    /// Default constructor; registers parameters via defaults_
    ConsensusMapMBRFilter();

    /// Destructor
    ~ConsensusMapMBRFilter() override = default;

    /**
      @brief Run the filter on a ConsensusMap.

      Mutates @p cmap by removing unreliable transferred FeatureHandles and
      recomputing the consensus position of affected ConsensusFeatures.

      @param[in,out] cmap ConsensusMap to be filtered
      @return Report containing per-pair stats and removal counts
    */
    Report filter(ConsensusMap& cmap);
  };

} // namespace OpenMS
