// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <string>
#include <vector>

namespace OpenMS
{

  /**
    @brief Normalizes the maps of a consensus map by a per-map offset that may vary with retention time.

    Runs of one experiment differ in their overall loading, and their ion signal can also drift along the
    gradient (spray stability, source contamination). A single factor per map removes the first, not the
    second. This algorithm removes both: every map gets a log2 offset @em o(t) that is estimated locally in
    retention time from features quantified across the maps, and each intensity is divided by @em 2^o(t).
    The recipe follows the run normalisation of andes (bigbio/andes, `quant::normalize`).

    - <b>Reference features</b> are the consensus features with an intensity in at least @p min_maps maps
      that pass the accession and description filters. Their intensities in all maps are reference cells.
    - <b>Backfitting.</b> All offsets start at 0. Each of @p rounds rounds
      1. sets the level of every reference feature to the mean over its maps of `log2(intensity) - o`,
      2. computes the residual `z = log2(intensity) - level` of every reference cell, and
      3. sets the offset of every reference cell to the median residual of the @p neighbors reference cells
         of the same map before it and the @p neighbors after it in retention time; the cell's own residual
         is left out.
    - <b>Other features</b> (fewer maps, filtered out) take the median residual of the @p neighbors
      reference cells of the map on each side of their retention time.
    - A map with fewer than @p min_reference reference cells is not normalized. When fewer than two maps
      qualify, nothing is normalized.

    The level of a feature is its mean over all of its maps, its own map included, so a map that drifts
    does not pull the offsets of the other maps. Because the level is the mean over the maps, the offsets
    of one feature average to roughly zero and the overall intensity scale is kept.

    The retention time of a cell is the retention time of its feature handle. For maps that were aligned
    before linking, this is the aligned retention time; any retention time scale that is monotonic within a
    map gives the same result. The intensity of the consensus feature itself is left unchanged, as in
    ConsensusMapNormalizerAlgorithmMedian.
  */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmRTLocal
  {
  public:
    /// Settings of the normalization
    struct Settings
    {
      Size neighbors = 50;     ///< reference cells on each side (in retention time) that estimate a cell's offset
      Size rounds = 5;         ///< backfitting rounds
      Size min_maps = 3;       ///< a consensus feature is a reference feature with intensities in at least this many maps
      Size min_reference = 20; ///< maps with fewer reference cells are not normalized
    };

    /**
      @brief Computes the log2 offset of every feature handle.

      @param[in] map The consensus map
      @param[in] settings Settings of the normalization
      @param[in] acc_filter Regular expression the protein accessions of a reference feature must (partially) match; empty: no filter
      @param[in] desc_filter Regular expression the protein descriptions of a reference feature must (partially) match; empty: no filter
      @return One vector per consensus feature (in map order) with one offset per feature handle (in the
              order of ConsensusFeature::getFeatures()). Handles without a positive intensity get 0.
    */
    static std::vector<std::vector<double>> computeOffsets(const ConsensusMap& map, const Settings& settings,
                                                           const std::string& acc_filter, const std::string& desc_filter);

    /**
      @brief Normalizes the maps of the consensus map: divides every feature handle's intensity by 2^offset.

      @param[in,out] map The consensus map
      @param[in] settings Settings of the normalization
      @param[in] acc_filter Regular expression the protein accessions of a reference feature must (partially) match; empty: no filter
      @param[in] desc_filter Regular expression the protein descriptions of a reference feature must (partially) match; empty: no filter
    */
    static void normalizeMaps(ConsensusMap& map, const Settings& settings, const std::string& acc_filter, const std::string& desc_filter);
  };

} // namespace OpenMS
