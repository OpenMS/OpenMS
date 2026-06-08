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
    @brief Intensity normalisation of a @ref ConsensusMap by averaging
           per-feature intensity ratios against a reference map.

    A stateless utility class -- @ref computeCorrelation derives one
    scalar normalisation factor per map, and @ref normalizeMaps applies
    those factors to every feature handle in the consensus map.

    The reference map is the one with the most features (largest
    @c size in the column headers). For each other map, the ratio
    @c reference_intensity / @c map_intensity is computed per feature
    that is present in both maps, then averaged after dropping ratios
    that fall outside the configurable @c ratio_threshold window
    (see @ref computeCorrelation).

    Copy construction and assignment are deliberately suppressed.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmThreshold
  {
private:
    /// Copy construction is deliberately suppressed.
    ConsensusMapNormalizerAlgorithmThreshold(const ConsensusMapNormalizerAlgorithmThreshold & copyin);

    /// Assignment is deliberately suppressed.
    ConsensusMapNormalizerAlgorithmThreshold & operator=(const ConsensusMapNormalizerAlgorithmThreshold & rhs);

public:
    /**
      @brief Default constructor.
    */
    ConsensusMapNormalizerAlgorithmThreshold();

    /**
      @brief Destructor.
    */
    virtual ~ConsensusMapNormalizerAlgorithmThreshold();

    /**
      @brief Compute one normalisation factor per map.

      The reference map is the column-header entry with the largest
      @c size in @p map. For every other map @c j and every consensus
      feature that has a non-zero intensity in both the reference and
      @c j, the ratio @c reference_intensity / @c intensity_j is
      collected; ratios outside the window
      @c (ratio_threshold, 1/ratio_threshold) are discarded as outliers.
      The mean of the surviving ratios is returned as the factor for
      map @c j.

      Features are first filtered by @p acc_filter (regex on protein
      accessions) and @p desc_filter (regex on protein descriptions);
      only features that pass both filters contribute to the ratios.

      @warning When no feature passes the filters for a given map (the
               surviving ratio list is empty), the function returns a
               vector of @c 1.0 of length @c number_of_maps and logs a
               warning -- i.e. the result is "no normalisation".

      @param[in] map             ConsensusMap providing the per-map feature handles.
      @param[in] ratio_threshold Lower bound of the kept-ratio window. Use e.g. @c 0.5 to keep ratios in @c (0.5, 2.0).
      @param[in] acc_filter      Regex applied to accession strings; empty means "match anything".
      @param[in] desc_filter     Regex applied to description strings; empty means "match anything".
      @return One normalisation factor per map (size equals @c map.getColumnHeaders().size()).
    */
    static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold, const String & acc_filter, const String & desc_filter);

    /**
      @brief Multiply every feature handle's intensity by its map's factor.

      Iterates every consensus feature in @p map and, for each
      sub-feature handle, multiplies its intensity by
      @p ratios[handle.getMapIndex()]. Progress is reported to the
      console via @ref ProgressLogger.

      @param[in,out] map    Consensus map whose feature-handle intensities are scaled in place.
      @param[in]     ratios Per-map factors (typically obtained from @ref computeCorrelation). Must have at least @c (max-map-index + 1) entries.
    */
    static void normalizeMaps(ConsensusMap & map, const std::vector<double> & ratios);
  };

} // namespace OpenMS

