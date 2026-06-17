// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
    @brief Quantile-normalisation of feature intensities across the maps of a @ref ConsensusMap.

    Stateless utility class. @ref normalizeMaps rescales every map's
    sorted intensity distribution to a shared reference distribution
    derived from the average of the per-map quantiles, so that all maps
    end up sharing the same intensity distribution while their internal
    intensity ranks are preserved.

    Copy construction and assignment are deliberately suppressed.

    @ingroup Quantitation
  */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmQuantile
  {
private:
    /// Copy construction is deliberately suppressed.
    ConsensusMapNormalizerAlgorithmQuantile(const ConsensusMapNormalizerAlgorithmQuantile & copyin);

    /// Assignment is deliberately suppressed.
    ConsensusMapNormalizerAlgorithmQuantile & operator=(const ConsensusMapNormalizerAlgorithmQuantile & rhs);

public:
    /// Default constructor.
    ConsensusMapNormalizerAlgorithmQuantile();

    /// Destructor.
    virtual ~ConsensusMapNormalizerAlgorithmQuantile();

    /**
      @brief Quantile-normalise the feature intensities of @p map in place.

      Builds a reference intensity distribution by sorting each map's
      feature intensities, resampling each sorted distribution to a
      common length (the size of the largest map), and averaging the
      resampled distributions per quantile. The intensities are then
      written back to the consensus features rank-preservingly: the
      feature with the @c k -th smallest intensity in map @c j gets the
      @c k -th smallest value of the resampled reference distribution
      for map @c j (after a second resampling step that shrinks the
      reference back to map @c j 's original length).

      @param[in,out] map Consensus map whose feature-handle intensities
                         are rewritten in place. Per-map column-header
                         indices @c 0 ... @c N-1 must each be present
                         (@ref extractIntensityVectors throws
                         otherwise).

      @throws Exception::ElementNotFound when a column-header index in
                                         @c 0 ... @c N-1 is missing
                                         from @p map (propagated from
                                         @ref extractIntensityVectors).
    */
    static void normalizeMaps(ConsensusMap & map);

    /**
      @brief Linearly resample @p data_in to @p n_resampling_points evenly-spaced positions.

      The output is cleared, resized to @p n_resampling_points and
      filled by linear interpolation between consecutive @p data_in
      entries:

      - @c data_out.front() takes @c data_in.front() and
        @c data_out.back() takes @c data_in.back().
      - Intermediate positions are interpolated between the two
        neighbouring source entries.
      - When @p n_resampling_points is @c 0 the output is just cleared
        and the function returns.

      @param[in]  data_in             Source vector; assumed sorted by
                                      the caller (@ref normalizeMaps
                                      sorts before calling).
      @param[out] data_out            Destination vector; cleared and
                                      resized.
      @param[in]  n_resampling_points Number of points to produce.
    */
    static void resample(const std::vector<double> & data_in, std::vector<double> & data_out, UInt n_resampling_points);

    /**
      @brief Collect the feature-handle intensities of @p map into one vector per source map.

      @p out_intensities is cleared and resized to the number of
      column-header entries; the per-map row receives the intensities
      of every feature handle whose @c map_index points at it (in the
      iteration order of the consensus features).

      @param[in]  map             Source @ref ConsensusMap.
      @param[out] out_intensities Per-map intensity vectors; index
                                  @c j contains the intensities of all
                                  feature handles whose map-index is
                                  @c j.

      @throws Exception::ElementNotFound when the column-header lookup
                                         for one of the expected map
                                         indices fails.
    */
    static void extractIntensityVectors(const ConsensusMap & map, std::vector<std::vector<double> > & out_intensities);

    /**
      @brief Push the per-map intensity vectors back into the corresponding feature handles.

      Iterates the consensus features of @p map in the same order as
      @ref extractIntensityVectors and writes
      @c feature_ints[map_index][k] into the @c k -th feature handle of
      that map index. The caller is responsible for keeping
      @p feature_ints in lock-step with the iteration order used by
      @ref extractIntensityVectors -- @ref normalizeMaps does this by
      mutating the same per-map vector in place.

      @param[in]     feature_ints Per-map intensity vectors (same shape
                                  as the output of
                                  @ref extractIntensityVectors).
      @param[in,out] map          Consensus map whose feature-handle
                                  intensities are updated in place.
    */
    static void setNormalizedIntensityValues(const std::vector<std::vector<double> > & feature_ints, ConsensusMap & map);
  };

} // namespace OpenMS

