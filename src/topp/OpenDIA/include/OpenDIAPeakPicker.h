// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <vector>

namespace OpenMS::OpenDIA
{
  /** @brief Settings of the consensus peak-group picker. */
  struct PeakPickerOptions
  {
    std::size_t smoothing_passes = 1;   ///< [1 2 1]/4 passes over the consensus trace; each widens a peak by about one point.
    std::size_t apex_search = 2;        ///< The apex is the smoothed maximum within this many points of the hint.
    double boundary_fraction = 0.05;    ///< A boundary stops at this fraction of the smoothed apex if no valley comes first.
    std::size_t max_half_span = 64;     ///< Points walked at most on each side of the apex.
    std::size_t fragment_apex_tolerance = 1; ///< A fragment agrees if its raw maximum inside the peak is this close to the apex.
  };

  /** @brief Why a peak boundary stopped. */
  enum class PeakBoundary
  {
    VALLEY, ///< The smoothed trace rises again (next peak or noise).
    FLOOR,  ///< The smoothed trace reached boundary_fraction of the apex.
    EDGE,   ///< End of the trace or max_half_span reached; the true boundary may lie further out.
  };

  /**
    @brief One picked peak group. Positions are point indices; fractional ones are linearly interpolated.
  */
  struct PickedPeak
  {
    bool found = false;           ///< False for an empty or all-zero trace.
    std::size_t apex = 0;         ///< Smoothed consensus maximum near the hint.
    std::size_t left = 0;         ///< First point inside the peak.
    std::size_t right = 0;        ///< Last point inside the peak.
    double half_left = 0.0;       ///< Left half-maximum crossing of the smoothed consensus.
    double half_right = 0.0;      ///< Right half-maximum crossing of the smoothed consensus.
    double apex_height = 0.0;     ///< Smoothed consensus value at the apex.
    std::size_t fragments_with_signal = 0; ///< Fragments with any intensity inside [left, right].
    std::size_t agreeing_fragments = 0;    ///< Fragments whose raw maximum in [left, right] is within the apex tolerance.
    PeakBoundary left_reason = PeakBoundary::EDGE;
    PeakBoundary right_reason = PeakBoundary::EDGE;
  };

  /**
    @brief Pick the peak group around @p hint in a set of co-eluting fragment traces.

    Deliberately simple and allocation-light (one consensus buffer), so it can run on every candidate window:
    the fragment traces are summed into a consensus trace, smoothed with @p options.smoothing_passes [1 2 1]/4
    passes, the apex is the smoothed maximum within @p options.apex_search points of @p hint, and each boundary
    walks outwards while the smoothed trace keeps falling and stays above boundary_fraction x apex. Half-maximum
    crossings are interpolated for a FWHM. No library intensities, labels or decoy status are used.

    Element (f, s) of the traces is @p data[f * fragment_stride + s * point_stride].

    @param[in] data Fragment intensities.
    @param[in] fragments Number of fragment traces.
    @param[in] points Number of points (spectra or RT bins) per trace.
    @param[in] fragment_stride Offset between fragments.
    @param[in] point_stride Offset between points.
    @param[in] hint Expected apex index, e.g. a detected candidate apex; clamped to the trace.
    @param[in] options Picker settings.
    @param[in,out] scratch Reused consensus buffer (resized as needed).
    @return The picked peak; found is false if the consensus is all zero.
  */
  PickedPeak pickPeakGroup(const float* data, std::size_t fragments, std::size_t points, std::size_t fragment_stride,
                           std::size_t point_stride, std::size_t hint, const PeakPickerOptions& options,
                           std::vector<double>& scratch);

  /// Linear interpolation of @p values (e.g. spectrum RTs) at a fractional index, clamped to the ends.
  double interpolateAt(const std::vector<double>& values, double index);

  /// Short name of a boundary reason ("valley", "floor", "edge").
  const char* peakBoundaryName(PeakBoundary reason);
}
