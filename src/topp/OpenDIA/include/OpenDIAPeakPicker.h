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

  /** @brief Settings of the multi-peak native-scan picker (pickXicPeaks). No RT constants: peaks come from the trace shape. */
  struct XicPickerOptions
  {
    std::size_t smoothing_passes = 1; ///< [1 2 1]/4 passes over the consensus before finding maxima and valleys.
    double valley_fraction = 0.5;     ///< Two maxima stay separate peaks only if the valley between them is at most this fraction of the lower one.
    double boundary_fraction = 0.05;  ///< A boundary also stops where the smoothed consensus falls to this fraction of the apex.
    std::size_t correlation_half_window = 2; ///< The co-elution window covers at least this many scans on each side of the apex.
    std::size_t max_peaks = 5;        ///< Peaks reported per trace, best co-elution first.
  };

  /** @brief One peak of a native-scan fragment trace set. Indices are scans of the trace. */
  struct XicPeak
  {
    std::size_t apex = 0;            ///< Raw consensus maximum of the peak.
    std::size_t left = 0;            ///< First scan inside the peak.
    std::size_t right = 0;           ///< Last scan inside the peak.
    double apex_height = 0.0;        ///< Raw consensus (sum of fragments) at the apex.
    double area = 0.0;               ///< Trapezoid area of the raw consensus over RT between the boundaries.
    double fwhm = 0.0;               ///< Full width at half maximum of the raw consensus, in RT units (interpolated).
    double coelution = 0.0;          ///< Sum over fragments of max(0, Pearson r) with the summed other fragments; 0..fragments.
    std::size_t support = 0;         ///< Fragments with signal at the apex scan or an adjacent scan.
    std::size_t signal_scans = 0;    ///< Scans in [left, right] with any signal.
  };

  /**
    @brief Find every elution peak of a set of fragment traces sampled at native scans and rank them by co-elution.

    The consensus (sum of the fragments) is smoothed, every local maximum starts a peak, and adjacent peaks are
    merged while the valley between them is above valley_fraction of the lower maximum (noise bumps on one peak).
    Boundaries are the remaining valleys, or where the smoothed consensus falls to boundary_fraction of the apex.
    Apex, area and FWHM use the raw consensus. Co-elution correlates each fragment with the sum of the others over
    the peak (at least correlation_half_window scans each side). Peaks are returned by co-elution, then apex height.
    Uses no library intensities, labels or decoy status.

    @param[in] data Intensities, element (scan s, fragment k) at data[s * fragments + k]; negative values count as 0.
    @param[in] scans Number of scans (trace points).
    @param[in] fragments Number of fragment traces.
    @param[in] rts Retention time of each scan, ascending (scans values).
    @param[in] options Picker settings.
    @param[out] peaks_found Number of peaks after merging, before the max_peaks cut.
    @param[in,out] scratch Reused work buffer.
    @return Up to max_peaks peaks, best first.
  */
  std::vector<XicPeak> pickXicPeaks(const float* data, std::size_t scans, std::size_t fragments, const double* rts,
                                    const XicPickerOptions& options, std::size_t& peaks_found, std::vector<double>& scratch);

  /// Linear interpolation of @p values (e.g. spectrum RTs) at a fractional index, clamped to the ends.
  double interpolateAt(const std::vector<double>& values, double index);

  /// Short name of a boundary reason ("valley", "floor", "edge").
  const char* peakBoundaryName(PeakBoundary reason);
}
