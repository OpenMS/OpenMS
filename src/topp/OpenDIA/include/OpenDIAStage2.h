// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include "OpenDIATraces.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// Feature order of stage2Features() (the evidence block reuses columns 1..7). Feature 0 is unused (zero).
  enum class Stage2Feature : std::size_t
  {
    CNN_Z = 0,             ///< Unused, zero (kept so the column indices stay stable).
    FRAGMENT_FRACTION,     ///< Share of library fragments with signal in the core window.
    FRAGMENT_COELUTION,    ///< Mean pairwise Pearson correlation of fragment traces in the core window.
    LIBRARY_CORRELATION,   ///< Pearson correlation of sqrt core fragment areas and sqrt library intensities.
    MS2_ABS_PPM,           ///< |intensity-weighted fragment ppm error - run offset|, capped.
    MS1_ABS_PPM,           ///< |intensity-weighted MS1 monoisotope ppm error - run offset|, capped.
    MS1_COELUTION,         ///< Pearson correlation of the MS1 monoisotope trace and the summed fragment trace.
    MS1_ISOTOPE,           ///< 1 if both the MS1 monoisotope and first isotope carry signal in the core window.
    COUNT
  };
  constexpr std::size_t STAGE2_FEATURES = static_cast<std::size_t>(Stage2Feature::COUNT);

  struct Stage2ExtractionOptions
  {
    double rt_half_width = 15.0;     ///< Tensor RT half width in seconds (TraceCollector).
    double core_half_width = 4.0;    ///< Half width around the apex used by every feature.
    double ms1_offset_ppm = 0.0;     ///< Run MS1 offset subtracted before taking the absolute error.
    double ms2_offset_ppm = 0.0;     ///< Run fragment offset subtracted before taking the absolute error.
    double ppm_cap = 20.0;           ///< Absolute errors are capped here; missing signal also gets the cap.
    double im_half_width = 0.03;     ///< Tensor mobility half width in 1/K0 (TraceCollector im_half_width; matched filter only).
  };

  /**
    @brief Extraction features 1..7 of one window (feature 0 is unused).

    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats, fragment-major as in the collector.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence from TraceCollector::enableEvidence().
    @param[in] options Window geometry and run mass offsets.
    @return All features with feature 0 set to zero. Undefined correlations are zero.
  */
  std::array<double, STAGE2_FEATURES> stage2Features(const float* tensor, const std::array<double, 6>& library_intensities,
                                                     const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                     const Stage2ExtractionOptions& options);

  /// Peak-group boundary features: the window is re-picked with pickPeakGroup and scored over [left, right].
  enum class Stage2BoundaryFeature : std::size_t
  {
    LOG_FWHM = 0,          ///< log1p of the interpolated full width at half maximum in seconds.
    LOG_WIDTH,             ///< log1p of the picked base width in seconds.
    BOUNDARY_COELUTION,    ///< Mean pairwise Pearson of fragment traces over the picked peak (0 if fewer than 3 scans).
    BOUNDARY_LIBRARY_CORRELATION, ///< Pearson of sqrt fragment areas over the picked peak and sqrt library intensities.
    APEX_AGREEMENT,        ///< Share of library fragments whose maximum in the peak is within one scan of the consensus apex.
    COUNT
  };
  constexpr std::size_t STAGE2_BOUNDARY_FEATURES = static_cast<std::size_t>(Stage2BoundaryFeature::COUNT);
  /// Boundary column names, in Stage2BoundaryFeature order.
  const std::array<std::string, STAGE2_BOUNDARY_FEATURES>& stage2BoundaryFeatureNames();
  /**
    @brief Boundary features of one window, identical for targets and decoys and free of training.

    The ion-mobility-summed fragment traces are restricted to the RT bins that received an isolating MS2 scan (bins without a
    scan are not zero signal), the peak group around the window apex (centre bin) is picked with pickPeakGroup defaults, and
    width, co-elution, library correlation and apex agreement are measured over the picked boundaries.

    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats, fragment-major as in the collector.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @param[in] options Window geometry (only rt_half_width is used).
    @return Features in Stage2BoundaryFeature order; all zero if fewer than 3 scan bins or no signal.
  */
  std::array<double, STAGE2_BOUNDARY_FEATURES> stage2BoundaryFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                      const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                                      const Stage2ExtractionOptions& options);

  /// Evidence block (PLAN_EXTRA): core-window mass error, MS1 and ion mobility.
  enum class Stage2ExtraFeature : std::size_t
  {
    FRAGMENT_FRACTION = 0, ///< Stage2Feature::FRAGMENT_FRACTION.
    FRAGMENT_COELUTION,    ///< Stage2Feature::FRAGMENT_COELUTION.
    LIBRARY_CORRELATION,   ///< Stage2Feature::LIBRARY_CORRELATION.
    MS2_ABS_PPM,           ///< Stage2Feature::MS2_ABS_PPM.
    MS1_ABS_PPM,           ///< Stage2Feature::MS1_ABS_PPM.
    MS1_COELUTION,         ///< Stage2Feature::MS1_COELUTION.
    MS1_ISOTOPE,           ///< Stage2Feature::MS1_ISOTOPE.
    MS1_LOG_AREA,          ///< log1p of the MS1 monoisotope intensity summed over the core bins.
    MS2_LOG_AREA,          ///< log1p of the fragment intensity summed over the core scan bins.
    IM_CENTROID_OFFSET,    ///< |intensity-weighted mobility bin centroid - centre| in bins (core window; 4 without signal).
    IM_COELUTION,          ///< Mean pairwise Pearson of the fragments' mobility profiles (core window).
    IM_CENTRAL_SHARE,      ///< Share of core fragment intensity in the central 4 of the 8 mobility bins.
    COUNT
  };
  constexpr std::size_t STAGE2_EXTRA_FEATURES = static_cast<std::size_t>(Stage2ExtraFeature::COUNT);
  const std::array<std::string, STAGE2_EXTRA_FEATURES>& stage2ExtraFeatureNames();
  /**
    @brief Evidence block of one window; identical for targets and decoys and free of training.
    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats ([6,32,8]).
    @param[in] options Window geometry and run mass offsets (as stage2Features()).
  */
  std::array<double, STAGE2_EXTRA_FEATURES> stage2ExtraFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                                const Stage2ExtractionOptions& options);

  /// Matched-filter block (PLAN_MT) from the existing 6 x 32 x 8 trace tensor.
  enum class Stage2MatchedFilterFeature : std::size_t
  {
    Z_MED = 0, ///< slog of the median per-fragment robust z at the centre bin.
    Z_TRIM,    ///< slog of the mean robust z without its largest value.
    N_Z2,      ///< Number of fragments with robust z >= 2.
    N_Z3,      ///< Number of fragments with robust z >= 3.
    LOO_MED,   ///< Median leave-one-out Pearson correlation over centre +/- 2 RT bins.
    LOO_N,     ///< Number of leave-one-out correlations > 0.5.
    COS,       ///< Cosine of the centre-bin fragment vector with sqrt library intensities.
    DOM,       ///< log of the largest observed share over its library share.
    LLR,       ///< slog of sum_k lw_k * log1p(max(z_k, 0)).
    COUNT
  };
  constexpr std::size_t STAGE2_MATCHED_FILTER_FEATURES = static_cast<std::size_t>(Stage2MatchedFilterFeature::COUNT);
  const std::array<std::string, STAGE2_MATCHED_FILTER_FEATURES>& stage2MatchedFilterFeatureNames();
  /**
    @brief Robust z of the 3-tap box sum of @p x at @p centre against the background outside [@p a, @p b).

    s = 3-tap box sum (zero padded). Background is s[0, max(a-1,0)) plus s[b+1, end); with fewer than 3 values it is
    s[0, max(centre-3,0)) plus s[centre+4, end). z = (s[centre] - med) / (1.4826 * MAD + sqrt(max(med,0)) + 1), 0 when
    the background is empty.
  */
  double stage2RobustZ(const std::vector<double>& x, std::size_t centre, std::size_t a, std::size_t b);
  /**
    @brief Compute the matched-filter (PLAN_MT) features from the trace tensor.

    Channels are non-padding fragments (fewer than 2 gives zeros). Each RT bin is weighted over mobility with a Gaussian
    (sd 0.008 1/K0, i.e. 0.008 / (2 * options.im_half_width / 8) bins) centred on the smoothed core (centre +/- 2 RT bins)
    mobility apex. Robust z per fragment uses the peak picked around the centre over MS2-scanned bins (fallback centre-1
    .. centre+2). Because the scorer standardises by mean and SD and the z values are heavy tailed, Z_MED, Z_TRIM and
    LLR are compressed with sign(x) * log1p(|x|).
  */
  std::array<double, STAGE2_MATCHED_FILTER_FEATURES> stage2MatchedFilterFeatures(const float* tensor,
                                                                                  const std::array<double, 6>& library_intensities,
                                                                                  const std::array<double, 6>& fragment_mz,
                                                                                  const TraceEvidence& evidence,
                                                                                  const Stage2ExtractionOptions& options);

  /// b/y ion evidence of the isolating MS2 spectrum nearest the window apex (all theoretical ions, not only the library's six).
  enum class Stage2ByIonFeature : std::size_t
  {
    BYION_B = 0,  ///< Number of matched singly charged b ions (TraceEvidence::byion_b).
    BYION_Y,      ///< Number of matched singly charged y ions (TraceEvidence::byion_y).
    COUNT
  };
  constexpr std::size_t STAGE2_BYION_FEATURES = static_cast<std::size_t>(Stage2ByIonFeature::COUNT);
  /// b/y ion column names, in Stage2ByIonFeature order.
  const std::array<std::string, STAGE2_BYION_FEATURES>& stage2ByIonFeatureNames();
  /**
    @brief b/y ion count features of one window; zero if no isolating MS2 spectrum fell inside the core window.
    @param[in] evidence Window evidence collected with TraceCollector::enableByIonCounts().
    @return Features in Stage2ByIonFeature order.
  */
  std::array<double, STAGE2_BYION_FEATURES> stage2ByIonFeatures(const TraceEvidence& evidence);

  /**
    @brief Peak-group evidence of one search candidate.

    These come from the candidate search, not from the trace tensor, and are computed identically for targets and decoys.
  */
  struct SearchPeakFeatures
  {
    std::string precursor_id;
    double rt_apex = 0.0;
    double transition_support = 0.0;   ///< Fragments that formed the peak group.
    double matched_scans = 0.0;        ///< Scans that contributed points.
    double peak_width_seconds = 0.0;
    double abs_rt_error_seconds = 0.0; ///< |apex - (calibrated) library RT|; 0 without a library RT.
  };
  /// Number of search columns.
  constexpr std::size_t STAGE2_SEARCH_FEATURES = 4;
  /// Search column names: transition_support, log_matched_scans, log_peak_width, abs_library_rt_error.
  const std::array<std::string, STAGE2_SEARCH_FEATURES>& stage2SearchFeatureNames();
  /// The four search columns of one candidate: support, log1p(scans), log1p(width), |RT error|.
  std::array<double, STAGE2_SEARCH_FEATURES> stage2SearchFeatures(const SearchPeakFeatures& peak);

}

namespace OpenMS::OpenDIA
{
  /// Number of charge-state sibling features (PLAN_SIB): s_n, s_co, s_best, s_best_d.
  constexpr std::size_t STAGE2_SIBLING_FEATURES = 4;

  /**
    @brief Charge-state sibling features per window (PLAN_SIB).

    Siblings of a window are the windows of other precursors in the same group. With smin the minimum score over kept
    windows (score > not_kept): s_n = group_precursors[group] - 1; over sibling windows j with |rt_j - rt_i| <= rt_tolerance,
    s_co = number of distinct precursors, s_best = max of (score_j if kept else smin), smin if none, and s_best_d = |rt_j - rt_i|
    of the window attaining s_best (nearest on ties), rt_tolerance if none.

    @param[in] group Group index per window (< group_precursors.size()).
    @param[in] precursor Precursor index per window.
    @param[in] rt Apex RT per window (seconds).
    @param[in] score Score per window, not_kept for non-kept windows.
    @param[in] group_precursors Library-wide number of precursors per group.
    @return Row-major n x STAGE2_SIBLING_FEATURES.
  */
  std::vector<float> stage2SiblingFeatures(const std::vector<std::uint32_t>& group, const std::vector<std::uint32_t>& precursor,
                                           const std::vector<double>& rt, const std::vector<float>& score,
                                           const std::vector<std::uint32_t>& group_precursors, float not_kept,
                                           double rt_tolerance = 3.0);
}
