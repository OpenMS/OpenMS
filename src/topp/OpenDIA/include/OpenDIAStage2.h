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
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// Feature order of the experimental second-stage score. Feature 0 is the main (initial ranking) score.
  enum class Stage2Feature : std::size_t
  {
    CNN_Z = 0,             ///< CNN logit standardised by its fold's decoy precursor winners (median, 1.4826 MAD).
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
  /// Column names, in Stage2Feature order.
  const std::array<std::string, STAGE2_FEATURES>& stage2FeatureNames();

  struct Stage2ExtractionOptions
  {
    double rt_half_width = 15.0;     ///< Tensor RT half width in seconds (TraceCollector).
    double core_half_width = 4.0;    ///< Half width around the apex used by every feature.
    double ms1_offset_ppm = 0.0;     ///< Run MS1 offset subtracted before taking the absolute error.
    double ms2_offset_ppm = 0.0;     ///< Run fragment offset subtracted before taking the absolute error.
    double ppm_cap = 20.0;           ///< Absolute errors are capped here; missing signal also gets the cap.
  };

  /**
    @brief Extraction features 1..7 of one window (feature 0, the CNN score, is filled by the caller).

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

  /**
    @brief Which feature groups enter the second-stage discriminant. Feature 0 is always the standardised CNN logit.

    Row order: cnn_z, then the seven extraction features, then the search features, then the isotope features, then the
    boundary features, each group only if selected.
  */
  struct Stage2FeatureSet
  {
    bool extraction = true;   ///< Stage2Feature columns 1..7 (tensor and MS1 evidence; the original -stage2 set).
    bool search = false;      ///< The four search peak-group features (needs the OpenDIA -features_out table).
    bool isotope = false;     ///< The MS1 isotope-envelope features (Stage2IsotopeFeature).
    bool boundary = false;    ///< Experimental peak-group boundary features (Stage2BoundaryFeature).
  };
  /**
    @brief Parse a '+'-separated list of groups, e.g. "extraction", "search+isotope"; "extraction_search" means "extraction+search".
    @throws std::invalid_argument On an unknown or repeated group name.
  */
  Stage2FeatureSet parseStage2FeatureSet(const std::string& name);
  /// Column names of a feature set, in row order.
  std::vector<std::string> stage2FeatureSetNames(const Stage2FeatureSet& set);

  /// MS1 isotope-envelope features, following the OpenSWATH MS1 isotope scores with three extra isotopes.
  enum class Stage2IsotopeFeature : std::size_t
  {
    ISOTOPE_CORRELATION = 0, ///< Pearson of core M+0..M+3 intensities and the averagine pattern at the precursor mass.
    ISOTOPE_OVERLAP,         ///< Core M-1/z intensity over M+0 intensity, capped at 10 (a larger peak before the monoisotope).
    ISOTOPE_COELUTION,       ///< Mean Pearson of the M+1 and M+2 traces with the M+0 trace over the core MS1 bins.
    COUNT
  };
  constexpr std::size_t STAGE2_ISOTOPE_FEATURES = static_cast<std::size_t>(Stage2IsotopeFeature::COUNT);
  /// Isotope column names, in Stage2IsotopeFeature order.
  const std::array<std::string, STAGE2_ISOTOPE_FEATURES>& stage2IsotopeFeatureNames();
  /**
    @brief Isotope-envelope features of one window from its MS1 evidence.

    @param[in] evidence Window evidence from TraceCollector::enableEvidence().
    @param[in] precursor_mz Library precursor m/z.
    @param[in] charge Precursor charge; 0 or less gives all-zero features.
    @param[in] options Window geometry (only the RT and core half widths are used).
    @return Features in Stage2IsotopeFeature order; undefined correlations are zero.
  */
  std::array<double, STAGE2_ISOTOPE_FEATURES> stage2IsotopeFeatures(const TraceEvidence& evidence, double precursor_mz, int charge,
                                                                    const Stage2ExtractionOptions& options);

  /// Experimental peak-group boundary features: the window is re-picked with pickPeakGroup and scored over [left, right].
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

  /**
    @brief Peak-group evidence of one search candidate, read from the OpenDIA -features_out table.

    These come from the candidate search, not from the CNN tensor, and are computed identically for targets and decoys.
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
  /// Number of search columns appended by the SEARCH sets.
  constexpr std::size_t STAGE2_SEARCH_FEATURES = 4;
  /// Search column names: transition_support, log_matched_scans, log_peak_width, abs_library_rt_error.
  const std::array<std::string, STAGE2_SEARCH_FEATURES>& stage2SearchFeatureNames();
  /// The four search columns of one candidate: support, log1p(scans), log1p(width), |RT error|.
  std::array<double, STAGE2_SEARCH_FEATURES> stage2SearchFeatures(const SearchPeakFeatures& peak);

  /**
    @brief Read the OpenDIA -features_out candidate table, keyed by candidate_row (the window sample_id).

    @param[in] path TSV with candidate_row, precursor_id, rt_apex, peak_width_seconds, transition_support,
      matched_scans, has_library_rt and apex_rt_error_seconds columns.
    @return Map from candidate_row to its peak-group features.
    @throws std::invalid_argument If the file is missing, a column is absent, a value is not finite or a row is duplicated.
  */
  std::unordered_map<std::string, SearchPeakFeatures> readSearchPeakFeatures(const std::string& path);

  /// One window entering the second-stage discriminant. Every row of one call must have the same feature count.
  struct Stage2Row
  {
    std::string precursor_id;
    bool is_decoy = false;
    int fold = 0;
    std::vector<double> features;
  };

  struct Stage2Options
  {
    int iterations = 3;            ///< Semi-supervised LDA rounds per fold.
    double initial_fdr = 0.15;     ///< Positive cut in round 1, ranked by feature 0 (PyProphet ss_initial_fdr).
    double iteration_fdr = 0.05;   ///< Positive cut in later rounds (PyProphet ss_iteration_fdr).
    double ridge = 1e-3;           ///< Added to the diagonal of the standardised within-class covariance.
    std::size_t min_positives = 10;
  };

  /// Fitted weights of one held-out fold, on standardised features.
  struct Stage2FoldReport
  {
    int fold = 0;
    std::size_t positives = 0;
    std::size_t negatives = 0;
    int iterations = 0;
    std::vector<double> weights;
  };

  /**
    @brief Cross-fitted semi-supervised linear discriminant (PyProphet-style) over the rows.

    For each fold k, the discriminant is trained only on rows of the other folds:
    rank by feature 0, take precursor winners, use target winners at decoy-count q <= initial_fdr as
    positives and every decoy winner as negatives, fit a ridge LDA on standardised features, rescore and
    repeat with iteration_fdr. The final weights score fold k. Scores of different folds are therefore
    on different scales; standardise per fold before pooling (aggregatePairedCNNQValues pooled).

    @param[in] rows Windows with folds 0..folds-1.
    @param[in] folds Number of folds.
    @param[in] options Iteration settings.
    @param[out] report One entry per fold.
    @return One score per row, in input order.
    @throws std::invalid_argument If folds, features (non-finite, or differing counts) or a training set are invalid (too few positives or negatives).
  */
  std::vector<double> crossfitStage2LDA(const std::vector<Stage2Row>& rows, int folds, const Stage2Options& options,
                                        std::vector<Stage2FoldReport>& report);
}
