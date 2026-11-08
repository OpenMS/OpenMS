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

  /// One window entering the second-stage discriminant.
  struct Stage2Row
  {
    std::string precursor_id;
    bool is_decoy = false;
    int fold = 0;
    std::array<double, STAGE2_FEATURES> features{};
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
    std::array<double, STAGE2_FEATURES> weights{};
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
    @throws std::invalid_argument If folds, features or a training set are invalid (too few positives or negatives).
  */
  std::vector<double> crossfitStage2LDA(const std::vector<Stage2Row>& rows, int folds, const Stage2Options& options,
                                        std::vector<Stage2FoldReport>& report);
}
