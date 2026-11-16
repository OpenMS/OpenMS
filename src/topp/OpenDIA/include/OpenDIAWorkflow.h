// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIACalibration.h>
#include <OpenDIAScorer.h>
#include <OpenDIASearch.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief The library of one OpenDIA analysis: precursors with their target-decoy pairing, peptide family, fold, the
    search fragments (up to 6) and the scoring fragments (up to 12).

    Decoys are named DECOY_<target TRAML_ID>. A precursor's family is its parent target's unmodified sequence; the fold
    is FNV-1a 64 of the family mod 3, so a target, its decoy and all charge states and modified forms of one peptide
    share a fold.
  */
  struct AnalysisLibrary
  {
    struct Precursor
    {
      std::string id;          ///< PQP precursor ID.
      std::string traml_id;
      std::string sequence;    ///< Modified sequence.
      std::string family;      ///< Parent target's unmodified sequence.
      std::string pair_key;    ///< TRAML_ID of the target (own for a target, without DECOY_ for a decoy).
      std::vector<std::string> proteins;
      int charge = 0;
      double mz = 0.0;
      double library_rt = 0.0;
      bool has_library_rt = false;
      bool decoy = false;
      std::uint8_t fold = 0;
      std::uint32_t pair = 0;  ///< Dense pair index.
      std::array<double, 12> scoring_mz{};        ///< Library m/z of the scoring fragments (0 = absent).
      std::array<double, 12> scoring_library{};   ///< Relative to the precursor's most intense scoring fragment.
      std::array<bool, 12> scoring_in_search{};   ///< Also one of the search fragments (same annotation).
      std::uint8_t scoring_count = 0;
    };

    OpenSwath::LightTargetedExperiment search;  ///< Search library (calibrated in place by calibrate()).
    std::vector<Precursor> precursors;
    std::unordered_map<std::string, std::uint32_t> index;  ///< PQP precursor ID -> index.
    std::vector<std::pair<std::string, std::string>> traml_ids;
    std::size_t pairs = 0;

    /**
      @brief Read the search library and, if given, the scoring library (same precursor IDs; up to 12 fragments each,
      in transition ID order). Without a scoring library the search fragments are the scoring fragments.
      @throws Exception::InvalidValue On a decoy without its target or a precursor missing from the scoring library.
    */
    static AnalysisLibrary read(const std::string& search_pqp, const std::string& scoring_pqp = "");
  };

  /// Settings of one run.
  struct RunOptions
  {
    CalibrationOptions calibration;
    std::size_t window_panel_targets = 300000;  ///< Windowed second calibration pass (0 disables).
    std::size_t max_peaks = 5;
    std::size_t min_transitions = 1;
    std::size_t population_min_support = 3;     ///< Pairs whose best candidate has fewer supporting fragments are dropped.
    double grid_half_width_s = 15.0;
    double grid_extraction_ppm = 17.5;
    double byion_ppm = 10.0;
    double byion_intensity_min = 300.0;
    double core_half_width_s = 4.0;
    double anchor_apex_s = 5.0;                 ///< An anchor labels the window of its precursor with apex within this.
    std::size_t decoy_seeds = 6000;
    std::uint64_t seed = 0;
    int threads = 1;
    bool im_array_as_vssc = false;
    std::size_t ensemble = 1;  ///< CrossFitScorer::Options::ensemble (PLAN_ENSEMBLE).
    /// Optional: reuse an OpenDIACalibrate JSON and anchor TSV instead of calibrating (reproduction of earlier runs).
    std::string calibration_json;
    std::string anchors_tsv;
  };

  /// One scored candidate window of the run (the population), with its features and result.
  struct RunWindow
  {
    std::uint32_t precursor = 0;  ///< AnalysisLibrary::precursors index.
    std::uint32_t candidate = 0;  ///< Search candidate row.
    SearchCandidate peak;
    float score = CrossFitScorer::NOT_KEPT;
    double q_value = 1.0;
  };

  struct RunResult
  {
    std::string run;
    CalibrationRT calibration;
    std::string calibration_json;  ///< Written calibration (empty when a calibration was reused).
    std::vector<RunWindow> windows;
    std::vector<std::string> feature_names;
    std::vector<float> features;  ///< windows x feature_names, row-major.
    std::vector<std::string> log;
  };

  /**
    @brief The OpenDIA analysis of one DIA run, from raw spectra to scored precursors.

    1. Calibration (OpenDIACalibration): library RT, fragment m/z and 1/K0 to run coordinates and tolerances,
       identification-free; with the windowed second pass when the whole-run pass finds too few anchors.
    2. Candidate search (StreamSearch) with the calibrated library: up to max_peaks peaks per precursor.
    3. Population: every candidate of every pair whose best candidate has >= population_min_support fragments.
    4. One pass over the MS2 spectra collects, per population window, the 6-fragment trace tensor and the apex b/y ion
       counts (TraceCollector) and the 12-fragment grid (FragmentGridCollector).
    5. Features (24): search (4), boundary (5), b/y ions (2), 12-fragment evidence (6), fragment deficit (3), shared-signal
       interference (4).
    6. CrossFitScorer: three folds, seeds = calibration anchor windows (positives), their decoys' best windows and a
       random sample of decoy windows (negatives); q-values by pooled paired target-decoy competition.

    Steps 1, 2 and 4 each read the spectra once (calibration may take up to three passes).
  */
  class RunAnalysis
  {
  public:
    using Log = std::function<void(const std::string&)>;
    RunAnalysis(AnalysisLibrary& library, RunOptions options, Log log = {});

    /// Analyse one run (raw, mzML, ...). The library's search transitions are recalibrated for each run.
    RunResult run(const std::string& path);

  private:
    AnalysisLibrary& library_;
    RunOptions options_;
    Log log_;
  };

  /// Feature names of RunAnalysis, in column order.
  const std::vector<std::string>& runFeatureNames();
}
