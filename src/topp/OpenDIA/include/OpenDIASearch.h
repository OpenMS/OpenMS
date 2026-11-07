// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cstddef>
#include <iosfwd>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace OpenMS::OpenDIA
{
  struct TransitionQuery
  {
    std::string_view precursor_id;
    double precursor_mz = 0.0;
    double fragment_mz = 0.0;
    double library_rt = 0.0;
    bool has_library_rt = false;
    /// Library precursor mobility in inverse reduced mobility (1/K0).
    double precursor_im = 0.0;
    bool has_precursor_im = false;
    bool is_decoy = false;
  };

  struct SearchOptions
  {
    double fragment_tolerance_da = 0.02;
    double ion_mobility_tolerance_vssc = 0.03;
    double rt_window_seconds = 0.0;
    double peak_gap_seconds = 30.0;
    std::size_t min_transitions = 1;
    bool split_local_peaks = false;
    double peak_valley_fraction = 0.1;
    std::size_t min_peak_scans = 3;
  };

  /**
    @brief Stream DIA spectra and report fragment-group peak candidates.

    MS2 spectra are packed into bounded m/z-sorted grid batches by isolation
    window lane. Candidate blocks scan each batch in order, updating their
    owned candidate state without retaining the run.
  */
  class StreamSearch final : public Interfaces::IMSDataConsumer
  {
  public:
    StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output,
                 std::ostream* candidate_features = nullptr);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;

    /// Flush active groups at end of input.
    void finish();

    /// Whether compatible raw and library mobility metadata activated filtering.
    bool ionMobilityFilteringUsed() const;

  private:
    struct Candidate
    {
      std::string id;
      double precursor_mz = 0.0;
      double library_rt = 0.0;
      bool has_library_rt = false;
      double precursor_im = 0.0;
      bool has_precursor_im = false;
      bool is_decoy = false;
      bool active = false;
      bool awaiting_peak = false;
      double reopen_threshold = 0.0;
      std::size_t below_peak_scans = 0;
      double start_rt = 0.0;
      double apex_rt = 0.0;
      double last_rt = 0.0;
      double last_intensity = 0.0;
      double apex_intensity = 0.0;
      double area = 0.0;
      std::size_t points = 0;
      std::size_t max_transition_support = 0;
    };

    struct Query
    {
      double fragment_mz = 0.0;
      std::size_t candidate = 0;
    };

    struct LaneQuery
    {
      std::size_t query_index = 0;
      std::size_t order = 0;
    };

    struct GridPoint
    {
      double mz = 0.0;
      double intensity = 0.0;
      float ion_mobility = -1.0f;
    };

    struct GridCell
    {
      double rt = 0.0;
      double scalar_im = 0.0;
      std::size_t lane = 0;
      std::size_t point_begin = 0;
      std::size_t point_end = 0;
      bool has_per_peak_im = false;
      bool has_scalar_im = false;
    };

    struct GridLane
    {
      std::vector<std::pair<double, double>> isolation_windows;
      std::vector<std::vector<LaneQuery>> candidate_blocks;
    };

    struct CompletedCandidate
    {
      std::size_t candidate_index = 0;
      std::size_t cell_order = 0;
      std::size_t query_order = 0;
      double start_rt = 0.0;
      double apex_rt = 0.0;
      double last_rt = 0.0;
      double area = 0.0;
      double apex_intensity = 0.0;
      std::size_t points = 0;
      std::size_t transition_support = 0;
    };

    static constexpr std::size_t CANDIDATES_PER_BLOCK = 256;
    static constexpr std::size_t MAX_BATCH_CELLS = 32;
    static constexpr std::size_t MAX_BATCH_POINTS = 1'000'000;

    void consumeSpectrum_(const MSSpectrum& spectrum);
    void prepareLaneQueries_(GridLane& lane);
    void flushBatch_();
    void updateCandidate_(std::size_t candidate, double rt, double intensity, std::size_t transition_support,
                          std::size_t cell_order, std::size_t query_order,
                          std::vector<CompletedCandidate>& completed);
    void writeCandidate_(std::size_t candidate_index, const Candidate& candidate);
    void writeCandidate_(const CompletedCandidate& candidate);
    void writeCandidateFeatures_(std::size_t candidate_index, double start_rt, double apex_rt, double end_rt,
                                double area, double apex_intensity, std::size_t points,
                                std::size_t transition_support);
    bool queryMatches_(const Query& query, double rt, double peak_im, bool has_peak_im,
                       const std::vector<std::pair<double, double>>& isolation_windows) const;

    SearchOptions options_;
    std::ostream& output_;
    std::ostream* candidate_features_ = nullptr;
    std::size_t candidate_feature_row_ = 0;
    std::vector<Candidate> candidates_;
    std::vector<Query> queries_;
    std::vector<std::size_t> candidate_query_offsets_;
    std::vector<std::size_t> candidate_query_indices_;
    std::vector<std::size_t> candidates_by_mz_;
    std::vector<GridCell> batch_cells_;
    std::vector<GridPoint> grid_points_;
    std::vector<GridLane> grid_lanes_;
    std::map<std::vector<std::pair<double, double>>, std::size_t> lane_indices_;
    std::vector<double> candidate_scan_intensity_;
    std::vector<std::size_t> candidate_scan_support_;
    std::vector<std::size_t> candidate_generation_;
    std::vector<std::size_t> candidate_first_match_order_;
    std::size_t generation_ = 0;
    bool finished_ = false;
    bool has_library_im_candidates_ = false;
    bool ion_mobility_filtering_used_ = false;
  };
}
