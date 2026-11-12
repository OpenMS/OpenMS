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

#include "OpenDIAPeakPicker.h"

#include <cstddef>
#include <cstdint>
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
    /// When > 0, the fragment half width is this many ppm of each fragment m/z, replacing fragment_tolerance_da.
    double fragment_tolerance_ppm = 0.0;
    double ion_mobility_tolerance_vssc = 0.03;
    double rt_window_seconds = 0.0;
    double peak_gap_seconds = 30.0;
    std::size_t min_transitions = 1;
    bool split_local_peaks = false;
    double peak_valley_fraction = 0.1;
    std::size_t min_peak_scans = 3;
    /// Report up to this many separated local apexes per peak group (1 keeps the single global apex).
    std::size_t max_apexes_per_group = 1;
    /// Minimum RT distance between reported apexes of one peak group.
    double min_apex_separation_seconds = 10.0;
    /// Treat a spectrum ion-mobility array with a non-CCS unit as 1/K0 (VSSC).
    /// For mzML files whose unitless "Ion Mobility" array actually stores 1/K0.
    bool ion_mobility_array_as_vssc = false;
    /// Native-scan peak picking: keep every fragment's intensity at every scan of the candidate's isolation window
    /// inside its RT window, then pick peaks from the traces (pickXicPeaks). Ignores peak_gap_seconds,
    /// split_local_peaks and the apex options. Needs rt_window_seconds > 0 (e.g. from a calibration).
    bool xic_peaks = false;
    /// Picker settings of the xic_peaks mode.
    XicPickerOptions xic;
  };

  /**
    @brief Stream DIA spectra and report fragment-group peak candidates.

    MS2 spectra are packed into bounded m/z-sorted grid batches by isolation
    window lane. Candidate blocks scan each batch in order, updating their
    owned candidate state without retaining the run.

    In the xic_peaks mode the same grid matching stores each fragment's matched intensity per scan for the
    candidates whose RT window is open; a candidate is closed once the stream has passed its RT window, its
    native-scan traces (zeros where nothing matched) are peak-picked, and its buffer is released.
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

    /**
      @brief Inspection side channel of the xic_peaks mode: write the native-scan traces of selected precursors.

      One row per (precursor, scan of the candidate's RT window, fragment): precursor_id, scan (global MS2 cell id),
      rt, fragment slot, query fragment m/z, matched intensity (0 where nothing matched) and the observed m/z of the
      matched point (NaN where nothing matched). The rows are exactly the matrix that the peak picker sees.
      Call before the first spectrum. Does not change the search output.

      @param[in] output  Stream that receives a header and the rows
      @param[in] precursor_ids  Precursor ids to dump; ids not in the library are ignored

      @throws Exception::InvalidParameter if the search is not in the xic_peaks mode
    */
    void setXicTraceDump(std::ostream& output, const std::vector<std::string>& precursor_ids);

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
      double mobility_error_sum = 0.0;
      std::size_t mobility_match_count = 0;
      // Local-apex tracking, used only when max_apexes_per_group > 1.
      bool rising = false;
      std::vector<std::pair<double, double>> apexes;
    };

    struct Query
    {
      double fragment_mz = 0.0;
      std::size_t candidate = 0;
      std::uint32_t slot = 0; ///< Fragment position within its candidate (m/z order), the xic trace column.
    };

    /// One matched fragment intensity of the xic_peaks mode.
    struct XicEntry
    {
      std::uint32_t cell = 0;
      std::uint32_t slot = 0;
      float intensity = 0.0f;
    };

    /// One picked peak of the xic_peaks mode, in RT units.
    struct XicRow
    {
      double start_rt = 0.0;
      double apex_rt = 0.0;
      double end_rt = 0.0;
      XicPeak peak;
      std::size_t rank = 0;
      std::size_t peaks_found = 0;
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
      std::uint32_t id = 0;
      std::size_t point_begin = 0;
      std::size_t point_end = 0;
      bool has_per_peak_im = false;
      bool has_scalar_im = false;
    };

    struct GridLane
    {
      std::vector<std::pair<double, double>> isolation_windows;
      std::vector<std::vector<LaneQuery>> candidate_blocks;
      std::vector<double> cell_rts;          ///< xic_peaks: RT of every cell of the lane, acquisition order.
      std::vector<std::uint32_t> cell_ids;   ///< xic_peaks: global id of every cell of the lane.
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
      double mobility_error_sum = 0.0;
      std::size_t mobility_match_count = 0;
      std::vector<std::pair<double, double>> apexes;
    };

    static constexpr std::size_t CANDIDATES_PER_BLOCK = 256;
    static constexpr std::size_t MAX_BATCH_CELLS = 32;
    static constexpr std::size_t MAX_BATCH_POINTS = 1'000'000;

    void consumeSpectrum_(const MSSpectrum& spectrum);
    void prepareLaneQueries_(GridLane& lane);
    void flushBatch_();
    void updateCandidate_(std::size_t candidate, double rt, double intensity, std::size_t transition_support,
                          double mobility_error_sum, std::size_t mobility_match_count,
                          std::size_t cell_order, std::size_t query_order,
                          std::vector<CompletedCandidate>& completed);
    void offerApex_(Candidate& candidate, double rt, double intensity) const;
    std::vector<std::pair<double, double>> closeApexes_(Candidate& candidate) const;
    void writeCandidateRows_(std::size_t candidate_index, double start_rt, double apex_rt, double end_rt,
                             double area, double apex_intensity, std::size_t points,
                             std::size_t transition_support, double mobility_error_sum,
                             std::size_t mobility_match_count,
                             const std::vector<std::pair<double, double>>& apexes);
    void writeCandidate_(std::size_t candidate_index, Candidate& candidate);
    void writeCandidate_(const CompletedCandidate& candidate);
    void writeCandidateFeatures_(std::size_t candidate_index, double start_rt, double apex_rt, double end_rt,
                                double area, double apex_intensity, std::size_t points,
                                std::size_t transition_support, double mobility_error_sum,
                                std::size_t mobility_match_count);
    void closeXicCandidates_(double before_rt);
    std::vector<XicRow> pickXicCandidate_(std::size_t candidate_index, std::vector<float>& matrix,
                                          std::vector<double>& scratch) const;
    void writeXicRow_(std::size_t candidate_index, const XicRow& row);
    /// Scans (global cell id, RT) of every lane isolating the candidate inside its RT window, by cell id.
    std::vector<std::pair<std::uint32_t, double>> xicCells_(std::size_t candidate_index) const;
    void dumpXicTraces_(std::size_t candidate_index);
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
    std::vector<double> candidate_scan_mobility_error_sum_;
    std::vector<std::size_t> candidate_scan_support_;
    std::vector<std::size_t> candidate_scan_mobility_match_count_;
    std::vector<std::size_t> candidate_generation_;
    std::vector<std::size_t> candidate_first_match_order_;
    std::size_t generation_ = 0;
    // xic_peaks state: matched entries per candidate, lanes containing each candidate, close schedule.
    std::vector<std::vector<XicEntry>> xic_entries_;
    std::vector<std::vector<std::uint32_t>> candidate_lanes_;
    std::vector<std::size_t> xic_close_order_;
    std::size_t xic_close_cursor_ = 0;
    std::uint32_t next_cell_id_ = 0;
    // Optional trace dump: selected candidates and the observed m/z of each of their matched entries.
    std::ostream* xic_traces_ = nullptr;
    std::vector<bool> xic_dump_selected_;
    std::vector<std::vector<double>> xic_dump_mz_;
    bool finished_ = false;
    bool has_library_im_candidates_ = false;
    bool ion_mobility_filtering_used_ = false;
  };
}
