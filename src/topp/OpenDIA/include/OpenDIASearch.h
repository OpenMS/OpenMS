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

  /// One reported candidate peak, as the feature table row with the same candidate_row (= index in candidates()).
  struct SearchCandidate
  {
    std::uint32_t precursor = 0;    ///< Index of the precursor in the StreamSearch library order (first appearance).
    double rt_start = 0.0;
    double rt_apex = 0.0;
    double rt_end = 0.0;
    double area = 0.0;
    double apex_intensity = 0.0;
    std::uint32_t transition_support = 0;
    std::uint32_t matched_scans = 0;
    double apex_rt_error_seconds = 0.0; ///< apex - (calibrated) library RT; 0 without a library RT.
    double fwhm_seconds = 0.0;
    double coelution = 0.0;
    std::uint16_t peak_rank = 0;
    std::uint16_t peaks_in_window = 0;
  };

  struct SearchOptions
  {
    double fragment_tolerance_da = 0.02;
    /// When > 0, the fragment half width is this many ppm of each fragment m/z, replacing fragment_tolerance_da.
    double fragment_tolerance_ppm = 0.0;
    double ion_mobility_tolerance_vssc = 0.03;
    /// Half width of each precursor's RT window around its (calibrated) library RT; must be > 0. Precursors
    /// without a library RT are searched over the whole run.
    double rt_window_seconds = 0.0;
    /// Peaks with fewer supporting fragments (signal at or next to the apex) are not reported.
    std::size_t min_transitions = 1;
    /// Treat a spectrum ion-mobility array with a non-CCS unit as 1/K0 (VSSC).
    /// For mzML files whose unitless "Ion Mobility" array actually stores 1/K0.
    bool ion_mobility_array_as_vssc = false;
    /// Peak picker settings.
    XicPickerOptions picker;
  };

  /**
    @brief Stream DIA spectra and report the elution peaks of every library precursor.

    MS2 spectra are packed into bounded m/z-sorted batches by isolation-window lane. Precursors are grouped into
    blocks of nearby library RT; for each scan only the blocks whose RT window contains it are matched, in parallel,
    and every open precursor stores the most intense match of each fragment. Once the stream has passed a precursor's
    RT window, its native-scan fragment traces (every scan that isolated it, zeros where nothing matched) are picked
    with pickXicPeaks, the peaks are written and the buffer is released. Nothing run-wide is retained.

    Output (@p output): precursor_id, precursor_mz, rt_start, rt_apex, rt_end, area_score, transition_support,
    matched_scans. Optional feature table (@p candidate_features, header written by the caller): candidate_row,
    precursor_id, precursor_mz, is_decoy, rt_start, rt_apex, rt_end, peak_width_seconds, area, apex_intensity,
    transition_support, matched_scans, has_library_rt, apex_rt_error_seconds, has_library_im, fwhm_seconds,
    coelution, peak_rank, peaks_in_window. Rows are written in order of RT-window end, independent of the thread count.
  */
  class StreamSearch final : public Interfaces::IMSDataConsumer
  {
  public:
    /// Header of the optional feature table, without a line end.
    static const char* featureHeader();

    /// @throws Exception::InvalidParameter for non-positive tolerances or RT window, or a library without fragments
    /// @throws Exception::InvalidValue if transitions of one precursor disagree on the decoy label
    StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output,
                 std::ostream* candidate_features = nullptr);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override;
    void setExpectedSize(std::size_t, std::size_t) override;
    void setExperimentalSettings(const ExperimentalSettings&) override;

    /// Pick and write every precursor still open at end of input.
    void finish();

    /// Also keep every reported candidate in memory (candidates()). Call before the first spectrum.
    void keepCandidates() { keep_candidates_ = true; }
    /// Reported candidates in feature-table order (candidate_row = index); empty unless keepCandidates() was called.
    const std::vector<SearchCandidate>& candidates() const { return kept_; }
    /// Precursor ids and decoy flags, indexed by SearchCandidate::precursor.
    std::size_t precursorCount() const { return candidates_.size(); }
    const std::string& precursorId(std::size_t p) const { return candidates_[p].id; }
    bool precursorIsDecoy(std::size_t p) const { return candidates_[p].is_decoy; }
    bool precursorHasLibraryRT(std::size_t p) const { return candidates_[p].has_library_rt; }

    /// Whether compatible raw and library mobility metadata activated filtering.
    bool ionMobilityFilteringUsed() const;

    /**
      @brief Inspection side channel: write the native-scan traces of selected precursors.

      One row per (precursor, scan of the candidate's RT window, fragment): precursor_id, scan (global MS2 cell id),
      rt, fragment slot, query fragment m/z, matched intensity (0 where nothing matched) and the observed m/z of the
      matched point (NaN where nothing matched). The rows are exactly the matrix that the peak picker sees.
      Call before the first spectrum. Does not change the search output.

      @param[in] output  Stream that receives a header and the rows
      @param[in] precursor_ids  Precursor ids to dump; ids not in the library are ignored
    */
    void setTraceDump(std::ostream& output, const std::vector<std::string>& precursor_ids);

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
    };

    struct Query
    {
      double fragment_mz = 0.0;
      std::size_t candidate = 0;
      std::uint32_t slot = 0; ///< Fragment position within its candidate (m/z order), the trace column.
    };

    /// A fragment query as matched inside a lane: everything the inner loop reads, contiguous.
    struct LaneQuery
    {
      double low_mz = 0.0;
      double high_mz = 0.0;
      double library_rt = 0.0;
      double precursor_im = 0.0;
      std::uint32_t candidate = 0;
      std::uint32_t slot = 0;
      bool has_library_rt = false;
      bool has_precursor_im = false;
    };

    /// One matched fragment intensity (the most intense point of the scan).
    struct Entry
    {
      std::uint32_t cell = 0;
      std::uint16_t slot = 0;
      float intensity = 0.0f;
    };

    /// A match found by a worker, merged serially into entries_.
    struct Match
    {
      std::uint32_t candidate = 0;
      Entry entry;
      double mz = 0.0; ///< Observed m/z, for the trace dump.
    };

    /// CANDIDATES_PER_CHUNK candidates of one lane, adjacent in library RT: the unit of parallel matching.
    struct Chunk
    {
      double first_rt = 0.0;  ///< Earliest RT at which a candidate of the chunk is open (conservative).
      double last_rt = 0.0;   ///< Latest such RT.
      std::vector<LaneQuery> queries; ///< low_mz ascending.
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
      /// Chunks of candidates with a library RT (by RT, so both open-range ends are nondecreasing), then the others.
      std::vector<Chunk> chunks;
      std::size_t timed_chunks = 0;
      std::vector<double> cell_rts;          ///< RT of every cell of the lane, acquisition order.
      std::vector<std::uint32_t> cell_ids;   ///< Global id of every cell of the lane.
    };

    static constexpr std::size_t CANDIDATES_PER_CHUNK = 256;
    static constexpr std::size_t MAX_BATCH_CELLS = 32;
    static constexpr std::size_t MAX_BATCH_POINTS = 1'000'000;

    void consumeSpectrum_(const MSSpectrum& spectrum);
    void prepareLane_(GridLane& lane, std::uint32_t lane_index);
    void flushBatch_();
    void matchChunk_(const GridCell& cell, const Chunk& chunk, std::vector<Match>& matches) const;
    /// Picks and writes close_order_[close_cursor_, end).
    void closeCandidates_(std::size_t end);
    /// Picks one candidate and formats its peak and feature rows (feature rows without the leading row number).
    void pickCandidate_(std::size_t candidate_index, std::vector<float>& matrix, std::vector<double>& scratch,
                        std::string& peak_rows, std::string& feature_rows, std::vector<SearchCandidate>* kept) const;
    /// Scans (global cell id, RT) of every lane isolating the candidate inside its RT window, by cell id.
    std::vector<std::pair<std::uint32_t, double>> candidateCells_(std::size_t candidate_index) const;
    void dumpTraces_(std::size_t candidate_index);
    double halfWidth_(double mz) const;

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
    std::vector<std::pair<std::size_t, std::size_t>> units_;  ///< (batch cell, lane chunk) work units of a batch.
    std::vector<std::vector<Match>> unit_matches_;
    std::vector<GridPoint> grid_points_;
    std::vector<GridLane> grid_lanes_;
    std::map<std::vector<std::pair<double, double>>, std::size_t> lane_indices_;
    std::vector<std::vector<Entry>> entries_;
    std::vector<std::vector<std::uint32_t>> candidate_lanes_;
    std::vector<std::size_t> close_order_;  ///< Candidates by RT-window end; no library RT last.
    std::size_t close_cursor_ = 0;
    std::uint32_t next_cell_id_ = 0;
    // Optional trace dump: selected candidates and the observed m/z of each of their matched entries.
    std::ostream* traces_ = nullptr;
    std::vector<bool> dump_selected_;
    std::vector<std::vector<double>> dump_mz_;
    bool finished_ = false;
    bool keep_candidates_ = false;
    std::vector<SearchCandidate> kept_;
    bool has_library_im_candidates_ = false;
    bool ion_mobility_filtering_used_ = false;
  };
}
