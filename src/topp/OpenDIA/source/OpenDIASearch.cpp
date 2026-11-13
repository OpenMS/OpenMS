// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIASearch.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace OpenMS::OpenDIA
{
  namespace
  {
    // Same text as std::ostream << std::setprecision(precision) << value (printf %.<precision>g), without the stream.
    void appendDouble(std::string& out, double value, int precision)
    {
      char buffer[64];
      const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::general, precision);
      out.append(buffer, result.ptr);
    }

    void appendSize(std::string& out, std::size_t value)
    {
      char buffer[24];
      const auto result = std::to_chars(buffer, buffer + sizeof(buffer), value);
      out.append(buffer, result.ptr);
    }

    struct PrecursorKeyHash
    {
      std::size_t operator()(const std::pair<std::string_view, double>& key) const
      {
        return std::hash<std::string_view>{}(key.first) ^ (std::hash<double>{}(key.second) * 0x9e3779b97f4a7c15ULL);
      }
    };
  }

  const char* StreamSearch::featureHeader()
  {
    return "candidate_row\tprecursor_id\tprecursor_mz\tis_decoy\trt_start\trt_apex\trt_end\t"
           "peak_width_seconds\tarea\tapex_intensity\ttransition_support\tmatched_scans\t"
           "has_library_rt\tapex_rt_error_seconds\thas_library_im\t"
           "fwhm_seconds\tcoelution\tpeak_rank\tpeaks_in_window";
  }

  StreamSearch::StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output,
                             std::ostream* candidate_features) :
    options_(options),
    output_(output),
    candidate_features_(candidate_features)
  {
    const XicPickerOptions& picker = options_.picker;
    if (!std::isfinite(options_.fragment_tolerance_da) || !(options_.fragment_tolerance_da > 0.0) ||
        !std::isfinite(options_.fragment_tolerance_ppm) || options_.fragment_tolerance_ppm < 0.0 ||
        !std::isfinite(options_.ion_mobility_tolerance_vssc) || !(options_.ion_mobility_tolerance_vssc > 0.0) ||
        !std::isfinite(options_.rt_window_seconds) || !(options_.rt_window_seconds > 0.0) || options_.min_transitions == 0 ||
        !std::isfinite(picker.valley_fraction) || !(picker.valley_fraction > 0.0) || picker.valley_fraction > 1.0 ||
        !std::isfinite(picker.boundary_fraction) || picker.boundary_fraction < 0.0 || !(picker.boundary_fraction < 1.0) ||
        picker.max_peaks == 0 || !std::isfinite(picker.boundary_half_widths) || picker.boundary_half_widths < 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Search tolerances, the RT window and the minimum transition count must be positive.");
    }

    // Flatten the library into compact candidate/query arrays; candidates are numbered by first appearance.
    std::unordered_map<std::pair<std::string_view, double>, std::size_t, PrecursorKeyHash> candidate_indices;
    candidate_indices.reserve(queries.size() / 4 + 1);
    candidates_.reserve(queries.size() / 4 + 1);
    queries_.reserve(queries.size());
    for (auto& query : queries)
    {
      if (query.precursor_id.empty() || !std::isfinite(query.precursor_mz) || !(query.precursor_mz > 0.0) ||
          !std::isfinite(query.fragment_mz) || !(query.fragment_mz > 0.0))
      {
        continue;
      }
      const auto [candidate_it, inserted] = candidate_indices.try_emplace({query.precursor_id, query.precursor_mz}, candidates_.size());
      const std::size_t candidate_index = candidate_it->second;
      if (inserted)
      {
        Candidate candidate;
        candidate.id = std::string(query.precursor_id);
        candidate.precursor_mz = query.precursor_mz;
        candidate.library_rt = query.library_rt;
        candidate.has_library_rt = query.has_library_rt;
        candidate.precursor_im = query.precursor_im;
        candidate.has_precursor_im = query.has_precursor_im && std::isfinite(query.precursor_im) && query.precursor_im > 0.0;
        candidate.is_decoy = query.is_decoy;
        has_library_im_candidates_ = has_library_im_candidates_ || candidate.has_precursor_im;
        candidates_.push_back(std::move(candidate));
      }
      else if (candidates_[candidate_index].is_decoy != query.is_decoy)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transitions grouped under one precursor have inconsistent target/decoy labels.",
                                      candidates_[candidate_index].id);
      }
      queries_.push_back({query.fragment_mz, candidate_index});
    }
    if (queries_.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "The transition library contains no usable fragment queries.");
    }

    // Each candidate's fragment queries, m/z ascending; the position is the trace column (slot).
    std::sort(queries_.begin(), queries_.end(), [](const Query& lhs, const Query& rhs) { return lhs.fragment_mz < rhs.fragment_mz; });
    candidate_query_offsets_.assign(candidates_.size() + 1, 0);
    for (const Query& query : queries_)
    {
      ++candidate_query_offsets_[query.candidate + 1];
    }
    std::partial_sum(candidate_query_offsets_.begin(), candidate_query_offsets_.end(), candidate_query_offsets_.begin());
    candidate_query_indices_.resize(queries_.size());
    std::vector<std::size_t> next_query(candidate_query_offsets_.begin(), candidate_query_offsets_.end() - 1);
    for (std::size_t query_index = 0; query_index < queries_.size(); ++query_index)
    {
      const std::size_t candidate_index = queries_[query_index].candidate;
      queries_[query_index].slot = static_cast<std::uint32_t>(next_query[candidate_index] - candidate_query_offsets_[candidate_index]);
      candidate_query_indices_[next_query[candidate_index]++] = query_index;
    }

    candidates_by_mz_.resize(candidates_.size());
    std::iota(candidates_by_mz_.begin(), candidates_by_mz_.end(), 0);
    std::sort(candidates_by_mz_.begin(), candidates_by_mz_.end(), [this](std::size_t lhs, std::size_t rhs)
    {
      return candidates_[lhs].precursor_mz < candidates_[rhs].precursor_mz;
    });

    batch_cells_.reserve(MAX_BATCH_CELLS);
    grid_points_.reserve(65536);
    entries_.resize(candidates_.size());
    candidate_lanes_.resize(candidates_.size());
    // Close schedule: a candidate is complete once the stream has passed the end of its RT window.
    close_order_.resize(candidates_.size());
    std::iota(close_order_.begin(), close_order_.end(), 0);
    std::stable_sort(close_order_.begin(), close_order_.end(), [this](std::size_t lhs, std::size_t rhs)
    {
      const auto end = [this](std::size_t c)
      {
        return candidates_[c].has_library_rt ? candidates_[c].library_rt + options_.rt_window_seconds
                                             : std::numeric_limits<double>::infinity();
      };
      return end(lhs) < end(rhs);
    });
  }

  void StreamSearch::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (!finished_)
    {
      consumeSpectrum_(spectrum);
    }
  }

  void StreamSearch::consumeChromatogram(MSChromatogram&)
  {
  }

  void StreamSearch::setExpectedSize(std::size_t, std::size_t)
  {
  }

  void StreamSearch::setExperimentalSettings(const ExperimentalSettings&)
  {
  }

  double StreamSearch::halfWidth_(double mz) const
  {
    return options_.fragment_tolerance_ppm > 0.0 ? mz * options_.fragment_tolerance_ppm * 1e-6 : options_.fragment_tolerance_da;
  }

  void StreamSearch::consumeSpectrum_(const MSSpectrum& spectrum)
  {
    // An empty MS2 scan carries no match but is still a native trace point (a zero).
    if (spectrum.getMSLevel() != 2 || spectrum.getPrecursors().empty())
    {
      return;
    }
    const double rt = spectrum.getRT();
    if (!std::isfinite(rt))
    {
      return;
    }

    std::vector<std::pair<double, double>> isolation_windows;
    for (const Precursor& precursor : spectrum.getPrecursors())
    {
      const double center = precursor.getMZ();
      const double lower_offset = precursor.getIsolationWindowLowerOffset();
      const double upper_offset = precursor.getIsolationWindowUpperOffset();
      if (!std::isfinite(center) || !(center > 0.0) || !std::isfinite(lower_offset) || lower_offset < 0.0 ||
          !std::isfinite(upper_offset) || upper_offset < 0.0 || lower_offset + upper_offset <= 0.0)
      {
        continue;
      }
      isolation_windows.emplace_back(center - lower_offset, center + upper_offset);
    }
    if (isolation_windows.empty())
    {
      return;
    }
    // Merge overlapping isolation windows so a candidate is counted once per scan.
    if (isolation_windows.size() > 1)
    {
      std::sort(isolation_windows.begin(), isolation_windows.end());
      std::size_t merged_count = 0;
      for (std::size_t i = 1; i < isolation_windows.size(); ++i)
      {
        if (isolation_windows[i].first <= isolation_windows[merged_count].second)
        {
          isolation_windows[merged_count].second = std::max(isolation_windows[merged_count].second, isolation_windows[i].second);
        }
        else
        {
          isolation_windows[++merged_count] = isolation_windows[i];
        }
      }
      isolation_windows.resize(merged_count + 1);
    }

    // Bruker TDF exposes one calibrated VSSC mobility value per peak. Other compatible OpenMS readers may expose a
    // scalar VSSC drift time per spectrum. Unsupported or missing units keep m/z-only matching.
    const DataArrays::FloatDataArray* im_values = nullptr;
    double scalar_im = 0.0;
    bool has_scalar_im = false;
    if (spectrum.containsIMData())
    {
      const auto [im_index, im_unit] = spectrum.getIMData();
      const auto& arrays = spectrum.getFloatDataArrays();
      if (im_index < arrays.size() && arrays[im_index].size() == spectrum.size() &&
          (im_unit == DriftTimeUnit::VSSC || (options_.ion_mobility_array_as_vssc && im_unit != DriftTimeUnit::CCS)))
      {
        im_values = &arrays[im_index];
      }
    }
    if (im_values == nullptr && spectrum.getDriftTime() != IMTypes::DRIFTTIME_NOT_SET &&
        spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC)
    {
      scalar_im = spectrum.getDriftTime();
      has_scalar_im = std::isfinite(scalar_im) && scalar_im > 0.0;
    }
    const bool has_raw_im = im_values != nullptr || has_scalar_im;
    ion_mobility_filtering_used_ = ion_mobility_filtering_used_ || (has_raw_im && has_library_im_candidates_);
    // TDF spectra commonly repeat one mobility value for every centroid: store it once per cell.
    bool uniform_per_peak_im = im_values != nullptr && !im_values->empty() && std::isfinite((*im_values)[0]) && (*im_values)[0] > 0.0f;
    if (uniform_per_peak_im)
    {
      const float first_im = (*im_values)[0];
      uniform_per_peak_im = std::all_of(im_values->begin(), im_values->end(), [first_im](float value) { return value == first_im; });
      if (uniform_per_peak_im)
      {
        scalar_im = first_im;
        has_scalar_im = true;
      }
    }

    const auto [lane_it, inserted] = lane_indices_.try_emplace(isolation_windows, grid_lanes_.size());
    const std::size_t lane_index = lane_it->second;
    if (inserted)
    {
      grid_lanes_.push_back({std::move(isolation_windows), {}, 0, {}, {}});
      prepareLane_(grid_lanes_.back(), static_cast<std::uint32_t>(lane_index));
    }

    GridCell cell;
    cell.rt = rt;
    cell.scalar_im = scalar_im;
    cell.lane = lane_index;
    cell.has_per_peak_im = im_values != nullptr && !uniform_per_peak_im;
    cell.has_scalar_im = has_scalar_im;
    cell.id = next_cell_id_++;
    grid_lanes_[lane_index].cell_rts.push_back(rt);
    grid_lanes_[lane_index].cell_ids.push_back(cell.id);

    // Pack m/z-sorted cells into a bounded batch so one worker team handles several scans.
    if (!batch_cells_.empty() && (batch_cells_.size() >= MAX_BATCH_CELLS || grid_points_.size() >= MAX_BATCH_POINTS))
    {
      flushBatch_();
    }
    cell.point_begin = grid_points_.size();
    grid_points_.reserve(grid_points_.size() + spectrum.size());
    for (std::size_t peak_index = 0; peak_index < spectrum.size(); ++peak_index)
    {
      const auto& peak = spectrum[peak_index];
      const double peak_mz = peak.getMZ();
      const double peak_intensity = peak.getIntensity();
      if (!std::isfinite(peak_mz) || peak_mz < 0.0 || !std::isfinite(peak_intensity))
      {
        continue;
      }
      const float peak_im = im_values != nullptr ? (*im_values)[peak_index] : -1.0f;
      grid_points_.push_back({peak_mz, peak_intensity, cell.has_per_peak_im ? peak_im : -1.0f});
    }
    cell.point_end = grid_points_.size();
    const auto points_begin = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_begin);
    const auto points_end = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_end);
    const auto by_mz = [](const GridPoint& lhs, const GridPoint& rhs) { return lhs.mz < rhs.mz; };
    if (!std::is_sorted(points_begin, points_end, by_mz))
    {
      // Stable, so equal-m/z points keep their spectrum order (ties of the most intense match).
      std::stable_sort(points_begin, points_end, by_mz);
    }
    batch_cells_.push_back(cell);
    if (batch_cells_.size() >= MAX_BATCH_CELLS || grid_points_.size() >= MAX_BATCH_POINTS)
    {
      flushBatch_();
    }
  }

  void StreamSearch::prepareLane_(GridLane& lane, std::uint32_t lane_index)
  {
    // The lane's candidates (precursor m/z inside an isolation window, inclusive), by library RT.
    std::vector<std::size_t> members;
    for (const auto& [low_mz, high_mz] : lane.isolation_windows)
    {
      auto candidate = std::lower_bound(candidates_by_mz_.begin(), candidates_by_mz_.end(), low_mz,
                                        [this](std::size_t index, double mz) { return candidates_[index].precursor_mz < mz; });
      for (; candidate != candidates_by_mz_.end() && candidates_[*candidate].precursor_mz <= high_mz; ++candidate)
      {
        members.push_back(*candidate);
        candidate_lanes_[*candidate].push_back(lane_index);
      }
    }
    std::stable_partition(members.begin(), members.end(), [this](std::size_t c) { return candidates_[c].has_library_rt; });
    const auto timed_end = std::find_if(members.begin(), members.end(), [this](std::size_t c) { return !candidates_[c].has_library_rt; });
    std::stable_sort(members.begin(), timed_end, [this](std::size_t lhs, std::size_t rhs) { return candidates_[lhs].library_rt < candidates_[rhs].library_rt; });
    const std::size_t timed = static_cast<std::size_t>(timed_end - members.begin());

    // Chunks of CANDIDATES_PER_CHUNK consecutive candidates; timed and untimed candidates never share a chunk.
    const double window = options_.rt_window_seconds;
    const auto addChunk = [&](std::size_t first, std::size_t last, bool is_timed)
    {
      Chunk chunk;
      if (is_timed)
      {
        // Conservative open range (1 ms margin against rounding); the exact window test is per query.
        chunk.first_rt = candidates_[members[first]].library_rt - window - 1e-3;
        chunk.last_rt = candidates_[members[last - 1]].library_rt + window + 1e-3;
      }
      else
      {
        chunk.first_rt = -std::numeric_limits<double>::infinity();
        chunk.last_rt = std::numeric_limits<double>::infinity();
      }
      for (std::size_t m = first; m < last; ++m)
      {
        const std::size_t c = members[m];
        const Candidate& candidate = candidates_[c];
        for (std::size_t q = candidate_query_offsets_[c]; q < candidate_query_offsets_[c + 1]; ++q)
        {
          const std::size_t query_index = candidate_query_indices_[q];
          const Query& query = queries_[query_index];
          LaneQuery lane_query;
          lane_query.low_mz = query.fragment_mz - halfWidth_(query.fragment_mz);
          lane_query.high_mz = query.fragment_mz + halfWidth_(query.fragment_mz);
          lane_query.library_rt = candidate.library_rt;
          lane_query.precursor_im = candidate.precursor_im;
          lane_query.candidate = static_cast<std::uint32_t>(c);
          lane_query.slot = query.slot;
          lane_query.has_library_rt = candidate.has_library_rt;
          lane_query.has_precursor_im = candidate.has_precursor_im;
          chunk.queries.push_back(lane_query);
        }
      }
      // The point cursor of a chunk only moves forward, so its queries are walked by lower bound.
      std::stable_sort(chunk.queries.begin(), chunk.queries.end(), [](const LaneQuery& lhs, const LaneQuery& rhs) { return lhs.low_mz < rhs.low_mz; });
      lane.chunks.push_back(std::move(chunk));
    };
    for (std::size_t first = 0; first < timed; first += CANDIDATES_PER_CHUNK)
    {
      addChunk(first, std::min(timed, first + CANDIDATES_PER_CHUNK), true);
    }
    lane.timed_chunks = lane.chunks.size();
    for (std::size_t first = timed; first < members.size(); first += CANDIDATES_PER_CHUNK)
    {
      addChunk(first, std::min(members.size(), first + CANDIDATES_PER_CHUNK), false);
    }
  }

  void StreamSearch::matchChunk_(const GridCell& cell, const Chunk& chunk, std::vector<Match>& matches) const
  {
    const double window = options_.rt_window_seconds;
    const double tolerance_im = options_.ion_mobility_tolerance_vssc;
    const bool has_peak_im = cell.has_per_peak_im || cell.has_scalar_im;
    const auto points_begin = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_begin);
    const auto points_end = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_end);
    auto cursor = std::lower_bound(points_begin, points_end, chunk.queries.front().low_mz,
                                   [](const GridPoint& point, double mz) { return point.mz < mz; });
    for (const LaneQuery& query : chunk.queries)
    {
      while (cursor != points_end && cursor->mz < query.low_mz)
      {
        ++cursor;
      }
      if (cursor == points_end)
      {
        break;
      }
      if (query.has_library_rt && std::abs(cell.rt - query.library_rt) > window)
      {
        continue;
      }
      // The most intense point inside the m/z (and, when both sides have it, 1/K0) window; the first one wins ties.
      double intensity = 0.0;
      double matched_mz = 0.0;
      bool matched = false;
      for (auto point = cursor; point != points_end && point->mz <= query.high_mz; ++point)
      {
        if (query.has_precursor_im && has_peak_im)
        {
          const double peak_im = cell.has_per_peak_im ? point->ion_mobility : cell.scalar_im;
          if (!std::isfinite(peak_im) || std::abs(peak_im - query.precursor_im) > tolerance_im)
          {
            continue;
          }
        }
        if (!matched || point->intensity > intensity)
        {
          intensity = point->intensity;
          matched_mz = point->mz;
        }
        matched = true;
      }
      if (matched)
      {
        matches.push_back({query.candidate, {cell.id, query.slot, static_cast<float>(intensity)}, matched_mz});
      }
    }
  }

  void StreamSearch::flushBatch_()
  {
    if (batch_cells_.empty())
    {
      return;
    }
    // Work units: (cell, chunk) pairs whose chunk has a candidate open at the cell's RT.
    units_.clear();
    for (std::size_t cell_order = 0; cell_order < batch_cells_.size(); ++cell_order)
    {
      const GridCell& cell = batch_cells_[cell_order];
      const GridLane& lane = grid_lanes_[cell.lane];
      if (cell.point_begin == cell.point_end)
      {
        continue;
      }
      // Timed chunks are ordered by library RT, so both ends of their open ranges are nondecreasing.
      const auto timed_end = lane.chunks.begin() + static_cast<std::ptrdiff_t>(lane.timed_chunks);
      auto chunk = std::lower_bound(lane.chunks.begin(), timed_end, cell.rt, [](const Chunk& c, double rt) { return c.last_rt < rt; });
      for (; chunk != timed_end && chunk->first_rt <= cell.rt; ++chunk)
      {
        units_.emplace_back(cell_order, static_cast<std::size_t>(chunk - lane.chunks.begin()));
      }
      for (std::size_t c = lane.timed_chunks; c < lane.chunks.size(); ++c)
      {
        units_.emplace_back(cell_order, c);
      }
    }
    if (unit_matches_.size() < units_.size())
    {
      unit_matches_.resize(units_.size());
    }
    const std::ptrdiff_t unit_count = static_cast<std::ptrdiff_t>(units_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) if(unit_count > 1)
#endif
    for (std::ptrdiff_t u = 0; u < unit_count; ++u)
    {
      const auto [cell_order, chunk] = units_[static_cast<std::size_t>(u)];
      std::vector<Match>& matches = unit_matches_[static_cast<std::size_t>(u)];
      matches.clear();
      const GridCell& cell = batch_cells_[cell_order];
      matchChunk_(cell, grid_lanes_[cell.lane].chunks[chunk], matches);
    }
    // Serial merge in unit order. A (cell, fragment) has at most one match, so the order within a candidate does
    // not change any result.
    for (std::size_t u = 0; u < units_.size(); ++u)
    {
      for (const Match& match : unit_matches_[u])
      {
        entries_[match.candidate].push_back(match.entry);
        if (traces_ != nullptr && dump_selected_[match.candidate])
        {
          dump_mz_[match.candidate].push_back(match.mz);
        }
      }
    }

    // Spectra arrive in RT order: no later cell can match a candidate whose RT window ended before this one.
    const double last_rt = batch_cells_.back().rt;
    batch_cells_.clear();
    grid_points_.clear();
    std::size_t end = close_cursor_;
    while (end < close_order_.size())
    {
      const Candidate& candidate = candidates_[close_order_[end]];
      if (!(candidate.has_library_rt && candidate.library_rt + options_.rt_window_seconds < last_rt))
      {
        break;
      }
      ++end;
    }
    closeCandidates_(end);
  }

  void StreamSearch::closeCandidates_(std::size_t end)
  {
    const std::size_t begin = close_cursor_;
    if (end <= begin)
    {
      return;
    }
    const std::size_t count = end - begin;
    std::vector<std::string> peak_rows(count), feature_rows(count);
#ifdef _OPENMP
#pragma omp parallel if(count >= 64)
#endif
    {
      std::vector<float> matrix;
      std::vector<double> scratch;
#ifdef _OPENMP
#pragma omp for schedule(dynamic, 32)
#endif
      for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(count); ++i)
      {
        const std::size_t candidate_index = close_order_[begin + static_cast<std::size_t>(i)];
        pickCandidate_(candidate_index, matrix, scratch, peak_rows[static_cast<std::size_t>(i)], feature_rows[static_cast<std::size_t>(i)]);
        if (traces_ == nullptr || !dump_selected_[candidate_index])
        {
          std::vector<Entry>().swap(entries_[candidate_index]);
        }
      }
    }
    std::string numbered;
    for (std::size_t i = 0; i < count; ++i)
    {
      output_ << peak_rows[i];
      if (candidate_features_ != nullptr)
      {
        // Feature rows are formatted without their running row number.
        numbered.clear();
        for (std::size_t start = 0; start < feature_rows[i].size();)
        {
          const std::size_t line_end = feature_rows[i].find('\n', start);
          appendSize(numbered, candidate_feature_row_++);
          numbered += '\t';
          numbered.append(feature_rows[i], start, line_end + 1 - start);
          start = line_end + 1;
        }
        *candidate_features_ << numbered;
      }
      const std::size_t candidate_index = close_order_[begin + i];
      if (traces_ != nullptr && dump_selected_[candidate_index])
      {
        dumpTraces_(candidate_index);
      }
    }
    close_cursor_ = end;
  }

  void StreamSearch::pickCandidate_(std::size_t candidate_index, std::vector<float>& matrix, std::vector<double>& scratch,
                                    std::string& peak_rows, std::string& feature_rows) const
  {
    const std::vector<Entry>& entries = entries_[candidate_index];
    if (entries.empty())
    {
      return;
    }
    // Every scan that isolated the precursor inside its RT window, zero filled where nothing matched.
    const std::vector<std::pair<std::uint32_t, double>> cells = candidateCells_(candidate_index);
    const std::size_t fragments = candidate_query_offsets_[candidate_index + 1] - candidate_query_offsets_[candidate_index];
    if (cells.empty() || fragments == 0)
    {
      return;
    }
    matrix.assign(cells.size() * fragments, 0.0f);
    std::vector<double> rts(cells.size());
    for (std::size_t s = 0; s < cells.size(); ++s)
    {
      rts[s] = cells[s].second;
    }
    for (const Entry& entry : entries)
    {
      const auto row = std::lower_bound(cells.begin(), cells.end(), entry.cell, [](const auto& cell, std::uint32_t id) { return cell.first < id; });
      if (row != cells.end() && row->first == entry.cell)
      {
        float& value = matrix[static_cast<std::size_t>(row - cells.begin()) * fragments + entry.slot];
        value = std::max(value, entry.intensity);
      }
    }
    std::size_t peaks_found = 0;
    const std::vector<XicPeak> peaks = pickXicPeaks(matrix.data(), cells.size(), fragments, rts.data(), options_.picker, peaks_found, scratch);
    const Candidate& candidate = candidates_[candidate_index];
    for (std::size_t rank = 0; rank < peaks.size(); ++rank)
    {
      const XicPeak& peak = peaks[rank];
      if (peak.support < options_.min_transitions)
      {
        continue;
      }
      const double start_rt = rts[peak.left], apex_rt = rts[peak.apex], end_rt = rts[peak.right];
      // Peak table, %.10g.
      peak_rows += candidate.id;
      for (const double value : {candidate.precursor_mz, start_rt, apex_rt, end_rt, peak.area > 0.0 ? peak.area : peak.apex_height})
      {
        peak_rows += '\t';
        appendDouble(peak_rows, value, 10);
      }
      peak_rows += '\t';
      appendSize(peak_rows, peak.support);
      peak_rows += '\t';
      appendSize(peak_rows, peak.signal_scans);
      peak_rows += '\n';
      if (candidate_features_ == nullptr)
      {
        continue;
      }
      // Feature table, %.17g, without the leading candidate_row.
      const auto number = [&feature_rows](double value) { feature_rows += '\t'; appendDouble(feature_rows, value, 17); };
      const auto count = [&feature_rows](std::size_t value) { feature_rows += '\t'; appendSize(feature_rows, value); };
      feature_rows += candidate.id;
      number(candidate.precursor_mz);
      count(candidate.is_decoy ? 1 : 0);
      number(start_rt);
      number(apex_rt);
      number(end_rt);
      number(end_rt - start_rt);
      number(peak.area);
      number(peak.apex_height);
      count(peak.support);
      count(peak.signal_scans);
      count(candidate.has_library_rt ? 1 : 0);
      number(candidate.has_library_rt ? apex_rt - candidate.library_rt : 0.0);
      count(candidate.has_precursor_im ? 1 : 0);
      number(peak.fwhm);
      number(peak.coelution);
      count(rank);
      count(peaks_found);
      feature_rows += '\n';
    }
  }

  std::vector<std::pair<std::uint32_t, double>> StreamSearch::candidateCells_(std::size_t candidate_index) const
  {
    const Candidate& candidate = candidates_[candidate_index];
    std::vector<std::pair<std::uint32_t, double>> cells;
    for (const std::uint32_t lane_index : candidate_lanes_[candidate_index])
    {
      const GridLane& lane = grid_lanes_[lane_index];
      auto first = lane.cell_rts.begin();
      auto last = lane.cell_rts.end();
      if (candidate.has_library_rt)
      {
        first = std::lower_bound(lane.cell_rts.begin(), lane.cell_rts.end(), candidate.library_rt - options_.rt_window_seconds);
        last = std::upper_bound(first, lane.cell_rts.end(), candidate.library_rt + options_.rt_window_seconds);
      }
      for (auto it = first; it != last; ++it)
      {
        cells.emplace_back(lane.cell_ids[static_cast<std::size_t>(it - lane.cell_rts.begin())], *it);
      }
    }
    std::sort(cells.begin(), cells.end());
    return cells;
  }

  void StreamSearch::setTraceDump(std::ostream& output, const std::vector<std::string>& precursor_ids)
  {
    std::map<std::string_view, bool> wanted;
    for (const std::string& id : precursor_ids)
    {
      wanted.emplace(id, true);
    }
    dump_selected_.assign(candidates_.size(), false);
    for (std::size_t c = 0; c < candidates_.size(); ++c)
    {
      dump_selected_[c] = wanted.contains(std::string_view(candidates_[c].id));
    }
    dump_mz_.assign(candidates_.size(), {});
    traces_ = &output;
    output << "precursor_id\tscan\trt\tslot\tfragment_mz\tintensity\tobserved_mz\n";
  }

  void StreamSearch::dumpTraces_(std::size_t candidate_index)
  {
    const std::vector<Entry>& entries = entries_[candidate_index];
    const std::vector<double>& observed = dump_mz_[candidate_index];
    const std::vector<std::pair<std::uint32_t, double>> cells = candidateCells_(candidate_index);
    const std::size_t offset = candidate_query_offsets_[candidate_index];
    const std::size_t fragments = candidate_query_offsets_[candidate_index + 1] - offset;
    // Same reduction as pickCandidate_: the most intense match of a (scan, fragment) wins.
    std::vector<float> intensity(cells.size() * fragments, 0.0f);
    std::vector<double> mz(cells.size() * fragments, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t e = 0; e < entries.size(); ++e)
    {
      const auto row = std::lower_bound(cells.begin(), cells.end(), entries[e].cell, [](const auto& cell, std::uint32_t id) { return cell.first < id; });
      if (row == cells.end() || row->first != entries[e].cell)
      {
        continue;
      }
      const std::size_t at = static_cast<std::size_t>(row - cells.begin()) * fragments + entries[e].slot;
      if (std::isnan(mz[at]) || entries[e].intensity > intensity[at])
      {
        intensity[at] = entries[e].intensity;
        mz[at] = e < observed.size() ? observed[e] : std::numeric_limits<double>::quiet_NaN();
      }
    }
    std::ostream& out = *traces_;
    const std::string& id = candidates_[candidate_index].id;
    for (std::size_t s = 0; s < cells.size(); ++s)
    {
      for (std::size_t f = 0; f < fragments; ++f)
      {
        const std::size_t at = s * fragments + f;
        out << id << '\t' << cells[s].first << '\t' << std::setprecision(10) << cells[s].second << '\t' << f << '\t'
            << queries_[candidate_query_indices_[offset + f]].fragment_mz << '\t' << intensity[at] << '\t';
        if (std::isnan(mz[at]))
        {
          out << "nan";
        }
        else
        {
          out << mz[at];
        }
        out << '\n';
      }
    }
    std::vector<Entry>().swap(entries_[candidate_index]);
    std::vector<double>().swap(dump_mz_[candidate_index]);
  }

  void StreamSearch::finish()
  {
    if (finished_)
    {
      return;
    }
    flushBatch_();
    // Every remaining candidate, including those without a library RT (open for the whole run).
    closeCandidates_(close_order_.size());
    finished_ = true;
  }

  bool StreamSearch::ionMobilityFilteringUsed() const
  {
    return ion_mobility_filtering_used_;
  }
}
