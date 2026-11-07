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
#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <ostream>
#include <string_view>
#include <tuple>
#include <utility>

namespace OpenMS::OpenDIA
{
  StreamSearch::StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output,
                             std::ostream* candidate_features) :
    options_(options),
    output_(output),
    candidate_features_(candidate_features)
  {
    if (!std::isfinite(options_.fragment_tolerance_da) || !(options_.fragment_tolerance_da > 0.0) ||
        !std::isfinite(options_.ion_mobility_tolerance_vssc) || !(options_.ion_mobility_tolerance_vssc > 0.0) ||
        !std::isfinite(options_.rt_window_seconds) || options_.rt_window_seconds < 0.0 ||
        !std::isfinite(options_.peak_gap_seconds) || !(options_.peak_gap_seconds > 0.0) || options_.min_transitions == 0 ||
        (options_.split_local_peaks &&
         (!std::isfinite(options_.peak_valley_fraction) || !(options_.peak_valley_fraction > 0.0) ||
          !(options_.peak_valley_fraction < 1.0) || options_.min_peak_scans < 2)) ||
        options_.max_apexes_per_group == 0 || !std::isfinite(options_.min_apex_separation_seconds) ||
        options_.min_apex_separation_seconds < 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Search tolerances and minimum transition count must be positive.");
    }

    std::map<std::pair<std::string_view, double>, std::size_t> candidate_indices;
    candidates_.reserve(queries.size());
    queries_.reserve(queries.size());

    // The library is already grouped in memory by OpenMS. Flatten it once into
    // compact candidate/query arrays, then sort the query index by fragment m/z.
    for (auto& query : queries)
    {
      if (query.precursor_id.empty() || !std::isfinite(query.precursor_mz) || !(query.precursor_mz > 0.0) ||
          !std::isfinite(query.fragment_mz) || !(query.fragment_mz > 0.0))
      {
        continue;
      }

      const auto key = std::make_pair(query.precursor_id, query.precursor_mz);
      auto [candidate_it, inserted] = candidate_indices.try_emplace(key, candidates_.size());
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

    std::sort(queries_.begin(), queries_.end(), [](const Query& lhs, const Query& rhs)
    {
      return lhs.fragment_mz < rhs.fragment_mz;
    });
    // Keep each candidate's fragment-query range compact and m/z sorted.
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
      candidate_query_indices_[next_query[queries_[query_index].candidate]++] = query_index;
    }

    candidates_by_mz_.resize(candidates_.size());
    std::iota(candidates_by_mz_.begin(), candidates_by_mz_.end(), 0);
    std::sort(candidates_by_mz_.begin(), candidates_by_mz_.end(), [this](std::size_t lhs, std::size_t rhs)
    {
      return candidates_[lhs].precursor_mz < candidates_[rhs].precursor_mz;
    });

    batch_cells_.reserve(MAX_BATCH_CELLS);
    grid_points_.reserve(65536);
    candidate_scan_intensity_.resize(candidates_.size());
    candidate_scan_mobility_error_sum_.resize(candidates_.size());
    candidate_scan_support_.resize(candidates_.size());
    candidate_scan_mobility_match_count_.resize(candidates_.size());
    candidate_generation_.resize(candidates_.size());
    candidate_first_match_order_.resize(candidates_.size());
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

  void StreamSearch::consumeSpectrum_(const MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() != 2 || spectrum.getPrecursors().empty() || spectrum.empty())
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
      if (!std::isfinite(center) || !(center > 0.0) ||
          !std::isfinite(lower_offset) || lower_offset < 0.0 ||
          !std::isfinite(upper_offset) || upper_offset < 0.0 ||
          lower_offset + upper_offset <= 0.0)
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
          isolation_windows[merged_count].second =
            std::max(isolation_windows[merged_count].second, isolation_windows[i].second);
        }
        else
        {
          isolation_windows[++merged_count] = isolation_windows[i];
        }
      }
      isolation_windows.resize(merged_count + 1);
    }

    // Bruker TDF exposes one calibrated VSSC mobility value per peak. Other
    // compatible OpenMS readers may expose a scalar VSSC drift time per spectrum.
    // Unsupported/missing units leave the existing m/z-only path unchanged.
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

    // TDF spectra commonly repeat one mobility value for every centroid. Store
    // that once per cell; keep a sidecar only when peak mobilities vary.
    bool uniform_per_peak_im = im_values != nullptr && !im_values->empty() &&
                               std::isfinite((*im_values)[0]) && (*im_values)[0] > 0.0f;
    if (uniform_per_peak_im)
    {
      const float first_im = (*im_values)[0];
      uniform_per_peak_im = std::all_of(im_values->begin(), im_values->end(),
                                        [first_im](float value) { return value == first_im; });
      if (uniform_per_peak_im)
      {
        scalar_im = first_im;
        has_scalar_im = true;
      }
    }

    auto [lane_it, inserted] = lane_indices_.try_emplace(isolation_windows, grid_lanes_.size());
    const std::size_t lane_index = lane_it->second;
    if (inserted)
    {
      grid_lanes_.push_back({isolation_windows, {}});
      prepareLaneQueries_(grid_lanes_.back());
    }

    GridCell cell;
    cell.rt = rt;
    cell.scalar_im = scalar_im;
    cell.lane = lane_index;
    cell.has_per_peak_im = im_values != nullptr && !uniform_per_peak_im;
    cell.has_scalar_im = has_scalar_im;

    // Pack sorted cells into a bounded batch so candidate blocks can reuse one
    // worker team across several scans without retaining the run.
    if (!batch_cells_.empty() &&
        (batch_cells_.size() >= MAX_BATCH_CELLS || grid_points_.size() >= MAX_BATCH_POINTS))
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
    auto cell_points_begin = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_begin);
    auto cell_points_end = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_end);
    std::sort(cell_points_begin, cell_points_end,
              [](const GridPoint& lhs, const GridPoint& rhs) { return lhs.mz < rhs.mz; });
    batch_cells_.push_back(cell);
    if (batch_cells_.size() >= MAX_BATCH_CELLS || grid_points_.size() >= MAX_BATCH_POINTS)
    {
      flushBatch_();
    }
  }

  void StreamSearch::prepareLaneQueries_(GridLane& lane)
  {
    std::vector<std::size_t> lane_query_indices;
    for (const auto& [low_mz, high_mz] : lane.isolation_windows)
    {
      auto candidate = std::lower_bound(candidates_by_mz_.begin(), candidates_by_mz_.end(), low_mz,
                                        [this](std::size_t index, double mz)
                                        { return candidates_[index].precursor_mz < mz; });
      for (; candidate != candidates_by_mz_.end() && candidates_[*candidate].precursor_mz <= high_mz; ++candidate)
      {
        const std::size_t candidate_index = *candidate;
        lane_query_indices.insert(lane_query_indices.end(),
                                  candidate_query_indices_.begin() + static_cast<std::ptrdiff_t>(candidate_query_offsets_[candidate_index]),
                                  candidate_query_indices_.begin() + static_cast<std::ptrdiff_t>(candidate_query_offsets_[candidate_index + 1]));
      }
    }
    std::sort(lane_query_indices.begin(), lane_query_indices.end(), [this](std::size_t lhs, std::size_t rhs)
    {
      return queries_[lhs].fragment_mz < queries_[rhs].fragment_mz;
    });

    const std::size_t block_count = (candidates_.size() + CANDIDATES_PER_BLOCK - 1) / CANDIDATES_PER_BLOCK;
    lane.candidate_blocks.resize(block_count);
    for (std::size_t order = 0; order < lane_query_indices.size(); ++order)
    {
      const std::size_t query_index = lane_query_indices[order];
      const std::size_t block = queries_[query_index].candidate / CANDIDATES_PER_BLOCK;
      lane.candidate_blocks[block].push_back({query_index, order});
    }
  }

  void StreamSearch::flushBatch_()
  {
    if (batch_cells_.empty())
    {
      return;
    }

    std::vector<std::size_t> cell_generations(batch_cells_.size());
    const std::size_t block_count = (candidates_.size() + CANDIDATES_PER_BLOCK - 1) / CANDIDATES_PER_BLOCK;
    std::vector<bool> active_blocks(block_count, false);
    std::size_t active_block_count = 0;
    std::size_t batch_query_work = 0;
    for (std::size_t cell_order = 0; cell_order < batch_cells_.size(); ++cell_order)
    {
      if (++generation_ == 0)
      {
        std::fill(candidate_generation_.begin(), candidate_generation_.end(), 0);
        generation_ = 1;
      }
      cell_generations[cell_order] = generation_;

      const GridLane& lane = grid_lanes_[batch_cells_[cell_order].lane];
      for (std::size_t block = 0; block < lane.candidate_blocks.size(); ++block)
      {
        const std::vector<LaneQuery>& block_queries = lane.candidate_blocks[block];
        if (block_queries.empty())
        {
          continue;
        }
        if (!active_blocks[block])
        {
          active_blocks[block] = true;
          ++active_block_count;
        }
        if (batch_query_work <= std::numeric_limits<std::size_t>::max() - block_queries.size())
        {
          batch_query_work += block_queries.size();
        }
        else
        {
          batch_query_work = std::numeric_limits<std::size_t>::max();
        }
      }
    }

    std::vector<std::vector<CompletedCandidate>> completed_by_block(block_count);
    const auto processBlock = [this, &cell_generations, &completed_by_block](std::size_t block)
    {
      std::vector<std::size_t> touched;
      touched.reserve(std::min(CANDIDATES_PER_BLOCK, candidates_.size()));
      std::vector<CompletedCandidate>& completed = completed_by_block[block];

      for (std::size_t cell_order = 0; cell_order < batch_cells_.size(); ++cell_order)
      {
        const GridCell& cell = batch_cells_[cell_order];
        const GridLane& lane = grid_lanes_[cell.lane];
        if (block >= lane.candidate_blocks.size())
        {
          continue;
        }
        const std::vector<LaneQuery>& block_queries = lane.candidate_blocks[block];
        if (block_queries.empty())
        {
          continue;
        }

        touched.clear();
        const std::size_t cell_generation = cell_generations[cell_order];
        const auto point_end = grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_end);
        const double first_low_mz = queries_[block_queries.front().query_index].fragment_mz - options_.fragment_tolerance_da;
        auto point_cursor = std::lower_bound(grid_points_.begin() + static_cast<std::ptrdiff_t>(cell.point_begin),
                                             point_end, first_low_mz,
                                             [](const GridPoint& point, double mz) { return point.mz < mz; });
        for (const LaneQuery& lane_query : block_queries)
        {
          const std::size_t query_index = lane_query.query_index;
          const Query& query = queries_[query_index];
          const double low_mz = query.fragment_mz - options_.fragment_tolerance_da;
          const double high_mz = query.fragment_mz + options_.fragment_tolerance_da;
          while (point_cursor != point_end && point_cursor->mz < low_mz)
          {
            ++point_cursor;
          }

          double intensity = 0.0;
          double mobility_error = 0.0;
          bool has_mobility_evidence = false;
          bool matched = false;
          for (auto point = point_cursor; point != point_end && point->mz <= high_mz; ++point)
          {
            const bool has_peak_im = cell.has_per_peak_im || cell.has_scalar_im;
            const double peak_im = cell.has_per_peak_im ? point->ion_mobility : cell.scalar_im;
            if (!queryMatches_(query, cell.rt, peak_im, has_peak_im, lane.isolation_windows))
            {
              continue;
            }
            if (!matched || point->intensity > intensity)
            {
              intensity = point->intensity;
              const Candidate& candidate = candidates_[query.candidate];
              has_mobility_evidence = candidate.has_precursor_im && has_peak_im && std::isfinite(peak_im);
              mobility_error = has_mobility_evidence ? std::abs(peak_im - candidate.precursor_im) : 0.0;
            }
            matched = true;
          }
          if (!matched)
          {
            continue;
          }

          const std::size_t candidate_index = query.candidate;
          if (candidate_generation_[candidate_index] != cell_generation)
          {
            candidate_generation_[candidate_index] = cell_generation;
            candidate_scan_intensity_[candidate_index] = 0.0;
            candidate_scan_mobility_error_sum_[candidate_index] = 0.0;
            candidate_scan_support_[candidate_index] = 0;
            candidate_scan_mobility_match_count_[candidate_index] = 0;
            candidate_first_match_order_[candidate_index] = lane_query.order;
            touched.push_back(candidate_index);
          }
          candidate_scan_intensity_[candidate_index] += intensity;
          ++candidate_scan_support_[candidate_index];
          if (has_mobility_evidence)
          {
            candidate_scan_mobility_error_sum_[candidate_index] += mobility_error;
            ++candidate_scan_mobility_match_count_[candidate_index];
          }
        }

        for (const std::size_t candidate_index : touched)
        {
          updateCandidate_(candidate_index, cell.rt, candidate_scan_intensity_[candidate_index],
                           candidate_scan_support_[candidate_index],
                           candidate_scan_mobility_error_sum_[candidate_index],
                           candidate_scan_mobility_match_count_[candidate_index], cell_order,
                           candidate_first_match_order_[candidate_index], completed);
        }
      }
    };

#ifdef _OPENMP
    const bool use_parallel = active_block_count > 1 && batch_query_work >= 4096 &&
                              (batch_cells_.size() > 1 || batch_query_work >= 32768);
#pragma omp parallel for schedule(dynamic, 8) if(use_parallel)
    for (std::ptrdiff_t block = 0; block < static_cast<std::ptrdiff_t>(block_count); ++block)
    {
      if (active_blocks[static_cast<std::size_t>(block)])
      {
        processBlock(static_cast<std::size_t>(block));
      }
    }
#else
    for (std::size_t block = 0; block < block_count; ++block)
    {
      if (active_blocks[block])
      {
        processBlock(block);
      }
    }
#endif

    std::vector<CompletedCandidate> completed;
    for (auto& block_events : completed_by_block)
    {
      completed.insert(completed.end(), std::make_move_iterator(block_events.begin()),
                       std::make_move_iterator(block_events.end()));
    }
    std::sort(completed.begin(), completed.end(), [](const CompletedCandidate& lhs, const CompletedCandidate& rhs)
    {
      return std::tie(lhs.cell_order, lhs.query_order) < std::tie(rhs.cell_order, rhs.query_order);
    });
    for (const CompletedCandidate& candidate : completed)
    {
      writeCandidate_(candidate);
    }

    batch_cells_.clear();
    grid_points_.clear();
  }

  bool StreamSearch::queryMatches_(const Query& query, double rt, double peak_im, bool has_peak_im,
                                  const std::vector<std::pair<double, double>>& isolation_windows) const
  {
    const Candidate& candidate = candidates_[query.candidate];
    const auto window = std::upper_bound(isolation_windows.begin(), isolation_windows.end(), candidate.precursor_mz,
                                         [](double mz, const auto& interval) { return mz < interval.first; });
    if (window == isolation_windows.begin() || std::prev(window)->second < candidate.precursor_mz)
    {
      return false;
    }
    if (options_.rt_window_seconds != 0.0 && candidate.has_library_rt &&
        std::abs(rt - candidate.library_rt) > options_.rt_window_seconds)
    {
      return false;
    }
    // OpenSWATH PQP precursor mobility values use inverse reduced mobility
    // (VSSC). Apply the extra coordinate only when both sides provide it.
    return !candidate.has_precursor_im || !has_peak_im ||
           (std::isfinite(peak_im) && std::abs(peak_im - candidate.precursor_im) <= options_.ion_mobility_tolerance_vssc);
  }

  void StreamSearch::updateCandidate_(std::size_t candidate_index, double rt, double intensity,
                                      std::size_t transition_support, double mobility_error_sum,
                                      std::size_t mobility_match_count, std::size_t cell_order,
                                      std::size_t query_order, std::vector<CompletedCandidate>& completed)
  {
    Candidate& candidate = candidates_[candidate_index];
    if (candidate.active && (rt < candidate.last_rt || rt - candidate.last_rt > options_.peak_gap_seconds))
    {
      completed.push_back({candidate_index, cell_order, query_order, candidate.start_rt, candidate.apex_rt,
                           candidate.last_rt, candidate.area, candidate.apex_intensity, candidate.points,
                           candidate.max_transition_support, candidate.mobility_error_sum,
                           candidate.mobility_match_count, closeApexes_(candidate)});
      candidate.active = false;
      candidate.below_peak_scans = 0;
    }

    if (!candidate.active)
    {
      if (candidate.awaiting_peak && intensity < candidate.reopen_threshold)
      {
        return;
      }
      candidate.awaiting_peak = false;
      candidate.active = true;
      candidate.start_rt = rt;
      candidate.apex_rt = rt;
      candidate.last_rt = rt;
      candidate.last_intensity = intensity;
      candidate.apex_intensity = intensity;
      candidate.area = 0.0;
      candidate.points = 1;
      candidate.max_transition_support = transition_support;
      candidate.mobility_error_sum = mobility_error_sum;
      candidate.mobility_match_count = mobility_match_count;
      candidate.rising = true;
      candidate.apexes.clear();
      return;
    }

    if (options_.split_local_peaks && candidate.apex_intensity > 0.0 &&
        intensity < candidate.apex_intensity * options_.peak_valley_fraction)
    {
      ++candidate.below_peak_scans;
      if (candidate.below_peak_scans >= 2)
      {
        completed.push_back({candidate_index, cell_order, query_order, candidate.start_rt, candidate.apex_rt,
                             candidate.last_rt, candidate.area, candidate.apex_intensity, candidate.points,
                             candidate.max_transition_support, candidate.mobility_error_sum,
                             candidate.mobility_match_count, closeApexes_(candidate)});
        candidate.active = false;
        candidate.awaiting_peak = true;
        // ponytail: fixed fractional recovery keeps state O(1) per candidate.
        // Use adaptive prominence only if held-out traces show lost weak peaks.
        candidate.reopen_threshold = candidate.apex_intensity * options_.peak_valley_fraction;
        candidate.below_peak_scans = 0;
        return;
      }
    }
    else
    {
      candidate.below_peak_scans = 0;
    }

    // A rise followed by a drop marks the previous scan as a local apex.
    if (options_.max_apexes_per_group > 1)
    {
      if (intensity > candidate.last_intensity)
      {
        candidate.rising = true;
      }
      else if (intensity < candidate.last_intensity && candidate.rising)
      {
        offerApex_(candidate, candidate.last_rt, candidate.last_intensity);
        candidate.rising = false;
      }
    }

    candidate.area += 0.5 * (candidate.last_intensity + intensity) * (rt - candidate.last_rt);
    candidate.last_rt = rt;
    candidate.last_intensity = intensity;
    ++candidate.points;
    candidate.max_transition_support = std::max(candidate.max_transition_support, transition_support);
    candidate.mobility_error_sum += mobility_error_sum;
    candidate.mobility_match_count += mobility_match_count;
    if (intensity > candidate.apex_intensity)
    {
      candidate.apex_intensity = intensity;
      candidate.apex_rt = rt;
    }
  }

  void StreamSearch::offerApex_(Candidate& candidate, double rt, double intensity) const
  {
    // Streaming non-maximum suppression: a new apex replaces every kept apex
    // within the separation only if it is stronger than all of them.
    auto& apexes = candidate.apexes;
    bool stronger_than_neighbours = true;
    for (const auto& [kept_rt, kept_intensity] : apexes)
    {
      if (std::abs(kept_rt - rt) < options_.min_apex_separation_seconds && kept_intensity >= intensity)
      {
        stronger_than_neighbours = false;
        break;
      }
    }
    if (!stronger_than_neighbours)
    {
      return;
    }
    apexes.erase(std::remove_if(apexes.begin(), apexes.end(), [this, rt](const auto& apex)
                                { return std::abs(apex.first - rt) < options_.min_apex_separation_seconds; }),
                 apexes.end());
    apexes.emplace_back(rt, intensity);
    if (apexes.size() > options_.max_apexes_per_group)
    {
      // Drop the weakest; the earliest wins exact ties for determinism.
      auto weakest = std::min_element(apexes.begin(), apexes.end(), [](const auto& lhs, const auto& rhs)
      {
        return lhs.second != rhs.second ? lhs.second < rhs.second : lhs.first > rhs.first;
      });
      apexes.erase(weakest);
    }
  }

  std::vector<std::pair<double, double>> StreamSearch::closeApexes_(Candidate& candidate) const
  {
    if (options_.max_apexes_per_group <= 1)
    {
      return {};
    }
    if (candidate.rising)
    {
      offerApex_(candidate, candidate.last_rt, candidate.last_intensity);
    }
    std::vector<std::pair<double, double>> apexes = std::move(candidate.apexes);
    candidate.apexes.clear();
    // Strongest apex first, then RT, so the first row matches the global apex intensity.
    std::sort(apexes.begin(), apexes.end(), [](const auto& lhs, const auto& rhs)
    {
      return lhs.second != rhs.second ? lhs.second > rhs.second : lhs.first < rhs.first;
    });
    return apexes;
  }

  void StreamSearch::writeCandidateRows_(std::size_t candidate_index, double start_rt, double apex_rt, double end_rt,
                                         double area, double apex_intensity, std::size_t points,
                                         std::size_t transition_support, double mobility_error_sum,
                                         std::size_t mobility_match_count,
                                         const std::vector<std::pair<double, double>>& apexes)
  {
    if (transition_support < options_.min_transitions ||
        (options_.split_local_peaks && points < options_.min_peak_scans))
    {
      return;
    }
    const Candidate& candidate = candidates_[candidate_index];
    const double score = area > 0.0 ? area : apex_intensity;
    const auto writeRow = [&](double row_apex_rt, double row_apex_intensity)
    {
      output_ << candidate.id << '\t' << std::setprecision(10) << candidate.precursor_mz << '\t'
              << start_rt << '\t' << row_apex_rt << '\t' << end_rt << '\t'
              << score << '\t' << transition_support << '\t' << points << '\n';
      writeCandidateFeatures_(candidate_index, start_rt, row_apex_rt, end_rt, area, row_apex_intensity, points,
                              transition_support, mobility_error_sum, mobility_match_count);
    };
    if (apexes.empty())
    {
      writeRow(apex_rt, apex_intensity);
      return;
    }
    // Each local apex shares the group's boundaries, area and support evidence.
    for (const auto& [row_apex_rt, row_apex_intensity] : apexes)
    {
      writeRow(row_apex_rt, row_apex_intensity);
    }
  }

  void StreamSearch::writeCandidate_(std::size_t candidate_index, Candidate& candidate)
  {
    const auto apexes = closeApexes_(candidate);
    writeCandidateRows_(candidate_index, candidate.start_rt, candidate.apex_rt, candidate.last_rt, candidate.area,
                        candidate.apex_intensity, candidate.points, candidate.max_transition_support,
                        candidate.mobility_error_sum, candidate.mobility_match_count, apexes);
  }

  void StreamSearch::writeCandidate_(const CompletedCandidate& completed)
  {
    writeCandidateRows_(completed.candidate_index, completed.start_rt, completed.apex_rt, completed.last_rt,
                        completed.area, completed.apex_intensity, completed.points, completed.transition_support,
                        completed.mobility_error_sum, completed.mobility_match_count, completed.apexes);
  }

  void StreamSearch::writeCandidateFeatures_(std::size_t candidate_index, double start_rt, double apex_rt,
                                             double end_rt, double area, double apex_intensity,
                                             std::size_t points, std::size_t transition_support,
                                             double mobility_error_sum, std::size_t mobility_match_count)
  {
    if (candidate_features_ == nullptr)
    {
      return;
    }
    const Candidate& candidate = candidates_[candidate_index];
    const double apex_rt_error = candidate.has_library_rt ? apex_rt - candidate.library_rt : 0.0;
    *candidate_features_ << std::setprecision(17) << candidate_feature_row_++ << '\t' << candidate.id << '\t'
                         << candidate.precursor_mz << '\t' << (candidate.is_decoy ? 1 : 0) << '\t'
                         << start_rt << '\t' << apex_rt << '\t' << end_rt << '\t'
                         << end_rt - start_rt << '\t' << area << '\t' << apex_intensity << '\t'
                         << transition_support << '\t' << points << '\t'
                         << (candidate.has_library_rt ? 1 : 0) << '\t' << apex_rt_error << '\t'
                         << (candidate.has_precursor_im ? 1 : 0) << '\t' << (mobility_match_count > 0 ? 1 : 0) << '\t'
                         << (mobility_match_count > 0 ? mobility_error_sum / mobility_match_count : 0.0) << '\n';
  }

  void StreamSearch::finish()
  {
    if (finished_)
    {
      return;
    }
    flushBatch_();
    for (std::size_t candidate_index = 0; candidate_index < candidates_.size(); ++candidate_index)
    {
      Candidate& candidate = candidates_[candidate_index];
      if (candidate.active)
      {
        writeCandidate_(candidate_index, candidate);
      }
    }
    finished_ = true;
  }

  bool StreamSearch::ionMobilityFilteringUsed() const
  {
    return ion_mobility_filtering_used_;
  }
}
