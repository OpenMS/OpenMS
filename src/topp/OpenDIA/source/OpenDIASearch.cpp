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
#include <utility>

namespace OpenMS::OpenDIA
{
  StreamSearch::StreamSearch(std::vector<TransitionQuery> queries, SearchOptions options, std::ostream& output) :
    options_(options),
    output_(output)
  {
    if (!std::isfinite(options_.fragment_tolerance_da) || !(options_.fragment_tolerance_da > 0.0) ||
        !std::isfinite(options_.ion_mobility_tolerance_vssc) || !(options_.ion_mobility_tolerance_vssc > 0.0) ||
        !std::isfinite(options_.rt_window_seconds) || options_.rt_window_seconds < 0.0 ||
        !std::isfinite(options_.peak_gap_seconds) || !(options_.peak_gap_seconds > 0.0) || options_.min_transitions == 0)
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
        has_library_im_candidates_ = has_library_im_candidates_ || candidate.has_precursor_im;
        candidates_.push_back(std::move(candidate));
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

    query_intensity_.resize(queries_.size());
    candidate_scan_intensity_.resize(candidates_.size());
    candidate_scan_support_.resize(candidates_.size());
    candidate_generation_.resize(candidates_.size());
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
      if (im_unit == DriftTimeUnit::VSSC && im_index < arrays.size() && arrays[im_index].size() == spectrum.size())
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

    // Sort only this spectrum's cell. Doubles preserve legacy m/z boundaries
    // and intensity integration while the reusable grid scratch stays bounded.
    grid_points_.clear();
    grid_points_.reserve(spectrum.size());
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
    std::sort(grid_points_.begin(), grid_points_.end(),
              [](const GridPoint& lhs, const GridPoint& rhs) { return lhs.mz < rhs.mz; });
    cell.point_begin = 0;
    cell.point_end = grid_points_.size();
    extractCell_(cell, grid_lanes_[lane_index]);
  }

  void StreamSearch::prepareLaneQueries_(GridLane& lane)
  {
    for (const auto& [low_mz, high_mz] : lane.isolation_windows)
    {
      auto candidate = std::lower_bound(candidates_by_mz_.begin(), candidates_by_mz_.end(), low_mz,
                                        [this](std::size_t index, double mz)
                                        { return candidates_[index].precursor_mz < mz; });
      for (; candidate != candidates_by_mz_.end() && candidates_[*candidate].precursor_mz <= high_mz; ++candidate)
      {
        const std::size_t candidate_index = *candidate;
        lane.query_indices.insert(lane.query_indices.end(),
                                  candidate_query_indices_.begin() + static_cast<std::ptrdiff_t>(candidate_query_offsets_[candidate_index]),
                                  candidate_query_indices_.begin() + static_cast<std::ptrdiff_t>(candidate_query_offsets_[candidate_index + 1]));
      }
    }
    std::sort(lane.query_indices.begin(), lane.query_indices.end(), [this](std::size_t lhs, std::size_t rhs)
    {
      return queries_[lhs].fragment_mz < queries_[rhs].fragment_mz;
    });
  }

  void StreamSearch::extractCell_(const GridCell& cell, const GridLane& lane)
  {
    if (++generation_ == 0)
    {
      std::fill(candidate_generation_.begin(), candidate_generation_.end(), 0);
      generation_ = 1;
    }
    query_touched_.clear();
    candidate_touched_.clear();

    std::size_t point_cursor = cell.point_begin;
    for (const std::size_t query_index : lane.query_indices)
    {
      const Query& query = queries_[query_index];
      const double low_mz = query.fragment_mz - options_.fragment_tolerance_da;
      const double high_mz = query.fragment_mz + options_.fragment_tolerance_da;
      while (point_cursor < cell.point_end && grid_points_[point_cursor].mz < low_mz)
      {
        ++point_cursor;
      }

      double intensity = 0.0;
      bool matched = false;
      for (std::size_t point = point_cursor; point < cell.point_end && grid_points_[point].mz <= high_mz; ++point)
      {
        const GridPoint& raw_point = grid_points_[point];
        const bool has_peak_im = cell.has_per_peak_im || cell.has_scalar_im;
        const double peak_im = cell.has_per_peak_im ? raw_point.ion_mobility : cell.scalar_im;
        if (!queryMatches_(query, cell.rt, peak_im, has_peak_im, lane.isolation_windows))
        {
          continue;
        }
        const double peak_intensity = raw_point.intensity;
        intensity = matched ? std::max(intensity, peak_intensity) : peak_intensity;
        matched = true;
      }
      if (matched)
      {
        query_intensity_[query_index] = intensity;
        query_touched_.push_back(query_index);
      }
    }

    // Sum each transition once per scan, regardless of how many centroids fall in tolerance.
    for (const std::size_t query_index : query_touched_)
    {
      const std::size_t candidate_index = queries_[query_index].candidate;
      if (candidate_generation_[candidate_index] != generation_)
      {
        candidate_generation_[candidate_index] = generation_;
        candidate_scan_intensity_[candidate_index] = 0.0;
        candidate_scan_support_[candidate_index] = 0;
        candidate_touched_.push_back(candidate_index);
      }
      candidate_scan_intensity_[candidate_index] += query_intensity_[query_index];
      ++candidate_scan_support_[candidate_index];
    }
    for (const std::size_t candidate_index : candidate_touched_)
    {
      updateCandidate_(candidate_index, cell.rt, candidate_scan_intensity_[candidate_index],
                       candidate_scan_support_[candidate_index]);
    }
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
                                      std::size_t transition_support)
  {
    Candidate& candidate = candidates_[candidate_index];
    if (candidate.active && (rt < candidate.last_rt || rt - candidate.last_rt > options_.peak_gap_seconds))
    {
      writeCandidate_(candidate);
      candidate.active = false;
    }

    if (!candidate.active)
    {
      candidate.active = true;
      candidate.start_rt = rt;
      candidate.apex_rt = rt;
      candidate.last_rt = rt;
      candidate.last_intensity = intensity;
      candidate.apex_intensity = intensity;
      candidate.area = 0.0;
      candidate.points = 1;
      candidate.max_transition_support = transition_support;
      return;
    }

    candidate.area += 0.5 * (candidate.last_intensity + intensity) * (rt - candidate.last_rt);
    candidate.last_rt = rt;
    candidate.last_intensity = intensity;
    ++candidate.points;
    candidate.max_transition_support = std::max(candidate.max_transition_support, transition_support);
    if (intensity > candidate.apex_intensity)
    {
      candidate.apex_intensity = intensity;
      candidate.apex_rt = rt;
    }
  }

  void StreamSearch::writeCandidate_(const Candidate& candidate)
  {
    if (candidate.max_transition_support < options_.min_transitions)
    {
      return;
    }
    const double score = candidate.area > 0.0 ? candidate.area : candidate.apex_intensity;
    output_ << candidate.id << '\t' << std::setprecision(10) << candidate.precursor_mz << '\t'
            << candidate.start_rt << '\t' << candidate.apex_rt << '\t' << candidate.last_rt << '\t'
            << score << '\t' << candidate.max_transition_support << '\t' << candidate.points << '\n';
  }

  void StreamSearch::finish()
  {
    if (finished_)
    {
      return;
    }
    for (const Candidate& candidate : candidates_)
    {
      if (candidate.active)
      {
        writeCandidate_(candidate);
      }
    }
    finished_ = true;
  }

  bool StreamSearch::ionMobilityFilteringUsed() const
  {
    return ion_mobility_filtering_used_;
  }
}
