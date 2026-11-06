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
#include <map>
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
    query_intensity_.resize(queries_.size());
    query_generation_.resize(queries_.size());
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

    isolation_windows_.clear();
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
      isolation_windows_.emplace_back(center - lower_offset, center + upper_offset);
    }
    if (isolation_windows_.empty())
    {
      return;
    }

    // Merge overlapping isolation windows so a candidate is counted once per scan.
    if (isolation_windows_.size() > 1)
    {
      std::sort(isolation_windows_.begin(), isolation_windows_.end());
      std::size_t merged_count = 0;
      for (std::size_t i = 1; i < isolation_windows_.size(); ++i)
      {
        if (isolation_windows_[i].first <= isolation_windows_[merged_count].second)
        {
          isolation_windows_[merged_count].second =
            std::max(isolation_windows_[merged_count].second, isolation_windows_[i].second);
        }
        else
        {
          isolation_windows_[++merged_count] = isolation_windows_[i];
        }
      }
      isolation_windows_.resize(merged_count + 1);
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

    if (++generation_ == 0)
    {
      std::fill(query_generation_.begin(), query_generation_.end(), 0);
      std::fill(candidate_generation_.begin(), candidate_generation_.end(), 0);
      generation_ = 1;
    }
    query_touched_.clear();
    candidate_touched_.clear();

    // Visit only library fragments in the tolerance interval around each raw peak.
    for (std::size_t peak_index = 0; peak_index < spectrum.size(); ++peak_index)
    {
      const auto& peak = spectrum[peak_index];
      const double peak_mz = peak.getMZ();
      const double peak_intensity = peak.getIntensity();
      if (!std::isfinite(peak_mz) || !std::isfinite(peak_intensity))
      {
        continue;
      }
      const double low_mz = peak_mz - options_.fragment_tolerance_da;
      const double high_mz = peak_mz + options_.fragment_tolerance_da;
      auto first = std::lower_bound(queries_.begin(), queries_.end(), low_mz,
                                    [](const Query& query, double mz) { return query.fragment_mz < mz; });
      const auto last = std::upper_bound(first, queries_.end(), high_mz,
                                         [](double mz, const Query& query) { return mz < query.fragment_mz; });
      const double peak_im = im_values != nullptr ? static_cast<double>((*im_values)[peak_index]) : scalar_im;
      const bool has_peak_im = im_values != nullptr || has_scalar_im;
      for (auto it = first; it != last; ++it)
      {
        const std::size_t query_index = static_cast<std::size_t>(it - queries_.begin());
        if (!queryMatches_(*it, rt, peak_im, has_peak_im))
        {
          continue;
        }
        if (query_generation_[query_index] != generation_)
        {
          query_generation_[query_index] = generation_;
          query_intensity_[query_index] = peak_intensity;
          query_touched_.push_back(query_index);
        }
        else
        {
          query_intensity_[query_index] = std::max(query_intensity_[query_index], peak_intensity);
        }
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
      updateCandidate_(candidate_index, rt, candidate_scan_intensity_[candidate_index],
                       candidate_scan_support_[candidate_index]);
    }
  }

  bool StreamSearch::queryMatches_(const Query& query, double rt, double peak_im, bool has_peak_im) const
  {
    const Candidate& candidate = candidates_[query.candidate];
    const auto window = std::upper_bound(isolation_windows_.begin(), isolation_windows_.end(), candidate.precursor_mz,
                                         [](double mz, const auto& interval) { return mz < interval.first; });
    if (window == isolation_windows_.begin() || std::prev(window)->second < candidate.precursor_mz)
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
