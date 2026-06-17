// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflowScheduler.h>

#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <limits>
#include <utility>

namespace OpenMS
{
  namespace
  {
    UInt64 saturatingMultiply(UInt64 lhs, UInt64 rhs)
    {
      if (lhs != 0ULL && rhs > std::numeric_limits<UInt64>::max() / lhs)
      {
        return std::numeric_limits<UInt64>::max();
      }
      return lhs * rhs;
    }

    UInt64 estimateTransitionChromatogramBytes(
      UInt64 spectra_count,
      Size avg_transitions_per_swath,
      UInt64 bytes_per_chromatogram_point)
    {
      if (avg_transitions_per_swath == 0 || bytes_per_chromatogram_point == 0ULL)
      {
        return 0ULL;
      }

      const UInt64 transition_points = saturatingMultiply(
        spectra_count, static_cast<UInt64>(avg_transitions_per_swath));
      return saturatingMultiply(transition_points, bytes_per_chromatogram_point);
    }
  }

  UInt64 OpenSwathWorkflowScheduler::estimateAvailableMemoryForScoring(const Options& options)
  {
    size_t available_kb = 0;
    if (!SysInfo::getFreeSystemMemory(available_kb))
    {
      return 0ULL;
    }

    const UInt64 available_kb_64 = static_cast<UInt64>(available_kb);
    UInt64 available_bytes = 0;
    if (available_kb_64 > std::numeric_limits<UInt64>::max() / 1024ULL)
    {
      available_bytes = std::numeric_limits<UInt64>::max();
    }
    else
    {
      available_bytes = available_kb_64 * 1024ULL;
    }

    const UInt64 os_reserved = 2ULL * 1024ULL * 1024ULL * 1024ULL;
    const UInt64 io_buffer = 500ULL * 1024ULL * 1024ULL;
    const UInt64 reserved_bytes = os_reserved + io_buffer + options.osw_buffer_bytes;
    if (available_bytes <= reserved_bytes)
    {
      return 0ULL;
    }

    const double clamped_fraction = std::clamp(options.memory_usage_fraction, 0.05, 0.95);
    return static_cast<UInt64>((available_bytes - reserved_bytes) * clamped_fraction);
  }

  OpenSwathWorkflowScheduler::ConcurrencyEstimate OpenSwathWorkflowScheduler::estimateConcurrency(
    const std::vector<OpenSwath::SwathMap>& swath_maps,
    const Options& options)
  {
    ConcurrencyEstimate estimate;
    UInt64 total_spectra = 0;

    for (const OpenSwath::SwathMap& swath_map : swath_maps)
    {
      if (!swath_map.ms1 && swath_map.sptr)
      {
        total_spectra += swath_map.sptr->getNrSpectra();
        ++estimate.non_ms1_swath_count;
      }
    }

    if (estimate.non_ms1_swath_count == 0)
    {
      return estimate;
    }

    estimate.memory_budget_bytes = estimateAvailableMemoryForScoring(options);
    estimate.avg_spectra_per_swath =
      static_cast<Size>(total_spectra / estimate.non_ms1_swath_count);
    estimate.estimated_bytes_per_swath =
      static_cast<UInt64>(estimate.avg_spectra_per_swath) * options.bytes_per_spectrum +
      options.per_swath_overhead_bytes +
      estimateTransitionChromatogramBytes(
        static_cast<UInt64>(estimate.avg_spectra_per_swath),
        options.avg_transitions_per_swath,
        options.bytes_per_chromatogram_point);

    Size max_concurrent_swaths = estimate.non_ms1_swath_count;
    if (estimate.estimated_bytes_per_swath > 0 &&
        estimate.memory_budget_bytes > estimate.estimated_bytes_per_swath)
    {
      max_concurrent_swaths =
        static_cast<Size>(estimate.memory_budget_bytes / estimate.estimated_bytes_per_swath);
    }

    if (options.max_concurrent_swaths > 0)
    {
      max_concurrent_swaths = std::min<Size>(
        max_concurrent_swaths,
        static_cast<Size>(options.max_concurrent_swaths));
    }
    else
    {
      max_concurrent_swaths = std::min<Size>(
        max_concurrent_swaths,
        std::max<Size>(1, options.scoring_threads));
    }

    estimate.max_concurrent_swaths = std::max<Size>(
      1,
      std::min<Size>(estimate.non_ms1_swath_count, max_concurrent_swaths));
    return estimate;
  }

  std::vector<OpenSwathWorkflowScheduler::Wave> OpenSwathWorkflowScheduler::planWaves(
    const std::vector<OpenSwath::SwathMap>& swath_maps,
    const Options& options,
    const ConcurrencyEstimate& estimate)
  {
    std::vector<Wave> waves;
    Wave current_wave;

    const Size max_concurrent_swaths = std::max<Size>(1, estimate.max_concurrent_swaths);
    const UInt64 memory_budget_bytes = estimate.memory_budget_bytes > 0 ?
      estimate.memory_budget_bytes :
      std::numeric_limits<UInt64>::max();

    for (Size swath_index = 0; swath_index < swath_maps.size(); ++swath_index)
    {
      const OpenSwath::SwathMap& swath_map = swath_maps[swath_index];
      if (swath_map.ms1)
      {
        continue;
      }

      const UInt64 spectra_count = swath_map.sptr ? swath_map.sptr->getNrSpectra() : estimate.avg_spectra_per_swath;
      const UInt64 swath_bytes =
        spectra_count * options.bytes_per_spectrum +
        options.per_swath_overhead_bytes +
        estimateTransitionChromatogramBytes(
          spectra_count,
          options.avg_transitions_per_swath,
          options.bytes_per_chromatogram_point);

      const bool wave_full_by_count = current_wave.swath_indices.size() >= max_concurrent_swaths;
      const bool wave_full_by_memory =
        !current_wave.swath_indices.empty() &&
        current_wave.estimated_bytes + swath_bytes > memory_budget_bytes;
      if (wave_full_by_count || wave_full_by_memory)
      {
        waves.push_back(std::move(current_wave));
        current_wave = Wave();
      }

      current_wave.swath_indices.push_back(swath_index);
      current_wave.estimated_bytes += swath_bytes;
    }

    if (!current_wave.swath_indices.empty())
    {
      waves.push_back(std::move(current_wave));
    }

    return waves;
  }

  Size OpenSwathWorkflowScheduler::chooseInnerBatchSize(
    Size total_compounds,
    Size active_swaths,
    Size scoring_threads,
    int user_inner_batch_size,
    const Options& options)
  {
    if (total_compounds == 0)
    {
      return 0;
    }

    if (user_inner_batch_size > 0)
    {
      return std::min<Size>(static_cast<Size>(user_inner_batch_size), total_compounds);
    }

    const Size workers = std::max<Size>(1, scoring_threads);
    const Size active_jobs = std::max<Size>(1, active_swaths);
    const Size jobs_per_thread = std::max<Size>(1, options.target_jobs_per_thread);
    const Size target_jobs = std::max<Size>(active_jobs, workers * jobs_per_thread);
    Size batch_size = (total_compounds + target_jobs - 1) / target_jobs;

    const Size min_batch = std::max<Size>(1, options.min_inner_batch_size);
    const Size max_batch = std::max<Size>(min_batch, options.max_inner_batch_size);
    batch_size = std::clamp(batch_size, min_batch, max_batch);

    return std::min<Size>(batch_size, total_compounds);
  }

  OpenSwathWorkflowScheduler::ConcurrencyLimiter::ConcurrencyLimiter(Size max_concurrent_swaths) :
    max_concurrent_swaths_(std::max<Size>(1, max_concurrent_swaths))
  {
  }

  void OpenSwathWorkflowScheduler::ConcurrencyLimiter::acquireSlot()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [this]() { return active_swaths_ < max_concurrent_swaths_; });
    ++active_swaths_;
  }

  void OpenSwathWorkflowScheduler::ConcurrencyLimiter::releaseSlot()
  {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (active_swaths_ > 0)
      {
        --active_swaths_;
      }
    }
    cv_.notify_one();
  }

  OpenSwathWorkflowScheduler::ScopedSlot::ScopedSlot(ConcurrencyLimiter* limiter) :
    limiter_(limiter)
  {
    if (limiter_ != nullptr)
    {
      limiter_->acquireSlot();
    }
  }

  OpenSwathWorkflowScheduler::ScopedSlot::~ScopedSlot()
  {
    if (limiter_ != nullptr)
    {
      limiter_->releaseSlot();
    }
  }
}
