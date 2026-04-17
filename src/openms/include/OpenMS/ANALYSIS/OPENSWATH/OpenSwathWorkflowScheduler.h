// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

#include <condition_variable>
#include <mutex>
#include <vector>

namespace OpenMS
{
  /**
    @brief Plans bounded OpenSwathWorkflow work waves.

    The scheduler keeps the policy for memory-aware SWATH concurrency outside
    the extraction and scoring implementation. It estimates how many SWATH maps
    can be resident at once, then groups non-MS1 maps into waves that can be
    extracted before their inner scoring jobs are distributed across worker
    threads.
  */
  class OPENMS_DLLAPI OpenSwathWorkflowScheduler
  {
  public:
    /// Memory and concurrency options used for wave planning.
    struct OPENMS_DLLAPI Options
    {
      /// OSW writer memory reserved outside the scheduler budget.
      UInt64 osw_buffer_bytes = 2ULL * 1024ULL * 1024ULL * 1024ULL;
      /// Fraction of remaining free memory the scheduler may reserve.
      double memory_usage_fraction = 0.90;
      /// Estimated memory contribution of one cached spectrum.
      UInt64 bytes_per_spectrum = 600ULL * 1024ULL;
      /// Fixed estimated overhead per loaded SWATH.
      UInt64 per_swath_overhead_bytes = 100ULL * 1024ULL * 1024ULL;
      /// User override for concurrent SWATHs; <= 0 enables memory-based planning.
      int max_concurrent_swaths = -1;
    };

    /// Derived scheduler limits for one extraction/scoring run.
    struct OPENMS_DLLAPI ConcurrencyEstimate
    {
      Size non_ms1_swath_count = 0;
      Size avg_spectra_per_swath = 0;
      Size max_concurrent_swaths = 1;
      UInt64 memory_budget_bytes = 0;
      UInt64 estimated_bytes_per_swath = 0;
    };

    /// A group of SWATH indices that may be resident at the same time.
    struct OPENMS_DLLAPI Wave
    {
      std::vector<Size> swath_indices;
      UInt64 estimated_bytes = 0;
    };

    /// Estimate memory available for scoring after OS and writer reserves.
    static UInt64 estimateAvailableMemoryForScoring(const Options& options);

    /// Estimate a safe number of concurrently resident non-MS1 SWATH maps.
    static ConcurrencyEstimate estimateConcurrency(
      const std::vector<OpenSwath::SwathMap>& swath_maps,
      const Options& options);

    /// Plan memory-bounded waves over all non-MS1 SWATH maps.
    static std::vector<Wave> planWaves(
      const std::vector<OpenSwath::SwathMap>& swath_maps,
      const Options& options,
      const ConcurrencyEstimate& estimate);

    /// Simple blocking limiter for paths that still stream one SWATH per worker.
    class OPENMS_DLLAPI ConcurrencyLimiter
    {
    public:
      explicit ConcurrencyLimiter(Size max_concurrent_swaths);
      void acquireSlot();
      void releaseSlot();

    private:
      Size max_concurrent_swaths_;
      Size active_swaths_ = 0;
      std::mutex mutex_;
      std::condition_variable cv_;
    };

    /// RAII wrapper for ConcurrencyLimiter slots.
    class OPENMS_DLLAPI ScopedSlot
    {
    public:
      explicit ScopedSlot(ConcurrencyLimiter* limiter);
      ~ScopedSlot();

      ScopedSlot(const ScopedSlot&) = delete;
      ScopedSlot& operator=(const ScopedSlot&) = delete;

    private:
      ConcurrencyLimiter* limiter_;
    };
  };
}
