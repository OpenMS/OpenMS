// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
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
    @brief Plans memory-bounded OpenSwathWorkflow work waves and exposes a small concurrency-limiter pair.

    The scheduler keeps the policy for memory-aware SWATH concurrency outside the extraction
    and scoring implementation. It estimates how many SWATH maps can be resident at the same
    time given the configured per-SWATH memory model, then groups the non-MS1 maps of an
    @c OpenSwath::SwathMap vector into waves that can be loaded as a unit before their inner
    scoring jobs are distributed across worker threads. The class is a namespace in disguise:
    all entry points are static, all nested types are aggregates, and the two nested classes
    (@ref ConcurrencyLimiter / @ref ScopedSlot) are a tiny standalone slot-limiting primitive
    for legacy streaming paths.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathWorkflowScheduler
  {
  public:
    /**
      @brief Memory and concurrency knobs used by wave planning and batch sizing.

      Default values target a multi-gigabyte memory budget on a typical OpenSwath workflow.
      All quantities are byte counts unless otherwise noted.
    */
    struct OPENMS_DLLAPI Options
    {
      /// OSW writer memory reserved outside the scheduler's budget. Default @c 2 GiB.
      UInt64 osw_buffer_bytes = 2ULL * 1024ULL * 1024ULL * 1024ULL;
      /// Fraction of the remaining free memory the scheduler may reserve. Internally clamped to @c [0.05, 0.95] by @ref estimateAvailableMemoryForScoring. Default @c 0.90.
      double memory_usage_fraction = 0.90;
      /// Estimated memory contribution of one cached spectrum. Default @c 600 KiB.
      UInt64 bytes_per_spectrum = 600ULL * 1024ULL;
      /// Fixed estimated overhead per loaded SWATH. Default @c 100 MiB.
      UInt64 per_swath_overhead_bytes = 100ULL * 1024ULL * 1024ULL;
      /// Average number of chromatograms extracted per non-MS1 SWATH. @c 0 disables transition-density accounting (the extra term is dropped from the per-SWATH byte estimate).
      Size avg_transitions_per_swath = 0;
      /// Estimated bytes retained per extracted chromatogram point. Default @c 64.
      UInt64 bytes_per_chromatogram_point = 64ULL;
      /// User override for concurrent SWATHs. Values @c > @c 0 cap the auto-estimate; @c <= @c 0 enables memory-based planning constrained only by @ref scoring_threads.
      int max_concurrent_swaths = -1;
      /// Number of scoring worker threads available. When @ref max_concurrent_swaths is @c <= @c 0, the auto SWATH concurrency will not exceed this.
      Size scoring_threads = 1;
      /// Smallest scoring batch size that may be picked automatically by @ref chooseInnerBatchSize. Default @c 2000.
      Size min_inner_batch_size = 2000;
      /// Largest scoring batch size that may be picked automatically by @ref chooseInnerBatchSize. Default @c 10000.
      Size max_inner_batch_size = 10000;
      /// Target number of queued scoring jobs per worker thread for the @ref chooseInnerBatchSize heuristic. Default @c 3.
      Size target_jobs_per_thread = 3;
    };

    /// Derived scheduler limits computed by @ref estimateConcurrency for one extraction/scoring run.
    struct OPENMS_DLLAPI ConcurrencyEstimate
    {
      Size non_ms1_swath_count = 0;        ///< Number of non-MS1 SWATH maps with a non-null data pointer found in the input vector.
      Size avg_spectra_per_swath = 0;      ///< Mean spectra per non-MS1 SWATH (integer division across @ref non_ms1_swath_count).
      Size max_concurrent_swaths = 1;      ///< Final cap on simultaneously-resident non-MS1 SWATH maps. Always at least @c 1.
      UInt64 memory_budget_bytes = 0;      ///< Memory budget returned by @ref estimateAvailableMemoryForScoring; @c 0 if the budget could not be determined (treated as unbounded by @ref planWaves).
      UInt64 estimated_bytes_per_swath = 0;///< Predicted memory contribution of one non-MS1 SWATH based on the @ref Options model.
    };

    /// A wave: a group of SWATH indices (into the original @c swath_maps vector) that may be resident at the same time. Populated by @ref planWaves.
    struct OPENMS_DLLAPI Wave
    {
      std::vector<Size> swath_indices; ///< Indices into the original @c swath_maps vector. Always non-MS1.
      UInt64 estimated_bytes = 0;      ///< Sum of the per-SWATH byte estimates (derived from @ref Options) for the indices above.
    };

    /**
      @brief Estimate memory available for scoring after subtracting OS, I/O, and OSW-writer reserves.

      Queries the OS via @c SysInfo::getFreeSystemMemory, subtracts a hardcoded @c 2 GiB
      OS reserve plus a @c 500 MiB I/O-buffer reserve plus @p options.osw_buffer_bytes, and
      returns the remainder multiplied by @p options.memory_usage_fraction (clamped to
      @c [0.05, 0.95]).

      Returns @c 0 when @c SysInfo::getFreeSystemMemory fails or when the reserves already
      exceed the reported free memory.
    */
    static UInt64 estimateAvailableMemoryForScoring(const Options& options);

    /**
      @brief Derive per-run concurrency limits from a SWATH-map vector and an @ref Options block.

      Counts non-MS1 SWATH maps that carry a non-null data pointer, computes the mean
      spectra count, and derives the per-SWATH byte estimate as
      @c avg_spectra @c * @c bytes_per_spectrum @c + @c per_swath_overhead @c + @c (transition-density term).
      The transition-density term is included only when @p options.avg_transitions_per_swath
      and @p options.bytes_per_chromatogram_point are both non-zero.

      The reported @ref ConcurrencyEstimate::max_concurrent_swaths is the minimum of:
        - the number of non-MS1 SWATHs,
        - the memory budget divided by the per-SWATH byte estimate (skipped when either side is zero),
        - @p options.max_concurrent_swaths when @c > @c 0, or @c max(1, scoring_threads) otherwise.
      The final value is clamped to @c [1, non_ms1_swath_count].

      When the input vector contains no non-MS1 SWATHs, the returned estimate is
      default-constructed (all zeros except @c max_concurrent_swaths which stays at @c 1).
    */
    static ConcurrencyEstimate estimateConcurrency(
      const std::vector<OpenSwath::SwathMap>& swath_maps,
      const Options& options);

    /**
      @brief Group non-MS1 SWATHs into waves that simultaneously honour the count and memory caps.

      Walks @p swath_maps in input order, skips MS1 maps, and assigns each remaining SWATH
      to the current wave. The wave is closed (and a new one started) when either of the
      following becomes true:
        - the wave already holds @p estimate.max_concurrent_swaths entries, or
        - adding this SWATH would push the running byte total past @p estimate.memory_budget_bytes
          (a budget of @c 0 is treated as unbounded).

      Per-SWATH byte cost uses the actual spectra count when @c sptr is non-null; otherwise
      falls back to @p estimate.avg_spectra_per_swath. Returns waves in the same order their
      contents appear in the input vector.
    */
    static std::vector<Wave> planWaves(
      const std::vector<OpenSwath::SwathMap>& swath_maps,
      const Options& options,
      const ConcurrencyEstimate& estimate);

    /**
      @brief Pick a scoring batch size that keeps the worker pool occupied.

      Behaviour:
        - @p total_compounds @c == @c 0 returns @c 0.
        - @p user_inner_batch_size @c > @c 0 short-circuits to @c min(user_inner_batch_size, total_compounds).
        - Otherwise the heuristic targets @c max(active_swaths, scoring_threads * target_jobs_per_thread)
          jobs in flight: @c batch_size @c = @c ceil(total_compounds / target_jobs), then
          clamped to @c [min_inner_batch_size, max_inner_batch_size] and finally capped at
          @p total_compounds.

      @p scoring_threads and @p active_swaths are each treated as @c max(1, value) inside
      the heuristic so callers may pass @c 0 without surprises.
    */
    static Size chooseInnerBatchSize(
      Size total_compounds,
      Size active_swaths,
      Size scoring_threads,
      int user_inner_batch_size,
      const Options& options);

    /**
      @brief Counting semaphore for paths that still stream one SWATH per worker.

      A small wrapper around a @c std::mutex + @c std::condition_variable pair that
      enforces an upper bound on concurrently-active slots. The constructor clamps the
      requested cap to @c max(1, max_concurrent_swaths) so callers don't accidentally
      configure an unreachable limit. Not copyable nor movable (inherited from the mutex
      and condition variable members).

      Pair with @ref ScopedSlot for RAII slot management instead of calling
      @ref acquireSlot / @ref releaseSlot directly.
    */
    class OPENMS_DLLAPI ConcurrencyLimiter
    {
    public:
      /// Construct with the upper bound on concurrent slots; values @c < @c 1 are clamped to @c 1.
      explicit ConcurrencyLimiter(Size max_concurrent_swaths);
      /// Block until a slot is free and reserve it. Wakes via the condition variable.
      void acquireSlot();
      /// Release a previously reserved slot and notify one waiter. No-op (safe) if the active count is already zero.
      void releaseSlot();

    private:
      Size max_concurrent_swaths_;
      Size active_swaths_ = 0;
      std::mutex mutex_;
      std::condition_variable cv_;
    };

    /**
      @brief RAII slot guard for @ref ConcurrencyLimiter — acquires on construction, releases on destruction.

      The constructor calls @ref ConcurrencyLimiter::acquireSlot on @p limiter (which may
      block); the destructor calls @ref ConcurrencyLimiter::releaseSlot. Passing @c nullptr
      makes both operations no-ops, so guard sites that may or may not be limited can be
      written without branching. Not copyable nor movable.
    */
    class OPENMS_DLLAPI ScopedSlot
    {
    public:
      /// Acquire a slot from @p limiter (blocking). If @p limiter is @c nullptr, the constructor is a no-op.
      explicit ScopedSlot(ConcurrencyLimiter* limiter);
      /// Release the slot back to the limiter (no-op if constructed with @c nullptr).
      ~ScopedSlot();

      ScopedSlot(const ScopedSlot&) = delete;
      ScopedSlot& operator=(const ScopedSlot&) = delete;

    private:
      ConcurrencyLimiter* limiter_;
    };
  };
}
