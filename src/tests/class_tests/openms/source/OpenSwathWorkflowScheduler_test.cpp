// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflowScheduler.h>
///////////////////////////

#include <algorithm>
#include <vector>

using namespace OpenMS;
using namespace std;

using Scheduler = OpenSwathWorkflowScheduler;

namespace
{
  // SWATH map vector: one leading MS1 map, then `n_ms2` non-MS1 maps, plus a second
  // MS1 map spliced into the middle (so MS1 exclusion is exercised mid-stream too).
  std::vector<OpenSwath::SwathMap> makeMaps(Size n_ms2)
  {
    std::vector<OpenSwath::SwathMap> maps;
    maps.emplace_back(0.0, 0.0, 0.0, /*is_ms1=*/true);   // index 0: MS1
    double lo = 400.0;
    for (Size i = 0; i < n_ms2; ++i)
    {
      if (i == n_ms2 / 2)
      {
        maps.emplace_back(0.0, 0.0, 0.0, /*is_ms1=*/true); // an MS1 map spliced into the middle
      }
      maps.emplace_back(lo, lo + 25.0, lo + 12.5, /*is_ms1=*/false);
      lo += 25.0;
    }
    return maps;
  }

  // The non-MS1 indices in input order — what a valid wave partition must cover exactly.
  std::vector<Size> nonMs1Indices(const std::vector<OpenSwath::SwathMap>& maps)
  {
    std::vector<Size> out;
    for (Size i = 0; i < maps.size(); ++i) { if (!maps[i].ms1) { out.push_back(i); } }
    return out;
  }
}

START_TEST(OpenSwathWorkflowScheduler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static std::vector<Wave> planWaves(const std::vector<OpenSwath::SwathMap>& swath_maps, const Options& options, const ConcurrencyEstimate& estimate)))
{
  Scheduler::Options opts; // defaults: avg_transitions_per_swath = 0 -> no transition term

  // Reusable check: planWaves must produce a valid partition of the non-MS1 maps --
  // every non-MS1 index appears exactly once, in input order, with no MS1 index, and
  // no wave exceeding the concurrency cap.
  auto check_partition = [&opts](const std::vector<OpenSwath::SwathMap>& maps,
                                 const Scheduler::ConcurrencyEstimate& est)
  {
    std::vector<Scheduler::Wave> waves = Scheduler::planWaves(maps, opts, est);
    const std::vector<Size> expected = nonMs1Indices(maps);
    const Size cap = std::max<Size>(1, est.max_concurrent_swaths);

    std::vector<Size> collected;
    for (const Scheduler::Wave& w : waves)
    {
      TEST_EQUAL(w.swath_indices.size() <= cap, true)   // count cap honoured
      TEST_EQUAL(w.swath_indices.empty(), false)        // no empty waves emitted
      for (Size idx : w.swath_indices)
      {
        TEST_EQUAL(maps[idx].ms1, false)                // MS1 maps never appear in a wave
        collected.push_back(idx);
      }
    }
    // covers every non-MS1 map exactly once, in input order
    TEST_EQUAL(collected.size(), expected.size())
    TEST_EQUAL(collected == expected, true)
    return waves;
  };

  // 7 non-MS1 maps (+ 2 MS1 maps). memory_budget = 0 -> unbounded, so only the count
  // cap closes waves.
  std::vector<OpenSwath::SwathMap> maps = makeMaps(7);
  TEST_EQUAL(nonMs1Indices(maps).size(), 7)

  {
    // max_concurrent = 1 -> one map per wave -> 7 waves
    Scheduler::ConcurrencyEstimate est;
    est.max_concurrent_swaths = 1;
    est.memory_budget_bytes = 0;
    est.avg_spectra_per_swath = 0;
    TEST_EQUAL(check_partition(maps, est).size(), 7)
  }
  {
    // max_concurrent = 3 -> waves of 3,3,1 -> 3 waves
    Scheduler::ConcurrencyEstimate est;
    est.max_concurrent_swaths = 3;
    est.memory_budget_bytes = 0;
    est.avg_spectra_per_swath = 0;
    TEST_EQUAL(check_partition(maps, est).size(), 3)
  }
  {
    // max_concurrent larger than the SWATH count -> a single wave with all 7
    Scheduler::ConcurrencyEstimate est;
    est.max_concurrent_swaths = 100;
    est.memory_budget_bytes = 0;
    est.avg_spectra_per_swath = 0;
    std::vector<Scheduler::Wave> waves = check_partition(maps, est);
    TEST_EQUAL(waves.size(), 1)
    TEST_EQUAL(waves[0].swath_indices.size(), 7)
  }
  {
    // Memory-bounded: with avg_spectra=0 each SWATH costs per_swath_overhead_bytes
    // (100 MiB by default). A 250 MiB budget fits 2 SWATHs per wave (300 MiB would
    // overflow), so 7 maps -> 4 waves -- and the partition must still be complete.
    Scheduler::ConcurrencyEstimate est;
    est.max_concurrent_swaths = 100;                      // count does not bind here
    est.memory_budget_bytes = 250ULL * 1024ULL * 1024ULL; // 250 MiB
    est.avg_spectra_per_swath = 0;
    std::vector<Scheduler::Wave> waves = check_partition(maps, est);
    TEST_EQUAL(waves.size(), 4)
    TEST_EQUAL(waves[0].swath_indices.size(), 2)
    TEST_EQUAL(waves[3].swath_indices.size(), 1)
  }
  {
    // No non-MS1 maps -> no waves.
    std::vector<OpenSwath::SwathMap> ms1_only;
    ms1_only.emplace_back(0.0, 0.0, 0.0, true);
    ms1_only.emplace_back(0.0, 0.0, 0.0, true);
    Scheduler::ConcurrencyEstimate est;
    est.max_concurrent_swaths = 4;
    TEST_EQUAL(Scheduler::planWaves(ms1_only, opts, est).size(), 0)

    // Empty input -> no waves.
    TEST_EQUAL(Scheduler::planWaves(std::vector<OpenSwath::SwathMap>{}, opts, est).size(), 0)
  }
}
END_SECTION

START_SECTION((static Size chooseInnerBatchSize(Size total_compounds, Size active_swaths, Size scoring_threads, int user_inner_batch_size, const Options& options)))
{
  Scheduler::Options opts; // defaults: target_jobs_per_thread=3, min=2000, max=10000

  // total_compounds == 0 -> 0 (nothing to batch)
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(0, 1, 4, -1, opts), 0)

  // explicit user batch size short-circuits to min(user, total)
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(10000, 1, 4, 500, opts), 500)
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(1000, 1, 4, 50000, opts), 1000)   // user capped at total

  // auto (user <= 0): target_jobs = max(active, threads*jobs_per_thread); batch =
  // ceil(total/target_jobs), clamped to [min,max], then capped at total.
  // total=60000, threads=4 -> target_jobs=12 -> ceil=5000 (within [2000,10000])
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(60000, 1, 4, -1, opts), 5000)
  // total=1000000 -> ceil=83334 -> clamped down to max 10000
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(1000000, 1, 4, -1, opts), 10000)
  // total=1000 -> ceil=84 -> clamped up to min 2000 -> then capped at total 1000
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(1000, 1, 4, -1, opts), 1000)
  // active_swaths dominates: target_jobs = max(50, 12) = 50 -> ceil(200000/50)=4000
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(200000, 50, 4, -1, opts), 4000)
  // zero threads / zero active are treated as 1: target_jobs = max(1, 1*3)=3 ->
  // ceil(30000/3)=10000
  TEST_EQUAL(Scheduler::chooseInnerBatchSize(30000, 0, 0, -1, opts), 10000)
}
END_SECTION

START_SECTION([EXTRA] ConcurrencyLimiter / ScopedSlot acquire and release within capacity)
{
  // A null-limiter ScopedSlot is a no-op (lets guard sites be written without branching).
  {
    Scheduler::ScopedSlot s(nullptr);
  }

  // Acquire/release through ScopedSlot RAII, never exceeding capacity, must not
  // deadlock. (Over-acquiring a full limiter blocks by design, so the sequence below
  // holds at most 2 slots at a time against a capacity-2 limiter.)
  Scheduler::ConcurrencyLimiter limiter(2);
  bool completed = false;
  {
    Scheduler::ScopedSlot s1(&limiter); // active = 1
    {
      Scheduler::ScopedSlot s2(&limiter); // active = 2 (full)
    }                                     // s2 releases -> active = 1
    Scheduler::ScopedSlot s3(&limiter);   // re-fills the freed slot -> active = 2
  }                                       // s1, s3 release -> active = 0
  completed = true;
  TEST_TRUE(completed) // reached here without blocking: acquire/release works

  // A limiter requested with capacity 0 is clamped to 1 internally; a single slot
  // must still be acquirable and releasable.
  Scheduler::ConcurrencyLimiter clamped(0);
  bool clamped_ok = false;
  {
    Scheduler::ScopedSlot s(&clamped);
  }
  clamped_ok = true;
  TEST_TRUE(clamped_ok)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
