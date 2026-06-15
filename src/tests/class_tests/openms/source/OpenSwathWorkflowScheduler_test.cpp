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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
