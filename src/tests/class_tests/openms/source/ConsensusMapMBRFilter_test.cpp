// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapMBRFilter.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  /// Build a Feature with a unique id, RT, m/z, intensity, optional ID (sequence + charge).
  Feature makeFeature(UInt64 unique_id, double rt, double mz, double intensity,
                      const String& sequence = "", Int charge = 0)
  {
    Feature f;
    f.setUniqueId(unique_id);
    f.setRT(rt);
    f.setMZ(mz);
    f.setIntensity(intensity);
    f.setCharge(charge);
    if (!sequence.empty())
    {
      PeptideIdentification pid;
      pid.setScoreType("score");
      pid.setHigherScoreBetter(true);
      PeptideHit hit;
      hit.setSequence(AASequence::fromString(sequence));
      hit.setCharge(charge);
      hit.setScore(1.0);
      pid.setHits({hit});
      f.setPeptideIdentifications({pid});
    }
    return f;
  }

  /// Make a 4-map column header set on a ConsensusMap.
  void setupHeaders(ConsensusMap& cmap, Size n_maps)
  {
    auto& headers = cmap.getColumnHeaders();
    for (Size i = 0; i < n_maps; ++i)
    {
      ConsensusMap::ColumnHeader h;
      h.filename = "map_" + String(i) + ".mzML";
      h.label = "label_" + String(i);
      h.size = 0;
      h.unique_id = i + 1;
      headers[i] = h;
    }
  }

  /// Build a ConsensusFeature from features in maps. Each entry is (map_index, Feature).
  /// Uses the ConsensusFeature::insert(map_index, BaseFeature) path so map_index gets
  /// stamped onto every PID exactly as the linker does.
  ConsensusFeature makeCF(const std::vector<std::pair<UInt64, Feature>>& items)
  {
    ConsensusFeature cf;
    cf.setUniqueId(items.front().second.getUniqueId() * 1000ULL + items.size());
    for (const auto& kv : items)
    {
      cf.insert(kv.first, kv.second);
    }
    cf.computeConsensus();
    return cf;
  }
}

START_TEST(ConsensusMapMBRFilter, "$Id$")

ConsensusMapMBRFilter* ptr = nullptr;
ConsensusMapMBRFilter* nullPointer = nullptr;
START_SECTION((ConsensusMapMBRFilter()))
{
  ptr = new ConsensusMapMBRFilter();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~ConsensusMapMBRFilter()))
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] CF with all handles ID-bearing is left untouched)
{
  ConsensusMap cmap;
  setupHeaders(cmap, 2);
  ConsensusFeature cf = makeCF({
    {0, makeFeature(1, 100.0, 500.0, 1e5, "PEPTIDEK", 2)},
    {1, makeFeature(2, 200.0, 500.0, 1e5, "PEPTIDEK", 2)} // both have IDs => no transfer to filter
  });
  cmap.push_back(cf);

  ConsensusMapMBRFilter f;
  auto rep = f.filter(cmap);
  TEST_EQUAL(rep.handles_removed, 0)
  TEST_EQUAL(cmap[0].size(), 2)
}
END_SECTION

START_SECTION((Report filter(ConsensusMap& cmap)))
{
  ConsensusMap cmap;
  setupHeaders(cmap, 4);

  // -------- 12 clean anchor CFs across all map pairs --------
  // Inject per-pair RT bias: rt(map_i) = base_rt + bias[i] + jitter
  // bias[0]=0.0, bias[1]=0.5, bias[2]=-0.3, bias[3]=0.7
  std::vector<double> bias = {0.0, 0.5, -0.3, 0.7};
  std::vector<double> jitters = {-0.4, -0.2, 0.0, 0.2, 0.4, -0.3, 0.1, -0.1, 0.3, 0.0, -0.2, 0.2};
  // distinct valid amino-acid sequences for the 12 anchor CFs
  const char* anchor_seqs[12] = {
    "PEPTIDEAK", "PEPTIDECK", "PEPTIDEDK", "PEPTIDEEK", "PEPTIDEFK", "PEPTIDEGK",
    "PEPTIDEHK", "PEPTIDEIK", "PEPTIDELK", "PEPTIDEMK", "PEPTIDENK", "PEPTIDEQK"};
  UInt64 uid = 100;
  for (Size i = 0; i < 12; ++i)
  {
    const double base_rt = 100.0 + 50.0 * i;
    const String seq = anchor_seqs[i];
    std::vector<std::pair<UInt64, Feature>> items;
    for (Size m = 0; m < 4; ++m)
    {
      items.emplace_back(m, makeFeature(uid++, base_rt + bias[m] + jitters[i], 500.0 + i, 1e5, seq, 2));
    }
    cmap.push_back(makeCF(items));
  }

  // -------- 5 CFs with consistent transfer (3 IDs + 1 ID-less, RT consistent) --------
  // The unidentified handle in map 3 has RT consistent with the per-pair bias.
  const char* consistent_seqs[5] = {"AAAAAAK", "AAAAACK", "AAAAADK", "AAAAAEK", "AAAAAFK"};
  for (Size i = 0; i < 5; ++i)
  {
    const double base_rt = 1000.0 + 50.0 * i;
    const String seq = consistent_seqs[i];
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(0, makeFeature(uid++, base_rt + bias[0], 600.0 + i, 1e5, seq, 2));
    items.emplace_back(1, makeFeature(uid++, base_rt + bias[1], 600.0 + i, 1e5, seq, 2));
    items.emplace_back(2, makeFeature(uid++, base_rt + bias[2], 600.0 + i, 1e5, seq, 2));
    // ID-less in map 3, but RT close to the expected bias
    items.emplace_back(3, makeFeature(uid++, base_rt + bias[3] + 0.1, 600.0 + i, 5e4));
    cmap.push_back(makeCF(items));
  }

  // -------- 5 CFs with outlier transfer (3 IDs + 1 ID-less far away) --------
  const char* outlier_seqs[5] = {"GGGGGGK", "GGGGGAK", "GGGGGCK", "GGGGGDK", "GGGGGEK"};
  for (Size i = 0; i < 5; ++i)
  {
    const double base_rt = 2000.0 + 50.0 * i;
    const String seq = outlier_seqs[i];
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(0, makeFeature(uid++, base_rt + bias[0], 700.0 + i, 1e5, seq, 2));
    items.emplace_back(1, makeFeature(uid++, base_rt + bias[1], 700.0 + i, 1e5, seq, 2));
    items.emplace_back(2, makeFeature(uid++, base_rt + bias[2], 700.0 + i, 1e5, seq, 2));
    // ID-less in map 3, but 30 seconds off the expected RT
    items.emplace_back(3, makeFeature(uid++, base_rt + bias[3] + 30.0, 700.0 + i, 5e4));
    cmap.push_back(makeCF(items));
  }

  // -------- Edge case: singleton CF (1 handle) --------
  cmap.push_back(makeCF({{0, makeFeature(uid++, 5000.0, 800.0, 1e5, "SOLO", 2)}}));

  // -------- Edge case: all-no-ID CF --------
  {
    std::vector<std::pair<UInt64, Feature>> items;
    for (Size m = 0; m < 4; ++m)
    {
      items.emplace_back(m, makeFeature(uid++, 6000.0 + m, 900.0, 1e5));
    }
    cmap.push_back(makeCF(items));
  }

  // -------- Edge case: all-IDs CF (no transfers to filter) --------
  {
    std::vector<std::pair<UInt64, Feature>> items;
    for (Size m = 0; m < 4; ++m)
    {
      items.emplace_back(m, makeFeature(uid++, 7000.0 + bias[m], 1000.0, 1e5, "ALLID", 2));
    }
    cmap.push_back(makeCF(items));
  }

  // ---- run filter ----
  ConsensusMapMBRFilter f;
  Param p = f.getParameters();
  p.setValue("k_sigma", 3.0);
  p.setValue("min_anchors_per_pair", 5);
  p.setValue("reject_fraction", 0.5);
  p.setValue("min_pair_sigma", 0.5);
  p.setValue("use_pooled_fallback", "true");
  f.setParameters(p);

  auto rep = f.filter(cmap);

  // Pair (0,1) anchors: 12 clean CFs + 5 consistent-transfer (both maps 0/1 have IDs)
  // + 5 outlier-transfer (also both have IDs) + 1 all-IDs CF = 23 anchors.
  // Median: rt(0) - rt(1) = -bias[1] = -0.5.
  const ConsensusMapMBRFilter::PairKey k01(0, 1);
  TEST_EQUAL(rep.pair_stats.count(k01), 1)
  const ConsensusMapMBRFilter::PairStats& s01 = rep.pair_stats[k01];
  TEST_EQUAL(s01.n, 23)
  TEST_REAL_SIMILAR(s01.median, -0.5)
  // sigma should be small but positive (floored at min_pair_sigma)
  const bool sigma_ok = (s01.mad_sigma >= 0.5);
  TEST_EQUAL(sigma_ok, true)

  // All 12 anchor CFs untouched
  for (Size i = 0; i < 12; ++i)
  {
    TEST_EQUAL(cmap[i].size(), 4)
  }
  // 5 consistent transfer CFs untouched
  for (Size i = 12; i < 17; ++i)
  {
    TEST_EQUAL(cmap[i].size(), 4)
  }
  // 5 outlier transfer CFs each lost their map-3 handle
  for (Size i = 17; i < 22; ++i)
  {
    TEST_EQUAL(cmap[i].size(), 3)
    // The remaining handles are exactly maps 0..2
    Size mask = 0;
    for (const auto& fh : cmap[i].getFeatures()) { mask |= (1u << fh.getMapIndex()); }
    const Size expected_mask = 0b0111u;
    TEST_EQUAL(mask, expected_mask)
    // PIDs are intact
    const bool pids_ok = (cmap[i].getPeptideIdentifications().size() >= 3);
    TEST_EQUAL(pids_ok, true)
  }
  // Singleton CF untouched
  TEST_EQUAL(cmap[22].size(), 1)
  // All-no-ID CF untouched
  TEST_EQUAL(cmap[23].size(), 4)
  // All-IDs CF untouched
  TEST_EQUAL(cmap[24].size(), 4)

  // 5 handles removed in total
  TEST_EQUAL(rep.handles_removed, 5)
}
END_SECTION

START_SECTION([EXTRA] anchors require matching charge)
{
  // Two CFs of the same modified sequence but different charges should NOT be paired as anchors.
  ConsensusMap cmap;
  setupHeaders(cmap, 2);
  UInt64 uid = 1;
  // 6 CFs where map 0 is charge 2 and map 1 is charge 3 — same sequence, different precursors.
  // No anchors should be collected because the (mod_seq, charge) keys differ.
  for (Size i = 0; i < 6; ++i)
  {
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(0, makeFeature(uid++, 100.0 + 50.0 * i, 500.0 + i, 1e5, "PEPTIDEAK", 2));
    items.emplace_back(1, makeFeature(uid++, 100.0 + 50.0 * i + 0.3, 500.0 + i, 1e5, "PEPTIDEAK", 3));
    cmap.push_back(makeCF(items));
  }
  ConsensusMapMBRFilter f;
  auto rep = f.filter(cmap);
  // pair (0,1) should not exist in the stats because no anchors matched
  TEST_EQUAL(rep.pair_stats.empty(), true)
  TEST_EQUAL(rep.handles_removed, 0)
}
END_SECTION

START_SECTION([EXTRA] PIDs with no hits or no map_index meta value are ignored)
{
  // Build a CF where one of the contributing features has a PID with empty hits.
  // The filter should not crash and should treat that handle as ID-less.
  ConsensusMap cmap;
  setupHeaders(cmap, 2);

  Feature f0;
  f0.setUniqueId(1);
  f0.setRT(100.0); f0.setMZ(500.0); f0.setIntensity(1e5); f0.setCharge(2);
  PeptideIdentification empty_pid; // no hits set
  empty_pid.setScoreType("score");
  empty_pid.setHigherScoreBetter(true);
  f0.setPeptideIdentifications({empty_pid});

  Feature f1 = makeFeature(2, 200.0, 500.0, 1e5, "PEPTIDEAK", 2);

  ConsensusFeature cf;
  cf.setUniqueId(42);
  cf.insert(0, f0);
  cf.insert(1, f1);
  cf.computeConsensus();
  cmap.push_back(cf);

  ConsensusMapMBRFilter f;
  auto rep = f.filter(cmap);
  // f0's empty PID should not produce an anchor; f1 has a real ID, f0 looks ID-less.
  // With only 1 ID-bearing handle and no anchors at all, no per-pair stats exist;
  // the pooled fallback is also empty, so f0 is left alone.
  TEST_EQUAL(rep.handles_removed, 0)
  TEST_EQUAL(cmap[0].size(), 2)
}
END_SECTION

START_SECTION([EXTRA] under-anchored pair uses pooled fallback)
{
  ConsensusMap cmap;
  setupHeaders(cmap, 3);

  // 10 clean anchors for pair (0,1) and (0,2) but only 1 for pair (1,2).
  const char* foo_seqs[10] = {
    "AAAFOOAK", "AAAFOOCK", "AAAFOODK", "AAAFOOEK", "AAAFOOFK",
    "AAAFOOGK", "AAAFOOHK", "AAAFOOIK", "AAAFOOLK", "AAAFOOMK"};
  const char* bar_seqs[10] = {
    "BARAAAAK", "BARAACAK", "BARAADAK", "BARAAEAK", "BARAAFAK",
    "BARAAGAK", "BARAAHAK", "BARAAIAK", "BARAALAK", "BARAAMAK"};
  UInt64 uid = 1;
  for (Size i = 0; i < 10; ++i)
  {
    const double base_rt = 100.0 + 50.0 * i;
    const String seq = foo_seqs[i];
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(0, makeFeature(uid++, base_rt + 0.0, 500.0 + i, 1e5, seq, 2));
    items.emplace_back(1, makeFeature(uid++, base_rt + 0.3, 500.0 + i, 1e5, seq, 2));
    cmap.push_back(makeCF(items));
  }
  for (Size i = 0; i < 10; ++i)
  {
    const double base_rt = 1000.0 + 50.0 * i;
    const String seq = bar_seqs[i];
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(0, makeFeature(uid++, base_rt + 0.0, 600.0 + i, 1e5, seq, 2));
    items.emplace_back(2, makeFeature(uid++, base_rt - 0.2, 600.0 + i, 1e5, seq, 2));
    cmap.push_back(makeCF(items));
  }
  // single anchor for pair (1,2)
  {
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(1, makeFeature(uid++, 2000.0, 700.0, 1e5, "SINGLE", 2));
    items.emplace_back(2, makeFeature(uid++, 2000.0 + 0.4, 700.0, 1e5, "SINGLE", 2));
    cmap.push_back(makeCF(items));
  }

  // CF with map 1 ID + map 2 ID-less, 0.4 sec offset (consistent with pooled distribution)
  {
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(1, makeFeature(uid++, 3000.0, 800.0, 1e5, "TRANSFERAK", 2));
    items.emplace_back(2, makeFeature(uid++, 3000.4, 800.0, 5e4));
    cmap.push_back(makeCF(items));
  }
  const Size idx_consistent = cmap.size() - 1;

  // CF with map 1 ID + map 2 ID-less, 50 sec offset (clear outlier under any pooled stat)
  {
    std::vector<std::pair<UInt64, Feature>> items;
    items.emplace_back(1, makeFeature(uid++, 4000.0, 900.0, 1e5, "TRANSFERCK", 2));
    items.emplace_back(2, makeFeature(uid++, 4050.0, 900.0, 5e4));
    cmap.push_back(makeCF(items));
  }
  const Size idx_outlier = cmap.size() - 1;

  ConsensusMapMBRFilter f;
  Param p = f.getParameters();
  p.setValue("k_sigma", 3.0);
  p.setValue("min_anchors_per_pair", 5);
  p.setValue("reject_fraction", 0.5);
  p.setValue("min_pair_sigma", 0.5);
  p.setValue("use_pooled_fallback", "true");
  f.setParameters(p);

  auto rep = f.filter(cmap);

  // Pooled distribution should be populated from pairs (0,1) and (0,2)
  const bool pooled_ok = (rep.pooled.n >= 20);
  TEST_EQUAL(pooled_ok, true)

  // Consistent transfer CF unchanged
  TEST_EQUAL(cmap[idx_consistent].size(), 2)
  // Outlier transfer CF lost its map-2 handle
  TEST_EQUAL(cmap[idx_outlier].size(), 1)
  TEST_EQUAL(rep.handles_removed, 1)
}
END_SECTION

END_TEST
