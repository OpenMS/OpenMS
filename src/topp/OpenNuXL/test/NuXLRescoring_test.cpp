// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include "NuXLRescoring.h"
///////////////////////////

#include <OpenMS/METADATA/PeptideIdentification.h>

#include <random>
#include <string>

using namespace OpenMS;
using namespace std;

namespace
{
  // identification of one spectrum with a single hit
  PeptideIdentification identification(int scan_index, double score, bool decoy, bool xl)
  {
    PeptideIdentification pi;
    pi.setMetaValue("scan_index", scan_index);
    pi.setHigherScoreBetter(true);
    PeptideHit ph;
    ph.setScore(score);
    ph.setMetaValue("target_decoy", decoy ? "decoy" : "target");
    ph.setMetaValue("NuXL:isXL", xl ? 1 : 0);
    pi.setHits({ph});
    return pi;
  }
}

START_TEST(NuXLRescoring, "$Id$")

START_SECTION(([ScoreGap] bestOtherThan is the best score of any other sequence for any order of the hits))
{
  const vector<string> sequences{"PEPTIDEK", "PEPTIDER", "SAMPLER", "ELVISK", "LIVESK", "AAAK"};
  std::mt19937 rng(11);
  Size mismatches = 0, queries = 0;
  for (int trial = 0; trial < 20000; ++trial)
  {
    // 1-8 hits over a few sequences, with repeated sequences and tied scores
    vector<pair<float, string_view>> hits;
    const int n = 1 + static_cast<int>(rng() % 8);
    for (int i = 0; i < n; ++i)
    {
      hits.emplace_back(static_cast<float>(rng() % 5), sequences[rng() % 4]);
    }
    std::shuffle(hits.begin(), hits.end(), rng);
    NuXLRescoring::ScoreGap gap;
    for (const auto& [score, sequence] : hits) { gap.add(score, sequence); }
    for (const string& query : sequences) // includes sequences without a hit
    {
      float expected = -std::numeric_limits<float>::infinity();
      for (const auto& [score, sequence] : hits)
      {
        if (sequence != query) expected = std::max(expected, score);
      }
      if (gap.bestOtherThan(query) != expected) ++mismatches;
      ++queries;
    }
  }
  TEST_EQUAL(queries, 120000)
  TEST_EQUAL(mismatches, 0)

  // a single sequence: no other sequence
  NuXLRescoring::ScoreGap single;
  single.add(3.0f, "PEPTIDEK");
  single.add(5.0f, "PEPTIDEK");
  TEST_EQUAL(std::isinf(single.bestOtherThan("PEPTIDEK")), true)
  TEST_REAL_SIMILAR(single.bestOtherThan("SAMPLER"), 5.0)
}
END_SECTION

START_SECTION((static std::array<std::vector<Size>, 2> confidentHitsPerHalf(const PeptideIdentificationList& peptide_ids, double fdr, std::array<std::array<Size, 2>, 2>& counts)))
{
  PeptideIdentificationList ids;
  // even scan indices (half 0), peptides: targets with scores 20 ... 1 and one decoy at 10.5. Below the decoy
  // (D + 1) / T = 2 / k, above it 1 / k, so every target has q = 2 / 20 = 0.1.
  for (int i = 0; i < 20; ++i) { ids.push_back(identification(2 * i, 20.0 - i, false, false)); }
  ids.push_back(identification(100, 10.5, true, false));
  // odd scan indices (half 1), cross-links: a decoy ties the best target and wins the tie
  ids.push_back(identification(1, 5.0, true, true));
  ids.push_back(identification(3, 5.0, false, true));
  ids.push_back(identification(5, 4.0, false, true));
  // skipped: no hits, no scan index
  ids.push_back(PeptideIdentification());
  PeptideIdentification no_scan = identification(7, 50.0, false, true);
  no_scan.removeMetaValue("scan_index");
  ids.push_back(no_scan);

  std::array<std::array<Size, 2>, 2> counts{};
  auto training = NuXLRescoring::confidentHitsPerHalf(ids, 0.1, counts);
  TEST_EQUAL(training[0].size(), 20)
  TEST_EQUAL(counts[0][0], 20)
  TEST_EQUAL(counts[0][1], 0)
  // half 1: decoy (q 2), then the targets at 5 (2 / 1) and 4 (2 / 2): q = 1 for both targets
  TEST_EQUAL(training[1].size(), 0)

  training = NuXLRescoring::confidentHitsPerHalf(ids, 0.099, counts);
  TEST_EQUAL(training[0].size(), 0)

  training = NuXLRescoring::confidentHitsPerHalf(ids, 1.0, counts);
  TEST_EQUAL(counts[1][1], 2)
  TEST_EQUAL(counts[1][0], 0)
  // no decoy is ever selected
  for (const auto& half : training)
  {
    for (Size id : half) { TEST_EQUAL(ids[id].getHits()[0].isDecoy(), false) }
  }
  // training hits keep the order of the scores
  TEST_EQUAL(training[0].front(), 0)
}
END_SECTION

START_SECTION(static void rescoreHits(PeptideIdentificationList& peptide_ids))
{
  PeptideIdentification pi;
  pi.setHigherScoreBetter(true);
  PeptideHit xl, peptide, bare;
  xl.setScore(10.0);
  xl.setMetaValue("NuXL:score", 10.0);
  xl.setMetaValue("NuXL:isXL", 1);
  xl.setMetaValue("NuXL:delta_score", 2.0);
  xl.setMetaValue("NuXL:ion_llr", -10.0);
  peptide.setScore(9.5);
  peptide.setMetaValue("NuXL:score", 9.5);
  peptide.setMetaValue("NuXL:isXL", 0);
  peptide.setMetaValue("NuXL:delta_score", 5.0); // not added for peptides
  peptide.setMetaValue("NuXL:ion_llr", 1.0);
  bare.setScore(1.0); // no meta values: unchanged
  pi.setHits({xl, peptide, bare});
  PeptideIdentificationList ids;
  ids.push_back(pi);

  NuXLRescoring::rescoreHits(ids);
  const vector<PeptideHit>& hits = ids[0].getHits();
  TEST_EQUAL(hits.size(), 3)
  // peptide: 9.5 + 0.3 * 1 = 9.8; cross-link: 10 + 2 - 0.3 * 10 = 9: the peptide is now the best hit
  TEST_REAL_SIMILAR(hits[0].getScore(), 9.8)
  TEST_EQUAL(static_cast<int>(hits[0].getMetaValue("NuXL:isXL")), 0)
  TEST_REAL_SIMILAR(hits[1].getScore(), 9.0)
  TEST_REAL_SIMILAR(hits[2].getScore(), 1.0)
  TEST_EQUAL(hits[0].getRank(), 0)
  TEST_EQUAL(hits[1].getRank(), 1)
  TEST_EQUAL(hits[2].getRank(), 2)
  // NuXL:score keeps the search score
  TEST_REAL_SIMILAR(static_cast<double>(hits[1].getMetaValue("NuXL:score")), 10.0)
}
END_SECTION

END_TEST
