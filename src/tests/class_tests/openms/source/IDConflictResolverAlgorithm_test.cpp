// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/ID/IDConflictResolverAlgorithm.h>

using namespace OpenMS;
using namespace std;

START_TEST(IDConflictResolverAlgorithm, "$Id$")

START_SECTION(resolveBetweenFeatures())
{
  FeatureMap map;
  Feature f1;
  Feature f2;
  Feature f3;
  Feature f4;

  PeptideHit hit;
  hit.setScore(23);
  hit.setSequence(AASequence::fromString("MORRISSEY"));
  PeptideIdentification id;
  id.insertHit(hit);
  PeptideIdentificationList ids;
  ids.push_back(id);
  
  PeptideHit hit2;
  hit2.setScore(23);
  hit2.setSequence(AASequence::fromString("M(Oxidation)ORRISSEY"));
  PeptideIdentification id2;
  id2.insertHit(hit2);
  PeptideIdentificationList ids2;
  ids2.push_back(id2);
  
  f1.setRT(1600.5);
  f1.setMZ(400.7);
  f1.setIntensity(1000.0);
  f1.setCharge(2);
  f1.setOverallQuality(1.0);
  f1.setPeptideIdentifications(ids);
  
  f2.setRT(1600.5);
  f2.setMZ(400.7);
  f2.setIntensity(10000.0);
  f2.setCharge(2);
  f2.setOverallQuality(1.0);
  f2.setPeptideIdentifications(ids);
  
  f3.setRT(1600.5);
  f3.setMZ(400.7);
  f3.setIntensity(1000.0);
  f3.setCharge(3);
  f3.setOverallQuality(1.0);
  f3.setPeptideIdentifications(ids);
  
  f4.setRT(1600.5);
  f4.setMZ(400.7);
  f4.setIntensity(1001.0);
  f4.setCharge(2);
  f4.setOverallQuality(1.0);
  f4.setPeptideIdentifications(ids2);
  
  map.push_back(f1);
  map.push_back(f2);
  
  IDConflictResolverAlgorithm::resolveBetweenFeatures(map);
  
  for (FeatureMap::ConstIterator it = map.begin(); it != map.end(); ++it)
  {
    
    if ((it->getIntensity() == 1000.0) && (it->getCharge() == 2))
    {
      // This identification was removed by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), true)
    }
      
    if ((it->getIntensity() == 10000.0) && (it->getCharge() == 2))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }
    
    if ((it->getIntensity() == 1000.0) && (it->getCharge() == 3))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }
    
    if ((it->getIntensity() == 1001.0) && (it->getCharge() == 2))
    {
      // This identification remains unchanged by the resolveBetweenFeatures() method.
      TEST_EQUAL(it->getPeptideIdentifications().empty(), false)
    }

  }
      
}
END_SECTION

START_SECTION(resolveAllHitRankAggregation())
{
  // Test rank aggregation on a ConsensusMap where the winner by rank aggregation
  // DIFFERS from the winner by best single-run score. This validates the key value
  // of the rank aggregation approach over simply picking the best-scoring ID.
  //
  // Feature with 3 IDs (simulating 3 replicates):
  //   ID1: SEQ_B score=0.99 (rank 0), SEQ_A score=0.5 (rank 1)
  //   ID2: SEQ_A score=0.8  (rank 0), SEQ_B score=0.1 (rank 1)
  //   ID3: SEQ_A score=0.7  (rank 0), SEQ_B score=0.05 (rank 1)
  //
  // best_score picks SEQ_B (single-run high score 0.99 from ID1).
  //
  // max_hits = 2, n_runs = 3
  //
  // Rank sums (rank 0-based, penalty = max_hits = 2 for missing runs):
  //   SEQ_A: 1 (ID1) + 0 (ID2) + 0 (ID3) = 1, found in 3/3 -> score = 1 - 1/(2*3) = 5/6 ≈ 0.833
  //   SEQ_B: 0 (ID1) + 1 (ID2) + 1 (ID3) = 2, found in 3/3 -> score = 1 - 2/(2*3) = 4/6 ≈ 0.667
  //
  // rank_aggregation picks SEQ_A (consistently ranks first in 2/3 replicates).
  // Best original score for SEQ_A is 0.8 from ID2, so ID2 is kept.

  ConsensusMap cmap;
  ConsensusFeature cf;

  AASequence seqA = AASequence::fromString("SEQA");
  AASequence seqB = AASequence::fromString("SEQB");

  // ID1: SEQ_B (best single-run score 0.99), SEQ_A second
  PeptideHit hitB1; hitB1.setScore(0.99); hitB1.setSequence(seqB);
  PeptideHit hitA1; hitA1.setScore(0.5);  hitA1.setSequence(seqA);
  PeptideIdentification id1;
  id1.setHigherScoreBetter(true);
  id1.setScoreType("score");
  id1.insertHit(hitB1);
  id1.insertHit(hitA1);

  // ID2: SEQ_A (rank 0), SEQ_B second
  PeptideHit hitA2; hitA2.setScore(0.8); hitA2.setSequence(seqA);
  PeptideHit hitB2; hitB2.setScore(0.1); hitB2.setSequence(seqB);
  PeptideIdentification id2;
  id2.setHigherScoreBetter(true);
  id2.setScoreType("score");
  id2.insertHit(hitA2);
  id2.insertHit(hitB2);

  // ID3: SEQ_A (rank 0), SEQ_B second
  PeptideHit hitA3; hitA3.setScore(0.7);  hitA3.setSequence(seqA);
  PeptideHit hitB3; hitB3.setScore(0.05); hitB3.setSequence(seqB);
  PeptideIdentification id3;
  id3.setHigherScoreBetter(true);
  id3.setScoreType("score");
  id3.insertHit(hitA3);
  id3.insertHit(hitB3);

  PeptideIdentificationList pep_ids;
  pep_ids.push_back(id1);
  pep_ids.push_back(id2);
  pep_ids.push_back(id3);
  cf.setPeptideIdentifications(pep_ids);

  cmap.push_back(cf);

  IDConflictResolverAlgorithm::resolveAllHitRankAggregation(cmap);

  // After resolution: each feature should have exactly 1 PeptideIdentification with 1 hit
  TEST_EQUAL(cmap.size(), 1)
  const PeptideIdentificationList& result_ids = cmap[0].getPeptideIdentifications();
  TEST_EQUAL(result_ids.size(), 1)
  TEST_EQUAL(result_ids[0].getHits().size(), 1)
  // SEQ_A should be selected as the winner by rank aggregation,
  // even though SEQ_B has the highest single-run score (0.99).
  TEST_EQUAL(result_ids[0].getHits()[0].getSequence(), seqA)
  // ID2 has the best original score (0.8) for SEQ_A among all IDs
  TEST_REAL_SIMILAR(result_ids[0].getHits()[0].getScore(), 0.8)

  // The other 2 IDs should have been moved to unassigned
  TEST_EQUAL(cmap.getUnassignedPeptideIdentifications().size(), 2)
}
END_SECTION

END_TEST
