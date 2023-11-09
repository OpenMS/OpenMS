// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Andreas Bertsch, Sven Nahnsen, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/ConsensusIDAlgorithmRanks.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ConsensusIDAlgorithmRanks, "$Id$")

/////////////////////////////////////////////////////////////

ConsensusIDAlgorithm* ptr = nullptr;
ConsensusIDAlgorithm* null_pointer = nullptr;
START_SECTION(ConsensusIDAlgorithmRanks())
{
  ptr = new ConsensusIDAlgorithmRanks();
  TEST_NOT_EQUAL(ptr, null_pointer);
}
END_SECTION

START_SECTION(~ConsensusIDAlgorithmRanks())
{
  delete(ptr);
}
END_SECTION

// create 3 ID runs:
PeptideIdentification temp;
temp.setScoreType("Posterior Error Probability");
temp.setHigherScoreBetter(false);
vector<PeptideIdentification> ids(3, temp);
vector<PeptideHit> hits;
// the first ID has 5 hits
hits.resize(5);
hits[0].setSequence(AASequence::fromString("A"));
hits[0].setScore(0.1);
hits[1].setSequence(AASequence::fromString("B"));
hits[1].setScore(0.2);
hits[2].setSequence(AASequence::fromString("C"));
hits[2].setScore(0.3);
hits[3].setSequence(AASequence::fromString("D"));
hits[3].setScore(0.4);
hits[4].setSequence(AASequence::fromString("E"));
hits[4].setScore(0.5);
ids[0].setHits(hits);
// the second ID has 3 hits
hits.resize(3);
hits[0].setSequence(AASequence::fromString("C"));
hits[0].setScore(0.2);
hits[1].setSequence(AASequence::fromString("A"));
hits[1].setScore(0.4);
hits[2].setSequence(AASequence::fromString("B"));
hits[2].setScore(0.6);
ids[1].setHits(hits);
// the third ID has 10 hits
hits.resize(10);
hits[0].setSequence(AASequence::fromString("F"));
hits[0].setScore(0.0);
hits[1].setSequence(AASequence::fromString("C"));
hits[1].setScore(0.1);
hits[2].setSequence(AASequence::fromString("G"));
hits[2].setScore(0.2);
hits[3].setSequence(AASequence::fromString("D"));
hits[3].setScore(0.3);
hits[4].setSequence(AASequence::fromString("B"));
hits[4].setScore(0.4);
hits[5].setSequence(AASequence::fromString("E"));
hits[5].setScore(0.5);
hits[6].setSequence(AASequence::fromString("H"));
hits[6].setScore(0.6);
hits[7].setSequence(AASequence::fromString("I"));
hits[7].setScore(0.7);
hits[8].setSequence(AASequence::fromString("J"));
hits[8].setScore(0.8);
hits[9].setSequence(AASequence::fromString("K"));
hits[9].setScore(0.9);
ids[2].setHits(hits);

START_SECTION(void apply(std::vector<PeptideIdentification>& ids))
{
  TOLERANCE_ABSOLUTE(0.01)

  ConsensusIDAlgorithmRanks consensus;
  // define parameters:
  Param param;
  param.setValue("filter:considered_hits", 5);
  consensus.setParameters(param);
  // apply:
  map<String,String> empty;
  vector<PeptideIdentification> f = ids;
  consensus.apply(f, empty);

  TEST_EQUAL(f.size(), 1);
  hits = f[0].getHits();
  TEST_EQUAL(hits.size(), 7);

  TEST_EQUAL(hits[0].getRank(), 1);
  TEST_EQUAL(hits[0].getSequence(), AASequence::fromString("C"));
  TEST_REAL_SIMILAR(hits[0].getScore(), 0.8);

  TEST_EQUAL(hits[1].getRank(), 2);
  TEST_EQUAL(hits[1].getSequence(), AASequence::fromString("A"));
  TEST_REAL_SIMILAR(hits[1].getScore(), 0.6);

  TEST_EQUAL(hits[2].getRank(), 3);
  TEST_EQUAL(hits[2].getSequence(), AASequence::fromString("B"));
  TEST_REAL_SIMILAR(hits[2].getScore(), 0.5333);

  TEST_EQUAL(hits[3].getRank(), 4);
  TEST_EQUAL(hits[3].getSequence(), AASequence::fromString("F"));
  TEST_REAL_SIMILAR(hits[3].getScore(), 0.33333);

  TEST_EQUAL(hits[4].getRank(), 5);
  TEST_EQUAL(hits[4].getSequence(), AASequence::fromString("D"));
  TEST_REAL_SIMILAR(hits[4].getScore(), 0.26666);

  TEST_EQUAL(hits[5].getRank(), 6);
  TEST_EQUAL(hits[5].getSequence(), AASequence::fromString("G"));
  TEST_REAL_SIMILAR(hits[5].getScore(), 0.2);

  TEST_EQUAL(hits[6].getRank(), 7);
  TEST_EQUAL(hits[6].getSequence(), AASequence::fromString("E"));
  TEST_REAL_SIMILAR(hits[6].getScore(), 0.06666);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
