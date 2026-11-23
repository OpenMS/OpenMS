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
#include <OpenMS/ANALYSIS/ID/ClassFDR.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  PeptideIdentification makeId(double score, bool decoy, const std::string& cls)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDE"));
    hit.setScore(score);
    hit.setTargetDecoyType(decoy ? PeptideHit::TargetDecoyType::DECOY : PeptideHit::TargetDecoyType::TARGET);
    hit.setMetaValue("cls", cls);
    PeptideIdentification id;
    id.setHigherScoreBetter(true);
    id.setScoreType("score");
    id.insertHit(hit);
    return id;
  }

  std::string classOf(const PeptideHit& hit)
  {
    return ClassFDR::classFromMetaValue(hit, "cls");
  }
}

START_TEST(ClassFDR, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::vector<double> storeyQValues(const std::vector<double>& p, double& pi0))
{
  double pi0 = 0.0;
  // two of four p-values above 0.5: pi0 = 2 / (0.5 * 4) = 1
  const std::vector<double> q = ClassFDR::storeyQValues({0.6, 0.01, 0.9, 0.02}, pi0);
  TEST_REAL_SIMILAR(pi0, 1.0)
  TEST_EQUAL(q.size(), 4)
  TEST_REAL_SIMILAR(q[0], 0.8)  // 0.6 * 4 / 3
  TEST_REAL_SIMILAR(q[1], 0.04) // min(0.01 * 4 / 1, 0.02 * 4 / 2)
  TEST_REAL_SIMILAR(q[2], 0.9)
  TEST_REAL_SIMILAR(q[3], 0.04)

  // one of four above 0.5: pi0 = 0.5
  ClassFDR::storeyQValues({0.6, 0.1, 0.2, 0.3}, pi0);
  TEST_REAL_SIMILAR(pi0, 0.5)
  TEST_EQUAL(ClassFDR::storeyQValues({}, pi0).empty(), true)
}
END_SECTION

START_SECTION(static std::vector<ClassSummary> annotateClassQValues(PeptideIdentificationList& ids, const std::function<std::string(const PeptideHit&)>& class_of, bool use_all_hits = false))
{
  PeptideIdentificationList ids;
  // class A: 100 targets above 30 decoys -> pooled null
  for (int i = 0; i < 100; ++i) ids.push_back(makeId(1000.0 - i, false, "A"));
  for (int i = 0; i < 30; ++i) ids.push_back(makeId(1.0 + i, true, "A"));
  // class B: 5 targets, 3 decoys (fewer than 20) -> class target-decoy competition
  const double b_targets[] = {50, 40, 30, 20, 10};
  const double b_decoys[] = {45, 35, 5};
  for (double s : b_targets) ids.push_back(makeId(s, false, "B"));
  for (double s : b_decoys) ids.push_back(makeId(s, true, "B"));

  const auto summaries = ClassFDR::annotateClassQValues(ids, classOf);
  TEST_EQUAL(summaries.size(), 2)
  TEST_EQUAL(summaries[0].name, "A")
  TEST_EQUAL(summaries[0].targets, 100)
  TEST_EQUAL(summaries[0].decoys, 30)
  TEST_EQUAL(summaries[0].method == ClassFDR::Method::POOLED_NULL, true)
  TEST_EQUAL(summaries[1].name, "B")
  TEST_EQUAL(summaries[1].targets, 5)
  TEST_EQUAL(summaries[1].decoys, 3)
  TEST_EQUAL(summaries[1].method == ClassFDR::Method::CLASS_TDC, true)

  // A: every target outscores all 33 pooled decoys, p = 1/34 < 0.5, so pi0 = 0 and q = 0
  TEST_REAL_SIMILAR(summaries[0].pi0, 0.0)
  TEST_REAL_SIMILAR(ids[0].getHits()[0].getMetaValue(ClassFDR::CLASS_Q_VALUE), 0.0)
  // B: T50 D45 T40 D35 T30 T20 T10 D5; (D + 1) / T is lowest at T10 (3 / 5), so all targets get 0.6
  for (Size i = 130; i < 135; ++i)
  {
    TEST_REAL_SIMILAR(ids[i].getHits()[0].getMetaValue(ClassFDR::CLASS_Q_VALUE), 0.6)
  }
  // decoys are not annotated
  TEST_EQUAL(ids[100].getHits()[0].metaValueExists(ClassFDR::CLASS_Q_VALUE), false)
  TEST_EQUAL(ids[135].getHits()[0].metaValueExists(ClassFDR::CLASS_Q_VALUE), false)

  // a single class: nothing annotated
  PeptideIdentificationList one;
  for (int i = 0; i < 10; ++i) one.push_back(makeId(100.0 - i, false, "A"));
  for (int i = 0; i < 10; ++i) one.push_back(makeId(1.0 + i, true, "A"));
  TEST_EQUAL(ClassFDR::annotateClassQValues(one, classOf).size(), 1)
  TEST_EQUAL(one[0].getHits()[0].metaValueExists(ClassFDR::CLASS_Q_VALUE), false)

  // a class with many decoys in the pool's top: the guard switches to class target-decoy competition
  PeptideIdentificationList heavy;
  for (int i = 0; i < 1000; ++i) heavy.push_back(makeId(2000.0 - i, false, "A"));
  for (int i = 0; i < 500; ++i) heavy.push_back(makeId(1.0 + i, true, "A"));
  for (int i = 0; i < 30; ++i) heavy.push_back(makeId(600.0 + i, false, "C"));
  for (int i = 0; i < 30; ++i) heavy.push_back(makeId(550.0 + i, true, "C"));
  const auto heavy_summaries = ClassFDR::annotateClassQValues(heavy, classOf);
  TEST_EQUAL(heavy_summaries.size(), 2)
  TEST_EQUAL(heavy_summaries[1].name, "C")
  TEST_EQUAL(heavy_summaries[1].method == ClassFDR::Method::CLASS_TDC, true)

  // only the best hit of an identification is used by default
  PeptideIdentificationList multi;
  for (int i = 0; i < 20; ++i)
  {
    PeptideIdentification id = makeId(100.0 - i, false, "A");
    PeptideHit second = id.getHits()[0];
    second.setScore(-1.0);
    second.setMetaValue("cls", "B");
    id.insertHit(second);
    multi.push_back(id);
  }
  TEST_EQUAL(ClassFDR::annotateClassQValues(multi, classOf).size(), 1)
  TEST_EQUAL(ClassFDR::annotateClassQValues(multi, classOf, true).size(), 2)
}
END_SECTION

START_SECTION(static Size raiseToClassQValues(PeptideIdentificationList& ids))
{
  PeptideIdentificationList ids;
  PeptideIdentification id;
  id.setScoreType("q-value");
  id.setHigherScoreBetter(false);
  PeptideHit low, high, none;
  low.setScore(0.01);
  low.setMetaValue(ClassFDR::CLASS_Q_VALUE, 0.05);
  high.setScore(0.2);
  high.setMetaValue(ClassFDR::CLASS_Q_VALUE, 0.1);
  none.setScore(0.3);
  id.setHits({low, high, none});
  ids.push_back(id);
  TEST_EQUAL(ClassFDR::raiseToClassQValues(ids), 1)
  TEST_REAL_SIMILAR(ids[0].getHits()[0].getScore(), 0.05)
  TEST_REAL_SIMILAR(ids[0].getHits()[1].getScore(), 0.2)
  TEST_REAL_SIMILAR(ids[0].getHits()[2].getScore(), 0.3)
}
END_SECTION

START_SECTION(static std::string classFromMetaValue(const PeptideHit& hit, const std::string& name))
{
  PeptideHit hit;
  TEST_EQUAL(ClassFDR::classFromMetaValue(hit, "cls"), "NA")
  hit.setMetaValue("cls", "isoform");
  TEST_EQUAL(ClassFDR::classFromMetaValue(hit, "cls"), "isoform")
}
END_SECTION

START_SECTION(static std::string classFromModifications(const PeptideHit& hit, const std::vector<std::string>& ignore))
{
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPC(Carbamidomethyl)TIDE"));
  TEST_EQUAL(ClassFDR::classFromModifications(hit, {"Carbamidomethyl (C)"}), "unmodified")
  TEST_EQUAL(ClassFDR::classFromModifications(hit, {}), "Carbamidomethyl")
  hit.setSequence(AASequence::fromString("PEPS(Phospho)M(Oxidation)T(Phospho)IDE"));
  TEST_EQUAL(ClassFDR::classFromModifications(hit, {}), "Oxidation+Phospho")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
