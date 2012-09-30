// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/MRMFeature.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MRMFeature, "$Id$")

/////////////////////////////////////////////////////////////

MRMFeature* ptr = 0;
MRMFeature* nullPointer = 0;

START_SECTION(MRMFeature())
{
	ptr = new MRMFeature();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeature())
{
  delete ptr;
}
END_SECTION

START_SECTION (const PGScoresType & getScores() const)
{
  // tested with set/add score
  NOT_TESTABLE
}
END_SECTION

START_SECTION (double getScore(const String & score_name))
{
  // tested with set/add score
  NOT_TESTABLE
}
END_SECTION

START_SECTION (Feature & getFeature(String key))
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getFeature("chromatogram1").getMetaValue("dummy"), 1)
}
END_SECTION

START_SECTION (void setScores(const PGScoresType & scores))
{
  MRMFeature::PGScoresType scores;
  MRMFeature mrmfeature;
  scores["score1"] = 1;
  scores["score2"] = 2;
  mrmfeature.setScores(scores);
  TEST_EQUAL(mrmfeature.getScore("score1"), 1)
  TEST_EQUAL(mrmfeature.getScore("score2"), 2)
  TEST_EQUAL(mrmfeature.getScores().at("score1"), 1)
  TEST_EQUAL(mrmfeature.getScores().at("score2"), 2)
}
END_SECTION

START_SECTION (void addScore(const String & score_name, double score))
{
  MRMFeature::PGScoresType scores;
  MRMFeature mrmfeature;
  mrmfeature.addScore("score1",1);
  mrmfeature.addScore("score2",2);
  TEST_EQUAL(mrmfeature.getScore("score1"), 1)
  TEST_EQUAL(mrmfeature.getScore("score2"), 2)
  TEST_EQUAL(mrmfeature.getScores().at("score1"), 1)
  TEST_EQUAL(mrmfeature.getScores().at("score2"), 2)
}
END_SECTION

START_SECTION (void addFeature(Feature & feature, String key))
{
  // tested in getFeature
  NOT_TESTABLE
}
END_SECTION

START_SECTION (const std::vector<Feature> & getFeatures() const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getFeatures().size(), 2)
}
END_SECTION

START_SECTION (void getFeatureIDs(std::vector<String> & result) const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addFeature(f1, "chromatogram1");
  mrmfeature.addFeature(f1, "chromatogram2");
  std::vector<String> result;
  mrmfeature.getFeatureIDs(result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result[0], "chromatogram1")
  TEST_EQUAL(result[1], "chromatogram2")
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
