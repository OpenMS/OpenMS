// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MRMFeature.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MRMFeature, "$Id$")

/////////////////////////////////////////////////////////////

MRMFeature* ptr = nullptr;
MRMFeature* nullPointer = nullptr;

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

START_SECTION(MRMFeature(const MRMFeature &rhs))
{
  MRMFeature tmp;
  tmp.setIntensity(100.0);
  tmp.addScore("testscore", 200);

  MRMFeature tmp2 (tmp);

  TEST_REAL_SIMILAR(tmp2.getMetaValue("testscore"), 200)
  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
}
END_SECTION

START_SECTION((MRMFeature(const MRMFeature&& source)))
{
#ifndef OPENMS_COMPILER_MSVC
  // Ensure that MRMFeature has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  // Note that MSVS does not support noexcept move constructors for STL
  // constructs such as std::map.
  TEST_EQUAL(noexcept(MRMFeature(std::declval<MRMFeature&&>())), true)
#endif
}
END_SECTION

START_SECTION(MRMFeature& operator=(const MRMFeature &rhs))
{
  MRMFeature tmp;
  tmp.setIntensity(100.0);
  tmp.addScore("testscore", 200);

  MRMFeature tmp2;
  tmp2 = tmp;

  TEST_REAL_SIMILAR(tmp2.getMetaValue("testscore"), 200)
  TEST_REAL_SIMILAR(tmp2.getIntensity(), 100.0)
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
  MRMFeature mrmfeature;
  OpenSwath_Scores scores;
  scores.library_sangle = 99;
  mrmfeature.setScores(scores);

  TEST_REAL_SIMILAR(scores.library_sangle, mrmfeature.getScores().library_sangle)
}
END_SECTION

START_SECTION (void addScore(const String & score_name, double score))
{
  MRMFeature mrmfeature;
  mrmfeature.addScore("score1",1);
  mrmfeature.addScore("score2",2);
  TEST_REAL_SIMILAR(mrmfeature.getMetaValue("score1"), 1)
  TEST_REAL_SIMILAR(mrmfeature.getMetaValue("score2"), 2)
}
END_SECTION

START_SECTION (void addFeature(Feature & feature, const String & key))
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

START_SECTION (void addPrecursorFeature(Feature & feature, const String & key))
{
  // Initially, there should be no feature present
  MRMFeature mrmfeature;
  {
    std::vector<String> result;
    mrmfeature.getPrecursorFeatureIDs(result);
    TEST_EQUAL(result.size(), 0)
  }

  // After adding a feature, there should be one feature present
  Feature f1;
  mrmfeature.addPrecursorFeature(f1, "precursor_chromatogram1");
  {
    std::vector<String> result;
    mrmfeature.getPrecursorFeatureIDs(result);
    TEST_EQUAL(result.size(), 1)
  }
}
END_SECTION

START_SECTION (void getPrecursorFeatureIDs(std::vector<String> & result) const)
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addPrecursorFeature(f1, "chromatogram1");
  mrmfeature.addPrecursorFeature(f1, "chromatogram2");
  std::vector<String> result;
  mrmfeature.getPrecursorFeatureIDs(result);
  TEST_EQUAL(result.size(), 2)
  TEST_EQUAL(result[0], "chromatogram1")
  TEST_EQUAL(result[1], "chromatogram2")
}
END_SECTION

START_SECTION (Feature & getPrecursorFeature(String key))
{
  MRMFeature mrmfeature;
  Feature f1;
  f1.setMetaValue("dummy", 1);
  Feature f2;
  mrmfeature.addPrecursorFeature(f1, "chromatogram1");
  mrmfeature.addPrecursorFeature(f1, "chromatogram2");
  TEST_EQUAL(mrmfeature.getPrecursorFeature("chromatogram1").getMetaValue("dummy"), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
