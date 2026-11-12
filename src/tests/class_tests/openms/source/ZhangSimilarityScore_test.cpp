// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <iostream>
#include <cmath>

#include <OpenMS/COMPARISON/ZhangSimilarityScore.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>

///////////////////////////

START_TEST(ZhangSimilarityScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ZhangSimilarityScore* ptr = nullptr;
ZhangSimilarityScore* nullPointer = nullptr;

START_SECTION(ZhangSimilarityScore())
	ptr = new ZhangSimilarityScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ZhangSimilarityScore())
	delete ptr;
END_SECTION

ptr = new ZhangSimilarityScore();

START_SECTION(ZhangSimilarityScore(const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(ZhangSimilarityScore& operator = (const ZhangSimilarityScore& source))
	ZhangSimilarityScore copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec) const)
	PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);

  double score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1.82682);
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const)
  PeakSpectrum s1, s2;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);
  normalizer.filterSpectrum(s2);

  TOLERANCE_ABSOLUTE(0.01)

  double score = (*ptr)(s1, s2);
  TEST_REAL_SIMILAR(score, 1.82682)

  s2.resize(100);

  score = (*ptr)(s1, s2);

  normalizer.filterSpectrum(s2);
  TEST_REAL_SIMILAR(score, 0.328749)
END_SECTION

START_SECTION([EXTRA] Gaussian factor uses the tolerance of each call)
{
  // The Gaussian weight depends on the tolerance; it must not stick to the tolerance of the first
  // call in the process (it was a function-local static).
  PeakSpectrum s1, s2;
  Peak1D p;
  p.setIntensity(1.0);
  p.setMZ(500.0);
  s1.push_back(p);
  p.setMZ(500.1);
  s2.push_back(p);

  // one peak pair with intensities 1: score = sqrt(erfc(diff / (tol * 3 * sqrt(2))))
  const auto expected = [](double tol) { return std::sqrt(std::erfc(0.1 / (tol * 3.0 * std::sqrt(2.0)))); };

  ZhangSimilarityScore narrow, wide;
  Param param(narrow.getParameters());
  param.setValue("use_gaussian_factor", "true");
  param.setValue("tolerance", 0.2);
  narrow.setParameters(param);
  param.setValue("tolerance", 1.0);
  wide.setParameters(param);

  TOLERANCE_ABSOLUTE(1e-6)
  TEST_REAL_SIMILAR(narrow(s1, s2), expected(0.2))
  TEST_REAL_SIMILAR(wide(s1, s2), expected(1.0))
  TEST_REAL_SIMILAR(narrow(s1, s2), expected(0.2))
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
