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

#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////

START_TEST(SpectrumAlignmentScore, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SpectrumAlignmentScore* ptr = nullptr;
SpectrumAlignmentScore* nullPointer = nullptr;

START_SECTION(SpectrumAlignmentScore())
	ptr = new SpectrumAlignmentScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
	TEST_REAL_SIMILAR(score, 1.48268)

	s2.resize(100);

	score = (*ptr)(s1, s2);

	normalizer.filterSpectrum(s2);
	TEST_REAL_SIMILAR(score, 3.82472)
END_SECTION

START_SECTION(virtual ~SpectrumAlignmentScore())
	delete ptr;
END_SECTION

ptr = new SpectrumAlignmentScore();

START_SECTION(SpectrumAlignmentScore(const SpectrumAlignmentScore &source))
	SpectrumAlignmentScore sas1;
	Param p(sas1.getParameters());
	p.setValue("tolerance", 0.2);
	sas1.setParameters(p);

	SpectrumAlignmentScore sas2(sas1);

	TEST_EQUAL(sas1.getName(), sas2.getName())
	TEST_EQUAL(sas1.getParameters(), sas2.getParameters())

END_SECTION


START_SECTION(SpectrumAlignmentScore& operator=(const SpectrumAlignmentScore &source))
  SpectrumAlignmentScore sas1;
  Param p(sas1.getParameters());
  p.setValue("tolerance", 0.2);
  sas1.setParameters(p);

  SpectrumAlignmentScore sas2;

	sas2 = sas1;

  TEST_EQUAL(sas1.getName(), sas2.getName())
  TEST_EQUAL(sas1.getParameters(), sas2.getParameters())
END_SECTION

START_SECTION(double operator()(const PeakSpectrum &spec) const)
	PeakSpectrum s1;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);

  Normalizer normalizer;
  Param p(normalizer.getParameters());
  p.setValue("method", "to_one");
  normalizer.setParameters(p);
  normalizer.filterSpectrum(s1);

	double score = (*ptr)(s1);
	TEST_REAL_SIMILAR(score, 1.48268);
	
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
