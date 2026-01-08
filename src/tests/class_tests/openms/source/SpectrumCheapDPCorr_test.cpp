// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/COMPARISON/SpectrumCheapDPCorr.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumCheapDPCorr, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumCheapDPCorr* e_ptr = nullptr;
SpectrumCheapDPCorr* e_nullPointer = nullptr;
START_SECTION(SpectrumCheapDPCorr())
	e_ptr = new SpectrumCheapDPCorr;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~SpectrumCheapDPCorr())
	delete e_ptr;
END_SECTION

e_ptr = new SpectrumCheapDPCorr();

START_SECTION(SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(SpectrumCheapDPCorr& operator = (const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests_2.dta"), spec2);

	double score = (*e_ptr)(spec1, spec2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(score, 10145.4)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
	SpectrumCheapDPCorr corr;
	score = corr(spec1, spec2);
	TEST_REAL_SIMILAR(score, 10145.4)

	score = corr(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
END_SECTION

START_SECTION(const PeakSpectrum& lastconsensus() const)
	TEST_EQUAL(e_ptr->lastconsensus().size(), 121)
END_SECTION

START_SECTION((Map<UInt, UInt> getPeakMap() const))
	TEST_EQUAL(e_ptr->getPeakMap().size(), 121)
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
  DTAFile dta_file;
  PeakSpectrum spec1;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

  double score = (*e_ptr)(spec1);

  TEST_REAL_SIMILAR(score, 12295.5)

END_SECTION

START_SECTION(void setFactor(double f))
	e_ptr->setFactor(0.3);

	TEST_EXCEPTION(Exception::OutOfRange, e_ptr->setFactor(1.1))
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
