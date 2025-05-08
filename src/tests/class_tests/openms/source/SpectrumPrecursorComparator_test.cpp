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
#include <OpenMS/COMPARISON/SpectrumPrecursorComparator.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;


START_TEST(SpectrumPrecursorComparator, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumPrecursorComparator* e_ptr = nullptr;
SpectrumPrecursorComparator* e_nullPointer = nullptr;
START_SECTION(SpectrumPrecursorComparator())
	e_ptr = new SpectrumPrecursorComparator;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~SpectrumPrecursorComparator())
	delete e_ptr;
END_SECTION

e_ptr = new SpectrumPrecursorComparator();

START_SECTION(SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source))
	SpectrumPrecursorComparator copy(*e_ptr);
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION(SpectrumPrecursorComparator& operator = (const SpectrumPrecursorComparator& source))
	SpectrumPrecursorComparator copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests_2.dta"), spec2);

	double score = (*e_ptr)(spec1, spec2);

	TEST_REAL_SIMILAR(score, 1.7685)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_SIMILAR(score, 2)
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

	TEST_REAL_SIMILAR((*e_ptr)(spec1), 2.0)

END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
