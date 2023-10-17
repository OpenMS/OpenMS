// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TICFilter, "$Id$")

/////////////////////////////////////////////////////////////

TICFilter* e_ptr = nullptr;
TICFilter* e_nullPointer = nullptr;
START_SECTION((TICFilter()))
	e_ptr = new TICFilter;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~TICFilter()))
	delete e_ptr;
END_SECTION

e_ptr = new TICFilter();

START_SECTION((TICFilter(const TICFilter& source)))
	TICFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters());
	TEST_EQUAL(copy.getName(), e_ptr->getName());
END_SECTION

START_SECTION((TICFilter& operator=(const TICFilter& source)))
	TICFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters());
	TEST_EQUAL(copy.getName(), e_ptr->getName());
END_SECTION

START_SECTION((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	double filter  = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 533.5)
END_SECTION

START_SECTION((static FilterFunctor* create()))
	NOT_TESTABLE
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "TICFilter")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
