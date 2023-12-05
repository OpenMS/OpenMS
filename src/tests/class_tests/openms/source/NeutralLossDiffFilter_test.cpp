// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NeutralLossDiffFilter, "$Id$")

/////////////////////////////////////////////////////////////

NeutralLossDiffFilter* e_ptr = nullptr;
NeutralLossDiffFilter* e_nullPointer = nullptr;

START_SECTION((NeutralLossDiffFilter()))
	e_ptr = new NeutralLossDiffFilter;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~NeutralLossDiffFilter()))
	delete e_ptr;
END_SECTION

e_ptr = new NeutralLossDiffFilter();

START_SECTION((NeutralLossDiffFilter(const NeutralLossDiffFilter& source)))
	NeutralLossDiffFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((NeutralLossDiffFilter& operator = (const NeutralLossDiffFilter& source)))
	NeutralLossDiffFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	double filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 429.0)

	Param p(e_ptr->getParameters());
	p.setValue("tolerance", 10.0);
	e_ptr->setParameters(p);
	filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 2482.5)
END_SECTION

START_SECTION((static FilterFunctor* create()))
	FilterFunctor* ff = NeutralLossDiffFilter::create();
	NeutralLossDiffFilter filter;
	TEST_EQUAL(ff->getParameters(), filter.getParameters())
	TEST_EQUAL(ff->getName(), filter.getName())
	delete ff;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "NeutralLossDiffFilter");
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
