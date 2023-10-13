// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;


START_TEST(ComplementFilter, "$Id$")

/////////////////////////////////////////////////////////////

ComplementFilter* e_ptr = nullptr;
ComplementFilter* e_nullPointer = nullptr;

START_SECTION((ComplementFilter()))
	e_ptr = new ComplementFilter;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~ComplementFilter()))
	delete e_ptr;
END_SECTION

e_ptr = new ComplementFilter();

START_SECTION((ComplementFilter(const ComplementFilter& source)))
	ComplementFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ComplementFilter& operator = (const ComplementFilter& source)))
	ComplementFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	double filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 37.0)

	Param p;
	p.setValue("tolerance", 2.0);
	e_ptr->setParameters(p);
	filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 132.5)
	
END_SECTION

START_SECTION((static FilterFunctor* create()))
	FilterFunctor* ff = ComplementFilter::create();
	ComplementFilter cf;
	TEST_EQUAL(ff->getParameters(), cf.getParameters())
	TEST_EQUAL(ff->getName(), cf.getName())
	delete ff;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "ComplementFilter")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
