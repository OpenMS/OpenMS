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

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ComplementMarker, "$Id$")

/////////////////////////////////////////////////////////////

ComplementMarker* e_ptr = nullptr;
ComplementMarker* e_nullPointer = nullptr;

START_SECTION((ComplementMarker()))
	e_ptr = new ComplementMarker;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~ComplementMarker()))
	delete e_ptr;
END_SECTION

e_ptr = new ComplementMarker();

START_SECTION((ComplementMarker(const ComplementMarker& source)))
	ComplementMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ComplementMarker& operator = (const ComplementMarker& source)))
	ComplementMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void apply(std::map<double, bool> marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);
	
	TEST_EQUAL(marked.size(), 0)

	Param p(e_ptr->getParameters());
	p.setValue("marks", 10);
	p.setValue("tolerance", 10.0);
	e_ptr->setParameters(p);
	marked.clear();
	e_ptr->apply(marked, spec);
	TEST_EQUAL(marked.size(), 0)

END_SECTION

START_SECTION((static PeakMarker* create()))
	PeakMarker* pm = ComplementMarker::create();
	ComplementMarker cm;
	TEST_EQUAL(pm->getParameters(), cm.getParameters())
	TEST_EQUAL(pm->getName(), cm.getName())
	delete pm;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "ComplementMarker")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
