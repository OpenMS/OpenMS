// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MarkerMower, "$Id$")

/////////////////////////////////////////////////////////////

MarkerMower* e_ptr = nullptr;
MarkerMower* e_nullPointer = nullptr;

START_SECTION((MarkerMower()))
	e_ptr = new MarkerMower;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~MarkerMower()))
	delete e_ptr;
END_SECTION

e_ptr = new MarkerMower();

START_SECTION((MarkerMower(const MarkerMower& source)))
	MarkerMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((MarkerMower& operator = (const MarkerMower& source)))
	MarkerMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("n", 10);
	e_ptr->setParameters(p);
	e_ptr->filterSpectrum(spec);
	TEST_EQUAL(spec.size(), 0)
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	delete e_ptr;
	e_ptr = new MarkerMower();
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("n", 10);
	e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 0)
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	delete e_ptr;
	e_ptr = new MarkerMower();
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  TEST_EQUAL(spec.size(), 121)
	
	Param p(e_ptr->getParameters());
	p.setValue("n", 10);
	e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 0)
END_SECTION
			
START_SECTION((void insertmarker(PeakMarker* peak_marker)))
	IsotopeMarker marker;
	e_ptr->insertmarker(&marker);
	// this function is not really testable, as there are no accessors
	NOT_TESTABLE
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
