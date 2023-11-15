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

#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IsotopeMarker, "$Id$")

/////////////////////////////////////////////////////////////

IsotopeMarker* e_ptr = nullptr;
IsotopeMarker* e_nullPointer = nullptr;

START_SECTION((IsotopeMarker()))
	e_ptr = new IsotopeMarker;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~IsotopeMarker()))
	delete e_ptr;
END_SECTION

e_ptr = new IsotopeMarker();

START_SECTION((IsotopeMarker(const IsotopeMarker& source)))
	IsotopeMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((IsotopeMarker& operator=(const IsotopeMarker& source)))
	IsotopeMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);
	
	TEST_EQUAL(marked.size(), 48)

	Param iso_param = e_ptr->getParameters();
	iso_param.setValue("marks", 2);
	e_ptr->setParameters(iso_param);
	marked.clear();
	e_ptr->apply(marked, spec);
	TEST_EQUAL(marked.size(), 17)
END_SECTION

START_SECTION((static PeakMarker* create()))
	PeakMarker* pm = IsotopeMarker::create();
	IsotopeMarker im;
	TEST_EQUAL(pm->getParameters(), im.getParameters())
	TEST_EQUAL(pm->getName(), im.getName())
	delete pm;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "IsotopeMarker")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
