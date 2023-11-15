// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Scaler, "$Id$")

/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.01)

Scaler* e_ptr = nullptr;
Scaler* e_nullPointer = nullptr;
START_SECTION((Scaler()))
	e_ptr = new Scaler;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~Scaler()))
	delete e_ptr;
END_SECTION

e_ptr = new Scaler();

START_SECTION((Scaler(const Scaler& source)))
	Scaler copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((Scaler& operator = (const Scaler& source)))
	Scaler copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 121)

	spec.sortByIntensity();
	TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
	TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
	TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
	
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.addSpectrum(spec);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 121)

  pm.begin()->sortByIntensity();
  TEST_REAL_SIMILAR(pm.begin()->begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getPosition()[0], 136.077)

END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 121)

  spec.sortByIntensity();
  TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
