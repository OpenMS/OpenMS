// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Volker Mosthaf $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <cmath>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SqrtMower, "$Id$")

/////////////////////////////////////////////////////////////

SqrtMower* e_ptr = nullptr;
SqrtMower* e_nullPointer = nullptr;
START_SECTION((SqrtMower()))
	e_ptr = new SqrtMower;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~SqrtMower()))
	delete e_ptr;
END_SECTION

e_ptr = new SqrtMower();

START_SECTION((SqrtMower(const SqrtMower& source)))
	SqrtMower copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION((SqrtMower& operator=(const SqrtMower& source)))
	SqrtMower copy;
	copy = *e_ptr;
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterSpectrum(spec);
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.addSpectrum(spec);

	TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterPeakMap(pm);
	TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterPeakSpectrum(spec);
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
