// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/PROCESSING/SCALING/SqrtScaler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <cmath>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SqrtScaler, "$Id$")

/////////////////////////////////////////////////////////////

SqrtScaler* e_ptr = nullptr;
SqrtScaler* e_nullPointer = nullptr;
START_SECTION((SqrtScaler()))
	e_ptr = new SqrtScaler;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~SqrtScaler()))
	delete e_ptr;
END_SECTION

e_ptr = new SqrtScaler();

START_SECTION((SqrtScaler(const SqrtScaler& source)))
	SqrtScaler copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION((SqrtScaler& operator=(const SqrtScaler& source)))
	SqrtScaler copy;
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
