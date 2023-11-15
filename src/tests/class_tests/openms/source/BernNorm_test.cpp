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

#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(BernNorm, "$Id$")

/////////////////////////////////////////////////////////////

BernNorm* e_ptr = nullptr;
BernNorm* e_nullPointer = nullptr;

START_SECTION((BernNorm()))
	e_ptr = new BernNorm;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~BernNorm()))
	delete e_ptr;
END_SECTION

e_ptr = new BernNorm();

START_SECTION((BernNorm(const BernNorm& source)))
	BernNorm copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((BernNorm& operator=(const BernNorm& source)))
	BernNorm copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	TEST_EQUAL(spec.size(), 121)

	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("C2", 2000.0);
	e_ptr->setParameters(p);
	e_ptr->filterSpectrum(spec);

	TEST_EQUAL(spec.size(), 28)

END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	delete e_ptr;
	e_ptr = new BernNorm();

  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	
	PeakMap pm;
	pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("C2", 2000.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 28)


END_SECTION
			
START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	delete e_ptr;
	e_ptr = new BernNorm();

	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  TEST_EQUAL(spec.size(), 121)

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("C2", 2000.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 28)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
