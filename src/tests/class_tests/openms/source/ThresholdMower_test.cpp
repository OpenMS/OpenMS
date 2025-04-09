// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ThresholdMower, "$Id$")

/////////////////////////////////////////////////////////////

ThresholdMower* e_ptr = nullptr;
ThresholdMower* e_nullPointer = nullptr;
START_SECTION((ThresholdMower()))
  e_ptr = new ThresholdMower;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~ThresholdMower()))
  delete e_ptr;
END_SECTION

e_ptr = new ThresholdMower();

START_SECTION((ThresholdMower(const ThresholdMower& source)))
  ThresholdMower copy(*e_ptr);
  TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
  TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ThresholdMower& operator=(const ThresholdMower& source)))
  ThresholdMower copy;
  copy = *e_ptr;
  TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
  TEST_EQUAL(copy.getName(), e_ptr->getName());
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  
  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterSpectrum(spec);
  TEST_EQUAL(spec.size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);

  e_ptr->filterSpectrum(spec);
  TEST_EQUAL(spec.size(), 14)
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  PeakMap pm;
  pm.addSpectrum(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 14)

END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("threshold", 1.0);
  e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 121)

  p.setValue("threshold", 10.0);
  e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 14)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
