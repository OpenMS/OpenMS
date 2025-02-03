// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(FAIMSHelper, "$Id$")

/////////////////////////////////////////////////////////////

FAIMSHelper* e_ptr = nullptr;
FAIMSHelper* e_nullPointer = nullptr;

START_SECTION((FAIMSHelper()))
	e_ptr = new FAIMSHelper;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~FAIMSHelper()))
  delete e_ptr;
END_SECTION

e_ptr = new FAIMSHelper();

START_SECTION((std::set<double> getCompensationVoltages(PeakMap& exp)))
  delete e_ptr;
  e_ptr = new FAIMSHelper();
  MzMLFile IM_file;
  PeakMap exp;
  IM_file.load(OPENMS_GET_TEST_DATA_PATH("IM_FAIMS_test.mzML"), exp);

  TEST_EQUAL(exp.getSpectra().size(), 19)

  std::set<double> CVs = e_ptr->getCompensationVoltages(exp);
  
  TEST_EQUAL(CVs.size(), 3)
  TEST_EQUAL(CVs.find(-65.0) == CVs.end(), 0)
  TEST_EQUAL(CVs.find(-55.0) == CVs.end(), 0)
  TEST_EQUAL(CVs.find(-45.0) == CVs.end(), 0)

END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
