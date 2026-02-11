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
#include <OpenMS/IONMOBILITY/IMTypes.h>

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

START_SECTION((std::set<double> getCompensationVoltages() detects FAIMS beyond first spectrum and ignores sentinel))
{
  PeakMap exp2;
  // spectrum without FAIMS (first)
  MSSpectrum s0;
  s0.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  s0.setDriftTime(12.3);

  // FAIMS spectrum with valid CV
  MSSpectrum s1;
  s1.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  s1.setDriftTime(-50.0);

  // FAIMS spectrum with sentinel (should be ignored)
  MSSpectrum s2;
  s2.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
  s2.setDriftTime(IMTypes::DRIFTTIME_NOT_SET);

  exp2.addSpectrum(s0);
  exp2.addSpectrum(s1);
  exp2.addSpectrum(s2);

  const std::set<double> cvs2 = FAIMSHelper::getCompensationVoltages(exp2);
  TEST_EQUAL(cvs2.size(), 1)
  TEST_TRUE(cvs2.find(-50.0) != cvs2.end())
}
END_SECTION

START_SECTION((std::set<double> getCompensationVoltages() returns empty for non-FAIMS))
{
  PeakMap exp3;
  MSSpectrum a;
  a.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  a.setDriftTime(1.0);
  MSSpectrum b;
  b.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  b.setDriftTime(2.0);
  exp3.addSpectrum(a);
  exp3.addSpectrum(b);

  const std::set<double> cvs3 = FAIMSHelper::getCompensationVoltages(exp3);
  TEST_EQUAL(cvs3.empty(), true)
}
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
