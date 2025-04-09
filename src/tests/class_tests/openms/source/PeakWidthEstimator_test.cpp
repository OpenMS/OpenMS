// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;

START_TEST(PeakWidthEstimator, "$Id$")

PeakMap exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"), exp);

PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);

std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;

PeakMap exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

PeakWidthEstimator* nullPointer = nullptr;
PeakWidthEstimator* ptr;

START_SECTION(PeakWidthEstimator(const PeakMap & exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > & boundaries))
{
  PeakWidthEstimator estimator(exp_picked, boundaries_exp_s);
  TEST_REAL_SIMILAR(estimator.getPeakWidth(365.3),0.00886469661896705);
  ptr = new PeakWidthEstimator(exp_picked, boundaries_exp_s);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
}
END_SECTION

PeakWidthEstimator estimator2(exp_picked, boundaries_exp_s);

START_SECTION(double getPeakWidth(double mz))
{
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(365.3),0.00886469661896705);
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(305.1),0.00886699447290451);    // outside m/z range
  TEST_REAL_SIMILAR(estimator2.getPeakWidth(405.1),0.01184458329884600);    // outside m/z range
}
END_SECTION

END_TEST

