// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/QC/TIC.h>

///////////////////////////

START_TEST(TIC, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

TIC* ptr = nullptr;
TIC* nullPointer = nullptr;
START_SECTION(TIC())
ptr = new TIC();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~TIC())
delete ptr;
END_SECTION

TIC tic;
START_SECTION(const String& getName() const override)
TEST_EQUAL(tic.getName(), "TIC")
END_SECTION

START_SECTION(Status requirements() const override)
TEST_EQUAL((tic.requirements() == QCBase::Status(QCBase::Requires::RAWMZML)), true);
END_SECTION

START_SECTION(void compute(const MSExperiment& exp, float bin_size))
// build an experiment with a few MS1 spectra to verify the computed TIC metrics
MSExperiment exp;
{
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(0.5);
  spec.emplace_back(100.0, 30.0f);
  spec.emplace_back(150.0, 70.0f);
  exp.addSpectrum(spec);
}
{
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(1.0);
  spec.emplace_back(200.0, 500.0f);
  spec.emplace_back(250.0, 2000.0f);
  exp.addSpectrum(spec);
}
{
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(1.5);
  spec.emplace_back(300.0, 5.0f);
  spec.emplace_back(350.0, 15.0f);
  exp.addSpectrum(spec);
}
{
  MSSpectrum spec;
  spec.setMSLevel(1);
  spec.setRT(2.0);
  spec.emplace_back(400.0, 1000.0f);
  spec.emplace_back(450.0, 3000.0f);
  exp.addSpectrum(spec);
}

TIC::Result result = tic.compute(exp, 0);

TEST_EQUAL(result.intensities.size(), 4);
TEST_EQUAL(result.retention_times.size(), 4);
TEST_EQUAL(result.relative_intensities.size(), 4);

TEST_EQUAL(result.intensities[0], 100);
TEST_EQUAL(result.intensities[1], 2500);
TEST_EQUAL(result.intensities[2], 20);
TEST_EQUAL(result.intensities[3], 4000);

TEST_REAL_SIMILAR(result.relative_intensities[0], 2.5);
TEST_REAL_SIMILAR(result.relative_intensities[1], 62.5);
TEST_REAL_SIMILAR(result.relative_intensities[2], 0.5);
TEST_REAL_SIMILAR(result.relative_intensities[3], 100.0);

TEST_REAL_SIMILAR(result.retention_times[0], 0.5);
TEST_REAL_SIMILAR(result.retention_times[1], 1.0);
TEST_REAL_SIMILAR(result.retention_times[2], 1.5);
TEST_REAL_SIMILAR(result.retention_times[3], 2.0);

TEST_EQUAL(result.area, 6620);
TEST_EQUAL(result.jump, 2);
TEST_EQUAL(result.fall, 1);

// empty experiment still yields an empty result
MSExperiment empty_exp;
TEST_EQUAL(tic.compute(empty_exp, 0) == TIC::Result(), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
