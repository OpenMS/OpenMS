// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexSatelliteProfile.h>

using namespace OpenMS;

START_TEST(MultiplexSatelliteProfile, "$Id$")

MultiplexSatelliteProfile* nullPointer = nullptr;
MultiplexSatelliteProfile* ptr;

START_SECTION(MultiplexSatelliteProfile(double rt, double mz, double intensity))
    MultiplexSatelliteProfile satellite(2565.3, 618.4, 1000.0);
    TEST_EQUAL(satellite.getMZ(), 618.4);
    ptr = new MultiplexSatelliteProfile(2565.3, 618.4, 1000.0);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(size_t getMZ())
  MultiplexSatelliteProfile satellite(2565.3, 618.4, 1000.0);
  TEST_REAL_SIMILAR(satellite.getMZ(), 618.4);
END_SECTION

START_SECTION(size_t getRT())
  MultiplexSatelliteProfile satellite(2565.3, 618.4, 1000.0);
  TEST_REAL_SIMILAR(satellite.getRT(), 2565.3);
END_SECTION

START_SECTION(size_t getIntensity())
  MultiplexSatelliteProfile satellite(2565.3, 618.4, 1000.0);
  TEST_REAL_SIMILAR(satellite.getIntensity(), 1000.0);
END_SECTION

END_TEST
