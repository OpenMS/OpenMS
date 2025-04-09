// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexFilteredPeak.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexSatelliteProfile.h>

using namespace OpenMS;

START_TEST(MultiplexFilteredPeak, "$Id$")

MultiplexFilteredPeak* nullPointer = nullptr;
MultiplexFilteredPeak* ptr;

START_SECTION(MultiplexFilteredPeak())
    MultiplexFilteredPeak peak(654.32, 2345.67, 24, 42);
    TEST_EQUAL(peak.getMZidx(), 24);
    ptr = new MultiplexFilteredPeak(654.32, 2345.67, 24, 42);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexFilteredPeak peak(654.32, 2345.67, 24, 42);
MultiplexSatelliteCentroided satellite_centroided(26, 44);
MultiplexSatelliteProfile satellite_profile(2346.67, 655.32, 1000.0);
peak.addSatellite(25, 43, 3);
peak.addSatellite(satellite_centroided, 3);
peak.addSatelliteProfile(2347.67, 656.32, 1010.0, 4);
peak.addSatelliteProfile(satellite_profile, 4);
size_t n;

START_SECTION(double getMZ())
  TEST_REAL_SIMILAR(peak.getMZ(), 654.32);
END_SECTION

START_SECTION(double getRT())
  TEST_REAL_SIMILAR(peak.getRT(), 2345.67);
END_SECTION

START_SECTION(size_t getMZidx())
  TEST_EQUAL(peak.getMZidx(), 24);
END_SECTION

START_SECTION(size_t getRTidx())
  TEST_EQUAL(peak.getRTidx(), 42);
END_SECTION

START_SECTION(void addSatellite(size_t rt_idx, size_t mz_idx, size_t pattern_idx))
  n = peak.getSatellites().size();
  peak.addSatellite(25, 43, 3);
  TEST_EQUAL(peak.getSatellites().size(), n+1);
END_SECTION

START_SECTION(void addSatellite(const MultiplexSatelliteCentroided& satellite, size_t pattern_idx))
  n = peak.getSatellites().size();
  MultiplexSatelliteCentroided satellite_centroided_temp(27, 45);
  peak.addSatellite(satellite_centroided_temp, 3);
  TEST_EQUAL(peak.getSatellites().size(), n+1);
END_SECTION

START_SECTION(void addSatelliteProfile(double rt, double mz, double intensity, size_t pattern_idx))
  n = peak.getSatellitesProfile().size();
  peak.addSatelliteProfile(2348.67, 657.32, 1020.0, 5);
  TEST_EQUAL(peak.getSatellitesProfile().size(), n+1);
END_SECTION

START_SECTION(void addSatelliteProfile(const MultiplexSatelliteProfile& satellite, size_t pattern_idx))
  n = peak.getSatellitesProfile().size();
  MultiplexSatelliteProfile satellite_profile_temp(2349.67, 658.32, 1030.0);
  peak.addSatelliteProfile(satellite_profile_temp, 6);
  TEST_EQUAL(peak.getSatellitesProfile().size(), n+1);
END_SECTION

START_SECTION(getSatellites())
  TEST_EQUAL(peak.getSatellites().size(), 4);
END_SECTION

START_SECTION(getSatellitesProfile())
  TEST_EQUAL(peak.getSatellitesProfile().size(), 4);
END_SECTION

START_SECTION(size_t size())
  TEST_EQUAL(peak.size(), 4);
END_SECTION

START_SECTION(size_t sizeProfile())
  TEST_EQUAL(peak.sizeProfile(), 4);
END_SECTION

END_TEST
