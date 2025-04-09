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
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

using namespace OpenMS;

START_TEST(MultiplexFilteredMSExperiment, "$Id$")

MultiplexFilteredMSExperiment* nullPointer = nullptr;
MultiplexFilteredMSExperiment* ptr;

START_SECTION(MultiplexFilteredMSExperiment())
    MultiplexFilteredMSExperiment exp;
    TEST_EQUAL(exp.size(), 0);
    ptr = new MultiplexFilteredMSExperiment();
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexFilteredMSExperiment exp;
MultiplexFilteredPeak peak(654.32, 2345.67, 24, 110);
exp.addPeak(peak);
size_t n;

START_SECTION(addPeak(const MultiplexFilteredPeak& peak))
  n = exp.size();
  MultiplexFilteredPeak peak_temp(655.32, 2346.67, 25, 111);
  exp.addPeak(peak_temp);
  TEST_EQUAL(exp.size(), n + 1);
END_SECTION

START_SECTION(MultiplexFilteredPeak getPeak(size_t i))
  MultiplexFilteredPeak peak = exp.getPeak(0);
  TEST_REAL_SIMILAR(peak.getMZ(), 654.32);
END_SECTION

START_SECTION(double getMZ(size_t i))
  TEST_REAL_SIMILAR(exp.getMZ(0), 654.32);
END_SECTION

START_SECTION(std::vector<double> getMZ())
  TEST_REAL_SIMILAR(exp.getMZ()[0], 654.32);
END_SECTION

START_SECTION(double getRT(size_t i))
  TEST_REAL_SIMILAR(exp.getRT(0), 2345.67);
END_SECTION

START_SECTION(std::vector<double> getRT())
  TEST_REAL_SIMILAR(exp.getRT()[0], 2345.67);
END_SECTION

START_SECTION(size_t size())
  TEST_EQUAL(exp.size(), 2);
END_SECTION

END_TEST
