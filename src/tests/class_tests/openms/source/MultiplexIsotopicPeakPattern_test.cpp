// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>
#include <OpenMS/FEATUREFINDER/MultiplexIsotopicPeakPattern.h>

using namespace OpenMS;

START_TEST(MultiplexIsotopicPeakPattern, "$Id$")

MultiplexDeltaMasses mass_shifts;
mass_shifts.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"no_label"));
mass_shifts.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(6.031817,"Arg6"));

MultiplexIsotopicPeakPattern* nullPointer = nullptr;
MultiplexIsotopicPeakPattern* ptr;

START_SECTION(MultiplexIsotopicPeakPattern(int c, int ppp, MultiplexDeltaMasses ms, int msi))
    MultiplexIsotopicPeakPattern pattern(2, 4, mass_shifts, 3);
    TEST_EQUAL(pattern.getCharge(), 2);
    ptr = new MultiplexIsotopicPeakPattern(2, 4, mass_shifts, 3);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexIsotopicPeakPattern pattern(2, 4, mass_shifts, 3);

START_SECTION(int getCharge() const)
  TEST_EQUAL(pattern.getCharge(), 2);
END_SECTION

START_SECTION(int getPeaksPerPeptide() const)
  TEST_EQUAL(pattern.getPeaksPerPeptide(), 4);
END_SECTION

START_SECTION(std::vector<double> getMassShifts() const)
  TEST_EQUAL(pattern.getMassShifts().getDeltaMasses()[0].delta_mass, 0);
  TEST_EQUAL(pattern.getMassShifts().getDeltaMasses()[1].delta_mass, 6.031817);
END_SECTION

START_SECTION(int getMassShiftIndex() const)
  TEST_EQUAL(pattern.getMassShiftIndex(), 3);
END_SECTION

START_SECTION(unsigned getMassShiftCount() const)
  TEST_EQUAL(pattern.getMassShiftCount(), 2);
END_SECTION

START_SECTION(double getMassShiftAt(int i) const)
  TEST_EQUAL(pattern.getMassShiftAt(0), 0);
  TEST_EQUAL(pattern.getMassShiftAt(1), 6.031817);
END_SECTION

/*START_SECTION(double getMZShiftAt(int i) const)
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(0), -0.501677);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(1), 0);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(2), 0.501677);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(3), 1.00335);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(4), 1.50503);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(5), 2.51423);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(6), 3.01591);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(7), 3.51759);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(8), 4.01926);
  TEST_REAL_SIMILAR(pattern.getMZShiftAt(9), 4.52094);
END_SECTION

START_SECTION(unsigned getMZShiftCount() const)
  TEST_EQUAL(pattern.getMZShiftCount(), 10);
END_SECTION*/

END_TEST
