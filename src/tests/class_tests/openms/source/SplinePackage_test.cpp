// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/PROCESSING/MISC/SplinePackage.h>

using namespace OpenMS;

START_TEST(SplinePackage, "$Id$")

std::vector<double> mz;
mz.push_back(413.8);
mz.push_back(413.9);
mz.push_back(414.0);
mz.push_back(414.1);
mz.push_back(414.2);
std::vector<double> intensity;
intensity.push_back(0.0);
intensity.push_back(100.2);
intensity.push_back(20.3);
intensity.push_back(2000.4);
intensity.push_back(4.3);

std::vector<double> mz1;
mz1.push_back(413.9);
std::vector<double> intensity1;
intensity1.push_back(100.2);

std::vector<double> mz2;
mz2.push_back(413.8);
mz2.push_back(413.9);
std::vector<double> intensity2;
intensity2.push_back(0.0);
intensity2.push_back(100.2);

SplinePackage sp1(mz, intensity);

SplinePackage* nullPointer = nullptr;

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity))
  SplinePackage* sp2 = new SplinePackage(mz, intensity);
  TEST_NOT_EQUAL(sp2, nullPointer)
  delete sp2;
END_SECTION

START_SECTION(getPosMin())
  TEST_EQUAL(sp1.getPosMin(), 413.8);
END_SECTION

START_SECTION(getPosMax())
  TEST_EQUAL(sp1.getPosMax(), 414.2);
END_SECTION

START_SECTION(getPosStepWidth())
  TEST_REAL_SIMILAR(sp1.getPosStepWidth(), 0.1);
END_SECTION

START_SECTION(isInPackage(double mz))
  TEST_EQUAL(sp1.isInPackage(414.05), true);
END_SECTION

START_SECTION(eval(double mz))
  TEST_REAL_SIMILAR(sp1.eval(414.05), 1134.08593750018);
END_SECTION

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity))
  TEST_EXCEPTION(Exception::IllegalArgument, SplinePackage(mz1, intensity1));
END_SECTION

START_SECTION(SplinePackage(std::vector<double> mz, std::vector<double> intensity))
  SplinePackage* sp4 = new SplinePackage(mz2, intensity2);
  TEST_NOT_EQUAL(sp4, nullPointer);
  TEST_REAL_SIMILAR((*sp4).eval(413.85), 50.1);
  delete sp4;
END_SECTION

END_TEST
