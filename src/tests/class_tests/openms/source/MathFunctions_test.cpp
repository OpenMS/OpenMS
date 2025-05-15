// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// /////////////////////////

#include <OpenMS/MATH/MathFunctions.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

// /////////////////////////

START_TEST(Distribution, "$Id$")

// ///////////////////////////////////////////////////////////

START_SECTION((std::pair<double, double> zoomIn(const double left, const double right, const float factor, const float align)))
{
  {
    auto r = zoomIn(10, 20, 0.5, 0);
    TEST_REAL_SIMILAR(r.first, 10)
    TEST_REAL_SIMILAR(r.second, 15)
  }
  {
    auto r = zoomIn(10, 20, 0.5, 1);
    TEST_REAL_SIMILAR(r.first, 15)
    TEST_REAL_SIMILAR(r.second, 20)
  }
  {
    auto r = zoomIn(10, 20, 0.5, 0.5);
    TEST_REAL_SIMILAR(r.first, 12.5)
    TEST_REAL_SIMILAR(r.second, 17.5)
  }
  {
    auto r = zoomIn(10, 20, 2, 0);
    TEST_REAL_SIMILAR(r.first, 10)
    TEST_REAL_SIMILAR(r.second, 30)
  }
  {
    auto r = zoomIn(10, 20, 2, 0.5);
    TEST_REAL_SIMILAR(r.first, 5)
    TEST_REAL_SIMILAR(r.second, 25)
  }
  {
    auto r = zoomIn(10, 20, 2, 1);
    TEST_REAL_SIMILAR(r.first, 00)
    TEST_REAL_SIMILAR(r.second, 20)
  }
  // test round trip
  {
    auto r = zoomIn(10, 20, 2, 1);
    auto r2 = zoomIn(r.first, r.second, 0.5, 1);
    TEST_REAL_SIMILAR(r2.first, 10)
    TEST_REAL_SIMILAR(r2.second, 20)
  }
  // test round trip
  {
    auto r = zoomIn(10, 20, 2, 0);
    auto r2 = zoomIn(r.first, r.second, 0.5, 0);
    TEST_REAL_SIMILAR(r2.first, 10)
    TEST_REAL_SIMILAR(r2.second, 20)
  }
}
END_SECTION

START_SECTION((ceilDecimal))
	TEST_REAL_SIMILAR(ceilDecimal(12345.671,-2),12345.68)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,2),12400.0)
END_SECTION

START_SECTION((roundDecimal))
	TEST_REAL_SIMILAR(roundDecimal(12345.671,-2),12345.67)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,2),12300.0)
END_SECTION

START_SECTION((intervalTransformation))
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,600.0),300.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.25,1.0,0.0,600.0),200.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,0.75,0.0,600.0),400.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,150.0,600.0),375.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,450.0),225.0)
END_SECTION 

START_SECTION((linear2log))
	TEST_REAL_SIMILAR(linear2log(0.0),0.0)
	TEST_REAL_SIMILAR(linear2log(9.0),1.0)
	TEST_REAL_SIMILAR(linear2log(99.0),2.0)
	TEST_REAL_SIMILAR(linear2log(999.0),3.0)
END_SECTION

START_SECTION((log2linear))
	TEST_REAL_SIMILAR(log2linear(0.0),0.0)
	TEST_REAL_SIMILAR(log2linear(1.0),9.0)
	TEST_REAL_SIMILAR(log2linear(2.0),99.0)
	TEST_REAL_SIMILAR(log2linear(3.0),999.0)
END_SECTION

START_SECTION((isOdd))
	TEST_EQUAL(isOdd(0),false)
	TEST_EQUAL(isOdd(1),true)
	TEST_EQUAL(isOdd(2),false)
	TEST_EQUAL(isOdd(3),true)
END_SECTION

START_SECTION((template <typename T> T round (T x)))
	float f_down=14.49f;		 // expected 14
	float f_up = 14.50f;		 // expected 15
	double d_up = -999.49;   // expected -999
	double d_down = -675.77; // expected -676
	TEST_REAL_SIMILAR(round(f_down), 14.0)
	TEST_REAL_SIMILAR(round(f_up), 15.0)
	TEST_REAL_SIMILAR(round(d_up), -999)
	TEST_REAL_SIMILAR(round(d_down), -676)
END_SECTION


START_SECTION((bool approximatelyEqual(double a, double b, double tol)))
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.1), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.01), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.001), true)
	TEST_EQUAL(approximatelyEqual(1.1, 1.1002, 0.0001), false)
END_SECTION

START_SECTION((template <typename T> T getPPM(T mz_obs, T mz_ref)))
  TEST_REAL_SIMILAR(getPPM(1001.0, 1000.0), 1000.0)  // == 1 / 1000 * 1e6
  TEST_REAL_SIMILAR(getPPM( 999.0, 1000.0), -1000.0)  // == -1 / 1000 * 1e6
END_SECTION

START_SECTION((template <typename T> T getPPMAbs(T mz_obs, T mz_ref)))
  TEST_REAL_SIMILAR(getPPMAbs(1001.0, 1000.0), 1000.0)  // == abs(1 / 1000 * 1e6)
  TEST_REAL_SIMILAR(getPPMAbs( 999.0, 1000.0), 1000.0)  // == abs(-1 / 1000 * 1e6)
END_SECTION

START_SECTION((pair<double, double> getTolWindow(double val, double tol, bool ppm)))
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, true).first, 999.99)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, true).second, 1000.0100001)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, false).first, 990)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, false).second, 1010)
  TEST_REAL_SIMILAR(getTolWindow(500, 5, true).first, 499.9975)
  TEST_REAL_SIMILAR(getTolWindow(500, 5, true).second, 500.0025000125)
END_SECTION

START_SECTION((Math::RandomShuffle::portable_random_shuffle(BeginIT, EndIT)))
  vector<Size> seq{1,2,3,4,5,6};
  RandomShuffler r{0};
  r.portable_random_shuffle(seq.begin(),seq.end());
  TEST_EQUAL(seq[0],4)
  TEST_EQUAL(seq[1],3)
  TEST_EQUAL(seq[2],2)
  TEST_EQUAL(seq[3],6)
  TEST_EQUAL(seq[4],5)
  TEST_EQUAL(seq[5],1)
END_SECTION

/////////////////////////////////////////////////////////////);
/////////////////////////////////////////////////////////////
END_TEST
