// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MathFunctions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace OpenMS::Math;

START_TEST(MathFunctions, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((double log_binomial_coef(unsigned n, unsigned k)))
{
  // Test normal cases
  TEST_REAL_SIMILAR(Math::log_binomial_coef(10, 5), 5.5294);
  TEST_REAL_SIMILAR(Math::log_binomial_coef(20, 10), 12.12679);
  
  // Test edge cases
  TEST_REAL_SIMILAR(Math::log_binomial_coef(5, 0), 0.0);
  TEST_REAL_SIMILAR(Math::log_binomial_coef(5, 5), 0.0);
  
  // Test symmetry property
  TEST_REAL_SIMILAR(Math::log_binomial_coef(10, 3), Math::log_binomial_coef(10, 7));
  
  // Test invalid input
  TEST_EXCEPTION(std::invalid_argument, Math::log_binomial_coef(5, 6));
}
END_SECTION

START_SECTION((double log_sum_exp(double x, double y)))
{
  // Test normal cases
  TEST_REAL_SIMILAR(Math::log_sum_exp(1.0, 2.0), 2.31326169);
  TEST_REAL_SIMILAR(Math::log_sum_exp(10.0, 10.0), 10.6931472);
  
  // Test with one very small value (should be close to the larger value)
  TEST_REAL_SIMILAR(Math::log_sum_exp(100.0, 0.0), 100.0);
  TEST_REAL_SIMILAR(Math::log_sum_exp(0.0, 100.0), 100.0);
  
  // Test with negative infinity
  double neg_inf = -std::numeric_limits<double>::infinity();
  TEST_REAL_SIMILAR(Math::log_sum_exp(neg_inf, 5.0), 5.0);
  TEST_REAL_SIMILAR(Math::log_sum_exp(5.0, neg_inf), 5.0);
}
END_SECTION
START_SECTION((ceilDecimal))
{
	TEST_REAL_SIMILAR(ceilDecimal(12345.671,-2),12345.68)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(ceilDecimal(12345.67,2),12400.0)
}  
END_SECTION

START_SECTION((roundDecimal))
{
	TEST_REAL_SIMILAR(roundDecimal(12345.671,-2),12345.67)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,-1),12345.7)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,0),12346.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,1),12350.0)
	TEST_REAL_SIMILAR(roundDecimal(12345.67,2),12300.0)
}
END_SECTION

START_SECTION((intervalTransformation))
{
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,600.0),300.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.25,1.0,0.0,600.0),200.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,0.75,0.0,600.0),400.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,150.0,600.0),375.0)
	TEST_REAL_SIMILAR(intervalTransformation(0.5,0.0,1.0,0.0,450.0),225.0)
}
END_SECTION 

START_SECTION((linear2log))
{
	TEST_REAL_SIMILAR(linear2log(0.0),0.0)
	TEST_REAL_SIMILAR(linear2log(9.0),1.0)
	TEST_REAL_SIMILAR(linear2log(99.0),2.0)
	TEST_REAL_SIMILAR(linear2log(999.0),3.0)
}  
END_SECTION

START_SECTION((log2linear))
{
	TEST_REAL_SIMILAR(log2linear(0.0),0.0)
	TEST_REAL_SIMILAR(log2linear(1.0),9.0)
	TEST_REAL_SIMILAR(log2linear(2.0),99.0)
	TEST_REAL_SIMILAR(log2linear(3.0),999.0)
}
END_SECTION

START_SECTION((isOdd))
{
	TEST_TRUE(!isOdd(0))
	TEST_TRUE(isOdd(1))
	TEST_TRUE(!isOdd(2))
	TEST_TRUE(isOdd(3))
}  
END_SECTION

START_SECTION((template <typename T> T round (T x)))
{
	float f_down=14.49f;		 // expected 14
	float f_up = 14.50f;		 // expected 15
	double d_up = -999.49;   // expected -999
	double d_down = -675.77; // expected -676
	TEST_REAL_SIMILAR(Math::round(f_down), 14.0)
  TEST_REAL_SIMILAR(Math::round(f_up), 15.0)
  TEST_REAL_SIMILAR(Math::round(d_up), -999)
	TEST_REAL_SIMILAR(Math::round(d_down), -676)
}
END_SECTION

START_SECTION(template<typename T> T roundTo(const T value, int digits))
{
  TEST_REAL_SIMILAR(roundTo(3.14159265, 2), 3.14)
  TEST_REAL_SIMILAR(roundTo(1234.9, -2), 1200)
  TEST_REAL_SIMILAR(roundTo(1234.9, 0), 1235)
  TEST_REAL_SIMILAR(roundTo(1234.9, -1), 1230)
  TEST_REAL_SIMILAR(roundTo(1234.9, -3), 1000)
}
END_SECTION

START_SECTION(template<typename T> double percentOf(T value, T total, int digits))
{
  TEST_REAL_SIMILAR(percentOf(1.0 / 3, 1.0, 2), 33.33)
  TEST_REAL_SIMILAR(percentOf(1.0 / 3, 1.0, 3), 33.333)
  TEST_REAL_SIMILAR(percentOf(1.0 / 3, 1.0, 4), 33.3333)
  
  TEST_REAL_SIMILAR(percentOf(166.6666, 1000.0, 1), 16.7)

  TEST_EXCEPTION(Exception::InvalidValue, percentOf(-1.0, 1000.0, 2))
  TEST_EXCEPTION(Exception::InvalidValue, percentOf(1.0, -1000.0, 2))
}
END_SECTION

START_SECTION((bool approximatelyEqual(double a, double b, double tol)))
{
	TEST_TRUE(approximatelyEqual(1.1, 1.1002, 0.1))
	TEST_TRUE(approximatelyEqual(1.1, 1.1002, 0.01))
	TEST_TRUE(approximatelyEqual(1.1, 1.1002, 0.001))
	TEST_FALSE(approximatelyEqual(1.1, 1.1002, 0.0001))
}
END_SECTION

START_SECTION((template <typename T> T getPPM(T mz_obs, T mz_ref)))
{
  TEST_REAL_SIMILAR(getPPM(1001.0, 1000.0), 1000.0)  // == 1 / 1000 * 1e6
  TEST_REAL_SIMILAR(getPPM( 999.0, 1000.0), -1000.0)  // == -1 / 1000 * 1e6
}
END_SECTION

START_SECTION((template <typename T> T getPPMAbs(T mz_obs, T mz_ref)))
{
  TEST_REAL_SIMILAR(getPPMAbs(1001.0, 1000.0), 1000.0)  // == abs(1 / 1000 * 1e6)
  TEST_REAL_SIMILAR(getPPMAbs( 999.0, 1000.0), 1000.0)  // == abs(-1 / 1000 * 1e6)
}
END_SECTION

START_SECTION((pair<double, double> getTolWindow(double val, double tol, bool ppm)))
{
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, true).first, 999.99)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, true).second, 1000.0100001)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, false).first, 990)
  TEST_REAL_SIMILAR(getTolWindow(1000, 10, false).second, 1010)
  TEST_REAL_SIMILAR(getTolWindow(500, 5, true).first, 499.9975)
  TEST_REAL_SIMILAR(getTolWindow(500, 5, true).second, 500.0025000125)
}  
END_SECTION

START_SECTION((double binomial_cdf_complement(unsigned N, unsigned n, double p)))
{
  // Test normal cases
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 5, 0.5), 0.623046875);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(20, 10, 0.4), 0.24466);
  
  // Test edge cases
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 0, 0.3), 1.0);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 10, 0.7), 0.121060821);
  
  // Test with p = 0 and p = 1
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 0, 0.0), 1.0);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 1, 0.0), 0.0);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 0, 1.0), 1.0);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(10, 10, 1.0), 1.0);
  
  // Test invalid inputs
  TEST_EXCEPTION(std::invalid_argument, Math::binomial_cdf_complement(10, 11, 0.5));
  TEST_EXCEPTION(std::invalid_argument, Math::binomial_cdf_complement(10, 5, -0.1));
  TEST_EXCEPTION(std::invalid_argument, Math::binomial_cdf_complement(10, 5, 1.1));
  
  // Test the function used in AScore
  // This is similar to the test in AScore_test.cpp for computeCumulativeScoreTest_
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(1, 1, 0.1), 0.1);
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(3, 1, 0.1), 0.271);
  
  // Test with larger values
  TEST_REAL_SIMILAR(Math::binomial_cdf_complement(100, 60, 0.5), 0.02844);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

