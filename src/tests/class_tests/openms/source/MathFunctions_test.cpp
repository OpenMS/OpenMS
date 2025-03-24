// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Benjamin Sachsenhofer $
// $Authors: Benjamin Sachsenhofer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MathFunctions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

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
