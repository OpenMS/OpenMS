// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ML/REGRESSION/LinearRegressionWithoutIntercept.h>

#include <cmath>
#include <vector>

///////////////////////////

START_TEST(LinearRegressionWithoutIntercept, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

LinearRegressionWithoutIntercept* ptr = nullptr;
LinearRegressionWithoutIntercept* null_ptr = nullptr;

START_SECTION((LinearRegressionWithoutIntercept()))
{
  ptr = new LinearRegressionWithoutIntercept();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // no data yet -> slope is NaN
  TEST_EQUAL(std::isnan(ptr->getSlope()), true)
}
END_SECTION

START_SECTION((~LinearRegressionWithoutIntercept()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void addData(double x, double y)))
{
  LinearRegressionWithoutIntercept lr;
  // a single observation is not enough -> NaN
  lr.addData(1.0, 2.0);
  TEST_EQUAL(std::isnan(lr.getSlope()), true)

  // perfectly collinear through the origin: y = 2x
  lr.addData(2.0, 4.0);
  lr.addData(3.0, 6.0);
  // slope = sum(x*y) / sum(x*x) = 28 / 14 = 2
  TEST_REAL_SIMILAR(lr.getSlope(), 2.0)
}
END_SECTION

START_SECTION((void addData(std::vector<double>& x, std::vector<double>& y)))
{
  LinearRegressionWithoutIntercept lr;
  std::vector<double> x{1.0, 2.0};
  std::vector<double> y{2.1, 3.9};
  lr.addData(x, y);
  // slope = (1*2.1 + 2*3.9) / (1*1 + 2*2) = 9.9 / 5 = 1.98
  TEST_REAL_SIMILAR(lr.getSlope(), 1.98)
}
END_SECTION

START_SECTION((double getSlope() const))
{
  // fewer than two observations -> NaN
  LinearRegressionWithoutIntercept lr0;
  TEST_EQUAL(std::isnan(lr0.getSlope()), true)

  LinearRegressionWithoutIntercept lr1;
  lr1.addData(5.0, 10.0);
  TEST_EQUAL(std::isnan(lr1.getSlope()), true)

  // two observations, scaled line y = 3x
  LinearRegressionWithoutIntercept lr2;
  lr2.addData(1.0, 3.0);
  lr2.addData(2.0, 6.0);
  TEST_REAL_SIMILAR(lr2.getSlope(), 3.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
