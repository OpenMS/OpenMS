// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Christian Ehrlich, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>
///////////////////////////

using namespace OpenMS;
using namespace std;
using namespace Math;

START_TEST(QuadraticRegression, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QuadraticRegression* ptr = nullptr;
QuadraticRegression* null_ptr = nullptr;
START_SECTION(QuadraticRegression())
{
	ptr = new QuadraticRegression();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QuadraticRegression())
{
	delete ptr;
}
END_SECTION

  // Create a test data set
vector<double> x_axis(10);
vector<double> y_axis(10);
vector<double> y_axis0(10);
vector<double> weight(10);
for (int i=0; i < 10; ++i)
{
  x_axis[i] = i;
  y_axis[i] = 5.5*i*i + 2*i + 4;
  y_axis0[i] = 5.5*i*i + 2*i; // no intercept
  weight[i]=1+i;
}

QuadraticRegression q_reg, q_reg2;

START_SECTION((template < typename Iterator > void computeRegression(Iterator x_begin, Iterator x_end, Iterator y_begin)))
{
  q_reg.computeRegression(x_axis.begin(), x_axis.end(), y_axis.begin());
  TEST_REAL_SIMILAR(q_reg.getA(), 4.0)
  TEST_REAL_SIMILAR(q_reg.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg.getChiSquared(), 0.0)

  q_reg2.computeRegression(x_axis.begin(), x_axis.end(), y_axis0.begin());
  TEST_REAL_SIMILAR(q_reg2.getA(), 0.0)
  TEST_REAL_SIMILAR(q_reg2.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg2.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg2.getChiSquared(), 0.0)
}
END_SECTION

START_SECTION((template < typename Iterator > void computeRegressionWeighted(Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
{
  q_reg.computeRegressionWeighted(x_axis.begin(), x_axis.end(), y_axis.begin(), weight.begin());
  TEST_REAL_SIMILAR(q_reg.getA(), 4.0)
  TEST_REAL_SIMILAR(q_reg.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg.getChiSquared(), 0.0)

  q_reg2.computeRegressionWeighted(x_axis.begin(), x_axis.end(), y_axis0.begin(), weight.begin());
  TEST_REAL_SIMILAR(q_reg2.getA(), 0.0)
  TEST_REAL_SIMILAR(q_reg2.getB(), 2.0)
  TEST_REAL_SIMILAR(q_reg2.getC(), 5.5)
  TEST_REAL_SIMILAR(q_reg2.getChiSquared(), 0.0)
}
END_SECTION

START_SECTION((double eval(double x) const ))
{
  double x = 100.0;
  TEST_REAL_SIMILAR(q_reg.eval(x), x*x*5.5 + x*2 + 4)
}
END_SECTION

START_SECTION(static double eval(double A, double B, double C, double x))
{
  double x = 100.0;
  TEST_REAL_SIMILAR(QuadraticRegression::eval(4.0, 2.0, 5.5, x), x*x*5.5 + x*2 + 4)
}
END_SECTION

START_SECTION((double getA() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getB() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getC() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((double getChiSquared() const ))
{
  NOT_TESTABLE // tested above
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



