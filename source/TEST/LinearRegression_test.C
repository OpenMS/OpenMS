// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

///////////////////////////

START_TEST(LinearRegression<Iterator>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

LinearRegression* ptr;
CHECK(LinearRegression())
  ptr = new LinearRegression;
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~LinearRegression()))
  delete ptr;
RESULT

// Create a test data set
vector<double> x_axis(10);
vector<double> y_axis(10);
vector<double> weight(10);
for (int i=0; i < 10; ++i)
{
  x_axis[i]=i;
  y_axis[i]=2*i+4;
  weight[i]=1;
}

LinearRegression lin_reg;

CHECK((bool computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)))
  lin_reg.computeRegression(0.95,x_axis.begin(),x_axis.end(),y_axis.begin());

  TEST_REAL_EQUAL(lin_reg.getSlope(),2.0)
  TEST_REAL_EQUAL(lin_reg.getIntercept(),4.0)
RESULT

CHECK((int computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
  lin_reg.computeRegressionWeighted(0.95,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());
RESULT

CHECK((DoubleReal getChiSquared() const))
  TEST_REAL_EQUAL(lin_reg.getChiSquared(),0)
RESULT

CHECK((DoubleReal getIntercept() const))
  TEST_REAL_EQUAL(lin_reg.getIntercept(),4.0)
RESULT

CHECK((DoubleReal getLower() const))
  TEST_REAL_EQUAL(lin_reg.getLower(),-2.0)
RESULT

CHECK((DoubleReal getUpper() const))
  TEST_REAL_EQUAL(lin_reg.getUpper(),-2.0)
RESULT

CHECK((DoubleReal getSlope() const))
  TEST_REAL_EQUAL(lin_reg.getSlope(),2.0)
RESULT

CHECK((DoubleReal getStandDevRes() const))
  TEST_REAL_EQUAL(lin_reg.getStandDevRes(),0.0)
RESULT

CHECK((DoubleReal getStandErrSlope() const))
  TEST_REAL_EQUAL(lin_reg.getStandErrSlope(),0.0)
RESULT

CHECK((DoubleReal getRSquared() const))
  TEST_REAL_EQUAL(lin_reg.getRSquared(),1.0)
RESULT

CHECK((DoubleReal getTValue() const))
  TEST_REAL_EQUAL(lin_reg.getTValue(),2.306)
RESULT

CHECK((DoubleReal getXIntercept() const))
  TEST_REAL_EQUAL(lin_reg.getXIntercept(),-2.0)
RESULT

CHECK((DoubleReal getRSD() const))
  TEST_REAL_EQUAL(lin_reg.getRSD(),0.0)
RESULT

CHECK((DoubleReal getMeanRes() const))
  TEST_REAL_EQUAL(lin_reg.getMeanRes(),0.0)
RESULT

//test with no intercept
for (int i=0; i < 10; ++i)
{
  y_axis[i]=2*i;
}
CHECK((template <typename Iterator> bool computeRegressionNoIntercept(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin);))
  lin_reg.computeRegressionNoIntercept(0.95,x_axis.begin(),x_axis.end(),y_axis.begin());

  TEST_REAL_EQUAL(lin_reg.getSlope(),2.0)
  TEST_REAL_EQUAL(lin_reg.getIntercept(),0.0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
