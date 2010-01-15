// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

///////////////////////////

START_TEST(LinearRegression<Iterator>, "$Id: LinearRegression_test.C 4776 2009-03-05 14:14:35Z groepl $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace Math;
using namespace std;

LinearRegression* ptr;
START_SECTION((LinearRegression()))
  ptr = new LinearRegression;
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~LinearRegression()))
  delete ptr;
END_SECTION

// Create a test data set
vector<double> x_axis(10);
vector<double> y_axis(10);
vector<double> weight(10);
for (int i=0; i < 10; ++i)
{
  x_axis[i]=i;
  y_axis[i]=2*i+4;
  weight[i]=1+i;
}

LinearRegression lin_reg;

START_SECTION((template < typename Iterator > void computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)))
  lin_reg.computeRegression(0.95,x_axis.begin(),x_axis.end(),y_axis.begin());
  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),4.0)
END_SECTION

START_SECTION((template < typename Iterator > void computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
  lin_reg.computeRegressionWeighted(0.95,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());
  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),4.0)
END_SECTION

START_SECTION((DoubleReal getChiSquared() const))
  TEST_REAL_SIMILAR(lin_reg.getChiSquared(),0)
END_SECTION

START_SECTION((DoubleReal getIntercept() const))
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),4.0)
END_SECTION

START_SECTION((DoubleReal getLower() const))
  TEST_REAL_SIMILAR(lin_reg.getLower(),-2.0)
END_SECTION

START_SECTION((DoubleReal getUpper() const))
  TEST_REAL_SIMILAR(lin_reg.getUpper(),-2.0)
END_SECTION

START_SECTION((DoubleReal getSlope() const))
  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
END_SECTION

START_SECTION((DoubleReal getStandDevRes() const))
  TEST_REAL_SIMILAR(lin_reg.getStandDevRes(),0.0)
END_SECTION

START_SECTION((DoubleReal getStandErrSlope() const))
  TEST_REAL_SIMILAR(lin_reg.getStandErrSlope(),0.0)
END_SECTION

START_SECTION((DoubleReal getRSquared() const))
  TEST_REAL_SIMILAR(lin_reg.getRSquared(),1.0)
END_SECTION

START_SECTION((DoubleReal getTValue() const))
  TEST_REAL_SIMILAR(lin_reg.getTValue(),2.306)
END_SECTION

START_SECTION((DoubleReal getXIntercept() const))
  TEST_REAL_SIMILAR(lin_reg.getXIntercept(),-2.0)
END_SECTION

START_SECTION((DoubleReal getRSD() const))
  TEST_REAL_SIMILAR(lin_reg.getRSD(),0.0)
END_SECTION

START_SECTION((DoubleReal getMeanRes() const))
  TEST_REAL_SIMILAR(lin_reg.getMeanRes(),0.0)
END_SECTION

//test with no intercept
for (int i=0; i < 10; ++i)
{
  y_axis[i]=2*i;
}

START_SECTION((template < typename Iterator > void computeRegressionNoIntercept(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)))

  lin_reg.computeRegressionNoIntercept(0.95,x_axis.begin(),x_axis.end(),y_axis.begin());

  TEST_REAL_SIMILAR(lin_reg.getSlope(),2.0)
  TEST_REAL_SIMILAR(lin_reg.getIntercept(),0.0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
