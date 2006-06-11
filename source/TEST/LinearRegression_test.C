// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: LinearRegression_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Eva Lange  $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

///////////////////////////

START_TEST(LinearRegression<Iterator>, "$Id: LinearRegression_test.C,v 1.3 2006/03/28 08:03:34 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

LinearRegression<vector<double>::const_iterator>* linreg_ptr;
CHECK(LinearRegression())
  linreg_ptr = new LinearRegression<vector<double>::const_iterator>;
  TEST_NOT_EQUAL(linreg_ptr, 0)
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

CHECK((void computeInterceptXAxis(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)))
  double ci=0.95;
  int error = linreg_ptr->computeInterceptXAxis(ci,x_axis.begin(),x_axis.end(),y_axis.begin());
  TEST_EQUAL(error,0)
RESULT

CHECK((void computeInterceptXAxisWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
  double ci=0.95;
  int error = linreg_ptr->computeInterceptXAxisWeighted(ci,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());
  TEST_EQUAL(error,0);
RESULT


CHECK(LinearRegression(LinearRegression const& linreg))
	double ci=0.95;
  LinearRegression<vector<double>::const_iterator> linreg(*linreg_ptr);
  int error = linreg.computeInterceptXAxisWeighted(ci,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());
  TEST_EQUAL(error,0);
RESULT


CHECK(const double& getChiSquared() const)
  TEST_REAL_EQUAL(linreg_ptr->getChiSquared(),0)
RESULT

CHECK(const double& getIntercept() const)
  TEST_REAL_EQUAL(linreg_ptr->getIntercept(),4.0)
RESULT

CHECK(const double& getLower() const)
  TEST_REAL_EQUAL(linreg_ptr->getLower(),-2.0)
RESULT

CHECK(const double& getUpper() const)
  TEST_REAL_EQUAL(linreg_ptr->getUpper(),-2.0)
RESULT

CHECK(const double& getSlope() const)
  TEST_REAL_EQUAL(linreg_ptr->getSlope(),2.0)
RESULT

CHECK(const double& getStandDevRes() const)
  TEST_REAL_EQUAL(linreg_ptr->getStandDevRes(),0.0)
RESULT

CHECK(const double& getStandErrSlope() const)
  TEST_REAL_EQUAL(linreg_ptr->getStandErrSlope(),0.0)
RESULT

CHECK(const double& getRSquared() const)
  TEST_REAL_EQUAL(linreg_ptr->getRSquared(),1.0)
RESULT

CHECK(const double& getTValue() const)
  TEST_REAL_EQUAL(linreg_ptr->getTValue(),2.306)
RESULT

CHECK(const double& getXIntercept() const)
  TEST_REAL_EQUAL(linreg_ptr->getXIntercept(),-2.0)
RESULT

CHECK(~LinearRegression())
  delete linreg_ptr;
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
