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

LinearRegression<vector<double>::const_iterator>* linreg_ptr;
CHECK((LinearRegression()))
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

CHECK((int computeRegression(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin)))
  double ci=0.95;
  int error = linreg_ptr->computeRegression(ci,x_axis.begin(),x_axis.end(),y_axis.begin());
  TEST_EQUAL(error,0)
RESULT

CHECK((int computeRegressionWeighted(double confidence_interval_P, Iterator x_begin, Iterator x_end, Iterator y_begin, Iterator w_begin)))
  double ci=0.95;
  int error = linreg_ptr->computeRegressionWeighted(ci,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());
  TEST_EQUAL(error,0);
RESULT


CHECK((LinearRegression( LinearRegression const & arg )))
	double ci=0.95;
  int error = linreg_ptr->computeRegressionWeighted(ci,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());

 	LinearRegression<vector<double>::const_iterator> linreg_copy(*linreg_ptr);

  TEST_REAL_EQUAL(error,linreg_copy.getStandErrSlope());
  TEST_REAL_EQUAL(linreg_ptr->getChiSquared(),linreg_copy.getChiSquared())
  TEST_REAL_EQUAL(linreg_ptr->getIntercept(),linreg_copy.getIntercept())
  TEST_REAL_EQUAL(linreg_ptr->getLower(),linreg_copy.getLower())
  TEST_REAL_EQUAL(linreg_ptr->getUpper(),linreg_copy.getUpper())
  TEST_REAL_EQUAL(linreg_ptr->getSlope(),linreg_copy.getSlope())
  TEST_REAL_EQUAL(linreg_ptr->getStandDevRes(),linreg_copy.getStandDevRes())
  TEST_REAL_EQUAL(linreg_ptr->getStandErrSlope(),linreg_copy.getStandErrSlope())
  TEST_REAL_EQUAL(linreg_ptr->getRSquared(),linreg_copy.getRSquared())
  TEST_REAL_EQUAL(linreg_ptr->getTValue(),linreg_copy.getTValue())
  TEST_REAL_EQUAL(linreg_ptr->getXIntercept(),linreg_copy.getXIntercept())
RESULT

CHECK((LinearRegression& operator=(LinearRegression const &arg)))
	double ci=0.95;
  int error = linreg_ptr->computeRegressionWeighted(ci,x_axis.begin(),x_axis.end(),y_axis.begin(),weight.begin());

 	LinearRegression<vector<double>::const_iterator> linreg_copy;
  linreg_copy = (*linreg_ptr);

  TEST_REAL_EQUAL(error,linreg_copy.getStandErrSlope());
  TEST_REAL_EQUAL(linreg_ptr->getChiSquared(),linreg_copy.getChiSquared())
  TEST_REAL_EQUAL(linreg_ptr->getIntercept(),linreg_copy.getIntercept())
  TEST_REAL_EQUAL(linreg_ptr->getLower(),linreg_copy.getLower())
  TEST_REAL_EQUAL(linreg_ptr->getUpper(),linreg_copy.getUpper())
  TEST_REAL_EQUAL(linreg_ptr->getSlope(),linreg_copy.getSlope())
  TEST_REAL_EQUAL(linreg_ptr->getStandDevRes(),linreg_copy.getStandDevRes())
  TEST_REAL_EQUAL(linreg_ptr->getStandErrSlope(),linreg_copy.getStandErrSlope())
  TEST_REAL_EQUAL(linreg_ptr->getRSquared(),linreg_copy.getRSquared())
  TEST_REAL_EQUAL(linreg_ptr->getTValue(),linreg_copy.getTValue())
  TEST_REAL_EQUAL(linreg_ptr->getXIntercept(),linreg_copy.getXIntercept())
RESULT



CHECK((DoubleReal getChiSquared() const))
  TEST_REAL_EQUAL(linreg_ptr->getChiSquared(),0)
RESULT

CHECK((DoubleReal getIntercept() const))
  TEST_REAL_EQUAL(linreg_ptr->getIntercept(),4.0)
RESULT

CHECK((DoubleReal getLower() const))
  TEST_REAL_EQUAL(linreg_ptr->getLower(),-2.0)
RESULT

CHECK((DoubleReal getUpper() const))
  TEST_REAL_EQUAL(linreg_ptr->getUpper(),-2.0)
RESULT

CHECK((DoubleReal getSlope() const))
  TEST_REAL_EQUAL(linreg_ptr->getSlope(),2.0)
RESULT

CHECK((DoubleReal getStandDevRes() const))
  TEST_REAL_EQUAL(linreg_ptr->getStandDevRes(),0.0)
RESULT

CHECK((DoubleReal getStandErrSlope() const))
  TEST_REAL_EQUAL(linreg_ptr->getStandErrSlope(),0.0)
RESULT

CHECK((DoubleReal getRSquared() const))
  TEST_REAL_EQUAL(linreg_ptr->getRSquared(),1.0)
RESULT

CHECK((DoubleReal getTValue() const))
  TEST_REAL_EQUAL(linreg_ptr->getTValue(),2.306)
RESULT

CHECK((DoubleReal getXIntercept() const))
  TEST_REAL_EQUAL(linreg_ptr->getXIntercept(),-2.0)
RESULT

CHECK((DoubleReal getRSD() const))
  TEST_REAL_EQUAL(linreg_ptr->getRSD(),0.0)
RESULT

CHECK((DoubleReal getMeanRes() const))
  TEST_REAL_EQUAL(linreg_ptr->getMeanRes(),0.0)
RESULT

CHECK((virtual ~LinearRegression()))
  delete linreg_ptr;
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
