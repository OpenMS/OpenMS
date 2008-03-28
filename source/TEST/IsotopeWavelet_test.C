// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>

using namespace OpenMS;
using namespace std;


START_TEST(IsotopeWavelet, "$Id$")


CHECK(getInstance)
	TEST_EQUAL(IsotopeWavelet::getInstance(), NULL)
RESULT


CHECK(getMaxCharge)
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 1)
RESULT


CHECK(setMaxCharge)
	IsotopeWavelet::setMaxCharge(3);
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 3)
RESULT


CHECK(setTableSteps) //includes the test for getTableSteps
	IsotopeWavelet::setTableSteps(0.0001);
	TEST_EQUAL(IsotopeWavelet::getTableSteps(), 0.0001)
RESULT


CHECK(getLambdaL)
	TEST_REAL_EQUAL(IsotopeWavelet::getLambdaL(1000), 0.69628)
RESULT


CHECK(getLambdaQ)
	TEST_REAL_EQUAL(IsotopeWavelet::getLambdaQ(1000), 0.685792)
RESULT


IsotopeWavelet* iw = NULL;
CHECK(init)
	iw = IsotopeWavelet::init (4000, 4);
	TEST_NOT_EQUAL(iw, NULL)
	TEST_EQUAL (IsotopeWavelet::getMaxCharge(), 4)
RESULT


UInt size=-1;
CHECK(getAveragine)
	IsotopeWavelet::getAveragine (1000, &size);
	TEST_EQUAL (size, 3)	 
RESULT


DoubleReal v=-1;
CHECK(getValueByMass) 
	PRECISION (1e-4)
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByMass (HALF_NEUTRON_MASS/(c+1.), 1000, c+1, 1);
		TEST_REAL_EQUAL(v, 0)
	};
RESULT

CHECK(getValueByLambda) 
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambda (iw->getLambdaQ(1000*(c+1)-(c+1)*PROTON_MASS), HALF_NEUTRON_MASS*(c+1)+1);
		PRECISION (1e-4)
		TEST_REAL_EQUAL(v, 0)
	};
RESULT


CHECK(~IsotopeWavelet)
	delete (iw);
	TEST_EQUAL (IsotopeWavelet::getInstance(), NULL)
	TEST_EQUAL (IsotopeWavelet::getMaxCharge(), 1)
RESULT



END_TEST
