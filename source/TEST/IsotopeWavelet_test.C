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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>

using namespace OpenMS;
using namespace std;


START_TEST(IsotopeWavelet, "$Id$")


START_SECTION((static IsotopeWavelet* getInstance()))
	TEST_EQUAL(IsotopeWavelet::getInstance(), 0)
END_SECTION


START_SECTION(static UInt getMaxCharge())
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 1)
END_SECTION

START_SECTION(static Int getGammaTableMaxIndex())
	TEST_EQUAL(IsotopeWavelet::getGammaTableMaxIndex(), -1)
END_SECTION

START_SECTION(static Int getExpTableMaxIndex())
	TEST_EQUAL(IsotopeWavelet::getExpTableMaxIndex(), -1)
END_SECTION

START_SECTION((static void setMaxCharge(const UInt max_charge))) 
	IsotopeWavelet::setMaxCharge(3);
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 3)
END_SECTION

START_SECTION((static DoubleReal getTableSteps()))
	TEST_NOT_EQUAL(IsotopeWavelet::getTableSteps(), 0)
END_SECTION

START_SECTION((static void setTableSteps(const DoubleReal table_steps))) 
	IsotopeWavelet::setTableSteps(0.0001);
	TEST_EQUAL(IsotopeWavelet::getTableSteps(), 0.0001)
END_SECTION

START_SECTION((static DoubleReal getInvTableSteps())) 
	IsotopeWavelet::getInvTableSteps();
	TEST_EQUAL(IsotopeWavelet::getInvTableSteps(), 10000)
END_SECTION

START_SECTION((static DoubleReal getLambdaL(const DoubleReal m)))
	TEST_REAL_SIMILAR(IsotopeWavelet::getLambdaL(1000), 0.69628)
END_SECTION


START_SECTION((static DoubleReal getLambdaQ(const DoubleReal m)))
	TEST_REAL_SIMILAR(IsotopeWavelet::getLambdaQ(1000), 0.685792)
END_SECTION


IsotopeWavelet* iw = 0;
START_SECTION((static IsotopeWavelet* init(const DoubleReal max_m, const UInt max_charge)))
	iw = IsotopeWavelet::init (4000, 4);
	TEST_NOT_EQUAL(iw, 0)
	TEST_EQUAL (IsotopeWavelet::getMaxCharge(), 4)
END_SECTION


UInt size=0;
START_SECTION((static const IsotopeDistribution::ContainerType& getAveragine (const DoubleReal m, UInt* size=NULL)))
	IsotopeWavelet::getAveragine (1000, &size);
	TEST_EQUAL (size, 3)	 
END_SECTION


DoubleReal v=-1;
START_SECTION((static DoubleReal getValueByMass (const DoubleReal t, const DoubleReal m, const UInt z, const Int mode=+1))) 
	TOLERANCE_ABSOLUTE (1e-4)
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByMass (HALF_NEUTRON_MASS/(c+1.), 1000, c+1, 1);
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION

START_SECTION((static DoubleReal getValueByLambda (const DoubleReal lambda, const DoubleReal tz1))) 
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambda (iw->getLambdaQ(1000*(c+1)-(c+1)*PROTON_MASS), HALF_NEUTRON_MASS*(c+1)+1);
		TOLERANCE_ABSOLUTE (1e-4)
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION

START_SECTION((static DoubleReal getValueByLambdaExtrapol (const DoubleReal lambda, const DoubleReal tz1))) 
	for (UInt c=0; c<iw->getMaxCharge(); ++c)
	{
		v=iw->getValueByLambdaExtrapol (iw->getLambdaQ(1000*(c+1)-(c+1)*PROTON_MASS), HALF_NEUTRON_MASS*(c+1)+1);
		TOLERANCE_ABSOLUTE (1e-4)
		TEST_REAL_SIMILAR(v, 0)
	};
END_SECTION


START_SECTION(static float myPow(float a, float b))
	TEST_EQUAL (trunc(IsotopeWavelet::myPow(1.1, 3)*10), 13);
END_SECTION


START_SECTION((virtual ~IsotopeWavelet()))
	delete (iw);
	TEST_EQUAL (IsotopeWavelet::getInstance(), 0)
	TEST_EQUAL (IsotopeWavelet::getMaxCharge(), 1)
END_SECTION



END_TEST
