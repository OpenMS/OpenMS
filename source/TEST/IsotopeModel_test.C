// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>


///////////////////////////

START_TEST(IsotopeModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
IsotopeModel* ptr = 0;
CHECK(IsotopeModel())
	ptr = new IsotopeModel();
  TEST_EQUAL(ptr->getName(), "IsotopeModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~IsotopeModel())
	delete ptr;
RESULT


CHECK(const String getName())
	TEST_EQUAL(IsotopeModel::getProductName(),"IsotopeModel")
	TEST_EQUAL(IsotopeModel().getName(),"IsotopeModel")
RESULT

// assignment operator
CHECK(IsotopeModel& operator = (const IsotopeModel& source))
	IsotopeModel im1;
	im1.setParam(670.5, 3, 0.8);

  IsotopeModel im2;
  im2 = im1;

  IsotopeModel im3;
	im3.setParam(670.5, 3, 0.8);

  im1 = IsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
RESULT

// copy ctor
CHECK(IsotopeModel(const IsotopeModel& source))
	IsotopeModel im1;
	im1.setParam(670.5, 3, 0.8);

	IsotopeModel im2(im1);
  IsotopeModel im3;
	im3.setParam(670.5, 3, 0.8);

  im1 = IsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
RESULT

CHECK([EXTRA] DefaultParamHandler::setParameters(...))
	PRECISION(0.001)
	IsotopeModel im1;
	im1.setParam(670.5, 3, 0.8);

	IsotopeModel im2;
	im2.setParameters(im1.getParameters());

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	im1.getSamples(dpa1);
	im2.getSamples(dpa2);

	PRECISION(0.00001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK(void setParam(CoordinateType, UnsignedInt, CoordinateType))

	// Isotope distribution with mean at mz=1000.62094 and charge=2
	IsotopeModel im1;
	im1.setInterpolationStep(0.01);
	im1.setScalingFactor(25.068);
	im1.setParam(1000.62094, 2, 0.1);

	double dist =  1.000495/2;
	double mono = 1000.035168;

	PRECISION(0.001)
	TEST_REAL_EQUAL(im1.getCenter(), mono)

	PRECISION(0.1)
	TEST_REAL_EQUAL(im1.getIntensity(mono), 33.6352)
	TEST_REAL_EQUAL(im1.getIntensity(mono+dist), 33.8107)
	TEST_REAL_EQUAL(im1.getIntensity(mono+2*dist), 19.858)
	TEST_REAL_EQUAL(im1.getIntensity(mono+3*dist), 8.57358)
	TEST_REAL_EQUAL(im1.getIntensity(mono+4*dist), 2.9971)
	TEST_REAL_EQUAL(im1.getIntensity(mono+5*dist), 0.892128)
	TEST_REAL_EQUAL(im1.getIntensity(mono+6*dist), 0.233331)

	im1.setInterpolationStep(0.1);
	im1.setSamples();

	TEST_REAL_EQUAL(im1.getIntensity(mono), 33.6352)
	TEST_REAL_EQUAL(im1.getIntensity(mono+dist), 33.8107)
	TEST_REAL_EQUAL(im1.getIntensity(mono+2*dist), 19.858)
	TEST_REAL_EQUAL(im1.getIntensity(mono+3*dist), 8.57358)
	TEST_REAL_EQUAL(im1.getIntensity(mono+4*dist), 2.9971)
	TEST_REAL_EQUAL(im1.getIntensity(mono+5*dist), 0.892128)
	TEST_REAL_EQUAL(im1.getIntensity(mono+6*dist), 0.233331)
RESULT


CHECK(void setOffset(double offset))
	PRECISION(0.001)
	IsotopeModel im1;
	im1.setParam(670.5, 3, 0.8);
	im1.setOffset( im1.getOffset()+1.8 );

	IsotopeModel im2;
	im2.setParam(672.3, 3, 0.8);

	TEST_EQUAL(im1.getParameters(), im2.getParameters())
	TEST_REAL_EQUAL(im1.getCenter(), im2.getCenter())
	TEST_REAL_EQUAL(im1.getCenter(), 671.909)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	im1.getSamples(dpa1);
	im2.getSamples(dpa2);

	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
