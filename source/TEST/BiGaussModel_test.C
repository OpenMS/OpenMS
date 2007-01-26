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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>


///////////////////////////

START_TEST(BiGaussModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using std::stringstream;

// default ctor
BiGaussModel* ptr = 0;
CHECK(BiGaussModel())
	ptr = new BiGaussModel();
  TEST_EQUAL(ptr->getName(), "BiGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~BiGaussModel())
	delete ptr;
RESULT


CHECK(const String getName())
	TEST_EQUAL(BiGaussModel::getProductName(),"BiGaussModel")
	TEST_EQUAL(BiGaussModel().getName(),"BiGaussModel")
RESULT

// assignment operator
CHECK(BiGaussModel& operator = (const BiGaussModel& source))
	BiGaussModel bgm1;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);
	bgm1.setParam(680.1, 2.0, 5.0, 678.9, 789.0);

  BiGaussModel bgm2;
  bgm2 = bgm1;

  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParam(680.1, 2.0, 5.0, 678.9, 789.0);

  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
RESULT

// copy ctor
CHECK(BiGaussModel(const BiGaussModel& source))
	BiGaussModel bgm1;
	BasicStatistics<>  stat;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);
	bgm1.setParam(680.1, 2.0, 5.0, 678.9, 789.0);

	BiGaussModel bgm2(bgm1);
  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParam(680.1, 2.0, 5.0, 678.9, 789.0);

  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
RESULT

CHECK([EXTRA] DefaultParamHandler::setParameters(...))
	PRECISION(0.001)
	BiGaussModel bgm1;
	bgm1.setParam(680.1, 2.0, 5.0, 678.9, 789.0);
	bgm1.setOffset(680.0);

	BiGaussModel bgm2;
	bgm2.setParameters(bgm1.getParameters());
	TEST_REAL_EQUAL(bgm1.getCenter(), 681.2)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	PRECISION(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK(void setParam(CoordinateType,CoordinateType,CoordinateType,CoordinateType,CoordinateType))
	//Bigauss1(x) = N(0,sigma1)*sigma1/area
	//Bigauss2(x) = N(0,sigma2)*sigma2/area
	//area = sigma1/2 + sigma2/2
	BiGaussModel bgm1;
	bgm1.setInterpolationStep(0.001);
	bgm1.setParam(0.0, 0.1, 1.0, -1.0, 4.0);
	TEST_REAL_EQUAL(bgm1.getCenter(), 0.0)

	PRECISION(0.001)
	TEST_REAL_EQUAL(bgm1.getIntensity(-0.1), 0.576626091);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.0), 0.606190343);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.1), 0.60316696);
	TEST_REAL_EQUAL(bgm1.getIntensity(1.0), 0.36767303);
	TEST_REAL_EQUAL(bgm1.getIntensity(2.0), 0.08203894);

	PRECISION(0.1)
	bgm1.setScalingFactor(10.0);
	bgm1.setSamples();
	TEST_REAL_EQUAL(bgm1.getIntensity(-0.1), 5.76626091);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.0), 6.06190343);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.1), 6.0316696);
	TEST_REAL_EQUAL(bgm1.getIntensity(1.0), 3.6767303);
	TEST_REAL_EQUAL(bgm1.getIntensity(2.0), 0.8203894);

	bgm1.setScalingFactor(1.0);
	bgm1.setInterpolationStep(0.2);
	bgm1.setSamples();

	TEST_REAL_EQUAL(bgm1.getIntensity(-0.1), 0.576626091);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.0), 0.606190343);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.1), 0.60316696);
	TEST_REAL_EQUAL(bgm1.getIntensity(1.0), 0.36767303);
	TEST_REAL_EQUAL(bgm1.getIntensity(2.0), 0.08203894);


	bgm1.setInterpolationStep(0.001);
	bgm1.setParam(0.0, 1.0, 1.0, -4.0, 4.0);
	TEST_REAL_EQUAL(bgm1.getCenter(), 0.0)

	PRECISION(0.001)
	TEST_REAL_EQUAL(bgm1.getIntensity(-1.0), 0.24197072);
	TEST_REAL_EQUAL(bgm1.getIntensity(0.0), 0.39894228);
	TEST_REAL_EQUAL(bgm1.getIntensity(1.0), 0.24197072);
	TEST_REAL_EQUAL(bgm1.getIntensity(2.0), 0.05399097);
	TEST_REAL_EQUAL(bgm1.getIntensity(-0.5), bgm1.getIntensity(0.5));
	TEST_REAL_EQUAL(bgm1.getIntensity(-2.0), bgm1.getIntensity(2.0));
RESULT


CHECK(void setOffset(double offset))
	BiGaussModel bgm1;
	bgm1.setParam(680.1, 2.0, 5.0, 678.9, 789.0);
	bgm1.setOffset(680.9);

	BiGaussModel bgm2;
	bgm2.setParam(682.1, 2.0, 5.0, 680.9, 791.0);

	TEST_EQUAL(bgm1.getParameters(), bgm2.getParameters())
	TEST_REAL_EQUAL(bgm1.getCenter(), bgm2.getCenter())
	TEST_REAL_EQUAL(bgm1.getCenter(), 682.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	PRECISION(0.001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

	bgm1.setParam(0.0, 0.81, 0.81, -4.0, 4.001);
	bgm1.setOffset(0.123);
	TEST_REAL_EQUAL(bgm1.getCenter(), 4.123)

	PRECISION(0.001)
	TEST_REAL_EQUAL(bgm1.getIntensity(4.123), 0.4432692);
	TEST_REAL_EQUAL(bgm1.getIntensity(4.223), bgm1.getIntensity(4.023));
	TEST_REAL_EQUAL(bgm1.getIntensity(3.123), bgm1.getIntensity(5.123));

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
