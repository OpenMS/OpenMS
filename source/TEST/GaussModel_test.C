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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>

///////////////////////////

START_TEST(GaussModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using std::stringstream;

// default ctor
GaussModel* ptr = 0;
CHECK(GaussModel())
	ptr = new GaussModel();
  TEST_EQUAL(ptr->getName(), "GaussModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~GaussModel())
	delete ptr;
RESULT


CHECK(const String getProductName())
	TEST_EQUAL(GaussModel::getProductName(),"GaussModel")
	TEST_EQUAL(GaussModel().getProductName(),"GaussModel")
RESULT

// assignment operator
CHECK(GaussModel& operator = (const GaussModel& source))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	gm1.setScalingFactor(10.0);
	gm1.setInterpolationStep(0.3);
	gm1.setParam(stat, 678.9, 789.0);

  GaussModel gm2;
  gm2 = gm1;

  GaussModel gm3;
	gm3.setScalingFactor(10.0);
	gm3.setInterpolationStep(0.3);
	gm3.setParam(stat, 678.9, 789.0);

  gm1 = GaussModel();
	TEST_EQUAL(gm3.getParameters(), gm2.getParameters())
RESULT

// copy ctor
CHECK(GaussModel(const GaussModel& source))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	gm1.setScalingFactor(10.0);
	gm1.setInterpolationStep(0.3);
	gm1.setParam(stat, 678.9, 789.0);

	GaussModel gm2(gm1);
  	GaussModel gm3;
	gm3.setScalingFactor(10.0);
	gm3.setInterpolationStep(0.3);
	gm3.setParam(stat, 678.9, 789.0);

  	gm1 = GaussModel();
	TEST_EQUAL(gm3.getParameters(), gm2.getParameters())
RESULT

CHECK([EXTRA] DefaultParmHandler::setParameters(...))
	PRECISION(0.001)
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(679.1);
	stat.setVariance(2.0);

	gm1.setScalingFactor(10.0);
	gm1.setParam(stat, 678.9, 680.9);
	gm1.setOffset(680.0);

	TEST_REAL_EQUAL(gm1.getCenter(), 680.2)

	GaussModel gm2;
	const Param p(gm1.getParameters());
	gm2.setParameters(p);

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	gm1.getSamples(dpa1);
	gm2.getSamples(dpa2);

	PRECISION(0.0000001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK(void setParam(const BasicStatistics&,CoordinateType,CoordinateType))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(0.0);
	stat.setVariance(1.0);
	gm1.setInterpolationStep(0.001);
	gm1.setParam(stat, -4.0, 4.0);

	TEST_REAL_EQUAL(gm1.getCenter(), 0.0)

	PRECISION(0.001)
	TEST_REAL_EQUAL(gm1.getIntensity(-1.0), 0.24197072);
	TEST_REAL_EQUAL(gm1.getIntensity(0.0), 0.39894228);
	TEST_REAL_EQUAL(gm1.getIntensity(1.0), 0.24197072);
	TEST_REAL_EQUAL(gm1.getIntensity(2.0), 0.05399097);

	gm1.setInterpolationStep(0.2);
	gm1.setSamples();

	TEST_REAL_EQUAL(gm1.getIntensity(-1.0), 0.24197072);
	TEST_REAL_EQUAL(gm1.getIntensity(0.0), 0.39894228);
	TEST_REAL_EQUAL(gm1.getIntensity(1.0), 0.24197072);
	TEST_REAL_EQUAL(gm1.getIntensity(2.0), 0.05399097);

	gm1.setScalingFactor(10.0);
	gm1.setSamples();

	TEST_REAL_EQUAL(gm1.getIntensity(-1.0), 2.4197072);
	TEST_REAL_EQUAL(gm1.getIntensity(0.0), 3.9894228);
	TEST_REAL_EQUAL(gm1.getIntensity(1.0), 2.4197072);
	TEST_REAL_EQUAL(gm1.getIntensity(2.0), 0.5399097);
RESULT


CHECK(void setOffset(double offset))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	gm1.setParam(stat, 678.9, 789.0);
	gm1.setOffset(680.9);

	GaussModel gm2;
	stat.setMean(682.1);
	stat.setVariance(2.0);
	gm2.setParam(stat, 680.9, 791.0);

	TEST_EQUAL(gm1.getParameters(), gm2.getParameters())
	TEST_REAL_EQUAL(gm1.getCenter(), gm2.getCenter())
	TEST_REAL_EQUAL(gm1.getCenter(), 682.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	gm1.getSamples(dpa1);
	gm2.getSamples(dpa2);

	PRECISION(0.001)
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
