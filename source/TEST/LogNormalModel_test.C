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
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalModel.h>

///////////////////////////

START_TEST(LogNormalModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
LogNormalModel* ptr = 0;
CHECK(LogNormalModel())
	ptr = new LogNormalModel();
  TEST_EQUAL(ptr->getName(), "LogNormalModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~LogNormalModel())
	delete ptr;
RESULT


CHECK(const String getName())
	TEST_EQUAL(LogNormalModel::getProductName(),"LogNormalModel")
	TEST_EQUAL(LogNormalModel().getName(),"LogNormalModel")
RESULT

// assignment operator
CHECK(LogNormalModel& operator = (const LogNormalModel& source))
	LogNormalModel logm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	logm1.setInterpolationStep(0.2);
	logm1.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 2.0, 678.9, 789.0);

	LogNormalModel logm2;
	logm2 = logm1;
	
	LogNormalModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 2.0, 678.9, 789.0);

  	logm1 = LogNormalModel();
	TEST_EQUAL(em3.getParameters(), logm2.getParameters())
RESULT

// copy ctor
CHECK(LogNormalModel(const LogNormalModel& source))
	LogNormalModel logm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	logm1.setInterpolationStep(0.2);
	logm1.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 2.0, 678.9, 789.0);

	LogNormalModel logm2(logm1);
  	LogNormalModel logm3;
	logm3.setInterpolationStep(0.2);
	logm3.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 2.0, 678.9, 789.0);

 	logm1 = LogNormalModel();
	TEST_EQUAL(logm3.getParameters(), logm2.getParameters())
RESULT

CHECK(void setParam(Param param))
	PRECISION(0.001)
	LogNormalModel logm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	logm1.setParam(stat, 1000000.0, 20.0, 3.0, 400.0, 2.0, 678.9, 700.0);
	
	LogNormalModel logm2;
	logm2.setParam(stat, 1000000.0, 20.0, 3.0, 400.0, 2.0, 678.9, 700.0);

	TEST_REAL_EQUAL(logm1.getCenter(), 680.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	logm1.getSamples(dpa1);
	logm2.getSamples(dpa2);

	PRECISION(0.1)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

RESULT

CHECK(void setParam(const Math::BasicStatistics&, CoordinateType, CoordinateType, CoordinateType, CoordinateType, CoordinateType, CoordinateType, CoordinateType))

	LogNormalModel logm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(0.0);
	stat.setVariance(0.1);
	logm1.setInterpolationStep(0.1);
	logm1.setParam(stat, 100.0, 5.0, 2.0, 3.0, 2.0, -1.0, 4.0);

	TEST_REAL_EQUAL(logm1.getCenter(), 0.0)

	PRECISION(0.01)
	TEST_REAL_EQUAL(logm1.getIntensity(0.0), 0.047651);
	TEST_REAL_EQUAL(logm1.getIntensity(1.0), 29.7819);
	TEST_REAL_EQUAL(logm1.getIntensity(2.0), 83.2322);
	TEST_REAL_EQUAL(logm1.getIntensity(3.0), 100.0);

	logm1.setInterpolationStep(0.2);
	logm1.setSamples();

	TEST_REAL_EQUAL(logm1.getIntensity(0.0), 0.047651);
	TEST_REAL_EQUAL(logm1.getIntensity(1.0), 29.7819);
	TEST_REAL_EQUAL(logm1.getIntensity(2.0), 83.2322);
	TEST_REAL_EQUAL(logm1.getIntensity(3.0), 100.0);


RESULT


CHECK(void setOffset(double offset))

	LogNormalModel logm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	logm1.setParam(stat, 1000000.0, 20.0, 3.0, 400.0, 2.0, 678.9, 700.0);
	logm1.setOffset(680.9);

	LogNormalModel logm2;
	logm2.setParam(stat, 1000000.0, 20.0, 3.0, 400.0, 2.0, 678.9, 700.0);
	logm2.setOffset(680.9);

	TEST_EQUAL(logm1.getParameters(), logm2.getParameters())
	TEST_REAL_EQUAL(logm1.getCenter(), logm2.getCenter())
	TEST_REAL_EQUAL(logm1.getCenter(), 682.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	logm1.getSamples(dpa1);
	logm2.getSamples(dpa2);

	PRECISION(0.1)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

RESULT

// checked by other check-methods 
// It is not necessarily to test the methods again.
CHECK(const CoordinateType getCenter() const)
RESULT

CHECK(static BaseModel<1>* create())
RESULT

CHECK(void setSamples())
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
