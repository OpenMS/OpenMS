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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgModel.h>

///////////////////////////

START_TEST(EmgModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
EmgModel* ptr = 0;
CHECK(EmgModel())
	ptr = new EmgModel();
  TEST_EQUAL(ptr->getName(), "EmgModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~EmgModel())
	delete ptr;
RESULT


CHECK(const String getName())
	TEST_EQUAL(EmgModel::getProductName(),"EmgModel")
	TEST_EQUAL(EmgModel().getName(),"EmgModel")
RESULT

// assignment operator
CHECK(EmgModel& operator = (const EmgModel& source))
	EmgModel em1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	em1.setInterpolationStep(0.2);
	em1.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);

	EmgModel em2;
	em2 = em1;
	
	EmgModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);

  	em1 = EmgModel();
	TEST_EQUAL(em3.getParameters(), em2.getParameters())
RESULT

// copy ctor
CHECK(EmgModel(const EmgModel& source))
	EmgModel em1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	em1.setInterpolationStep(0.2);
	em1.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);

	EmgModel em2(em1);
  	EmgModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);

 	em1 = EmgModel();
	TEST_EQUAL(em3.getParameters(), em2.getParameters())
RESULT

CHECK(void setParam(Param param))
	PRECISION(0.001)
	EmgModel em1;
	Math::BasicStatistics<>  stat;
	stat.setMean(679.1);
	stat.setVariance(2.0);

	em1.setParam(stat, 100000.0, 5.0, 5.0, 1200.0, 678.9, 680.9);
	em1.setOffset(680.0);

	TEST_REAL_EQUAL(em1.getCenter(), 680.2)

	EmgModel em2;
	em2.setParameters(em1.getParameters());

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	em1.getSamples(dpa1);
	em2.getSamples(dpa2);

	PRECISION(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

RESULT

CHECK(void setParam(const Math::BasicStatistics&, CoordinateType, CoordinateType, CoordinateType, CoordinateType, CoordinateType, CoordinateType))
	EmgModel em1;
	Math::BasicStatistics<>  stat;
	stat.setMean(0.0);
	stat.setVariance(0.1);
	em1.setInterpolationStep(0.1);
	em1.setParam(stat, 10.0, 1.0, 2.0, 3.0, -1.0, 4.0);

	TEST_REAL_EQUAL(em1.getCenter(), 0.0)

	PRECISION(0.01)
	TEST_REAL_EQUAL(em1.getIntensity(-1.0), 0.0497198);
	TEST_REAL_EQUAL(em1.getIntensity(0.0), 0.164882);
	TEST_REAL_EQUAL(em1.getIntensity(1.0), 0.54166);
	TEST_REAL_EQUAL(em1.getIntensity(2.0), 1.69364);

	em1.setInterpolationStep(0.2);
	em1.setSamples();

	TEST_REAL_EQUAL(em1.getIntensity(-1.0), 0.0497198);
	TEST_REAL_EQUAL(em1.getIntensity(0.0), 0.164882);
	TEST_REAL_EQUAL(em1.getIntensity(1.0), 0.54166);
	TEST_REAL_EQUAL(em1.getIntensity(2.0), 1.69364);

RESULT


CHECK(void setOffset(double offset))

	EmgModel em1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	em1.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);
	em1.setOffset(680.9);

	EmgModel em2;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	em2.setParam(stat, 100000.0, 5.0, 5.0, 725.0, 678.9, 789.0);
	em2.setOffset(680.9);

	TEST_EQUAL(em1.getParameters(), em2.getParameters())
	TEST_REAL_EQUAL(em1.getCenter(), em2.getCenter())
	TEST_REAL_EQUAL(em1.getCenter(), 682.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	em1.getSamples(dpa1);
	em2.getSamples(dpa2);

	PRECISION(0.01)
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
