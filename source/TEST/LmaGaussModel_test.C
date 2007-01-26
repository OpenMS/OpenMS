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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussModel.h>

///////////////////////////

START_TEST(LmaGaussModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
LmaGaussModel* ptr = 0;
CHECK(LmaGaussModel())
	ptr = new LmaGaussModel();
  TEST_EQUAL(ptr->getName(), "LmaGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~LmaGaussModel())
	delete ptr;
RESULT


CHECK(const String getName())
	TEST_EQUAL(LmaGaussModel::getProductName(),"LmaGaussModel")
	TEST_EQUAL(LmaGaussModel().getName(),"LmaGaussModel")
RESULT

// assignment operator
CHECK(LmaGaussModel& operator = (const LmaGaussModel& source))
	LmaGaussModel lm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	lm1.setInterpolationStep(0.3);
	lm1.setParam(stat, 1000000.0, 2.0, 680.0, 678.9, 789.0);

  LmaGaussModel lm2;
  lm2 = lm1;

  LmaGaussModel lm3;
	lm3.setInterpolationStep(0.3);
	lm3.setParam(stat, 1000000.0, 2.0, 680.0, 678.9, 789.0);

	TEST_EQUAL(lm3.getParameters(), lm2.getParameters())
RESULT

// copy ctor
CHECK(LmaGaussModel(const LmaGaussModel& source))
	LmaGaussModel lm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	lm1.setInterpolationStep(0.3);
	lm1.setParam(stat, 10.0, 2.0, 680.0, 678.9, 789.0);

	LmaGaussModel lm2(lm1);
  LmaGaussModel lm3;
	lm3.setInterpolationStep(0.3);
	lm3.setParam(stat, 10.0, 2.0, 680.0, 678.9, 789.0);

	TEST_EQUAL(lm3.getParameters(), lm2.getParameters())
RESULT

CHECK(void setParam(Param param))
	PRECISION(0.001)
	LmaGaussModel lm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(679.1);
	stat.setVariance(2.0);

	lm1.setParam(stat, 10.0, 2.0, 700.0, 678.9, 680.9);
	lm1.setOffset(680.0);

	TEST_REAL_EQUAL(lm1.getCenter(), 680.2)

	LmaGaussModel lm2;
	lm2.setParameters(lm1.getParameters());

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	lm1.getSamples(dpa1);
	lm2.getSamples(dpa2);

	PRECISION(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK(void setParam(const Math::BasicStatistics&,CoordinateType,CoordinateType,CoordinateType,CoordinateType,CoordinateType))
	LmaGaussModel lm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(0.0);
	stat.setVariance(0.1);
	lm1.setInterpolationStep(0.001);
	lm1.setParam(stat, 1.0, 2.0, 3.0, -1.0, 4.0);

	TEST_REAL_EQUAL(lm1.getCenter(), 0.0)

	PRECISION(0.001)
	TEST_REAL_EQUAL(lm1.getIntensity(-1.0), 0.0269955);
	TEST_REAL_EQUAL(lm1.getIntensity(0.0), 0.0647588);
	TEST_REAL_EQUAL(lm1.getIntensity(1.0), 0.120985);
	TEST_REAL_EQUAL(lm1.getIntensity(2.0), 0.176033);

	lm1.setInterpolationStep(0.2);
	lm1.setSamples();

	TEST_REAL_EQUAL(lm1.getIntensity(-1.0), 0.0269955);
	TEST_REAL_EQUAL(lm1.getIntensity(0.0), 0.0647588);
	TEST_REAL_EQUAL(lm1.getIntensity(1.0), 0.120985);
	TEST_REAL_EQUAL(lm1.getIntensity(2.0), 0.176033);

	PRECISION(0.1)
	lm1.setParam(stat, 10.0, 2.0, 3.0, -1.0, 4.0);
	lm1.setSamples();
	
	TEST_REAL_EQUAL(lm1.getIntensity(-1.0), 0.269955);
	TEST_REAL_EQUAL(lm1.getIntensity(0.0), 0.647588);
	TEST_REAL_EQUAL(lm1.getIntensity(1.0), 1.20985);
	TEST_REAL_EQUAL(lm1.getIntensity(2.0), 1.76033);


RESULT


CHECK(void setOffset(double offset))

	LmaGaussModel lm1;
	Math::BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	lm1.setParam(stat, 10.0, 2.0, 700.0, 678.9, 789.0);
	lm1.setOffset(680.9);

	LmaGaussModel lm2;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	lm2.setParam(stat, 10.0, 2.0, 700.0, 678.9, 789.0);
	lm2.setOffset(680.9);

	TEST_EQUAL(lm1.getParameters(), lm2.getParameters())
	TEST_REAL_EQUAL(lm1.getCenter(), lm2.getCenter())
	TEST_REAL_EQUAL(lm1.getCenter(), 682.1)

	DPeakArray<1> dpa1;
	DPeakArray<1> dpa2;
	lm1.getSamples(dpa1);
	lm2.getSamples(dpa2);

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
