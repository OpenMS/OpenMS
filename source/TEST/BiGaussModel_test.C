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
CHECK((BiGaussModel()))
	ptr = new BiGaussModel();
        TEST_EQUAL(ptr->getName(), "BiGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~BiGaussModel()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(BiGaussModel::getProductName(),"BiGaussModel")
	TEST_EQUAL(BiGaussModel().getName(),"BiGaussModel")
RESULT

CHECK( static BaseModel<1>* create() )
	BaseModel<1>* ptr = BiGaussModel::create();
	TEST_EQUAL(ptr->getName(), "BiGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// assignment operator
CHECK((virtual BiGaussModel& operator=(const BiGaussModel &source)))
	BiGaussModel bgm1;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);

  BiGaussModel bgm2;
  bgm2 = bgm1;

  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParameters(tmp);

  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
RESULT

// copy ctor
CHECK((BiGaussModel(const BiGaussModel& source)))
	BiGaussModel bgm1;
	BasicStatistics<>  stat;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);

	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);

	BiGaussModel bgm2(bgm1);
  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParameters(tmp);
  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
RESULT

CHECK([EXTRA] DefaultParamHandler::setParameters(...))
	PRECISION(0.001)
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.0);

	BiGaussModel bgm2;
	bgm2.setParameters(bgm1.getParameters());
	TEST_REAL_EQUAL(bgm1.getCenter(), 681.2)

	DPeakArray<DPeak<1> > dpa1;
	DPeakArray<DPeak<1> > dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	PRECISION(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (UInt i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK((void setOffset(CoordinateType offset)))
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.9);

	BiGaussModel bgm2;
	tmp.setValue("bounding_box:min", 680.9);
	tmp.setValue("bounding_box:max", 791.0);
	tmp.setValue("statistics:mean", 682.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm2.setParameters(tmp);

	TEST_EQUAL(bgm1.getParameters(), bgm2.getParameters())
	TEST_REAL_EQUAL(bgm1.getCenter(), bgm2.getCenter())
	TEST_REAL_EQUAL(bgm1.getCenter(), 682.1)

	DPeakArray<DPeak<1> > dpa1;
	DPeakArray<DPeak<1> > dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	PRECISION(0.001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (UInt i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

	tmp.setValue("bounding_box:min", -4.0);
	tmp.setValue("bounding_box:max", 4.001);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance1",  0.81);
	tmp.setValue("statistics:variance2", 0.81 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(0.123);
	TEST_REAL_EQUAL(bgm1.getCenter(), 4.123)

	PRECISION(0.001)
	TEST_REAL_EQUAL(bgm1.getIntensity(4.123), 0.4432692);
	TEST_REAL_EQUAL(bgm1.getIntensity(4.223), bgm1.getIntensity(4.023));
	TEST_REAL_EQUAL(bgm1.getIntensity(3.123), bgm1.getIntensity(5.123));

RESULT

CHECK( CoordinateType getCenter() const )
	
	// already test above, but just for the sake of it
	PRECISION(0.001)
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.0);
	TEST_REAL_EQUAL(bgm1.getCenter(), 681.2)

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
