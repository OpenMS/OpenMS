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
CHECK((EmgModel()))
	ptr = new EmgModel();
  	TEST_EQUAL(ptr->getName(), "EmgModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~EmgModel()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(EmgModel::getProductName(),"EmgModel")
	TEST_EQUAL(EmgModel().getName(),"EmgModel")
RESULT

CHECK((static BaseModel<1>* create()))
	BaseModel<1>* ptr = EmgModel::create();
	TEST_EQUAL(ptr->getName(), "EmgModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// assignment operator
CHECK((virtual EmgModel& operator=(const EmgModel &source)))
	EmgModel em1;
	em1.setInterpolationStep(0.2);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height",100000.0);
	tmp.setValue("emg:width",5.0);
	tmp.setValue("emg:symmetry",5.0);
	tmp.setValue("emg:retention",725.0);
	em1.setParameters(tmp);

	EmgModel em2;
	em2 = em1;
	
	EmgModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParameters(tmp);
  	TEST_EQUAL(em3.getParameters(), em2.getParameters())
RESULT

// copy ctor
CHECK((EmgModel(const EmgModel& source)))
	EmgModel em1;
	em1.setInterpolationStep(0.2);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height",100000.0);
	tmp.setValue("emg:width",5.0);
	tmp.setValue("emg:symmetry",5.0);
	tmp.setValue("emg:retention",725.0);
	em1.setParameters(tmp);

	EmgModel em2(em1);
  	EmgModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParameters(tmp);

 	em1 = EmgModel();
	TEST_EQUAL(em3.getParameters(), em2.getParameters())
RESULT


CHECK([EXTRA] DefaultParamHandler::setParameters(...))

	PRECISION(0.001)
	EmgModel em1;

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 680.9);
	tmp.setValue("statistics:mean", 679.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height", 100000.0);
	tmp.setValue("emg:width", 5.0);
	tmp.setValue("emg:symmetry", 5.0);
	tmp.setValue("emg:retention", 1200.0);
	em1.setParameters(tmp);
	em1.setOffset(680.0);

	TEST_REAL_EQUAL(em1.getCenter(), 680.2)

	EmgModel em3;
	em3.setParameters(em1.getParameters());

	DPeakArray<DPeak<1> > dpa1;
	DPeakArray<DPeak<1> > dpa2;
	em1.getSamples(dpa1);
	em3.getSamples(dpa2);

	PRECISION(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (UInt i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

	EmgModel em2;
	em2.setInterpolationStep(0.1);

	tmp.setValue("bounding_box:min", -1.0);
	tmp.setValue("bounding_box:max", 4.0);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance",  0.1);
	tmp.setValue("emg:height", 10.0);
	tmp.setValue("emg:width", 1.0);
	tmp.setValue("emg:symmetry", 2.0);
	tmp.setValue("emg:retention", 3.0);
	em2.setParameters(tmp);

	TEST_REAL_EQUAL(em2.getCenter(), 0.0)

	PRECISION(0.01)
	TEST_REAL_EQUAL(em2.getIntensity(-1.0), 0.0497198);
	TEST_REAL_EQUAL(em2.getIntensity(0.0), 0.164882);
	TEST_REAL_EQUAL(em2.getIntensity(1.0), 0.54166);
	TEST_REAL_EQUAL(em2.getIntensity(2.0), 1.69364);

	em2.setInterpolationStep(0.2);
	em2.setSamples();

	TEST_REAL_EQUAL(em2.getIntensity(-1.0), 0.0497198);
	TEST_REAL_EQUAL(em2.getIntensity(0.0), 0.164882);
	TEST_REAL_EQUAL(em2.getIntensity(1.0), 0.54166);
	TEST_REAL_EQUAL(em2.getIntensity(2.0), 1.69364);


	// checked small values of parameter symmetry
	tmp.setValue("bounding_box:min", 0.0);
	tmp.setValue("bounding_box:max", 10.0);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance",  0.1);
	tmp.setValue("emg:height", 10.0);
	tmp.setValue("emg:width", 6.0);
	tmp.setValue("emg:symmetry", 1.0);
	tmp.setValue("emg:retention", 3.0);
	em2.setParameters(tmp);

	TEST_REAL_EQUAL(em2.getIntensity(2.0), 747203);

	tmp.setValue("emg:symmetry", 0.1);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(em2.getIntensity(2.0)))

	tmp.setValue("emg:symmetry", 0.16);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(em2.getIntensity(2.0)))

	tmp.setValue("emg:symmetry", 0.17);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(!em2.getIntensity(2.0)))

	tmp.setValue("emg:symmetry", 0.2);
	em2.setParameters(tmp);
	ABORT_IF(!std::isinf(em2.getIntensity(2.0)))

RESULT

CHECK((void setOffset(CoordinateType offset)))
	EmgModel em1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height", 100000.0);
	tmp.setValue("emg:width", 5.0);
	tmp.setValue("emg:symmetry", 5.0);
	tmp.setValue("emg:retention", 725.0);
	em1.setParameters(tmp);
	em1.setOffset(680.9);

	EmgModel em2;
	em2.setParameters(tmp);
	em2.setOffset(680.9);

	TEST_EQUAL(em1.getParameters(), em2.getParameters())
	TEST_REAL_EQUAL(em1.getCenter(), em2.getCenter())
	TEST_REAL_EQUAL(em1.getCenter(), 682.1)

	DPeakArray<DPeak<1> > dpa1;
	DPeakArray<DPeak<1> > dpa2;
	em1.getSamples(dpa1);
	em2.getSamples(dpa2);

	PRECISION(0.01)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (UInt i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
RESULT

CHECK((CoordinateType getCenter() const))
	PRECISION(0.001)
	EmgModel em1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height",  100000.0);
	tmp.setValue("emg:width",  5.0);
	tmp.setValue("emg:symmetry",  5.0);
	tmp.setValue("emg:retention",  725.0);
	em1.setParameters(tmp);
	em1.setOffset(680.0);
	TEST_REAL_EQUAL(em1.getCenter(), 681.2)

RESULT

CHECK((void setSamples()))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
