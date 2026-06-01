// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FEATUREFINDER/EmgModel.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <boost/math/special_functions/fpclassify.hpp>


///////////////////////////

START_TEST(EmgModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
EmgModel* ptr = nullptr;
EmgModel* nullPointer = nullptr;
START_SECTION((EmgModel()))
	ptr = new EmgModel();
  	TEST_EQUAL(ptr->getName(), "EmgModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~EmgModel()))
	delete ptr;
END_SECTION

// assignment operator
START_SECTION((virtual EmgModel& operator=(const EmgModel &source)))
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
END_SECTION

// copy ctor
START_SECTION((EmgModel(const EmgModel& source)))
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
END_SECTION


START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))

	TOLERANCE_ABSOLUTE(0.001)
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

	TEST_REAL_SIMILAR(em1.getCenter(), 680.2)

	EmgModel em3;
	em3.setParameters(em1.getParameters());

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	em1.getSamples(dpa1);
	em3.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
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

	TEST_REAL_SIMILAR(em2.getCenter(), 0.0)

	TOLERANCE_ABSOLUTE(0.01)
	TEST_REAL_SIMILAR(em2.getIntensity(-1.0), 0.0497198);
	TEST_REAL_SIMILAR(em2.getIntensity(0.0), 0.164882);
	TEST_REAL_SIMILAR(em2.getIntensity(1.0), 0.54166);
	TEST_REAL_SIMILAR(em2.getIntensity(2.0), 1.69364);

	em2.setInterpolationStep(0.2);
	em2.setSamples();

	TEST_REAL_SIMILAR(em2.getIntensity(-1.0), 0.0497198);
	TEST_REAL_SIMILAR(em2.getIntensity(0.0), 0.164882);
	TEST_REAL_SIMILAR(em2.getIntensity(1.0), 0.54166);
	TEST_REAL_SIMILAR(em2.getIntensity(2.0), 1.69364);


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

	TEST_REAL_SIMILAR(em2.getIntensity(2.0), 747203);

	tmp.setValue("emg:symmetry", 0.1);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(em2.getIntensity(2.0)))

	tmp.setValue("emg:symmetry", 0.16);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(em2.getIntensity(2.0)))

	tmp.setValue("emg:symmetry", 0.17);
	em2.setParameters(tmp);
	ABORT_IF(std::isinf(float(!em2.getIntensity(2.0))))

END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))
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
	TEST_REAL_SIMILAR(em1.getCenter(), em2.getCenter())
	TEST_REAL_SIMILAR(em1.getCenter(), 682.1)

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	em1.getSamples(dpa1);
	em2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.01)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION((CoordinateType getCenter() const))
	TOLERANCE_ABSOLUTE(0.001)
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
	TEST_REAL_SIMILAR(em1.getCenter(), 681.2)

END_SECTION

START_SECTION((void setSamples()))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
