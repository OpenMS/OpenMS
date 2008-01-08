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
CHECK((LogNormalModel()))
	ptr = new LogNormalModel();
  TEST_EQUAL(ptr->getName(), "LogNormalModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK((virtual ~LogNormalModel()))
	delete ptr;
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(LogNormalModel::getProductName(),"LogNormalModel")
	TEST_EQUAL(LogNormalModel().getName(),"LogNormalModel")
RESULT

CHECK((static BaseModel<1>* create()))
	BaseModel<1>* ptr = LogNormalModel::create();
	TEST_EQUAL(ptr->getName(), "LogNormalModel")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// assignment operator
CHECK((virtual LogNormalModel& operator=(const LogNormalModel &source)))
	LogNormalModel logm1;
	logm1.setInterpolationStep(0.2);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lognormal:height",  100000.0);
	tmp.setValue("lognormal:width",  5.0);
	tmp.setValue("lognormal:symmetry",  5.0);
	tmp.setValue("lognormal:retention",  725.0);
	tmp.setValue("lognormal:r",  2.0);
	logm1.setParameters(tmp);

	LogNormalModel logm2;
	logm2 = logm1;
	
	LogNormalModel logm3;
	logm3.setInterpolationStep(0.2);
	logm3.setParameters(tmp);

  	logm1 = LogNormalModel();
	TEST_EQUAL(logm3.getParameters(), logm2.getParameters())
RESULT

// copy ctor
CHECK((LogNormalModel(const LogNormalModel& source)))
	LogNormalModel logm1;
	logm1.setInterpolationStep(0.2);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lognormal:height",  100000.0);
	tmp.setValue("lognormal:width",  5.0);
	tmp.setValue("lognormal:symmetry",  5.0);
	tmp.setValue("lognormal:retention",  725.0);
	tmp.setValue("lognormal:r",  2.0);
	logm1.setParameters(tmp);

	LogNormalModel logm2(logm1);
  	LogNormalModel logm3;
	logm3.setInterpolationStep(0.2);
	logm3.setParameters(tmp);

 	logm1 = LogNormalModel();
	TEST_EQUAL(logm3.getParameters(), logm2.getParameters())
RESULT

CHECK([EXTRA] DefaultParamHandler::setParameters(...))
	
	LogNormalModel logm1;
	logm1.setInterpolationStep(0.1);
	
	Param tmp;
	tmp.setValue("bounding_box:min", -1.0);
	tmp.setValue("bounding_box:max", 4.0);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance",  0.1);
	tmp.setValue("lognormal:height",  100.0);
	tmp.setValue("lognormal:width",  5.0);
	tmp.setValue("lognormal:symmetry",  2.0);
	tmp.setValue("lognormal:retention",  3.0);
	tmp.setValue("lognormal:r",  2.0);
	logm1.setParameters(tmp);

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

	// symmetry cannot be 1, because the log(1)=0 => division by zero
	tmp.setValue("lognormal:symmetry",  1.0);
	logm1.setParameters(tmp);
	ABORT_IF(std::isnan(logm1.getIntensity(1.0)))

	// symmetry cannot be 0, cause division by zero
	tmp.setValue("lognormal:symmetry",  0.0);
	logm1.setParameters(tmp);
	ABORT_IF(std::isnan(logm1.getIntensity(1.0)))

	// small values for the parameter symmetry are valid
	tmp.setValue("lognormal:symmetry",  1.001);
	logm1.setParameters(tmp);
	ABORT_IF(!std::isnan(logm1.getIntensity(1.0)))
	ABORT_IF(!std::isinf(logm1.getIntensity(1.0)))

	tmp.setValue("lognormal:symmetry",  0.998);
	logm1.setParameters(tmp);
	ABORT_IF(!std::isinf(logm1.getIntensity(1.0)))
	ABORT_IF(!std::isnan(logm1.getIntensity(1.0)))
	
	tmp.setValue("lognormal:symmetry",  0.001);
	logm1.setParameters(tmp);
	ABORT_IF(!std::isinf(logm1.getIntensity(1.0)))
	
	tmp.setValue("lognormal:symmetry",  -0.001);
	logm1.setParameters(tmp);
	ABORT_IF(!std::isinf(logm1.getIntensity(1.0)))

RESULT

CHECK((void setOffset(CoordinateType offset)))
	LogNormalModel logm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 700.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lognormal:height",  1000000.0);
	tmp.setValue("lognormal:width",  20.0);
	tmp.setValue("lognormal:symmetry",  3.0);
	tmp.setValue("lognormal:retention",  400.0);
	tmp.setValue("lognormal:r",  2.0);
	
	logm1.setParameters(tmp);
	logm1.setOffset(680.9);

	LogNormalModel logm2;
	logm2.setParameters(tmp);
	logm2.setOffset(680.9);

	TEST_EQUAL(logm1.getParameters(), logm2.getParameters())
	TEST_REAL_EQUAL(logm1.getCenter(), logm2.getCenter())
	TEST_REAL_EQUAL(logm1.getCenter(), 682.1)

	DPeakArray<DPeak<1> > dpa1;
	DPeakArray<DPeak<1> > dpa2;
	logm1.getSamples(dpa1);
	logm2.getSamples(dpa2);

	PRECISION(0.1)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (UInt i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_EQUAL(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_EQUAL(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

RESULT

CHECK((void setSamples()))
	// already test above
RESULT

CHECK((CoordinateType getCenter() const))
	// already test above, but just for the sake of it
	PRECISION(0.001)
	LogNormalModel logm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lognormal:height",  100000.0);
	tmp.setValue("lognormal:width",  5.0);
	tmp.setValue("lognormal:symmetry",  5.0);
	tmp.setValue("lognormal:retention",  725.0);
	tmp.setValue("lognormal:r",  2.0);
	logm1.setParameters(tmp);
	logm1.setOffset(680.0);
	TEST_REAL_EQUAL(logm1.getCenter(), 681.2)

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
