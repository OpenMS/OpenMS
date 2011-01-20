// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
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
START_SECTION((LmaGaussModel()))
	ptr = new LmaGaussModel();
  TEST_EQUAL(ptr->getName(), "LmaGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

// destructor
START_SECTION((virtual ~LmaGaussModel()))
	delete ptr;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(LmaGaussModel::getProductName(),"LmaGaussModel")
	TEST_EQUAL(LmaGaussModel().getName(),"LmaGaussModel")
END_SECTION

START_SECTION((static BaseModel<1>* create()))
	BaseModel<1>* ptr = LmaGaussModel::create();
	TEST_EQUAL(ptr->getName(), "LmaGaussModel")
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

// assignment operator
START_SECTION((virtual LmaGaussModel& operator=(const LmaGaussModel &source)))
	LmaGaussModel lm1;
	lm1.setInterpolationStep(0.3);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lma:scale_factor", 1000000.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 680.0);
	lm1.setParameters(tmp);

	LmaGaussModel lm2;
  	lm2 = lm1;

  	LmaGaussModel lm3;
	lm3.setInterpolationStep(0.3);
	lm3.setParameters(tmp);
	
	TEST_EQUAL(lm3.getParameters(), lm2.getParameters())
END_SECTION

// copy ctor
START_SECTION((LmaGaussModel(const LmaGaussModel& source)))
	LmaGaussModel lm1;
	lm1.setInterpolationStep(0.3);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lma:scale_factor", 10.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 680.0);
	lm1.setParameters(tmp);

	LmaGaussModel lm2(lm1);
  	LmaGaussModel lm3;
	lm3.setInterpolationStep(0.3);
	lm3.setParameters(tmp);

	TEST_EQUAL(lm3.getParameters(), lm2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	TOLERANCE_ABSOLUTE(0.001)
	LmaGaussModel lm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 680.9);
	tmp.setValue("statistics:mean", 679.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lma:scale_factor", 10.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 700.0);
	lm1.setParameters(tmp);
	lm1.setOffset(680.0);

	TEST_REAL_SIMILAR(lm1.getCenter(), 680.2)

	LmaGaussModel lm2;
	lm2.setParameters(lm1.getParameters());

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	lm1.getSamples(dpa1);
	lm2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	LmaGaussModel lm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", -1.0);
	tmp.setValue("bounding_box:max", 4.0);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance",  0.1);
	tmp.setValue("lma:scale_factor", 1.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 3.0);
	lm1.setParameters(tmp);

	TEST_REAL_SIMILAR(lm1.getCenter(), 0.0)

	TOLERANCE_ABSOLUTE(0.001)
	TEST_REAL_SIMILAR(lm1.getIntensity(-1.0), 0.0269955);
	TEST_REAL_SIMILAR(lm1.getIntensity(0.0), 0.0647588);
	TEST_REAL_SIMILAR(lm1.getIntensity(1.0), 0.120985);
	TEST_REAL_SIMILAR(lm1.getIntensity(2.0), 0.176033);

	lm1.setInterpolationStep(0.2);
	lm1.setSamples();

	TEST_REAL_SIMILAR(lm1.getIntensity(-1.0), 0.0269955);
	TEST_REAL_SIMILAR(lm1.getIntensity(0.0), 0.0647588);
	TEST_REAL_SIMILAR(lm1.getIntensity(1.0), 0.120985);
	TEST_REAL_SIMILAR(lm1.getIntensity(2.0), 0.176033);

	TOLERANCE_ABSOLUTE(0.1)
	tmp.setValue("lma:scale_factor", 10.0);
	lm1.setParameters(tmp);
	lm1.setSamples();
	
	TEST_REAL_SIMILAR(lm1.getIntensity(-1.0), 0.269955);
	TEST_REAL_SIMILAR(lm1.getIntensity(0.0), 0.647588);
	TEST_REAL_SIMILAR(lm1.getIntensity(1.0), 1.20985);
	TEST_REAL_SIMILAR(lm1.getIntensity(2.0), 1.76033);
END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))

	LmaGaussModel lm1;

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lma:scale_factor", 10.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 700.0);
	lm1.setParameters(tmp);
	lm1.setOffset(680.9);

	LmaGaussModel lm2;
	lm2.setParameters(tmp);
	lm2.setOffset(680.9);

	TEST_EQUAL(lm1.getParameters(), lm2.getParameters())
	TEST_REAL_SIMILAR(lm1.getCenter(), lm2.getCenter())
	TEST_REAL_SIMILAR(lm1.getCenter(), 682.1)

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	lm1.getSamples(dpa1);
	lm2.getSamples(dpa2);

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
	// already test above, but just for the sake of it
	TOLERANCE_ABSOLUTE(0.001)
	LmaGaussModel lm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("lma:scale_factor", 10.0);
	tmp.setValue("lma:standard_deviation", 2.0);
	tmp.setValue("lma:expected_value", 700.0);	
	lm1.setParameters(tmp);
	lm1.setOffset(680.0);

	TEST_REAL_SIMILAR(lm1.getCenter(), 681.2)

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
