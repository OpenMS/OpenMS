// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche  $
// $Authors: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SIMULATION/ElutionModel.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ElutionModel, "$Id: ElutionModel_test.C 5383 2009-06-18 13:07:26Z cbielow $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ElutionModel* ptr = 0;
START_SECTION(ElutionModel())
{
	ptr = new ElutionModel();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~ElutionModel())
{
	delete ptr;
}
END_SECTION

START_SECTION((ElutionModel(const ElutionModel &source)))
{
	ElutionModel em1;
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

	ElutionModel em2(em1);
  ElutionModel em3;
	em3.setInterpolationStep(0.2);
	em3.setParameters(tmp);

 	em1 = ElutionModel();
	TEST_EQUAL(em3.getParameters(), em2.getParameters())
	TEST_EQUAL(em3==em2,true)
}
END_SECTION

START_SECTION((virtual ElutionModel& operator=(const ElutionModel &source)))
{
	ElutionModel em1;
	em1.setInterpolationStep(0.2);

	Param tmp;
	tmp.setValue("bounding_box:min", 678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 686.0 );
	tmp.setValue("statistics:variance",  2.0);
	tmp.setValue("emg:height",100000.0);
	tmp.setValue("emg:width",5.0);
	tmp.setValue("emg:symmetry",5.0);
	tmp.setValue("emg:retention",725.0);
	em1.setParameters(tmp);

	ElutionModel em2;
	em2 = em1;
	
	ElutionModel em3;
	em3.setInterpolationStep(0.9);
	em3.setParameters(tmp);
  TEST_EQUAL(em3.getParameters(), em2.getParameters())
	TEST_EQUAL(em3==em2,true)
}
END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))
{
	ElutionModel em1;
	
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

	ElutionModel em2;
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
}
END_SECTION

START_SECTION((void setSamples()))
{ 
	// dummy subtest only since already tested above
	TEST_EQUAL(1,1)
}
END_SECTION

START_SECTION((CoordinateType getCenter() const ))
{
	TOLERANCE_ABSOLUTE(0.001)
	ElutionModel em1;
	
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
}
END_SECTION

START_SECTION((static BaseModel<1>* create()))
{
	BaseModel<1>* ptr = ElutionModel::create();
	TEST_EQUAL(ptr->getName(), "ElutionModel")
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(ElutionModel::getProductName(),"ElutionModel")
	TEST_EQUAL(ElutionModel().getName(),"ElutionModel")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



