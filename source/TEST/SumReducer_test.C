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
// $Maintainer:Cornelia Friedle $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////////

#include <OpenMS/FILTERING/DATAREDUCTION/SumReducer.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FileHandler.h>
using namespace OpenMS;

START_TEST(SumReducer, "$Id$")

SumReducer * ptr = 0;
CHECK((SumReducer()))
		ptr  = new SumReducer();
    TEST_NOT_EQUAL(ptr,0)
RESULT

CHECK((~SumReducer()))
			delete ptr;
RESULT

CHECK((static const String getName()))
	SumReducer s;
  TEST_EQUAL(s.getName(),"SumReducer")	
RESULT

CHECK((void applyReduction(const ExperimentType& in, ExperimentType& out )))
	DTA2DFile dta;
  MSExperiment<> in;
  MSExperiment<> out;
  FileHandler().loadExperiment("data/SumReducer_test.dta2d",in);
  Param param;	
  SumReducer sumreducer;
  param.setValue("Rangeperstep",0.5);
  sumreducer.setParameters(param);
  in.updateRanges();
  sumreducer.applyReduction(in,out);
  TEST_EQUAL(in.size(),3)
	TEST_EQUAL(out.size(),3)
	TEST_EQUAL(out[0].size(),2)
	TEST_EQUAL(out[1].size(),2)
	TEST_EQUAL(out[2].size(),2)
	TEST_REAL_EQUAL(out[0].getContainer()[0].getPosition()[0],5)
	TEST_REAL_EQUAL(out[0].getContainer()[1].getPosition()[0],10)
	TEST_REAL_EQUAL(out[1].getContainer()[0].getPosition()[0],6)
	TEST_REAL_EQUAL(out[1].getContainer()[1].getPosition()[0],11)
	TEST_REAL_EQUAL(out[2].getContainer()[0].getPosition()[0],12)
	TEST_REAL_EQUAL(out[2].getContainer()[1].getPosition()[0],35)
	TEST_EQUAL(out[0].getContainer()[0].getIntensity(),15)
	TEST_EQUAL(out[0].getContainer()[1].getIntensity(),40)
	TEST_REAL_EQUAL(out[1].getContainer()[0].getIntensity(),21)
	TEST_REAL_EQUAL(out[1].getContainer()[1].getIntensity(),45)
	TEST_REAL_EQUAL(out[2].getContainer()[0].getIntensity(),93)
	TEST_REAL_EQUAL(out[2].getContainer()[1].getIntensity(),40)

// 	TEST_EQUAL(in.size(),3)
// 	TEST_EQUAL(out.size(),3)
// 	TEST_EQUAL(out[0].size(),5)
// 	TEST_EQUAL(out[1].size(),4)
// 	TEST_EQUAL(out[2].size(),5)
// 	TEST_REAL_EQUAL(out[0].getRetentionTime(),1)
// 	TEST_REAL_EQUAL(out[1].getRetentionTime(),2)
// 	TEST_REAL_EQUAL(out[2].getRetentionTime(),3)
	
// 	TEST_REAL_EQUAL(out[0].getContainer()[0].getPosition()[0],2)
// 	TEST_REAL_EQUAL(out[0].getContainer()[1].getPosition()[0],4)
// 	TEST_REAL_EQUAL(out[0].getContainer()[2].getPosition()[0],6)
// 	TEST_REAL_EQUAL(out[0].getContainer()[3].getPosition()[0],8)
// 	TEST_REAL_EQUAL(out[0].getContainer()[4].getPosition()[0],10)
	
// 	TEST_REAL_EQUAL(out[1].getContainer()[0].getPosition()[0],3)
// 	TEST_REAL_EQUAL(out[1].getContainer()[1].getPosition()[0],6)
// 	TEST_REAL_EQUAL(out[1].getContainer()[2].getPosition()[0],9)
// 	TEST_REAL_EQUAL(out[1].getContainer()[3].getPosition()[0],11)
	
// 	TEST_REAL_EQUAL(out[2].getContainer()[0].getPosition()[0],7)
// 	TEST_REAL_EQUAL(out[2].getContainer()[1].getPosition()[0],12)
// 	TEST_REAL_EQUAL(out[2].getContainer()[2].getPosition()[0],20)
// 	TEST_REAL_EQUAL(out[2].getContainer()[3].getPosition()[0],31)
// 	TEST_REAL_EQUAL(out[2].getContainer()[4].getPosition()[0],35)
	
// 	TEST_REAL_EQUAL(out[0].getContainer()[0].getIntensity(),3)
// 	TEST_REAL_EQUAL(out[0].getContainer()[1].getIntensity(),7)
// 	TEST_REAL_EQUAL(out[0].getContainer()[2].getIntensity(),11)
// 	TEST_REAL_EQUAL(out[0].getContainer()[3].getIntensity(),15)
// 	TEST_REAL_EQUAL(out[0].getContainer()[4].getIntensity(),19)

// 	TEST_REAL_EQUAL(out[1].getContainer()[0].getIntensity(),6)
// 	TEST_REAL_EQUAL(out[1].getContainer()[1].getIntensity(),15)
// 	TEST_REAL_EQUAL(out[1].getContainer()[2].getIntensity(),24)
// 	TEST_REAL_EQUAL(out[1].getContainer()[3].getIntensity(),21)
	
// 	TEST_REAL_EQUAL(out[2].getContainer()[0].getIntensity(),28)
// 	TEST_REAL_EQUAL(out[2].getContainer()[1].getIntensity(),53)
// 	TEST_REAL_EQUAL(out[2].getContainer()[2].getIntensity(),18)
// 	TEST_REAL_EQUAL(out[2].getContainer()[3].getIntensity(),24)
// 	TEST_REAL_EQUAL(out[2].getContainer()[4].getIntensity(),10)
		 
RESULT

CHECK(static DataReducer* create())
	DataReducer* ptr2 = SumReducer::create();
	TEST_EQUAL("SumReducer",ptr2->getName());
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
