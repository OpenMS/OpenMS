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
// $Maintainer:Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////////

#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FileHandler.h>
using namespace OpenMS;

START_TEST(MaxReducer, "$Id$")

MaxReducer * ptr = 0;

CHECK((MaxReducer()))
 	ptr  = new MaxReducer();
  TEST_NOT_EQUAL(ptr,0)
RESULT

CHECK((~MaxReducer()))
			delete ptr;
RESULT

CHECK((static const String getProductName()))
	 	MaxReducer m;
    TEST_EQUAL(m.getName(),"max_reducer")	
RESULT

CHECK((void applyReduction(const ExperimentType& in, ExperimentType& out )))
	
  DTA2DFile dta;
  MSExperiment<> in;
  MSExperiment<> out;
  FileHandler().loadExperiment("data/MaxReducer_test.dta2d",in);
  Param param;	
  MaxReducer maxreducer;
  param.setValue("peaks_per_step",3);
  maxreducer.setParameters(param);
  in.updateRanges();
  maxreducer.applyReduction(in,out);
	TEST_EQUAL(out.size(),3)
	TEST_EQUAL(out[0].size(),4)
	TEST_EQUAL(out[1].size(),4)
 	TEST_EQUAL(out[2].size(),4)
	TEST_REAL_EQUAL(out[0].getRT(),1)
 	TEST_REAL_EQUAL(out[1].getRT(),2)
 	TEST_REAL_EQUAL(out[2].getRT(),3)	
	TEST_REAL_EQUAL(out[0].getContainer()[0].getPosition()[0],3)
 	TEST_REAL_EQUAL(out[0].getContainer()[1].getPosition()[0],6)
 	TEST_REAL_EQUAL(out[0].getContainer()[2].getPosition()[0],9)
 	TEST_REAL_EQUAL(out[0].getContainer()[3].getPosition()[0],10)

  TEST_REAL_EQUAL(out[1].getContainer()[0].getPosition()[0],3)
 	TEST_REAL_EQUAL(out[1].getContainer()[1].getPosition()[0],6)
 	TEST_REAL_EQUAL(out[1].getContainer()[2].getPosition()[0],9)
 	TEST_REAL_EQUAL(out[1].getContainer()[3].getPosition()[0],11)
  TEST_REAL_EQUAL(out[2].getContainer()[0].getPosition()[0],3)
 	TEST_REAL_EQUAL(out[2].getContainer()[1].getPosition()[0],6)
 	TEST_REAL_EQUAL(out[2].getContainer()[2].getPosition()[0],9)
 	TEST_REAL_EQUAL(out[2].getContainer()[3].getPosition()[0],12)

// 	TEST_EQUAL(out.size(),3)
// 	TEST_EQUAL(out[0].size(),5)
// 	TEST_EQUAL(out[1].size(),6)
// 	TEST_EQUAL(out[2].size(),6)
// 	TEST_REAL_EQUAL(out[0].getRT(),1)
// 	TEST_REAL_EQUAL(out[1].getRT(),2)
// 	TEST_REAL_EQUAL(out[2].getRT(),3)
	
// 	TEST_REAL_EQUAL(out[0].getContainer()[0].getPosition()[0],2)
// 	TEST_REAL_EQUAL(out[0].getContainer()[1].getPosition()[0],4)
// 	TEST_REAL_EQUAL(out[0].getContainer()[2].getPosition()[0],6)
// 	TEST_REAL_EQUAL(out[0].getContainer()[3].getPosition()[0],8)
// 	TEST_REAL_EQUAL(out[0].getContainer()[4].getPosition()[0],10)
	
// 	TEST_REAL_EQUAL(out[1].getContainer()[0].getPosition()[0],2)
// 	TEST_REAL_EQUAL(out[1].getContainer()[1].getPosition()[0],4)
// 	TEST_REAL_EQUAL(out[1].getContainer()[2].getPosition()[0],6)
// 	TEST_REAL_EQUAL(out[1].getContainer()[3].getPosition()[0],8)
// 	TEST_REAL_EQUAL(out[1].getContainer()[4].getPosition()[0],10)
// 	TEST_REAL_EQUAL(out[1].getContainer()[5].getPosition()[0],11)
	
// 	TEST_REAL_EQUAL(out[2].getContainer()[0].getPosition()[0],2)
// 	TEST_REAL_EQUAL(out[2].getContainer()[1].getPosition()[0],4)
// 	TEST_REAL_EQUAL(out[2].getContainer()[2].getPosition()[0],6)
// 	TEST_REAL_EQUAL(out[2].getContainer()[3].getPosition()[0],8)
// 	TEST_REAL_EQUAL(out[2].getContainer()[4].getPosition()[0],10)
// 	TEST_REAL_EQUAL(out[2].getContainer()[5].getPosition()[0],12)
	
// 	TEST_REAL_EQUAL(out[0].getContainer()[0].getIntensity(),2)
// 	TEST_REAL_EQUAL(out[0].getContainer()[1].getIntensity(),4)
// 	TEST_REAL_EQUAL(out[0].getContainer()[2].getIntensity(),6)
// 	TEST_REAL_EQUAL(out[0].getContainer()[3].getIntensity(),8)
// 	TEST_REAL_EQUAL(out[0].getContainer()[4].getIntensity(),10)
	
// 	TEST_REAL_EQUAL(out[1].getContainer()[0].getIntensity(),2)
// 	TEST_REAL_EQUAL(out[1].getContainer()[1].getIntensity(),4)
// 	TEST_REAL_EQUAL(out[1].getContainer()[2].getIntensity(),6)
// 	TEST_REAL_EQUAL(out[1].getContainer()[3].getIntensity(),8)
// 	TEST_REAL_EQUAL(out[1].getContainer()[4].getIntensity(),10)
// 	TEST_REAL_EQUAL(out[1].getContainer()[5].getIntensity(),11)
	
// 	TEST_REAL_EQUAL(out[2].getContainer()[0].getIntensity(),2)
// 	TEST_REAL_EQUAL(out[2].getContainer()[1].getIntensity(),4)
// 	TEST_REAL_EQUAL(out[2].getContainer()[2].getIntensity(),6)
// 	TEST_REAL_EQUAL(out[2].getContainer()[3].getIntensity(),8)
// 	TEST_REAL_EQUAL(out[2].getContainer()[4].getIntensity(),10)
// 	TEST_REAL_EQUAL(out[2].getContainer()[5].getIntensity(),12)
	
RESULT

CHECK(static DataReducer* create())
	DataReducer* ptr2 = MaxReducer::create();
	TEST_EQUAL("max_reducer",ptr2->getName());
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
