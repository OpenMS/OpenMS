// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FILTERING/DATAREDUCTION/MaxReducer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/FileHandler.h>
using namespace OpenMS;

START_TEST(MaxReducer, "$Id: MaxReducer_test.C 463  $")

CHECK((virtual void applyReduction(const  MSExperiment<>& , MSExperiment<>&  )))
	
  DTA2DFile dta;
  MSExperiment<> in;
  MSExperiment<> out;
  FileHandler().loadExperiment("data/MaxReducer_test.dta2d",in);
  in.setName("MaxReducer_test.dta2d");
  Param param;	
  MaxReducer maxreducer;
  param.setValue("Ratio",20);
  maxreducer.setParameter(param);
  in.updateRanges();
  maxreducer.applyReduction(in,out);

   TEST_EQUAL(out[0].size(),5)
   TEST_EQUAL(out[1].size(),5)
   TEST_EQUAL(out[2].size(),5)
	 TEST_EQUAL(out[0].getRetentionTime(),1)
   TEST_EQUAL(out[1].getRetentionTime(),2)
   TEST_EQUAL(out[2].getRetentionTime(),3)

	 TEST_EQUAL(out[0].getContainer()[0].getPosition()[0],2)
   TEST_EQUAL(out[0].getContainer()[1].getPosition()[0],4)
   TEST_EQUAL(out[0].getContainer()[2].getPosition()[0],6)
	 TEST_EQUAL(out[0].getContainer()[3].getPosition()[0],8)
	 TEST_EQUAL(out[0].getContainer()[4].getPosition()[0],10)
   
   TEST_EQUAL(out[1].getContainer()[0].getPosition()[0],2)
   TEST_EQUAL(out[1].getContainer()[1].getPosition()[0],4)
   TEST_EQUAL(out[1].getContainer()[2].getPosition()[0],6)
	 TEST_EQUAL(out[1].getContainer()[3].getPosition()[0],8)
	 TEST_EQUAL(out[1].getContainer()[4].getPosition()[0],10)
   
   TEST_EQUAL(out[2].getContainer()[0].getPosition()[0],2)
   TEST_EQUAL(out[2].getContainer()[1].getPosition()[0],4)
   TEST_EQUAL(out[2].getContainer()[2].getPosition()[0],6)
	 TEST_EQUAL(out[2].getContainer()[3].getPosition()[0],8)
	 TEST_EQUAL(out[2].getContainer()[4].getPosition()[0],10)

	 TEST_EQUAL(out[0].getContainer()[0].getIntensity(),2)
   TEST_EQUAL(out[0].getContainer()[1].getIntensity(),4)
   TEST_EQUAL(out[0].getContainer()[2].getIntensity(),6)
	 TEST_EQUAL(out[0].getContainer()[3].getIntensity(),8)
	 TEST_EQUAL(out[0].getContainer()[4].getIntensity(),10)
   
   TEST_EQUAL(out[1].getContainer()[0].getIntensity(),2)
   TEST_EQUAL(out[1].getContainer()[1].getIntensity(),4)
   TEST_EQUAL(out[1].getContainer()[2].getIntensity(),6)
	 TEST_EQUAL(out[1].getContainer()[3].getIntensity(),8)
	 TEST_EQUAL(out[1].getContainer()[4].getIntensity(),10)
   
   TEST_EQUAL(out[2].getContainer()[0].getIntensity(),2)
   TEST_EQUAL(out[2].getContainer()[1].getIntensity(),4)
   TEST_EQUAL(out[2].getContainer()[2].getIntensity(),6)
	 TEST_EQUAL(out[2].getContainer()[3].getIntensity(),8)
	 TEST_EQUAL(out[2].getContainer()[4].getIntensity(),10)


RESULT

MaxReducer * ptr = 0;
CHECK(static const  String getName())
	 	ptr  = new MaxReducer();
    TEST_EQUAL(ptr->getName(),"MaxReducer")	
RESULT

CHECK(MaxReducer())
 	ptr  = new MaxReducer();
  TEST_NOT_EQUAL(ptr,0)
RESULT

CHECK(~MaxReducer())
			delete ptr;
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
