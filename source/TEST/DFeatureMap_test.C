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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DFeature.h>
#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DFeatureMap<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DFeatureMap<2>* pl_ptr = 0;
CHECK((DFeatureMap()))
	pl_ptr = new DFeatureMap<2>();
	TEST_NOT_EQUAL(pl_ptr, 0)

	TEST_EQUAL(pl_ptr->getMin(), DFeatureMap<2>::PositionType::max)
	TEST_EQUAL(pl_ptr->getMax(), DFeatureMap<2>::PositionType::min_negative)
	TEST_REAL_EQUAL(pl_ptr->getMinInt(), numeric_limits<DFeatureMap<2>::IntensityType>::max())
	TEST_REAL_EQUAL(pl_ptr->getMaxInt(), -numeric_limits<DFeatureMap<2>::IntensityType>::max())
RESULT

CHECK((~DFeatureMap()))
	delete pl_ptr;
RESULT

DFeature<2> feature1;
feature1.getPosition()[0] = 2.0;
feature1.getPosition()[1] = 3.0;
feature1.getIntensity() = 1.0;

DFeature<2> feature2;
feature2.getPosition()[0] = 0.0;
feature2.getPosition()[1] = 2.5;
feature2.getIntensity() = 0.5;

DFeature<2> feature3;
feature3.getPosition()[0] = 10.5;
feature3.getPosition()[1] = 0.0;
feature3.getIntensity() = 0.01;

CHECK((updateRanges_(this->begin(), this->end())))
  DFeatureMap<2> s;
  s.push_back(feature1);
  s.push_back(feature2);
  s.push_back(feature3);
  
  s.updateRanges();
  s.updateRanges(); //second time to check the initialization
   
  TEST_REAL_EQUAL(s.getMaxInt(),1.0)
  TEST_REAL_EQUAL(s.getMinInt(),0.01)
  TEST_REAL_EQUAL(s.getMax()[0],10.5)
  TEST_REAL_EQUAL(s.getMax()[1],3.0)
  TEST_REAL_EQUAL(s.getMin()[0],0.0)
  TEST_REAL_EQUAL(s.getMin()[1],0.0)
RESULT

CHECK((DFeatureMap(const DFeatureMap& map)))
	DFeatureMap<2> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
		
	DFeatureMap<2> map2(map1);
	
	TEST_EQUAL(map2.size(),3);
  TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
  TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
RESULT

CHECK((DFeatureMap& operator = (const DFeatureMap& rhs)))
	DFeatureMap<2> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
	
	//assignment
	DFeatureMap<2> map2;
	map2 = map1;
	
	TEST_EQUAL(map2.size(),3);
    TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
    TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
  	
    //assignment of empty object
     map2 = DFeatureMap<2>();
	
	TEST_EQUAL(map2.size(),0);
	TEST_REAL_EQUAL(map2.getMinInt(), numeric_limits<DFeatureMap<2>::IntensityType>::max())
	TEST_REAL_EQUAL(map2.getMaxInt(), -numeric_limits<DFeatureMap<2>::IntensityType>::max())
  TEST_EQUAL(map2.getType(),ExperimentalSettings::UNKNOWN)
RESULT

CHECK((bool operator == (const DFeatureMap& rhs) const))
	DFeatureMap<2> empty,edit;
	
	TEST_EQUAL(empty==edit, true);
	
	edit.setType(ExperimentalSettings::MS);
	TEST_EQUAL(empty==edit, false);
	
	edit = empty;
	edit.push_back(feature1);
	TEST_EQUAL(empty==edit, false);

	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty==edit, false);
RESULT

CHECK((bool operator != (const DFeatureMap& rhs) const))
	DFeatureMap<2> empty,edit;
	
	TEST_EQUAL(empty!=edit, false);
	
	edit.setType(ExperimentalSettings::MS);
	TEST_EQUAL(empty!=edit, true);
	
	edit = empty;
	edit.push_back(feature1);
	TEST_EQUAL(empty!=edit, true);

	edit = empty;
	edit.push_back(feature1);
	edit.push_back(feature2);
	edit.updateRanges();
	edit.clear();
	TEST_EQUAL(empty!=edit, true);
RESULT


CHECK((void sortByIntensity()))
	
	DFeatureMap<2> to_be_sorted;
	
	DFeature<2> f1;
	f1.setIntensity(10);
	to_be_sorted.push_back(f1);
	
	DFeature<2> f2;
	f2.setIntensity(5);
	to_be_sorted.push_back(f2);
	
	DFeature<2> f3;
	f3.setIntensity(3);
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByIntensity();
	
	TEST_EQUAL(to_be_sorted[0].getIntensity(),3);
	TEST_EQUAL(to_be_sorted[1].getIntensity(),5);
	TEST_EQUAL(to_be_sorted[2].getIntensity(),10);
	
RESULT

CHECK((void sortByPosition()))
	
	DFeatureMap<2> to_be_sorted;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 10;
	to_be_sorted.push_back(f1);
	
	DFeature<2> f2;
	f2.getPosition()[0] = 5;
	to_be_sorted.push_back(f2);
	
	DFeature<2> f3;
	f3.getPosition()[0] = 3;
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByPosition();
	
	TEST_EQUAL(to_be_sorted[0].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],5);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],10);
	
RESULT

CHECK((void sortByNthPosition(UnsignedInt i) throw(Exception::NotImplemented)))
	
	DFeatureMap<2> to_be_sorted;
	
	DFeature<2> f1;
	f1.getPosition()[0] = 10;
	f1.getPosition()[1] = 25;
	to_be_sorted.push_back(f1);
	
	DFeature<2> f2;
	f2.getPosition()[0] = 5;
	f2.getPosition()[1] = 15;
	to_be_sorted.push_back(f2);
	
	DFeature<2> f3;
	f3.getPosition()[0] = 3;
	f3.getPosition()[1] = 10;
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByNthPosition(0);
	
	TEST_EQUAL(to_be_sorted[0].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],5);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],10);
	
	to_be_sorted.sortByNthPosition(1);
	
	TEST_EQUAL(to_be_sorted[0].getPosition()[1],10);
	TEST_EQUAL(to_be_sorted[1].getPosition()[1],15);
	TEST_EQUAL(to_be_sorted[2].getPosition()[1],25);
	
RESULT

CHECK(void sortByOverallQuality() )
	

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
