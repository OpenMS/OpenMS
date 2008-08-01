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
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <string>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(FeatureMap<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


FeatureMap<>* pl_ptr = 0;
CHECK((FeatureMap()))
	pl_ptr = new FeatureMap<>();
	TEST_NOT_EQUAL(pl_ptr, 0)

	TEST_EQUAL(pl_ptr->getMin(), FeatureMap<>::PositionType::max)
	TEST_EQUAL(pl_ptr->getMax(), FeatureMap<>::PositionType::min_negative)
	TEST_REAL_EQUAL(pl_ptr->getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_EQUAL(pl_ptr->getMaxInt(), -numeric_limits<DoubleReal>::max())
RESULT

CHECK((~FeatureMap()))
	delete pl_ptr;
RESULT

Feature feature1;
feature1.getPosition()[0] = 2.0;
feature1.getPosition()[1] = 3.0;
feature1.setIntensity(1.0);

Feature feature2;
feature2.getPosition()[0] = 0.0;
feature2.getPosition()[1] = 2.5;
feature2.setIntensity(0.5);

Feature feature3;
feature3.getPosition()[0] = 10.5;
feature3.getPosition()[1] = 0.0;
feature3.setIntensity(0.01);

//feature with convex hulls
Feature feature4;
feature4.getPosition()[0] = 5.25;
feature4.getPosition()[1] = 1.5;
feature4.setIntensity(0.5);
std::vector< ConvexHull2D > hulls(1);
hulls[0].addPoint(DPosition<2>(-1.0,2.0));
hulls[0].addPoint(DPosition<2>(4.0,1.2));
hulls[0].addPoint(DPosition<2>(5.0,3.123));
feature4.setConvexHulls(hulls);

CHECK( void updateRanges() )
	//test without convex hulls
  FeatureMap<> s;
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
  
  //test with convex hull
  s.push_back(feature4);
  s.updateRanges();
  TEST_REAL_EQUAL(s.getMaxInt(),1.0)
  TEST_REAL_EQUAL(s.getMinInt(),0.01)
  TEST_REAL_EQUAL(s.getMax()[0],10.5)
  TEST_REAL_EQUAL(s.getMax()[1],3.123)
  TEST_REAL_EQUAL(s.getMin()[0],-1.0)
  TEST_REAL_EQUAL(s.getMin()[1],0.0)
	
RESULT

CHECK((FeatureMap(const FeatureMap& map)))
	FeatureMap<> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
		
	FeatureMap<> map2(map1);
	
	TEST_EQUAL(map2.size(),3);
  TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
  TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
RESULT

CHECK((FeatureMap& operator = (const FeatureMap& rhs)))
	FeatureMap<> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
	
	//assignment
	FeatureMap<> map2;
	map2 = map1;
	
	TEST_EQUAL(map2.size(),3);
    TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
    TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
  	
    //assignment of empty object
     map2 = FeatureMap<>();
	
	TEST_EQUAL(map2.size(),0);
	TEST_REAL_EQUAL(map2.getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_EQUAL(map2.getMaxInt(), -numeric_limits<DoubleReal>::max())
  TEST_EQUAL(map2.getType(),ExperimentalSettings::UNKNOWN)
RESULT

CHECK((bool operator == (const FeatureMap& rhs) const))
	FeatureMap<> empty,edit;
	
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

CHECK((bool operator != (const FeatureMap& rhs) const))
	FeatureMap<> empty,edit;
	
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
	
	FeatureMap<> to_be_sorted;
	
	Feature f1;
	f1.setIntensity(10);
	to_be_sorted.push_back(f1);
	
	Feature f2;
	f2.setIntensity(5);
	to_be_sorted.push_back(f2);
	
	Feature f3;
	f3.setIntensity(3);
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByIntensity();
	
	TEST_EQUAL(to_be_sorted[0].getIntensity(),3);
	TEST_EQUAL(to_be_sorted[1].getIntensity(),5);
	TEST_EQUAL(to_be_sorted[2].getIntensity(),10);
	
RESULT

CHECK((void sortByPosition()))
	
	FeatureMap<> to_be_sorted;
	
	Feature f1;
	f1.getPosition()[0] = 10;
	to_be_sorted.push_back(f1);
	
	Feature f2;
	f2.getPosition()[0] = 5;
	to_be_sorted.push_back(f2);
	
	Feature f3;
	f3.getPosition()[0] = 3;
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByPosition();
	
	TEST_EQUAL(to_be_sorted[0].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],5);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],10);
	
RESULT

CHECK((void sortByNthPosition(UInt i) throw(Exception::NotImplemented)))
	
	FeatureMap<> to_be_sorted;
	
	Feature f1;
	f1.getPosition()[0] = 10;
	f1.getPosition()[1] = 25;
	to_be_sorted.push_back(f1);
	
	Feature f2;
	f2.getPosition()[0] = 5;
	f2.getPosition()[1] = 15;
	to_be_sorted.push_back(f2);
	
	Feature f3;
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

CHECK(void swap(FeatureMap& from))
	FeatureMap<> fm1, fm2;
	fm1.setComment("stupid comment");
	fm1.push_back(feature1);
	fm1.push_back(feature2);
	fm1.updateRanges();
	
	fm1.swap(fm2);
	
	TEST_EQUAL(fm1.getComment(),"")
	TEST_EQUAL(fm1.size(),0)
	TEST_REAL_EQUAL(fm1.getMinInt(),DRange<1>().min()[0])

	TEST_EQUAL(fm2.getComment(),"stupid comment")
	TEST_EQUAL(fm2.size(),2)
	TEST_REAL_EQUAL(fm2.getMinInt(),0.5)
	
RESULT

CHECK(void sortByOverallQuality() )
	
	FeatureMap<> to_be_sorted;
	
	Feature f1;
	f1.getPosition()[0] = 1;
	f1.getPosition()[1] = 1;
	f1.setOverallQuality(10);
	to_be_sorted.push_back(f1);
	
	Feature f2;
	f2.getPosition()[0] = 2;
	f2.getPosition()[1] = 2;
	f2.setOverallQuality(30);
	to_be_sorted.push_back(f2);
	
	Feature f3;
	f3.getPosition()[0] = 3;
	f3.getPosition()[1] = 3;
	f3.setOverallQuality(20);
	to_be_sorted.push_back(f3);
	
	to_be_sorted.sortByOverallQuality();
	
	TEST_EQUAL(to_be_sorted[0].getPosition()[0],1);
	TEST_EQUAL(to_be_sorted[1].getPosition()[0],3);
	TEST_EQUAL(to_be_sorted[2].getPosition()[0],2);
	
	TEST_EQUAL(to_be_sorted[0].getOverallQuality(),10);
	TEST_EQUAL(to_be_sorted[1].getOverallQuality(),20);
	TEST_EQUAL(to_be_sorted[2].getOverallQuality(),30);

RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
