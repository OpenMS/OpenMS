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
CHECK(DFeatureMap<2>())
	pl_ptr = new DFeatureMap<2>();
	TEST_NOT_EQUAL(pl_ptr, 0)

	TEST_EQUAL(pl_ptr->getMin(), DFeatureMap<2>::PositionType::max)
	TEST_EQUAL(pl_ptr->getMax(), DFeatureMap<2>::PositionType::min_negative)
	TEST_REAL_EQUAL(pl_ptr->getMinInt(), numeric_limits<DFeatureMap<2>::IntensityType>::max())
	TEST_REAL_EQUAL(pl_ptr->getMaxInt(), -numeric_limits<DFeatureMap<2>::IntensityType>::max())
RESULT

CHECK(~DFeatureMap<2>())
	delete pl_ptr;
RESULT

CHECK(const String& getName() const)
	DFeatureMap<2> tmp;
	TEST_EQUAL(tmp.getName(), "")
RESULT

CHECK(void setName(const String& name))
	DFeatureMap<2> tmp;
	tmp.setName("TEST");
	TEST_EQUAL(tmp.getName(), "TEST")
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

CHECK(void updateRanges())
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

CHECK(DFeatureMap<2>(const DFeatureMap& p))
	DFeatureMap<2> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.setName("test");
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
		
	DFeatureMap<2> map2(map1);
	
	TEST_EQUAL(map2.size(),3);
	TEST_EQUAL(map2.getName(),"test");
  TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
  TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
RESULT

CHECK(DFeatureMap& operator = (const DFeatureMap& rhs))
	DFeatureMap<2> map1;
	map1.push_back(feature1);
	map1.push_back(feature2);
	map1.push_back(feature3);
	map1.setName("test");
	map1.updateRanges();
	map1.setType(ExperimentalSettings::MS);
	
	//assignment
	DFeatureMap<2> map2;
	map2 = map1;
	
	TEST_EQUAL(map2.size(),3);
	TEST_EQUAL(map2.getName(),"test");
  TEST_REAL_EQUAL(map2.getMaxInt(),1.0)
  TEST_EQUAL(map2.getType(),ExperimentalSettings::MS)
  	
  //assignment of empty object
  map2 = DFeatureMap<2>();
	
	TEST_EQUAL(map2.size(),0);
	TEST_EQUAL(map2.getName(),"");
	TEST_REAL_EQUAL(map2.getMinInt(), numeric_limits<DFeatureMap<2>::IntensityType>::max())
	TEST_REAL_EQUAL(map2.getMaxInt(), -numeric_limits<DFeatureMap<2>::IntensityType>::max())
  TEST_EQUAL(map2.getType(),ExperimentalSettings::UNKNOWN)
RESULT

CHECK(bool operator == (const DFeatureMap& rhs) const)
	DFeatureMap<2> empty,edit;
	
	TEST_EQUAL(empty==edit, true);
	
	edit.setType(ExperimentalSettings::MS);
	TEST_EQUAL(empty==edit, false);
	
	edit = empty;
	edit.setName("TEST");
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

CHECK(bool operator != (const DFeatureMap& rhs) const)
	DFeatureMap<2> empty,edit;
	
	TEST_EQUAL(empty!=edit, false);
	
	edit.setType(ExperimentalSettings::MS);
	TEST_EQUAL(empty!=edit, true);
	
	edit = empty;
	edit.setName("TEST");
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
