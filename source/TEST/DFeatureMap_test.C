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
// $Id: DFeatureMap_test.C,v 1.6 2006/04/18 15:22:54 ole_st Exp $
// $Author: ole_st $
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

START_TEST(DFeatureMap<D>, "$Id: DFeatureMap_test.C,v 1.6 2006/04/18 15:22:54 ole_st Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DFeatureMap<2>* pl_ptr = 0;
CHECK(DFeatureMap<2>())
	pl_ptr = new DFeatureMap<2>();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK(~DFeatureMap<2>())
	delete pl_ptr;
RESULT

CHECK(DFeatureMap<2>(const DFeatureMap& p))
	DFeatureMap<2> map1;
	DFeature<2> feature;
	DPosition<2> pos;
	pos[0] = 1.0;
	pos[1] = 2.0;
	feature.getIntensity() = 1.0;
	feature.getPosition() = pos;
	map1.push_back(feature);
	TEST_EQUAL(map1.size(), 1)
	feature.getIntensity() = 3.0;
	map1.push_back(feature);
		
	DFeatureMap<2> map2(map1);
	TEST_EQUAL(map2.size(), 2)
	
	DFeature<2> feat1 = map2.back();
	map2.pop_back();
	TEST_EQUAL(feat1.getIntensity(),3.0);
	
	DFeature<2> feat2 = map2.back();
	map2.pop_back();
	TEST_EQUAL(feat2.getIntensity(),1.0);
	TEST_EQUAL(feat2.getPosition()[0],1.0);
	TEST_EQUAL(feat2.getPosition()[1],2.0);
	
RESULT

DFeatureMap<2> map;

CHECK( empty() const)
	TEST_EQUAL(map.empty(), true)
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

CHECK( size() const)
	TEST_EQUAL(map.size(), 0)
	
	map.push_back(feature1);
	TEST_EQUAL(map.size(), 1)

	map.push_back(feature2);
	TEST_EQUAL(map.size(), 2)

	map.push_back(feature3);
	TEST_EQUAL(map.size(), 3)
RESULT

CHECK( empty() const)
	TEST_EQUAL(map.empty(), false)
RESULT

CHECK(DFeatureMap& operator = (const DFeatureMap& rhs))
	DFeatureMap<2> edit,empty;
	
	// assignment of an empty object
	DFeatureMap<2>::FeatureType f;
	f.getPosition()[0] = 1.0;
	edit.push_back(f);
	edit.getSample().setName("TEST");
	edit = empty;
	TEST_REAL_EQUAL(edit.size(),0)
	TEST_EQUAL(edit.getSample().getName(),empty.getSample().getName())	
	
	// normal assignment
	edit.push_back(f);
	edit.getSample().setName("TEST2");
	empty = edit;
	TEST_REAL_EQUAL(edit[0].getPosition()[0],empty[0].getPosition()[0])
	TEST_EQUAL(edit.getSample().getName(),empty.getSample().getName())
RESULT

CHECK(bool operator == (const DFeatureMap& rhs) const)
	DFeatureMap<1> empty,edit;
	
	TEST_EQUAL(empty==edit, true);
	
	edit.getSample().setName("TEST");
	TEST_EQUAL(empty==edit, false);
	
	edit = empty;
	DFeatureMap<1>::FeatureType f;
	f.getPosition()[0] = 1.0;
	edit.push_back(f);
	TEST_EQUAL(empty==edit, false);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
