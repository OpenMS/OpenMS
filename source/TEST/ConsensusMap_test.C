// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMap* ptr = 0;
CHECK((ConsensusMap()))
	ptr = new ConsensusMap();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusMap()))
	delete ptr;
RESULT

CHECK((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.setFileName(0,"blub");
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy;
  cons_map_copy = cons_map;
  
  TEST_EQUAL(cons_map_copy.getFileNames()[0] == "blub", true)
RESULT

CHECK((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.setFileName(0,"blub");
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy(cons_map);
  
  TEST_EQUAL(cons_map_copy.getFileNames()[0] == "blub", true)
RESULT

CHECK((ConsensusMap(Base::size_type n)))
  ConsensusMap cons_map(5);
  
  TEST_REAL_EQUAL(cons_map.size(),5)
RESULT

CHECK((const Map<UInt,String>& getFileNames() const ))
  ConsensusMap cons_map;
  
  TEST_REAL_EQUAL(cons_map.getFileNames().size(),0)
RESULT

CHECK((void setFileName(UInt index, const String &name)))
  ConsensusMap cons_map;
  cons_map.setFileName(0,"blub");
  
  TEST_REAL_EQUAL(cons_map.getFileNames().size(),1)
  TEST_EQUAL(cons_map.getFileNames()[0] == "blub", true)
RESULT

CHECK((void merge(ConsensusMap& new_map)))
  ConsensusMap cons_map;
  ConsensusMap cons_map_2;
  vector<FeatureMap<>* > map_vector(4);
  Feature feat_1;
  feat_1.setRT(1);
  feat_1.setMZ(4);
  feat_1.setIntensity(23);
  Feature feat_2;
  feat_2.setRT(1.5);
  feat_2.setMZ(5);
  feat_2.setIntensity(23);
  Feature feat_3;
  feat_3.setRT(1.2);
  feat_3.setMZ(4.5);
  feat_3.setIntensity(23);
  Feature feat_4;
  feat_4.setRT(2.2);
  feat_4.setMZ(4.8);
  feat_4.setIntensity(23);
  FeatureMap<> feat_map_1, feat_map_2,feat_map_3,feat_map_4;
  feat_map_1.push_back(feat_1);
  feat_map_2.push_back(feat_2);
  feat_map_3.push_back(feat_3);
  feat_map_4.push_back(feat_4);
  map_vector[0] = &feat_map_1;
  map_vector[1] = &feat_map_2;
  map_vector[2] = &feat_map_3;
  map_vector[3] = &feat_map_4;
  
  ConsensusFeature cons_1(0,0,feat_1);
  cons_1.insert(1,0,feat_2);
  ConsensusFeature cons_2(2,0,feat_3);
  cons_2.insert(3,0,feat_4);
  cons_map.push_back(cons_1);
  cons_map.push_back(cons_2);
   
  TEST_REAL_EQUAL(cons_map.size(),2)
  cons_map.merge(cons_map_2);
  TEST_REAL_EQUAL(cons_map_2.size(),1)
RESULT

CHECK((bool isValid() const))
	ConsensusMap cm;
	//empty map
	TEST_EQUAL(cm.isValid(),true)
	//one, valid feature
	ConsensusFeature f;
	f.insert(1,1,Feature());
	cm.push_back(f);
	cm.setFileName(1,"bla");
	TEST_EQUAL(cm.isValid(),true)
	//one valid and one invalid feature
	f.insert(2,1,Feature());
	cm.push_back(f);
	TEST_EQUAL(cm.isValid(),false)
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



