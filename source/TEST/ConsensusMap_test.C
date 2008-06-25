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
// $Maintainer: Clemens Groepl $
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
	TEST_REAL_EQUAL(ptr->getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_EQUAL(ptr->getMaxInt(), -numeric_limits<DoubleReal>::max())
RESULT

CHECK((~ConsensusMap()))
	delete ptr;
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

Feature feature4;
feature4.getPosition()[0] = 5.25;
feature4.getPosition()[1] = 1.5;
feature4.setIntensity(0.5);

CHECK(void updateRanges())
  ConsensusMap map;
	ConsensusFeature f;
	f.setIntensity(1.0);
	f.setRT(2.0);
	f.setMZ(3.0);
	f.insert(1,1,feature1);
	map.push_back(f);
  
  map.updateRanges();
  TEST_REAL_EQUAL(map.getMaxInt(),1.0)
  TEST_REAL_EQUAL(map.getMinInt(),1.0)
  TEST_REAL_EQUAL(map.getMax()[0],2.0)
  TEST_REAL_EQUAL(map.getMax()[1],3.0)
  TEST_REAL_EQUAL(map.getMin()[0],2.0)
  TEST_REAL_EQUAL(map.getMin()[1],3.0)
  
  //second time to check the initialization
  map.updateRanges();
   
  TEST_REAL_EQUAL(map.getMaxInt(),1.0)
  TEST_REAL_EQUAL(map.getMinInt(),1.0)
  TEST_REAL_EQUAL(map.getMax()[0],2.0)
  TEST_REAL_EQUAL(map.getMax()[1],3.0)
  TEST_REAL_EQUAL(map.getMin()[0],2.0)
  TEST_REAL_EQUAL(map.getMin()[1],3.0)
  
  //two points
	f.insert(1,2,feature2);
	map.push_back(f);
	map.updateRanges();
	
  TEST_REAL_EQUAL(map.getMaxInt(),1.0)
  TEST_REAL_EQUAL(map.getMinInt(),0.5)
  TEST_REAL_EQUAL(map.getMax()[0],2.0)
  TEST_REAL_EQUAL(map.getMax()[1],3.0)
  TEST_REAL_EQUAL(map.getMin()[0],0.0)
  TEST_REAL_EQUAL(map.getMin()[1],2.5)
  
	//four points
	f.insert(1,3,feature3);
	f.insert(1,4,feature4);
	map.push_back(f);
	map.updateRanges();
	
  TEST_REAL_EQUAL(map.getMaxInt(),1.0)
  TEST_REAL_EQUAL(map.getMinInt(),0.01)
  TEST_REAL_EQUAL(map.getMax()[0],10.5)
  TEST_REAL_EQUAL(map.getMax()[1],3.0)
  TEST_REAL_EQUAL(map.getMin()[0],0.0)
  TEST_REAL_EQUAL(map.getMin()[1],0.0)
  	
RESULT

CHECK((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.getFileDescriptions()[0].filename = "blub";
  cons_map.getFileDescriptions()[0].size = 47;
  cons_map.getFileDescriptions()[0].label = "label";
	cons_map.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy;
  cons_map_copy = cons_map;
  
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
RESULT

CHECK((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap cons_map;
  cons_map.getFileDescriptions()[0].filename = "blub";
  cons_map.getFileDescriptions()[0].size = 47;
  cons_map.getFileDescriptions()[0].label = "label";
	cons_map.getFileDescriptions()[0].setMetaValue("meta",String("meta"));
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap cons_map_copy(cons_map);
  
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].filename == "blub", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].label == "label", true)
  TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].size == 47, true)
	TEST_EQUAL(cons_map_copy.getFileDescriptions()[0].getMetaValue("meta") == "meta", true)
RESULT

CHECK((ConsensusMap(Base::size_type n)))
  ConsensusMap cons_map(5);
  
  TEST_REAL_EQUAL(cons_map.size(),5)
RESULT

CHECK((const FileDescriptions& getFileDescriptions() const ))
  ConsensusMap cons_map;
  
  TEST_REAL_EQUAL(cons_map.getFileDescriptions().size(),0)
RESULT

CHECK((bool isValid() const))
	ConsensusMap cm;
	//empty map
	TEST_EQUAL(cm.isValid(),true)
	//one, valid feature
	ConsensusFeature f;
	f.insert(1,1,Feature());
	cm.push_back(f);
  cm.getFileDescriptions()[1].filename = "bla";
	cm.getFileDescriptions()[1].size = 5;
	TEST_EQUAL(cm.isValid(),true)
	//two valid features
	f.insert(1,4,Feature());
	cm.push_back(f);
	TEST_EQUAL(cm.isValid(),true)
	//two valid and one invalid feature (map index)
	f.insert(2,1,Feature());
	cm.push_back(f);
	TEST_EQUAL(cm.isValid(),false)
	//one invalid feature (element index)
	cm.clear();
	ConsensusFeature f2;
	f2.insert(2,1,Feature());
	cm.push_back(f2);
	
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



