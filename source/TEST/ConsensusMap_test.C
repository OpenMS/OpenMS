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
RESULT

CHECK((~ConsensusMap()))
	delete ptr;
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

CHECK((const Map<UInt,FileDescription>& getFileDescriptions() const ))
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



