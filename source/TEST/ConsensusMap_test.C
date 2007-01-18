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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusMap.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMap<ConsensusFeature < FeatureMap > >* ptr = 0;
CHECK((ConsensusMap()))
	ptr = new ConsensusMap<ConsensusFeature < FeatureMap > >();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusMap()))
	delete ptr;
RESULT

CHECK((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFilenames(name_vector);
  vector<FeatureMap> map_vector(4);
  cons_map.setMapVector(map_vector);
  
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map_copy;
  cons_map_copy = cons_map;
  
  TEST_EQUAL(cons_map_copy.getFilenames()[0] == "blub", true)
  TEST_REAL_EQUAL(cons_map_copy.getMapVector().size(),4)
RESULT

CHECK((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFilenames(name_vector);
  vector<FeatureMap> map_vector(4);
  cons_map.setMapVector(map_vector);
  
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map_copy(cons_map);
  
  TEST_EQUAL(cons_map_copy.getFilenames()[0] == "blub", true)
  TEST_REAL_EQUAL(cons_map_copy.getMapVector().size(),4)
RESULT

CHECK((ConsensusMap(typename Base::size_type n)))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map(5);
  
  TEST_REAL_EQUAL(cons_map.size(),5)
RESULT

// CHECK((ConsensusMap(typename Base::size_type n, const ConsensusElementType& element)))
//   DPosition<2> pos(1,2);
//   ConsensusFeature<> cons(pos,200);
//   DFeature<2> feat;
//   feat.setPosition(pos);
//   feat.setIntensity(200);
//   
//   IndexTuple<> ind(1,3,feat);
//   cons.insert(ind); 
//   ConsensusMap<ConsensusFeature < FeatureMap > > cons_map(2, cons);
//   
//   TEST_REAL_EQUAL(cons_map.size(),2)
//   TEST_REAL_EQUAL(cons_map[0].getPosition()[0],1)
//   TEST_REAL_EQUAL(cons_map[0].getPosition()[1],2)
//   TEST_REAL_EQUAL(cons_map[0].getIntensity(),200)
//   TEST_EQUAL(cons_map[0].getPositionRange() == cons.getPositionRange(), true)
//   TEST_EQUAL(cons_map[0].getIntensityRange() == cons.getIntensityRange(), true)
//   TEST_REAL_EQUAL((cons_map[0].begin())->getMapIndex(),1)
//   TEST_REAL_EQUAL((cons_map[0].begin())->getElementIndex(),3)
//   TEST_REAL_EQUAL((cons_map[0].begin())->getElement().getIntensity(),200)
//   TEST_REAL_EQUAL(cons_map[1].getPosition()[0],1)
//   TEST_REAL_EQUAL(cons_map[1].getPosition()[1],2)
//   TEST_REAL_EQUAL(cons_map[1].getIntensity(),200)
//   TEST_EQUAL(cons_map[1].getPositionRange() == cons.getPositionRange(), true)
//   TEST_EQUAL(cons_map[1].getIntensityRange() == cons.getIntensityRange(), true)
//   TEST_REAL_EQUAL((cons_map[1].begin())->getMapIndex(),1)
//   TEST_REAL_EQUAL((cons_map[1].begin())->getElementIndex(),3)
//   TEST_REAL_EQUAL((cons_map[1].begin())->getElement().getIntensity(),200)
// RESULT

CHECK((const std::vector< String >& getFilenames() const))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  
  TEST_REAL_EQUAL(cons_map.getFilenames().size(),0)
RESULT

CHECK((const std::vector< typename ConsensusElementType::ElementContainerType >& getMapVector() const))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  
  TEST_REAL_EQUAL(cons_map.getMapVector().size(),0)
RESULT

CHECK((std::vector< String >& getFilenames()))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  cons_map.getFilenames().resize(1);
  String name="blub";
  cons_map.getFilenames()[0] = name;
  
  TEST_REAL_EQUAL(cons_map.getFilenames().size(),1)
  TEST_EQUAL(cons_map.getFilenames()[0] == "blub", true)
RESULT

CHECK((std::vector< typename ConsensusElementType::ElementContainerType >& getMapVector()))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  cons_map.getMapVector().resize(1);
  FeatureMap feat_map;
  Feature feat;
  feat.getPosition()[0] = 1;
  feat.getPosition()[1] = 4;
  feat.getIntensity() = 23;
  feat_map.push_back(feat);
  cons_map.getMapVector()[0] = feat_map;
  
  TEST_REAL_EQUAL(cons_map.getMapVector().size(),1)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getPosition()[0],1)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getPosition()[1],4)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getIntensity(),23)
RESULT

CHECK((void setFilenames(const std::vector < String >& filenames)))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFilenames(name_vector);
  
  TEST_REAL_EQUAL(cons_map.getFilenames().size(),1)
  TEST_EQUAL(cons_map.getFilenames()[0] == "blub", true)
RESULT

CHECK((void setMapVector(const std::vector < typename ConsensusElementType::ElementContainerType >& map_vector)))
  ConsensusMap<ConsensusFeature < FeatureMap > > cons_map;
  vector<FeatureMap> map_vector(1);
  FeatureMap feat_map;
  Feature feat;
  feat.getPosition()[0] = 1;
  feat.getPosition()[1] = 4;
  feat.getIntensity() = 23;
  feat_map.push_back(feat);
  map_vector[0] = feat_map;
  cons_map.setMapVector(map_vector);
  
  TEST_REAL_EQUAL(cons_map.getMapVector().size(),1)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getPosition()[0],1)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getPosition()[1],4)
  TEST_REAL_EQUAL((cons_map.getMapVector()[0])[0].getIntensity(),23)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



