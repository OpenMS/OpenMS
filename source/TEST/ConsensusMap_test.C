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
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusMap< ConsensusFeature< FeatureMap<> > >* ptr = 0;
CHECK((ConsensusMap()))
	ptr = new ConsensusMap<ConsensusFeature < FeatureMap<> > >();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusMap()))
	delete ptr;
RESULT

CHECK((ConsensusMap& operator = (const ConsensusMap& source)))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFileNames(name_vector);
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map_copy;
  cons_map_copy = cons_map;
  
  TEST_EQUAL(cons_map_copy.getFileNames()[0] == "blub", true)
RESULT

CHECK((ConsensusMap(const ConsensusMap& source)))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFileNames(name_vector);
  vector< FeatureMap<>* > map_vector(4);
  
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map_copy(cons_map);
  
  TEST_EQUAL(cons_map_copy.getFileNames()[0] == "blub", true)
RESULT

CHECK((ConsensusMap(typename Base::size_type n)))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map(5);
  
  TEST_REAL_EQUAL(cons_map.size(),5)
RESULT

CHECK((const std::vector< String >& getFileNames() const))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  
  TEST_REAL_EQUAL(cons_map.getFileNames().size(),0)
RESULT

CHECK((std::vector< String >& getFileNames()))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  cons_map.getFileNames().resize(1);
  String name="blub";
  cons_map.getFileNames()[0] = name;
  
  TEST_REAL_EQUAL(cons_map.getFileNames().size(),1)
  TEST_EQUAL(cons_map.getFileNames()[0] == "blub", true)
RESULT

CHECK((void setFileNames(const std::vector < String >& filenames)))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  String name="blub";
  vector<String> name_vector(1,name);
  cons_map.setFileNames(name_vector);
  
  TEST_REAL_EQUAL(cons_map.getFileNames().size(),1)
  TEST_EQUAL(cons_map.getFileNames()[0] == "blub", true)
RESULT

CHECK((void merge(ConsensusMap<ConsensusElementT>& new_map)))
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map_2;
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
  
  ConsensusFeature<FeatureMap<> > cons_1(0,0,feat_1,1,0,feat_2);
  ConsensusFeature<FeatureMap<> > cons_2(2,0,feat_3,3,0,feat_4);
  cons_map.push_back(cons_1);
  cons_map.push_back(cons_2);
   
  TEST_REAL_EQUAL(cons_map.size(),2)
  cons_map.merge(cons_map_2);
  TEST_REAL_EQUAL(cons_map_2.size(),1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



