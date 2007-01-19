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
#include <OpenMS/FORMAT/DFeatureMapFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/StarAlignment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef ConsensusFeature<FeatureMap> ConsensusFeatureType;
typedef DPosition < 2, KernelTraits > PositionType;

START_TEST(StarAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StarAlignment<ConsensusFeatureType>* ptr = 0;
CHECK(StarAlignment())
	ptr = new StarAlignment<ConsensusFeatureType>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~StarAlignment())
	delete ptr;
RESULT

CHECK(StarAlignment& operator = (StarAlignment source))
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  Param param;
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  alignment.setParam(param);
   
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  alignment.setElementMapVector(map_vector);
  
  String name="blub";
  vector<String> name_vector(1,name);
  alignment.setFileNames(name_vector);
  alignment.setMapType("feature_map");
  alignment.setReferenceMapIndex(0);
  
  StarAlignment< ConsensusFeature<FeatureMap> > alignment_copy;
  alignment_copy = alignment;

  TEST_EQUAL(alignment.getParam() == alignment_copy.getParam(),true)
  TEST_EQUAL(alignment_copy.getElementMapVector().size() == 1, true)
  TEST_EQUAL(alignment_copy.getFileNames().size() == 1, true)
  TEST_EQUAL((alignment_copy.getFileNames())[0] == "blub", true)
  TEST_EQUAL(alignment_copy.getMapType() == "feature_map", true)
RESULT

CHECK(StarAlignment(const StarAlignment& source))
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  Param param;
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  alignment.setParam(param);
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  alignment.setElementMapVector(map_vector);
  
  String name="blub";
  vector<String> name_vector(1,name);
  alignment.setFileNames(name_vector);
  alignment.setMapType("feature_map");
  alignment.setReferenceMapIndex(0);
    
  StarAlignment< ConsensusFeature<FeatureMap> > alignment_copy(alignment);

  TEST_EQUAL(alignment.getParam() == alignment_copy.getParam(),true)
  TEST_EQUAL(alignment_copy.getElementMapVector().size() == 1, true)
  TEST_EQUAL(alignment_copy.getFileNames().size() == 1, true)
  TEST_EQUAL((alignment_copy.getFileNames())[0] == "blub", true)
  TEST_EQUAL(alignment_copy.getMapType() == "feature_map", true)
  TEST_REAL_EQUAL(alignment_copy.getReferenceMapIndex(),0)
RESULT

CHECK(String getAlignmentTree() const)
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  vector<FeatureMap*> map_vector;
  FeatureMap map1;
  FeatureMap map2;
  FeatureMap map3;
  FeatureMap map4;
  FeatureMap map5;
  FeatureMap map6;
  map_vector.push_back(&map1);
  map_vector.push_back(&map2);
  map_vector.push_back(&map3);
  map_vector.push_back(&map4);
  map_vector.push_back(&map5);
  map_vector.push_back(&map6);
  alignment.setElementMapVector(map_vector);
  alignment.setReferenceMapIndex(3);
  
  TEST_EQUAL(alignment.getAlignmentTree() == "((3:0,0:1):0,(3:0,1:2):0,(3:0,2:3):0,(3:0,4:5):0,(3:0,5:6):0)", true)
RESULT

CHECK(const UnsignedInt& getReferenceMapIndex() const)
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  
  TEST_REAL_EQUAL(alignment.getReferenceMapIndex(),0)
RESULT

CHECK(void run() throw(Exception::InvalidValue))
  DFeatureMapFile feature_file;
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  vector<String> name_vector(2);
  String file_name = "data/MapAlignmentFeatureMap1.xml";
  FeatureMap map1;
  feature_file.load(file_name,map1);
  name_vector[0] = file_name;
  file_name = "data/MapAlignmentFeatureMap2.xml";
  FeatureMap map2;
  feature_file.load(file_name,map2);
  name_vector[1] = file_name;
  vector<FeatureMap*> map_vector(2);
  map_vector[0] = &map1;
  map_vector[1] = &map2;
  alignment.setElementMapVector(map_vector);
  alignment.setFileNames(name_vector);
  
  Param param;
  param.setValue("map_type","feature_map");
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  param.setValue("matching_algorithm:superimposer:type","poseclustering_affine");
  param.setValue("matching_algorithm:pairfinder:type","delaunay");
  alignment.setParam(param);
  alignment.run();
  
  PRECISION(0.01)
  ConsensusFeature<FeatureMap> cons_feature = alignment.getFinalConsensusMap()[0];
  TEST_REAL_EQUAL(cons_feature.getPosition()[0],1273.27)  
  TEST_REAL_EQUAL(cons_feature.getPosition()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  ConsensusFeature<FeatureMap>::Group::const_iterator it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getElement().getPosition()[0],1273.27)  
  TEST_REAL_EQUAL(it->getElement().getPosition()[1],904.47)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),3.12539e+07)
    
  cons_feature = alignment.getFinalConsensusMap()[5];
  TEST_REAL_EQUAL(cons_feature.getPosition()[0],1194.82)  
  TEST_REAL_EQUAL(cons_feature.getPosition()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1194.82)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],777.101)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],1.78215e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],1.78215e+07)
  it = cons_feature.begin();
  TEST_REAL_EQUAL(it->getElement().getPosition()[0],1194.82)  
  TEST_REAL_EQUAL(it->getElement().getPosition()[1],777.101)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),1.78215e+07)
  ++it;
  TEST_REAL_EQUAL(it->getElement().getPosition()[0],2401.64)  
  TEST_REAL_EQUAL(it->getElement().getPosition()[1],777.201)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),1.78215e+07)
RESULT

CHECK(void setReferenceMapIndex(UnsignedInt index) throw(Exception::InvalidValue))
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  vector<FeatureMap*> map_vector;
  FeatureMap map1;
  FeatureMap map2;
  map_vector.push_back(&map1);
  map_vector.push_back(&map2);
  alignment.setElementMapVector(map_vector);
  alignment.setReferenceMapIndex(2);
  
  TEST_REAL_EQUAL(alignment.getReferenceMapIndex(),2)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



