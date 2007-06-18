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
#include <OpenMS/FORMAT/FeatureXMLFile.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/StarAlignment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef ConsensusFeature<FeatureMap<> > ConsensusFeatureType;
typedef DPosition <2> PositionType;

START_TEST(StarAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

StarAlignment<ConsensusFeatureType>* ptr = 0;
CHECK((StarAlignment()))
	ptr = new StarAlignment<ConsensusFeatureType>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~StarAlignment()))
	delete ptr;
RESULT

CHECK((virtual String getAlignmentTree() const))
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  vector<FeatureMap<>*> map_vector;
  FeatureMap<>map1;
  FeatureMap<>map2;
  FeatureMap<>map3;
  FeatureMap<>map4;
  FeatureMap<>map5;
  FeatureMap<>map6;
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

CHECK((UInt getReferenceMapIndex() const))
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  
  TEST_REAL_EQUAL(alignment.getReferenceMapIndex(),0)
RESULT

CHECK((virtual void run() throw (Exception::InvalidValue)))
  FeatureXMLFile feature_file;
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  vector<String> name_vector(2);
  String file_name = "data/MapAlignmentFeatureMap1.xml";
  FeatureMap<>map1;
  feature_file.load(file_name,map1);
  name_vector[0] = file_name;
  file_name = "data/MapAlignmentFeatureMap2.xml";
  FeatureMap<>map2;
  feature_file.load(file_name,map2);
  name_vector[1] = file_name;
  vector<FeatureMap<>*> map_vector(2);
  map_vector[0] = &map1;
  map_vector[1] = &map2;
  alignment.setElementMapVector(map_vector);
  alignment.setFileNames(name_vector);
  
  Param param;
  param.setValue("map_type","feature_map");
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  param.setValue("matching_algorithm:superimposer:type","poseclustering_affine");
  param.setValue("matching_algorithm:pairfinder:type","DelaunayPairFinder");
  alignment.setParameters(param);
  alignment.run();
  
  PRECISION(0.01)
  ConsensusFeature<FeatureMap<> > cons_feature = alignment.getFinalConsensusMap()[0];
  TEST_REAL_EQUAL(cons_feature.getPosition()[0],1273.27)  
  TEST_REAL_EQUAL(cons_feature.getPosition()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensity(),3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[0],1273.27)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().min()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getPositionRange().max()[1],904.47)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().min()[0],3.12539e+07)
  TEST_REAL_EQUAL(cons_feature.getIntensityRange().max()[0],3.12539e+07)
  ConsensusFeature<FeatureMap<> >::Group::const_iterator it = cons_feature.begin();
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

CHECK((void setReferenceMapIndex(UInt index) throw(Exception::InvalidValue)))
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  vector<FeatureMap<>*> map_vector;
  FeatureMap<>map1;
  FeatureMap<>map2;
  map_vector.push_back(&map1);
  map_vector.push_back(&map2);
  alignment.setElementMapVector(map_vector);
  alignment.setReferenceMapIndex(2);
  
  TEST_REAL_EQUAL(alignment.getReferenceMapIndex(),2)
RESULT

CHECK((void merge()))
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  ConsensusMap<ConsensusFeature < FeatureMap<> > > cons_map;
  vector<FeatureMap<>* > map_vector(3);
  Feature feat_1;
  feat_1.setRT(1);
  feat_1.setMZ(4);
  feat_1.setIntensity(23);
  Feature feat_2;
  feat_2.setRT(1.5);
  feat_2.setMZ(4);
  feat_2.setIntensity(23);
  Feature feat_3;
  feat_3.setRT(1.2);
  feat_3.setMZ(4);
  feat_3.setIntensity(23);
  FeatureMap<> feat_map_1, feat_map_2,feat_map_3;
  feat_map_1.push_back(feat_1);
  feat_map_2.push_back(feat_2);
  feat_map_3.push_back(feat_3);
  map_vector[0] = &feat_map_1;
  map_vector[1] = &feat_map_2;
  map_vector[2] = &feat_map_3;
  cons_map.setMapVector(map_vector);
  
  ConsensusFeature<FeatureMap<> > cons_1(0,0,feat_1,1,0,feat_2);
  ConsensusFeature<FeatureMap<> > cons_2(2,0,feat_3);
  cons_map.push_back(cons_1);
  cons_map.push_back(cons_2);
  
  alignment.setFinalConsensusMap(cons_map);
  
   
  TEST_REAL_EQUAL(cons_map.size(),2)
  cons_map.merge();
  TEST_REAL_EQUAL(cons_map.size(),1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



