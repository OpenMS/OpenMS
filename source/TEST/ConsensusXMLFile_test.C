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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusXMLFile* ptr = 0;
CHECK((ConsensusXMLFile()))
	ptr = new ConsensusXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~ConsensusXMLFile()))
	delete ptr;
RESULT

CHECK((template<typename AlignmentT> void store(const String& filename, const AlignmentT& alignment) const throw(Exception::UnableToCreateFile)))
  std::string tmp_filename;
  ConsensusMap< ConsensusFeature<FeatureMap<> > > cons_map;
  ConsensusXMLFile cons_file;
  FeatureMap<> feat_map_1;
  FeatureMap<> feat_map_2;
  std::vector< FeatureMap<>* > feature_maps(2);
  feature_maps[0] = &feat_map_1;
  feature_maps[1] = &feat_map_2;
  cons_map.setMapVector(feature_maps);
  
  cons_file.load("data/ConsensusXMLFile.xml",cons_map);
  LinearMapping trafo_rt(0.5,-5.99959);
  LinearMapping trafo_mz(0.999999,-0.0990517);
  BaseMapping* bm_rt = &trafo_rt;
  BaseMapping* bm_mz = &trafo_mz;
  std::vector<BaseMapping*> mapping(2);
  mapping[0] = bm_rt;
  mapping[1] = bm_mz;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  grid[0].setMappings(mapping);
  std::vector< Grid > grid_vector(2);
  grid_vector[1] = grid; 
  
  StarAlignment< ConsensusFeature<FeatureMap<> > > alignment;
  Param param;
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  alignment.setParameters(param);
  alignment.setTransformationVector(grid_vector);
  alignment.setFinalConsensusMap(cons_map);
  alignment.setFileNames(cons_map.getFilenames());
  alignment.setMapType("feature_map");
  alignment.setReferenceMapIndex(0);
  alignment.setElementMapVector(cons_map.getMapVector());
    
  NEW_TMP_FILE(tmp_filename);
  cons_file.store(tmp_filename,alignment);
  PRECISION(0.01)
  TEST_FILE(tmp_filename.c_str(),"data/ConsensusXMLFile.xml");
  TEST_EQUAL(cons_file.isValid(tmp_filename),true);
RESULT

CHECK((template <typename ElementT> void load(const String &filename, ConsensusMap< ElementT > &map, bool load_element_maps=true) throw (Exception::FileNotFound, Exception::ParseError)))
  ConsensusMap< ConsensusFeature<FeatureMap<> > > cons_map;
  ConsensusXMLFile cons_file;
  FeatureMap<> feat_map_1;
  FeatureMap<> feat_map_2;
  std::vector< FeatureMap<>* > feature_maps(2);
  feature_maps[0] = &feat_map_1;
  feature_maps[1] = &feat_map_2;
  cons_map.setMapVector(feature_maps);
  cons_file.load("data/ConsensusXMLFile.xml", cons_map);
  TEST_EQUAL(cons_map.getFilenames()[0] == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(cons_map.getFilenames()[1] == "data/MapAlignmentFeatureMap2.xml", true)

  ConsensusFeature<FeatureMap<> > cons_feature = cons_map[0];
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
    
  cons_feature = cons_map[5];
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

CHECK(static bool isValid(const String& filename))
	//tested above
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



