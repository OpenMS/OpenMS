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
CHECK(ConsensusXMLFile())
	ptr = new ConsensusXMLFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ConsensusXMLFile())
	delete ptr;
RESULT

CHECK((template<typename AlignmentT> void store(const String& filename, const AlignmentT& alignment) const throw(Exception::UnableToCreateFile)))
  std::string tmp_filename;
  ConsensusMap< ConsensusFeature<FeatureMap> > cons_map;
  ConsensusXMLFile cons_file;
 
  cons_file.load("data/ConsensusXMLFile.xml",cons_map);
  DLinearMapping<1> trafo_rt(0.5,-5.99959);
  DLinearMapping<1> trafo_mz(0.999999,-0.0990517);
  DBaseMapping<1>* bm_rt = &trafo_rt;
  DBaseMapping<1>* bm_mz = &trafo_mz;
  std::vector<DBaseMapping<1>*> mapping(2);
  mapping[0] = bm_rt;
  mapping[1] = bm_mz;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  grid[0].setMappings(mapping);
  std::vector< DGrid<2> > grid_vector(2);
  grid_vector[1] = grid; 
  
  StarAlignment< ConsensusFeature<FeatureMap> > alignment;
  Param param;
  param.setValue("matching_algorithm:type","poseclustering_pairwise");
  alignment.setParam(param);
  alignment.setTransformationVector(grid_vector);
  alignment.setFinalConsensusMap(cons_map);
  alignment.setFileNames(cons_map.getFilenames());
  alignment.setMapType("feature_map");
  alignment.setReferenceMapIndex(0);
  vector<FeatureMap> map_vector = cons_map.getMapVector();
  unsigned int n = map_vector.size();
  vector<FeatureMap*> map_pointer_vector(n);
  for (unsigned int i = 0; i < n; ++i)
  {
    map_pointer_vector[i] = &(map_vector[i]); 
  } 
  alignment.setElementMapVector(map_pointer_vector);
  
  NEW_TMP_FILE(tmp_filename);
  cons_file.store(tmp_filename,alignment);
  PRECISION(0.01)
  TEST_FILE(tmp_filename.c_str(),"data/ConsensusXMLFile.xml");
RESULT

CHECK((template<typename ElementT> void load(const String& filename, ConsensusMap<ElementT>& map) throw(Exception::FileNotFound, Exception::ParseError)))
  ConsensusMap< ConsensusFeature<FeatureMap> > cons_map;
  ConsensusXMLFile cons_file;
  cons_file.load("data/ConsensusXMLFile.xml", cons_map);
  TEST_EQUAL(cons_map.getFilenames()[0] == "data/MapAlignmentFeatureMap1.xml", true)
  TEST_EQUAL(cons_map.getFilenames()[1] == "data/MapAlignmentFeatureMap2.xml", true)

  ConsensusFeature<FeatureMap> cons_feature = cons_map[0];
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



