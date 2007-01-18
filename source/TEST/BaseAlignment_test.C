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
#include <OpenMS/KERNEL/ConsensusFeature.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseAlignment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef ConsensusFeature<FeatureMap> ConsensusFeatureType;
typedef DPosition < 2, KernelTraits > PositionType;

class TestAlignment : public BaseAlignment<ConsensusFeatureType>
{
  public:
    TestAlignment() : BaseAlignment<ConsensusFeatureType>()
    {}
    TestAlignment(const TestAlignment& bpf) : BaseAlignment<ConsensusFeatureType>(bpf)
    {}
    TestAlignment& operator=(const TestAlignment& bpf)
    {
      BaseAlignment<ConsensusFeatureType>::operator=(bpf);
      return *this;
    }
    virtual void run() throw (Exception::InvalidValue)
    {}
    virtual String getAlignmentTree() const
    {
      return String();
    }
};

START_TEST(BaseAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestAlignment* ptr = 0;
CHECK((BaseAlignment()))
  ptr = new TestAlignment();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~BaseAlignment()))
  delete ptr;
RESULT

CHECK((BaseAlignment& operator = (const BaseAlignment& source)))
  TestAlignment al;
  Param param;
  param.setValue("consensus_algorithm","delaunay");
  al.setParam(param);
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  al.setElementMapVector(map_vector);
  String name="blub";
  vector<String> name_vector(1,name);
  al.setFileNames(name_vector);
  al.setMapType("feature_map");
  DLinearMapping<1> trafo_rt(0.5,-5.99959);
  DLinearMapping<1> trafo_mz(0.999999,-0.0990517);
  DBaseMapping<1>* bm_rt = &trafo_rt;
  DBaseMapping<1>* bm_mz = &trafo_mz;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  std::vector<DBaseMapping<1>*> mapping(2);
  mapping[0] = bm_rt;
  mapping[1] = bm_mz;
  grid[0].setMappings(mapping);
  std::vector< DGrid<2> > grid_vector(2);
  grid_vector[1] = grid; 
  al.setTransformationVector(grid_vector);

  TestAlignment al_copy;
  al_copy = al;

  TEST_EQUAL(al.getTransformationVector() == al_copy.getTransformationVector(),true)
  TEST_EQUAL(al.getParam() == al_copy.getParam(),true)
  TEST_EQUAL(al_copy.getElementMapVector().size() == 1, true)
  TEST_EQUAL(al_copy.getFileNames().size() == 1, true)
  TEST_EQUAL((al_copy.getFileNames())[0] == "blub", true)
  TEST_EQUAL(al_copy.getMapType() == "feature_map", true)
RESULT

CHECK((BaseAlignment(const BaseAlignment& source)))
  TestAlignment al;
  Param param;
  param.setValue("consensus_algorithm","delaunay");
  al.setParam(param);
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  al.setElementMapVector(map_vector);
  String name="blub";
  vector<String> name_vector(1,name);
  al.setFileNames(name_vector);
  al.setMapType("feature_map");
  DLinearMapping<1> trafo_rt(0.5,-5.99959);
  DLinearMapping<1> trafo_mz(0.999999,-0.0990517);
  DBaseMapping<1>* bm_rt = &trafo_rt;
  DBaseMapping<1>* bm_mz = &trafo_mz;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  std::vector<DBaseMapping<1>*> mapping(2);
  mapping[0] = bm_rt;
  mapping[1] = bm_mz;
  grid[0].setMappings(mapping);
  std::vector< DGrid<2> > grid_vector(2);
  grid_vector[1] = grid; 
  al.setTransformationVector(grid_vector);

  TestAlignment al_copy(al);

  TEST_EQUAL(al.getTransformationVector() == al_copy.getTransformationVector(),true)
  TEST_EQUAL(al.getParam() == al_copy.getParam(),true)
  TEST_EQUAL(al_copy.getElementMapVector().size() == 1, true)
  TEST_EQUAL(al_copy.getFileNames().size() == 1, true)
  TEST_EQUAL((al_copy.getFileNames())[0] == "blub", true)
  TEST_EQUAL(al_copy.getMapType() == "feature_map", true)
RESULT

CHECK((String getAlignmentTree() const))

RESULT

CHECK((const Param& getParam() const))
  TestAlignment al;
  Param param;
  param.setValue("consensus_algorithm","delaunay");
  al.setParam(param);
  
  TEST_EQUAL(al.getParam() == param, true)
RESULT

CHECK((const String& getMapType() const))
  TestAlignment al;
  al.setMapType("feature_map");

  TEST_EQUAL(al.getMapType() == "feature_map", true)
RESULT

CHECK((const std::vector< ConsensusElementType >& getFinalConsensusMap() const))
  TestAlignment al;
  
  TEST_EQUAL(al.getFinalConsensusMap().size() == 0,true)
RESULT

CHECK((const std::vector< ElementContainerType* >& getElementMapVector() const))
  TestAlignment al;
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  al.setElementMapVector(map_vector);

  TEST_EQUAL(al.getElementMapVector().size() == 1, true)
  TEST_EQUAL((al.getElementMapVector()[0]) == &map, true)
RESULT

CHECK((const std::vector< String >& getFileNames() const))
  TestAlignment al;
  String name="blub";
  vector<String> name_vector(1,name);
  al.setFileNames(name_vector);

  TEST_EQUAL(al.getFileNames().size() == 1, true)
  TEST_EQUAL((al.getFileNames())[0] == "blub", true)
RESULT

CHECK((const std::vector< GridType >& getTransformationVector() const))
  TestAlignment alignment;
  
  TEST_EQUAL(alignment.getTransformationVector().size() == 0, true)
RESULT

CHECK((void run() throw(Exception::InvalidValue)))

RESULT

CHECK((void setElementMapVector(const std::vector< ElementContainerType* >& element_map_vector)))
  TestAlignment al;
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  al.setElementMapVector(map_vector);

  TEST_EQUAL(al.getElementMapVector().size() == 1, true)
  TEST_EQUAL(al.getElementMapVector()[0] == &map, true)
RESULT

CHECK((std::vector< ElementContainerType* >& getElementMapVector()))
  TestAlignment al;
  vector<FeatureMap*> map_vector;
  FeatureMap map;
  map_vector.push_back(&map);
  al.getElementMapVector() = map_vector;

  TEST_EQUAL(al.getElementMapVector().size() == 1, true)
  TEST_EQUAL(al.getElementMapVector()[0] == &map, true)
RESULT

CHECK((void setFileNames(const std::vector< String >& file_names)))
  TestAlignment al;
  String name="blub";
  vector<String> name_vector(1,name);
  al.setFileNames(name_vector);

  TEST_EQUAL(al.getFileNames().size() == 1, true)
  TEST_EQUAL((al.getFileNames())[0] == "blub", true)
RESULT

CHECK((void setMapType(const String& map_type)))
  TestAlignment al;
  al.setMapType("peak_map");

  TEST_EQUAL(al.getMapType() == "peak_map", true)
RESULT

CHECK((void setParam(const Param& param)))
  TestAlignment al;
  Param param;
  param.setValue("consensus_algorithm","delaunay");
  al.setParam(param);
  
  TEST_EQUAL(al.getParam() == param,true)
RESULT

CHECK((void setFinalConsensusMap(const std::vector < ConsensusElementType >& final_consensus_map)))
  TestAlignment al;
  std::vector<ConsensusFeatureType> cons_map(4);
  al.setFinalConsensusMap(cons_map);
  
  TEST_EQUAL(al.getFinalConsensusMap().size() == 4,true)
RESULT

CHECK((void setTransformationVector(const std::vector< GridType >& transformations)))
  TestAlignment alignment;
  std::vector< DGrid<2> > grid_vector(2);
  alignment.setTransformationVector(grid_vector);

  TEST_EQUAL(alignment.getTransformationVector().size() == 2, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



