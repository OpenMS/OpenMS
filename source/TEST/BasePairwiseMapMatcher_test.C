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

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairwiseMapMatcher.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DFeature<2, KernelTraits> FeatureType;
typedef DFeatureMap<2, FeatureType> FeatureMapType;
typedef DFeaturePair < 2, FeatureType > FeaturePairType;
typedef DFeaturePairVector < 2, FeatureType > FeaturePairVectorType;
typedef DGrid<2> GridType;
typedef DPosition < 2, KernelTraits > PositionType;

class TestPairwiseMapMatcher : public BasePairwiseMapMatcher<FeatureMapType>
{
  public:
  TestPairwiseMapMatcher() : BasePairwiseMapMatcher<FeatureMapType>(){}
  TestPairwiseMapMatcher(const TestPairwiseMapMatcher& bpf) : BasePairwiseMapMatcher<FeatureMapType>(bpf){}
  TestPairwiseMapMatcher& operator=(const TestPairwiseMapMatcher& bpf)
  {
     BasePairwiseMapMatcher<FeatureMapType>::operator=(bpf);
     return *this;
  }
  virtual void run() 
  {
  }
};

START_TEST(BasePairwiseMapMatcher, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairwiseMapMatcher* ptr = 0;
CHECK(BasePairwiseMapMatcher())
	ptr = new TestPairwiseMapMatcher();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~BasePairwiseMapMatcher())
	delete ptr;
RESULT

CHECK(BasePairwiseMapMatcher& operator = (const BasePairwiseMapMatcher& source))
  Param param;
  param.setValue("bla",3);
  FeatureMapType first;
  FeatureMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setParam(param);
  bpmm.setFeatureMap(0,first);
  bpmm.setFeatureMap(1,second);
  
  TestPairwiseMapMatcher bpmm_copy;
  bpmm_copy = bpmm;
  
  TEST_EQUAL(bpmm.getParam() == bpmm_copy.getParam(),true)
  TEST_EQUAL(&(bpmm.getFeatureMap(0)) == &(bpmm_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bpmm.getFeatureMap(1)) == &(bpmm_copy.getFeatureMap(1)),true)
RESULT

CHECK(BasePairwiseMapMatcher(const BasePairwiseMapMatcher& source))
  Param param;
  param.setValue("bla",3);
  FeatureMapType first;
  FeatureMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setParam(param);
  bpmm.setFeatureMap(0,first);
  bpmm.setFeatureMap(1,second);
  
  TestPairwiseMapMatcher bpmm_copy(bpmm);
  
  TEST_EQUAL(bpmm.getParam() == bpmm_copy.getParam(),true)
  TEST_EQUAL(&(bpmm.getFeatureMap(0)) == &(bpmm_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bpmm.getFeatureMap(1)) == &(bpmm_copy.getFeatureMap(1)),true)
RESULT

CHECK(const FeaturePairVectorType& getFeaturePairs() const)
  TestPairwiseMapMatcher bpmm;
  
  TEST_EQUAL(bpmm.getFeaturePairs().size() == 0,true)
RESULT


CHECK(const GridType& getGrid() const)
  GridType grid;
  TestPairwiseMapMatcher bpmm;
    
  TEST_EQUAL(bpmm.getGrid() == grid,true) 
RESULT

CHECK(const Param& getParam() const)
  Param param;
  param.setValue("bla",3);
  TestPairwiseMapMatcher bpmm;
  bpmm.setParam(param);
  const TestPairwiseMapMatcher bpmm_copy(bpmm);
 
  TEST_EQUAL(bpmm_copy.getParam() == param,true)
RESULT

CHECK(const PointMapType& getFeatureMap(Size index) const)
  FeatureMapType first;
  FeatureMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setFeatureMap(0,first);
  bpmm.setFeatureMap(1,second);
  
  TEST_EQUAL(&(bpmm.getFeatureMap(0)) == &first,true)
  TEST_EQUAL(&(bpmm.getFeatureMap(1)) == &second,true)
RESULT

CHECK(void initGridTransformation(const PointMapType& scene_map))
  FeatureMapType scene;
  FeatureType feat1;
  FeatureType feat2;
  PositionType pos1(0,0);
  PositionType pos2(2,3);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat1.setIntensity(300);
  scene.push_back(feat1);
  scene.push_back(feat2);
  
  TestPairwiseMapMatcher bpmm;
  bpmm.initGridTransformation(scene);
  
  TEST_EQUAL(bpmm.getGrid().size() == 1,true)   
RESULT

CHECK(UnsignedInt getNumberBuckets(Size index) const)
  TestPairwiseMapMatcher bpmm;
    
  TEST_EQUAL(bpmm.getNumberBuckets(0) == 1,true)
  TEST_EQUAL(bpmm.getNumberBuckets(1) == 1,true)
RESULT

CHECK((void setNumberBuckets(Size const index, UnsignedInt number)))
  TestPairwiseMapMatcher bpmm;
  bpmm.setNumberBuckets(0,3);
  bpmm.setNumberBuckets(1,4);
    
  TEST_EQUAL(bpmm.getNumberBuckets(0) == 3,true)
  TEST_EQUAL(bpmm.getNumberBuckets(1) == 4,true)
RESULT 

CHECK(void run())
  // ???
RESULT

CHECK(static void registerChildren())
  // ???
RESULT

CHECK((void setFeatureMap(Size const index, const PointMapType& feature_map)))
  FeatureMapType first;
  TestPairwiseMapMatcher bpmm;
  bpmm.setFeatureMap(0,first);
   
  TEST_EQUAL(&(bpmm.getFeatureMap(0)) == &first,true)
RESULT

CHECK(void setParam(const Param& param))
  Param param;
  TestPairwiseMapMatcher bpmm;
  param.setValue("bla",3);
  bpmm.setParam(param);
  const TestPairwiseMapMatcher bpmm_copy(bpmm);
 
  TEST_EQUAL(bpmm_copy.getParam() == param,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



