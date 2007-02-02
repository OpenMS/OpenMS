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

typedef DFeature<2, KernelTraits> ElementType;
typedef DFeatureMap<2, ElementType> ElementMapType;
typedef DFeaturePair < 2, ElementType > ElementPairType;
typedef DFeaturePairVector < 2, ElementType > ElementPairVectorType;
typedef DGrid<2> GridType;
typedef DPosition < 2, KernelTraits > PositionType;

class TestPairwiseMapMatcher : public BasePairwiseMapMatcher<ElementMapType>
{
  public:
  TestPairwiseMapMatcher() 
  	: BasePairwiseMapMatcher<ElementMapType>()
	{ 
		defaultsToParam_();
	}
  TestPairwiseMapMatcher(const TestPairwiseMapMatcher& bpf) 
  	: BasePairwiseMapMatcher<ElementMapType>(bpf)
  {
  	updateMembers_();
  }
  TestPairwiseMapMatcher& operator=(const TestPairwiseMapMatcher& bpf)
  {
     BasePairwiseMapMatcher<ElementMapType>::operator=(bpf);
     
     updateMembers_();
     
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
CHECK((BasePairwiseMapMatcher()))
	ptr = new TestPairwiseMapMatcher();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~BasePairwiseMapMatcher()))
	delete ptr;
RESULT

CHECK((BasePairwiseMapMatcher& operator = (const BasePairwiseMapMatcher& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setParameters(param);
  bpmm.setElementMap(0,first);
  bpmm.setElementMap(1,second);
  
  TestPairwiseMapMatcher bpmm_copy;
  bpmm_copy = bpmm;
  
  TEST_EQUAL(bpmm.getParameters() == bpmm_copy.getParameters(),true)
  TEST_EQUAL(&(bpmm.getElementMap(0)) == &(bpmm_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bpmm.getElementMap(1)) == &(bpmm_copy.getElementMap(1)),true)
RESULT

CHECK((BasePairwiseMapMatcher(const BasePairwiseMapMatcher& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setParameters(param);
  bpmm.setElementMap(0,first);
  bpmm.setElementMap(1,second);
  
  TestPairwiseMapMatcher bpmm_copy(bpmm);
  
  TEST_EQUAL(bpmm.getParameters() == bpmm_copy.getParameters(),true)
  TEST_EQUAL(&(bpmm.getElementMap(0)) == &(bpmm_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bpmm.getElementMap(1)) == &(bpmm_copy.getElementMap(1)),true)
RESULT

CHECK((const ElementPairVectorType& getElementPairs() const))
  TestPairwiseMapMatcher bpmm;
  
  TEST_EQUAL(bpmm.getElementPairs().size() == 0,true)
RESULT


CHECK((const GridType& getGrid() const))
  GridType grid;
  TestPairwiseMapMatcher bpmm;
    
  TEST_EQUAL(bpmm.getGrid() == grid,true) 
RESULT

CHECK((const Param& getParameters() const))
  Param param;
  param.setValue("number_buckets:MZ",1);
  param.setValue("number_buckets:RT",1);
  
  TestPairwiseMapMatcher bpmm;
  TEST_EQUAL(bpmm.getParameters() == param, true)
RESULT

CHECK((void setParameters(const Param& param)))
  Param param;
  param.setValue("number_buckets:MZ",3);
  param.setValue("number_buckets:RT",1);
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setParameters(param);
  const TestPairwiseMapMatcher bpmm_copy(bpmm);
 
  TEST_EQUAL(bpmm_copy.getParameters() == param, true)
RESULT

CHECK((const PointMapType& getElementMap(Size index) const))
  ElementMapType first;
  ElementMapType second;
  
  TestPairwiseMapMatcher bpmm;
  bpmm.setElementMap(0,first);
  bpmm.setElementMap(1,second);
  
  TEST_EQUAL(&(bpmm.getElementMap(0)) == &first,true)
  TEST_EQUAL(&(bpmm.getElementMap(1)) == &second,true)
RESULT

CHECK((void initGridTransformation(const PointMapType& scene_map)))
  ElementMapType scene;
  ElementType feat1;
  ElementType feat2;
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

CHECK((UnsignedInt getNumberBuckets(Size index) const))
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

CHECK((void run()))
 
RESULT

CHECK((static void registerChildren()))
  
RESULT

CHECK((void setElementMap(Size const index, const PointMapType& element_map)))
  ElementMapType first;
  TestPairwiseMapMatcher bpmm;
  bpmm.setElementMap(0,first);
   
  TEST_EQUAL(&(bpmm.getElementMap(0)) == &first,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



