// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DLinearMapping< 1, KernelTraits > TransformationType;
typedef DFeature<2, KernelTraits> FeatureType;
typedef DFeatureMap<2, FeatureType> FeatureMapType;
typedef DFeaturePair < 2, FeatureType > FeaturePairType;
typedef DFeaturePairVector < 2, FeatureType > FeaturePairVectorType;

class TestPairFinder : public BasePairFinder<FeatureMapType>
{
  public:
	TestPairFinder() : BasePairFinder<FeatureMapType>(){}
	TestPairFinder(const TestPairFinder& bpf) : BasePairFinder<FeatureMapType>(bpf){}
	TestPairFinder& operator=(const TestPairFinder& bpf)
	{
		 BasePairFinder<FeatureMapType>::operator=(bpf);
		 return *this;
	}
	virtual void run() 
	{
	}

};

START_TEST(BasePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairFinder* ptr = 0;
CHECK(BasePairFinder())
	ptr = new TestPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~BasePairFinder())
	delete ptr;
RESULT

CHECK(BasePairFinder& operator = (BasePairFinder source))
  Param param;
  param.setValue("bla",3);
  FeatureMapType first;
  FeatureMapType second;
  FeaturePairVectorType pairs;
  
  TestPairFinder bpf;
  bpf.setParam(param);
  bpf.setFeatureMap(0,first);
  bpf.setFeatureMap(1,second);
  bpf.setFeaturePairs(pairs);
  
  TestPairFinder bpf_copy;
  bpf_copy = bpf;
  
  TEST_EQUAL(bpf.getParam() == bpf_copy.getParam(),true)
  TEST_EQUAL(&(bpf.getFeatureMap(0)) == &(bpf_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bpf.getFeatureMap(1)) == &(bpf_copy.getFeatureMap(1)),true)
	TEST_EQUAL(&(bpf.getFeaturePairs()) == &(bpf_copy.getFeaturePairs()),true)
RESULT

CHECK(BasePairFinder(const BasePairFinder& source))
  Param param;
  param.setValue("bla",3);
  FeatureMapType first;
  FeatureMapType second;
  FeaturePairVectorType pairs;
  
  TestPairFinder bpf;
  bpf.setParam(param);
  bpf.setFeatureMap(0,first);
  bpf.setFeatureMap(1,second);
  bpf.setFeaturePairs(pairs);
  
  TestPairFinder bpf_copy(bpf);
  
  TEST_EQUAL(bpf.getParam() == bpf_copy.getParam(),true)
  TEST_EQUAL(&(bpf.getFeatureMap(0)) == &(bpf_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bpf.getFeatureMap(1)) == &(bpf_copy.getFeatureMap(1)),true)
	TEST_EQUAL(&(bpf.getFeaturePairs()) == &(bpf_copy.getFeaturePairs()),true)
RESULT

CHECK(Param& getParam())
  Param param;
  param.setValue("bla",3);
  TestPairFinder bpf;
  bpf.getParam()=param;
  TEST_EQUAL(bpf.getParam() == param,true)
RESULT

CHECK(const FeaturePairVectorType getFeaturePairs() const)
	FeaturePairVectorType pairs;
	TestPairFinder bpf;
	bpf.setFeaturePairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  
  TEST_EQUAL(&(bpf_copy.getFeaturePairs()) == &pairs,true)
RESULT

CHECK(const Param& getParam() const)
  Param param;
  param.setValue("bla",3);
  TestPairFinder bpf;
  bpf.setParam(param);
  const TestPairFinder bpf_copy(bpf);
 
  TEST_EQUAL(bpf_copy.getParam() == param,true)
RESULT

CHECK(const PointMapType& getFeatureMap(Size index) const)
  FeatureMapType map;
  TestPairFinder bpf;
  bpf.setFeatureMap(0,map);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getFeatureMap(0)) == &map,true)
RESULT

CHECK(static void registerChildren())
  // ???
RESULT

CHECK(void run())
  // ???
RESULT

CHECK((void setFeatureMap(Size const index, const PointMapType& feature_map)))
  FeatureMapType map;
  TestPairFinder bpf;
  bpf.setFeatureMap(0,map);
  TEST_EQUAL(&(bpf.getFeatureMap(0)) == &map,true)
RESULT

CHECK(void setFeaturePairs(FeaturePairVectorType& feature_pairs))
  FeaturePairVectorType pairs;
  TestPairFinder bpf;
  bpf.setFeaturePairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getFeaturePairs()) == &pairs,true)
RESULT

CHECK(void setParam(const Param& param))
  Param param;
  TestPairFinder bpf;
  param.setValue("bla",3);
  bpf.setParam(param);
  const TestPairFinder bpf_copy(bpf);
 
  TEST_EQUAL(bpf_copy.getParam() == param,true)
RESULT

CHECK((void setTransformation(Size dim, const TransformationType& trafo)))
  TransformationType trafo(1.,2.);
  TestPairFinder bpf;
  bpf.setTransformation(0,trafo);
    
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(const TransformationType& getTransformation(Size dim) const)
  TransformationType trafo(1.,2.);
  TestPairFinder bpf;
  bpf.setTransformation(0,trafo);
  const TestPairFinder bpf_copy(bpf);
  
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



