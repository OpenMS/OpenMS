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
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DLinearMapping< 1, KernelTraits > TransformationType;
typedef DFeature<2, KernelTraits> ElementType;
typedef DFeatureMap<2, ElementType> ElementMapType;
typedef DFeaturePair < 2, ElementType > ElementPairType;
typedef DFeaturePairVector < 2, ElementType > ElementPairVectorType;

class TestPairFinder 
	: public BasePairFinder<ElementMapType>
{
  public:
	TestPairFinder() 
		: BasePairFinder<ElementMapType>()
	{
		check_defaults_ = false; 
	}
	
	TestPairFinder(const TestPairFinder& bpf) 
	: BasePairFinder<ElementMapType>(bpf)
	{	
	}
	
	TestPairFinder& operator=(const TestPairFinder& bpf)
	{
		if (&bpf==this) return *this;
		
		BasePairFinder<ElementMapType>::operator=(bpf);
		
		return *this;
	}
	
	virtual void findElementPairs() 
	{
		
	}

};

START_TEST(BasePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairFinder* ptr = 0;
CHECK((BasePairFinder()))
	ptr = new TestPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~BasePairFinder()))
	delete ptr;
RESULT

CHECK((BasePairFinder& operator = (const BasePairFinder& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  ElementPairVectorType pairs;
  
  TestPairFinder bpf;
  bpf.setParameters(param);
  bpf.setElementMap(0,first);
  bpf.setElementMap(1,second);
  bpf.setElementPairs(pairs);
  
  TestPairFinder bpf_copy;
  bpf_copy = bpf;
  
  TEST_EQUAL(bpf.getParameters() == bpf_copy.getParameters(),true)
  TEST_EQUAL(&(bpf.getElementMap(0)) == &(bpf_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bpf.getElementMap(1)) == &(bpf_copy.getElementMap(1)),true)
	TEST_EQUAL(&(bpf.getElementPairs()) == &(bpf_copy.getElementPairs()),true)
RESULT

CHECK((BasePairFinder(const BasePairFinder& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  ElementPairVectorType pairs;
  
  TestPairFinder bpf;
  bpf.setParameters(param);
  bpf.setElementMap(0,first);
  bpf.setElementMap(1,second);
  bpf.setElementPairs(pairs);
  
  TestPairFinder bpf_copy(bpf);
  
  TEST_EQUAL(bpf.getParameters() == bpf_copy.getParameters(),true)
  TEST_EQUAL(&(bpf.getElementMap(0)) == &(bpf_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bpf.getElementMap(1)) == &(bpf_copy.getElementMap(1)),true)
	TEST_EQUAL(&(bpf.getElementPairs()) == &(bpf_copy.getElementPairs()),true)
RESULT

CHECK((const ElementPairVectorType& getElementPairs() const))
	ElementPairVectorType pairs;
	TestPairFinder bpf;
	bpf.setElementPairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  
  TEST_EQUAL(&(bpf_copy.getElementPairs()) == &pairs,true)
RESULT

CHECK((const Param& getParameters() const))
  Param param;
  param.setValue("bla",3);
  TestPairFinder bpf;
  bpf.setParameters(param);
  const TestPairFinder bpf_copy(bpf);
 
  TEST_EQUAL(bpf_copy.getParameters() == param,true)
RESULT

CHECK((const PointMapType& getElementMap(Size index) const))
  ElementMapType map;
  TestPairFinder bpf;
  bpf.setElementMap(0,map);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementMap(0)) == &map,true)
RESULT

CHECK((static void registerChildren()))
  
RESULT

CHECK((void findElementPairs()))
  
RESULT

CHECK((void setElementMap(Size const index, const PointMapType& element_map)))
  ElementMapType map;
  TestPairFinder bpf;
  bpf.setElementMap(0,map);
  TEST_EQUAL(&(bpf.getElementMap(0)) == &map,true)
RESULT

CHECK((void setElementPairs(ElementPairVectorType& element_pairs)))
  ElementPairVectorType pairs;
  TestPairFinder bpf;
  bpf.setElementPairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementPairs()) == &pairs,true)
RESULT

CHECK((void setParameters(const Param& param)))
  Param param;
  TestPairFinder bpf;
  param.setValue("bla",3);
  bpf.setParameters(param);
  const TestPairFinder bpf_copy(bpf);
 
  TEST_EQUAL(bpf_copy.getParameters() == param,true)
RESULT

CHECK((void setTransformation(Size dim, const TransformationType& trafo)))
  TransformationType trafo(1.,2.);
  TestPairFinder bpf;
  bpf.setTransformation(0,trafo);
    
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK((const TransformationType& getTransformation(Size dim) const))
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



