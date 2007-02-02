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
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DLinearMapping< 1, KernelTraits > TransformationType;
typedef DFeature<2, KernelTraits> ElementType;
typedef DFeatureMap<2, ElementType> ElementMapType;

class TestSuperimposer : public BaseSuperimposer<ElementMapType>
{
  public:
	TestSuperimposer() 
		: BaseSuperimposer<ElementMapType>()
	{ 
		check_defaults_ = false; 
	}
	TestSuperimposer(const TestSuperimposer& bpf) 
		: BaseSuperimposer<ElementMapType>(bpf)
	{
	}
	
	TestSuperimposer& operator=(const TestSuperimposer& bpf)
	{
		 BaseSuperimposer<ElementMapType>::operator=(bpf);
		 return *this;
	}
	virtual void run() 
	{
	}
};

START_TEST(BaseSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSuperimposer* ptr = 0;
CHECK(BaseSuperimposer())
	ptr = new TestSuperimposer();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~BaseSuperimposer())
	delete ptr;
RESULT

CHECK(BaseSuperimposer& operator = (BaseSuperimposer source))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  TransformationType trafo(1.,2.);
 
  TestSuperimposer bsi;
  bsi.setParameters(param);
  bsi.setElementMap(0,first);
  bsi.setElementMap(1,second);
  bsi.setTransformation(0,trafo);
  
  TestSuperimposer bsi_copy;
  bsi_copy = bsi;
  
  TEST_EQUAL(bsi.getParameters() == bsi_copy.getParameters(),true)
  TEST_EQUAL(&(bsi.getElementMap(0)) == &(bsi_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bsi.getElementMap(1)) == &(bsi_copy.getElementMap(1)),true)
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(BaseSuperimposer(const BaseSuperimposer& source))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  TransformationType trafo(1.,2.);
 
  TestSuperimposer bsi;
  bsi.setParameters(param);
  bsi.setElementMap(0,first);
  bsi.setElementMap(1,second);
  bsi.setTransformation(0,trafo);
  
  TestSuperimposer bsi_copy(bsi);
  
  TEST_EQUAL(bsi.getParameters() == bsi_copy.getParameters(),true)
  TEST_EQUAL(&(bsi.getElementMap(0)) == &(bsi_copy.getElementMap(0)),true)
  TEST_EQUAL(&(bsi.getElementMap(1)) == &(bsi_copy.getElementMap(1)),true)
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(const Param& getParameters() const)
  Param param;
  param.setValue("bla",3);
  TestSuperimposer bsi;
  bsi.setParameters(param);
  const TestSuperimposer bsi_copy(bsi);
  TEST_EQUAL(bsi_copy.getParameters() == param,true)
RESULT

CHECK(const PointMapType& getElementMap(Size index))
  
RESULT

CHECK(const PointMapType& getElementMap(Size index) const)
  ElementMapType map;
  TestSuperimposer bpf;
  bpf.setElementMap(0,map);
  const TestSuperimposer bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementMap(0)) == &map,true)
RESULT

CHECK(const TransformationType& getTransformation(Size dim) const)
  TransformationType trafo(1.,2.);
  TestSuperimposer bsi;
  bsi.setTransformation(0,trafo);
  const TestSuperimposer bsi_copy(bsi);
  
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(void run())
  
RESULT

CHECK((void setElementMap(Size const index, const PointMapType& Element_map)))
  ElementMapType map;
  TestSuperimposer bsi;
  bsi.setElementMap(0,map);
  TEST_EQUAL(&(bsi.getElementMap(0)) == &map,true)
RESULT

CHECK(void setParameters(const Param& param))
  Param param;
  TestSuperimposer bsi;
  param.setValue("bla",3);
  bsi.setParameters(param);
  const TestSuperimposer bsi_copy(bsi);
  TEST_EQUAL(bsi_copy.getParameters() == param,true)
RESULT

CHECK((void setTransformation(Size dim, const TransformationType& trafo)))
  TransformationType trafo(1.,2.);
  TestSuperimposer bsi;
  bsi.setTransformation(0,trafo);
    
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



