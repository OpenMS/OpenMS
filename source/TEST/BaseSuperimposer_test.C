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
// $Maintainer: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DLinearMapping< 1, KernelTraits > TransformationType;
typedef DFeature<2, KernelTraits> FeatureType;
typedef DFeatureMap<2, FeatureType> FeatureMapType;

class TestSuperimposer : public BaseSuperimposer<FeatureMapType>
{
  public:
	TestSuperimposer() : BaseSuperimposer<FeatureMapType>(){}
	TestSuperimposer(const TestSuperimposer& bpf) : BaseSuperimposer<FeatureMapType>(bpf){}
	TestSuperimposer& operator=(const TestSuperimposer& bpf)
	{
		 BaseSuperimposer<FeatureMapType>::operator=(bpf);
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
  FeatureMapType first;
  FeatureMapType second;
  TransformationType trafo(1.,2.);
 
  TestSuperimposer bsi;
  bsi.setParam(param);
  bsi.setFeatureMap(0,first);
  bsi.setFeatureMap(1,second);
  bsi.setTransformation(0,trafo);
  
  TestSuperimposer bsi_copy;
  bsi_copy = bsi;
  
  TEST_EQUAL(bsi.getParam() == bsi_copy.getParam(),true)
  TEST_EQUAL(&(bsi.getFeatureMap(0)) == &(bsi_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bsi.getFeatureMap(1)) == &(bsi_copy.getFeatureMap(1)),true)
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(BaseSuperimposer(const BaseSuperimposer& source))
  Param param;
  param.setValue("bla",3);
  FeatureMapType first;
  FeatureMapType second;
  TransformationType trafo(1.,2.);
 
  TestSuperimposer bsi;
  bsi.setParam(param);
  bsi.setFeatureMap(0,first);
  bsi.setFeatureMap(1,second);
  bsi.setTransformation(0,trafo);
  
  TestSuperimposer bsi_copy(bsi);
  
  TEST_EQUAL(bsi.getParam() == bsi_copy.getParam(),true)
  TEST_EQUAL(&(bsi.getFeatureMap(0)) == &(bsi_copy.getFeatureMap(0)),true)
  TEST_EQUAL(&(bsi.getFeatureMap(1)) == &(bsi_copy.getFeatureMap(1)),true)
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK(Param& getParam())
  Param param;
  param.setValue("bla",3);
  TestSuperimposer bpf;
  bpf.getParam()=param;
  TEST_EQUAL(bpf.getParam() == param,true)
RESULT

CHECK(const Param& getParam() const)
  Param param;
  param.setValue("bla",3);
  TestSuperimposer bsi;
  bsi.setParam(param);
  const TestSuperimposer bsi_copy(bsi);
  TEST_EQUAL(bsi_copy.getParam() == param,true)
RESULT

CHECK(const PointMapType& getFeatureMap(Size index))
  // ???
RESULT

CHECK(const PointMapType& getFeatureMap(Size index) const)
  FeatureMapType map;
  TestSuperimposer bpf;
  bpf.setFeatureMap(0,map);
  const TestSuperimposer bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getFeatureMap(0)) == &map,true)
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
  // ???
RESULT

CHECK((void setFeatureMap(Size const index, const PointMapType& feature_map)))
  FeatureMapType map;
  TestSuperimposer bsi;
  bsi.setFeatureMap(0,map);
  TEST_EQUAL(&(bsi.getFeatureMap(0)) == &map,true)
RESULT

CHECK(void setParam(const Param& param))
  Param param;
  TestSuperimposer bsi;
  param.setValue("bla",3);
  bsi.setParam(param);
  const TestSuperimposer bsi_copy(bsi);
  TEST_EQUAL(bsi_copy.getParam() == param,true)
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



