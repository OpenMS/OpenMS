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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseQuality.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ProductModel.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>


///////////////////////////

START_TEST(BaseQuality, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;

class TestQuality : public BaseQuality
{
  public:
	TestQuality()
		: BaseQuality()
	{
		setName(TestQuality::getProductName());
		check_defaults_ = false;
	}

	double evaluate(const IndexSet& /*set*/, const BaseModel<2>& /*model*/)
	{
		return 5.0;
	}
  
  double evaluate(const IndexSet& /*set*/, const BaseModel<1>& /*model*/, UnsignedInt /*dim*/)
	{
		return 3.0;
	}
	
	static const String getProductName()
	{ 
		return "TestQuality"; 
	}
};

// default ctor
TestQuality* ptr = 0;
CHECK(TestQuality())
	ptr = new TestQuality();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~TestQuality())
	delete ptr;
RESULT

// assignment operator
CHECK(TestQuality& operator = (const TestQuality& source))
	TestQuality tm1;
  TestQuality tm2;
  tm2 = tm1;

  TestQuality tm3;

  tm1 = TestQuality();
	TEST_EQUAL(tm3,tm2)
RESULT

// copy constructor
CHECK(TestQuality(const TestQuality& source))
	TestQuality fp1;	
 
  TestQuality fp2(fp1);

  TestQuality fp3;
 
  fp1 = TestQuality();
	TEST_EQUAL(fp2, fp3)
RESULT

CHECK(double evaluate(const IndexSet& /*set*/, const BaseModel<2>& /*model*/))
	TestQuality ft;
	FeaFiModule::IndexSet  inds;
	inds.insert(std::make_pair(7,7));
	
	// compile model
	GaussModel * gm1 = new GaussModel();
	GaussModel * gm2 = new GaussModel();
		
	BasicStatistics<>  stat;
	stat.setMean(2.5);
	stat.setVariance(3.0);
	
	gm1->setScalingFactor(5.0);
	gm1->setInterpolationStep(0.3);
	gm1->setParam(stat,1,5);
	
	gm2->setScalingFactor(5.0);
	gm2->setInterpolationStep(0.3);
	gm2->setParam(stat,1,5);
	
	ProductModel<2> pm1;
	pm1.setModel(0,gm1);
	pm1.setModel(1,gm2);
	
	TEST_REAL_EQUAL(ft.evaluate(inds,pm1),5.0);
RESULT

CHECK(double evaluate(const IndexSet& /*set*/, const BaseModel<1>& /*model*/, UnsignedInt dim))
	TestQuality ft;
	FeaFiModule::IndexSet  inds;
	inds.insert(std::make_pair(7,7));
	GaussModel bm;
	TEST_REAL_EQUAL(ft.evaluate(inds,bm,1),3.0);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
