// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

typedef FeatureMap< Feature> ElementMapType;
typedef ElementPair < Feature > ElementPairType;
typedef vector < Feature > ElementPairVectorType;

class TestPairwiseMapMatcher : public BasePairwiseMapMatcher<ElementMapType>
{
  public:
  TestPairwiseMapMatcher() 
  	: BasePairwiseMapMatcher<ElementMapType>()
	{ 
		defaultsToParam_();
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

CHECK((const ElementPairVectorType& getElementPairs() const))
  TestPairwiseMapMatcher bpmm;
  
  TEST_EQUAL(bpmm.getElementPairs().size() == 0,true)
RESULT


CHECK((const Grid& getGrid() const))
  Grid grid;
  TestPairwiseMapMatcher bpmm;
    
  TEST_EQUAL(bpmm.getGrid() == grid,true) 
RESULT

CHECK((const PointMapType& getElementMap(UInt index) const))
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
  Feature feat1;
  Feature feat2;
  Feature::PositionType pos1(0,0);
  Feature::PositionType pos2(2,3);
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

CHECK((UInt getNumberBuckets(UInt index) const))
  TestPairwiseMapMatcher bpmm;
    
  TEST_EQUAL(bpmm.getNumberBuckets(0) == 1,true)
  TEST_EQUAL(bpmm.getNumberBuckets(1) == 1,true)
RESULT

CHECK((void setNumberBuckets(UInt dim, UInt number)))
  TestPairwiseMapMatcher bpmm;
  bpmm.setNumberBuckets(0,3);
  bpmm.setNumberBuckets(1,4);
    
  TEST_EQUAL(bpmm.getNumberBuckets(0) == 3,true)
  TEST_EQUAL(bpmm.getNumberBuckets(1) == 4,true)
RESULT 

CHECK((void run()))
	NOT_TESTABLE
RESULT

CHECK((void setElementMap(UInt const index, const PointMapType& element_map)))
  ElementMapType first;
  TestPairwiseMapMatcher bpmm;
  bpmm.setElementMap(0,first);
   
  TEST_EQUAL(&(bpmm.getElementMap(0)) == &first,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



