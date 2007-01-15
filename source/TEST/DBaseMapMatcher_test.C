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

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapMatcher.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestMapMatcher : public DBaseMapMatcher<Feature>
{
  public:
  TestMapMatcher() : DBaseMapMatcher<Feature>(){}
  TestMapMatcher(const TestMapMatcher& tm) : DBaseMapMatcher<Feature>(tm){}
  TestMapMatcher& operator=(const TestMapMatcher& tm)
  {
     DBaseMapMatcher<Feature>::operator=(tm);
     return *this;
  }
    
   virtual void estimateTransform() {}
};


START_TEST(DBaseMapMatcher, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestMapMatcher* ptr = 0;
CHECK(DBaseMapMatcher())
	ptr = new TestMapMatcher();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~DBaseMapMatcher())
	delete ptr;
RESULT

CHECK(DBaseMapMatcher& operator = (const DBaseMapMatcher& source))
  TestMapMatcher tmm;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy;
  tmm_copy = tmm;
  
  TEST_EQUAL(tmm_copy.getGrid() == grid,true)
  TEST_REAL_EQUAL(tmm_copy.getFeaturePairs().size(),0)
  TEST_REAL_EQUAL(tmm_copy.getMinQuality(),0.2)
RESULT

CHECK(DBaseMapMatcher(const DBaseMapMatcher& source))
  TestMapMatcher tmm;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy(tmm);
  
  TEST_EQUAL(tmm_copy.getGrid() == grid,true)
  TEST_REAL_EQUAL(tmm_copy.getFeaturePairs().size(),0)
  TEST_REAL_EQUAL(tmm_copy.getMinQuality(),0.2)
RESULT

CHECK(FeaturePairVector& getFeaturePairs())
  TestMapMatcher tmm;
  DFeaturePairVector< 2, Feature > pairs;  
  Feature feat1;
  Feature feat2;
  DFeaturePair<2,Feature> pair(feat1,feat2);
  pairs.push_back(pair);
  tmm.getFeaturePairs() = pairs;
  
  TEST_EQUAL(tmm.getFeaturePairs() == pairs,true)
RESULT

CHECK(Grid& getGrid())
  TestMapMatcher tmm;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  tmm.getGrid() = grid;
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK(QualityType& getMinQuality())
  TestMapMatcher tmm;
  tmm.getMinQuality() = 0.2;
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),0.2)
RESULT

CHECK(bool operator == (const DBaseMapMatcher& rhs))
  TestMapMatcher tmm;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy(tmm);
  
  TEST_EQUAL(tmm_copy == tmm,true)
RESULT

CHECK(const FeaturePairVector& getFeaturePairs() const)
  TestMapMatcher tmm;
  
  TEST_REAL_EQUAL(tmm.getFeaturePairs().size(),0)
RESULT

CHECK(const Grid& getGrid() const)
  TestMapMatcher tmm;
  DGrid<2> grid;
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK(const QualityType& getMinQuality() const)
  TestMapMatcher tmm;
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),-1)
RESULT

CHECK(void estimateTransform())

RESULT

CHECK(void setFeaturePairs(const FeaturePairVector& plist))
  TestMapMatcher tmm;
  DFeaturePairVector< 2, Feature > pairs;  
  Feature feat1;
  Feature feat2;
  DFeaturePair<2,Feature> pair(feat1,feat2);
  pairs.push_back(pair);
  tmm.setFeaturePairs(pairs);
  
  TEST_EQUAL(tmm.getFeaturePairs() == pairs,true)
RESULT

CHECK(void setGrid(const Grid& g))
  TestMapMatcher tmm;
  DGrid<2> grid;
  grid.push_back(DGridCell<2>(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK(void setMinQuality(const QualityType& qu))
  TestMapMatcher tmm;
  tmm.setMinQuality(0.1);
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),0.1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



