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
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapMatcher.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestMapMatcher : public BaseMapMatcher<Feature>
{
  public:
  TestMapMatcher() : BaseMapMatcher<Feature>(){}
  TestMapMatcher(const TestMapMatcher& tm) : BaseMapMatcher<Feature>(tm){}
  TestMapMatcher& operator=(const TestMapMatcher& tm)
  {
     BaseMapMatcher<Feature>::operator=(tm);
     return *this;
  }
    
   virtual void estimateTransform() {}
};


START_TEST(BaseMapMatcher, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestMapMatcher* ptr = 0;
CHECK((BaseMapMatcher()))
	ptr = new TestMapMatcher();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~BaseMapMatcher()))
	delete ptr;
RESULT

CHECK((BaseMapMatcher& operator = (const BaseMapMatcher& source)))
  TestMapMatcher tmm;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy;
  tmm_copy = tmm;
  
  TEST_EQUAL(tmm_copy.getGrid() == grid,true)
  TEST_REAL_EQUAL(tmm_copy.getElementPairs().size(),0)
  TEST_REAL_EQUAL(tmm_copy.getMinQuality(),0.2)
RESULT

CHECK((BaseMapMatcher(const BaseMapMatcher& source)))
  TestMapMatcher tmm;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy(tmm);
  
  TEST_EQUAL(tmm_copy.getGrid() == grid,true)
  TEST_REAL_EQUAL(tmm_copy.getElementPairs().size(),0)
  TEST_REAL_EQUAL(tmm_copy.getMinQuality(),0.2)
RESULT

CHECK((ElementPairVector& getElementPairs()))
  TestMapMatcher tmm;
  vector < ElementPair< Feature > > pairs;  
  Feature feat1;
  Feature feat2;
  ElementPair< Feature > pair(feat1,feat2);
  pairs.push_back(pair);
  tmm.getElementPairs() = pairs;
  
  TEST_EQUAL(tmm.getElementPairs() == pairs,true)
RESULT

CHECK((Grid& getGrid()))
  TestMapMatcher tmm;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  tmm.getGrid() = grid;
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK((QualityType& getMinQuality()))
  TestMapMatcher tmm;
  tmm.getMinQuality() = 0.2;
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),0.2)
RESULT

CHECK((bool operator == (const BaseMapMatcher& rhs)))
  TestMapMatcher tmm;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  tmm.setMinQuality(0.2);
  
  TestMapMatcher tmm_copy(tmm);
  
  TEST_EQUAL(tmm_copy == tmm,true)
RESULT

CHECK((const ElementPairVector& getElementPairs() const))
  TestMapMatcher tmm;
  
  TEST_REAL_EQUAL(tmm.getElementPairs().size(),0)
RESULT

CHECK((const Grid& getGrid() const))
  TestMapMatcher tmm;
  Grid grid;
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK((QualityType getMinQuality() const))
  TestMapMatcher tmm;
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),-1)
RESULT

CHECK((virtual void estimateTransform()=0))

RESULT

CHECK((void setElementPairs(const ElementPairVector &plist)))
  TestMapMatcher tmm;
  vector < ElementPair< Feature > > pairs;    
  Feature feat1;
  Feature feat2;
  ElementPair< Feature > pair(feat1,feat2);
  pairs.push_back(pair);
  tmm.setElementPairs(pairs);
  
  TEST_EQUAL(tmm.getElementPairs() == pairs,true)
RESULT

CHECK((void setGrid(const Grid& g)))
  TestMapMatcher tmm;
  Grid grid;
  grid.push_back(GridCell(1816,603.449,3108.3,1002.35));
  tmm.setGrid(grid);
  
  TEST_EQUAL(tmm.getGrid() == grid,true)
RESULT

CHECK((void setMinQuality(QualityType qu)))
  TestMapMatcher tmm;
  tmm.setMinQuality(0.1);
  
  TEST_REAL_EQUAL(tmm.getMinQuality(),0.1)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



