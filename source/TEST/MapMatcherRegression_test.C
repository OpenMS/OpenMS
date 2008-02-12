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


#include<OpenMS/CONCEPT/ClassTest.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/MapMatcherRegression.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/ElementPair.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>

#include<OpenMS/KERNEL/FeatureMap.h>
#include<OpenMS/KERNEL/Feature.h>
#include<string>

///////////////////////////

using namespace std;
using namespace OpenMS;

typedef ElementPair< Feature > ElementPairType;
typedef vector< ElementPairType > ElementPairVector;
typedef BaseMapping MappingType;
typedef std::vector<MappingType*> MappingVector;

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(MapMatcherRegression<ElementT>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapMatcherRegression<>* pl_ptr = 0;
CHECK((MapMatcherRegression()))
	pl_ptr = new MapMatcherRegression<>();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK((~MapMatcherRegression()))
	delete pl_ptr;
RESULT

CHECK((bool operator == (const MapMatcherRegression& rhs)))
	
	// first pair
	ElementPairType pair1;
	Feature feat1, feat2;
	
  feat1.setMZ(1.0);
  feat1.setRT(2.0);
  feat2.setMZ(2.0);
  feat2.setRT(5.0);
  
  pair1.setFirst(feat1);
  pair1.setSecond(feat2);
  pair1.setQuality(5);
  
  // second pair
  ElementPairType pair2;
  Feature feat3, feat4;
  
  feat3.setMZ(2.0);
  feat3.setRT(4.0);
  feat4.setMZ(4.0);
  feat4.setRT(9.0);
	
	pair2.setFirst(feat3);
	pair2.setSecond(feat4);
	pair2.setQuality(5);

	ElementPairVector pairs;
	pairs.push_back(pair1);
	pairs.push_back(pair2);
	
	GridCell  cell1(0,0,20,20);
		
	Grid the_grid;
	the_grid.push_back(cell1);
	
	MapMatcherRegression<> mmatcher;
	mmatcher.setGrid(the_grid);
	mmatcher.setElementPairs(pairs);
	
	MapMatcherRegression<> mmatcher2;
	mmatcher2.setGrid(the_grid);
	mmatcher2.setElementPairs(pairs);
	
	TEST_EQUAL(mmatcher==mmatcher2,true)
	
	
RESULT


CHECK((void estimateTransform()))

	// first pair
	ElementPairType pair1;
	Feature feat1, feat2;
	
  feat1.setMZ(1.0);
  feat1.setRT(2.0);
  feat2.setMZ(2.0);
  feat2.setRT(5.0);
  
  pair1.setFirst(feat1);
  pair1.setSecond(feat2);
  pair1.setQuality(5);
  
  // second pair
  ElementPairType pair2;
  Feature feat3, feat4;
  
  feat3.setMZ(2.0);
  feat3.setRT(4.0);
  feat4.setMZ(4.0);
  feat4.setRT(9.0);;;
	
	pair2.setFirst(feat3);
	pair2.setSecond(feat4);
	pair2.setQuality(5);
	
	// third pair
	ElementPairType pair3;
	Feature feat5, feat6;
	
  feat5.setMZ(3.0);
  feat5.setRT(6.0);
  feat6.setMZ(6.0);
  feat6.setRT(13.0);
	
	pair3.setFirst(feat5);
	pair3.setSecond(feat6);
	pair3.setQuality(5);
	
	ElementPairVector pairs;
	pairs.push_back(pair1);
	pairs.push_back(pair2);
	pairs.push_back(pair3);
	
	GridCell  cell1(0,0,20,20);
		
	Grid the_grid;
	the_grid.push_back(cell1);
	
	MapMatcherRegression<> mmatcher;
	mmatcher.setGrid(the_grid);
	mmatcher.setElementPairs(pairs);
	
	mmatcher.estimateTransform();
	
	Grid grid2 = mmatcher.getGrid();
	Grid::const_iterator cit = grid2.begin();
	MappingVector mvec1 = cit->getMappings();
	
	// we expect two mappings, one for each dimension
	TEST_EQUAL(mvec1.size(),2);
	
	// now let's see how these mappings look like
	//PRECISION(0.05)
	LinearMapping* lmap1 = dynamic_cast<LinearMapping* >(mvec1.at(0));
	TEST_REAL_EQUAL(lmap1->getSlope(),2)
	TEST_REAL_EQUAL(lmap1->getIntercept(),1)
	
	LinearMapping* lmap2 = dynamic_cast<LinearMapping* >(mvec1.at(1));
	TEST_REAL_EQUAL(lmap2->getSlope(),2)
	TEST_REAL_EQUAL(lmap2->getIntercept(),0)
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
