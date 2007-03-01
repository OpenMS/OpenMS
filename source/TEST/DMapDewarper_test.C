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


#include<OpenMS/CONCEPT/ClassTest.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

#include<OpenMS/KERNEL/FeatureMap.h>
#include<OpenMS/KERNEL/Feature.h>
#include<string>

///////////////////////////

using namespace std;
using namespace OpenMS;

typedef DGrid<2> Grid;
typedef DFeaturePairVector<2> FeaturePairVector;
typedef DBaseMapping<1> MappingType;
typedef std::vector<MappingType*> MappingVector;

enum DimensionId
		{
			RT = DimensionDescription < LCMS_Tag >::RT,
			MZ = DimensionDescription < LCMS_Tag >::MZ
		};	

///////////////////////////

/////////////////////////////////////////////////////////////

START_TEST(DMapDewarper<MapT>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


DMapDewarper<>* pl_ptr = 0;
CHECK((DMapDewarper()))
	pl_ptr = new DMapDewarper<>();
	TEST_NOT_EQUAL(pl_ptr, 0)
RESULT

CHECK((~DMapDewarper()))
	delete pl_ptr;
RESULT

CHECK((DMapDewarper(const DMapDewarper& source)))
	
	Feature feat1, feat2, feat3;
	
	feat1.setMZ(1.0);
	feat1.setRT(2.0);
	feat2.setMZ(2.0);
	feat2.setRT(5.0);
	feat3.setMZ(2.0);
	feat3.setRT(4.0);
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	FeatureMap<> feat_map;
	feat_map.push_back(feat1);
	feat_map.push_back(feat2);
	feat_map.push_back(feat3);
	
	DMapDewarper<> dewarper1;
	dewarper1.setMap(feat_map);
	dewarper1.setGrid(the_grid);
	
	DMapDewarper<> dewarper2(dewarper1);
	
	TEST_EQUAL(dewarper1.getMap()==dewarper1.getMap(),true);
	TEST_EQUAL(dewarper1.getGrid()==dewarper1.getGrid(),true);
				
RESULT

CHECK((DMapDewarper& operator = (const DMapDewarper& source)))

	Feature feat1, feat2, feat3;
	
	feat1.setMZ(1.0);
	feat1.setRT(2.0);
	feat2.setMZ(2.0);
	feat2.setRT(5.0);
	feat3.setMZ(2.0);
	feat3.setRT(4.0);
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	FeatureMap<> feat_map;
	feat_map.push_back(feat1);
	feat_map.push_back(feat2);
	feat_map.push_back(feat3);
	
	DMapDewarper<> dewarper1;
	dewarper1.setMap(feat_map);
	dewarper1.setGrid(the_grid);
	
	DMapDewarper<> dewarper2 = dewarper1;
	
	TEST_EQUAL(dewarper1.getMap()==dewarper1.getMap(),true);
	TEST_EQUAL(dewarper1.getGrid()==dewarper1.getGrid(),true);
				
RESULT


CHECK((bool operator == (const DMapDewarper& rhs)))
	
	Feature feat1, feat2, feat3;
	
	feat1.setMZ(1.0);
	feat1.setRT(2.0);
	feat2.setMZ(2.0);
	feat2.setRT(5.0);
	feat3.setMZ(2.0);
	feat3.setRT(4.0);
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	FeatureMap<> feat_map;
	feat_map.push_back(feat1);
	feat_map.push_back(feat2);
	feat_map.push_back(feat3);
	
	DMapDewarper<> dewarper1;
	dewarper1.setMap(feat_map);
	dewarper1.setGrid(the_grid);
	
	DMapDewarper<> dewarper2;
	dewarper2.setMap(feat_map);
	dewarper2.setGrid(the_grid);
	
	TEST_EQUAL(dewarper1==dewarper2,true);
	
RESULT


CHECK((void dewarp()))
	
	// first pair
	DFeaturePair<2> pair1;
	Feature feat1, feat2;
	
	feat1.setMZ(1.0);
	feat1.setRT(2.0);
	feat2.setMZ(2.0);
	feat2.setRT(5.0);
	
	pair1.setFirst(feat1);
	pair1.setSecond(feat2);
	pair1.setQuality(5);
	
	// second pair
	DFeaturePair<2> pair2;
	Feature feat3, feat4;
	
	feat3.setMZ(2.0);
	feat3.setRT(4.0);
	feat4.setMZ(4.0);
	feat4.setRT(9.0);
	
	pair2.setFirst(feat3);
	pair2.setSecond(feat4);
	pair2.setQuality(5);
	
	// third pair
	DFeaturePair<2> pair3;
	Feature feat5, feat6;
	
	feat5.setMZ(3.0);
	feat5.setRT(6.0);
	feat6.setMZ(6.0);
	feat6.setRT(13.0);
	
	pair3.setFirst(feat5);
	pair3.setSecond(feat6);
	pair3.setQuality(5);
	
	FeaturePairVector pairs;
	pairs.push_back(pair1);
	pairs.push_back(pair2);
	pairs.push_back(pair3);
	
	DGridCell<2> cell1(0,0,20,20);
		
	Grid the_grid;
	the_grid.push_back(cell1);
	
	DMapMatcherRegression<> mmatcher;
	mmatcher.setGrid(the_grid);
	mmatcher.setFeaturePairs(pairs);
	
	// estimate mappings
	mmatcher.estimateTransform();
		
	// now we apply these mappings and check the results
	FeatureMap<> feat_map;
	feat_map.push_back(feat1);
	feat_map.push_back(feat3);
	feat_map.push_back(feat5);
	
	Grid grid2 = mmatcher.getGrid();
	DMapDewarper<> dewarper;
	dewarper.setMap(feat_map);
	dewarper.setGrid(grid2);
	dewarper.dewarp();
	
	FeatureMap<> dewarped = dewarper.getMap();
	FeatureMap<>::const_iterator map_iter = dewarped.begin();
	TEST_REAL_EQUAL(map_iter->getRT(),5.0);
	TEST_REAL_EQUAL(map_iter->getMZ(),2.0);
	
	map_iter++;
	TEST_REAL_EQUAL(map_iter->getRT(),9.0);
	TEST_REAL_EQUAL(map_iter->getMZ(),4.0);
	
	map_iter++;
	TEST_REAL_EQUAL(map_iter->getRT(),13.0);
	TEST_REAL_EQUAL(map_iter->getMZ(),6.0);	
		
RESULT

CHECK((Grid& getGrid()))
	Grid agrid;
	DMapDewarper<> dewarper;
	dewarper.setGrid(agrid);
	
	TEST_EQUAL( agrid == dewarper.getGrid(),true )
	
RESULT

CHECK((const Grid& getGrid() const))
	Grid agrid;
	DMapDewarper<> dewarper;
	dewarper.setGrid(agrid);	
	const Grid agrid2 = dewarper.getGrid();
	
	TEST_EQUAL( agrid==agrid2, true )
	
RESULT

CHECK((void setGrid(Grid& g)))
	Grid agrid;
	DMapDewarper<> dewarper;
	dewarper.setGrid(agrid);	
	Grid agrid2 = dewarper.getGrid();
	
	TEST_EQUAL(agrid == agrid2,true)
RESULT

CHECK((MapType& getMap()))
	FeatureMap<> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);	
		
	TEST_EQUAL( map==dewarper.getMap(), true )
	
RESULT

CHECK((void setMap(MapType& elem)))
	FeatureMap<> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);	
		
	TEST_EQUAL( map==dewarper.getMap(), true )
	
RESULT

CHECK((const MapType& getMap() const))
	FeatureMap<> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);		
	const  FeatureMap<> map2 = dewarper.getMap();
	
	TEST_EQUAL( map==map2, true )
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
