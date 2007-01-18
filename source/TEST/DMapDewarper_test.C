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

#include<OpenMS/KERNEL/DFeatureMap.h>
#include<OpenMS/KERNEL/DFeature.h>
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
	
	DFeature<2> feat1, feat2, feat3;
	
	feat1.getPosition()[MZ] = 1.0;
	feat1.getPosition()[RT] = 2.0;
	feat2.getPosition()[MZ] = 2.0;
	feat2.getPosition()[RT] = 5.0;
	feat3.getPosition()[MZ] = 2.0;
	feat3.getPosition()[RT] = 4.0;
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	DFeatureMap<2> feat_map;
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

	DFeature<2> feat1, feat2, feat3;
	
	feat1.getPosition()[MZ] = 1.0;
	feat1.getPosition()[RT] = 2.0;
	feat2.getPosition()[MZ] = 2.0;
	feat2.getPosition()[RT] = 5.0;
	feat3.getPosition()[MZ] = 2.0;
	feat3.getPosition()[RT] = 4.0;
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	DFeatureMap<2> feat_map;
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
	
	DFeature<2> feat1, feat2, feat3;
	
	feat1.getPosition()[MZ] = 1.0;
	feat1.getPosition()[RT] = 2.0;
	feat2.getPosition()[MZ] = 2.0;
	feat2.getPosition()[RT] = 5.0;
	feat3.getPosition()[MZ] = 2.0;
	feat3.getPosition()[RT] = 4.0;
	
	DGridCell<2> cell1(0,0,20,20);
	Grid the_grid;
	the_grid.push_back(cell1);
		
	DFeatureMap<2> feat_map;
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
	DFeature<2> feat1, feat2;
	
	feat1.getPosition()[MZ] = 1.0;
	feat1.getPosition()[RT] = 2.0;
	feat2.getPosition()[MZ] =  2.0;
	feat2.getPosition()[RT] = 5.0;
	
	pair1.setFirst(feat1);
	pair1.setSecond(feat2);
	pair1.setQuality(5);
	
	// second pair
	DFeaturePair<2> pair2;
	DFeature<2> feat3, feat4;
	
	feat3.getPosition()[MZ] = 2.0;
	feat3.getPosition()[RT] = 4.0;
	feat4.getPosition()[MZ] = 4.0;
	feat4.getPosition()[RT] = 9.0;
	
	pair2.setFirst(feat3);
	pair2.setSecond(feat4);
	pair2.setQuality(5);
	
	// third pair
	DFeaturePair<2> pair3;
	DFeature<2> feat5, feat6;
	
	feat5.getPosition()[MZ] = 3.0;
	feat5.getPosition()[RT] = 6.0;
	feat6.getPosition()[MZ] = 6.0;
	feat6.getPosition()[RT] = 13.0;
	
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
	DFeatureMap<2> feat_map;
	feat_map.push_back(feat1);
	feat_map.push_back(feat3);
	feat_map.push_back(feat5);
	
	Grid grid2 = mmatcher.getGrid();
	DMapDewarper<> dewarper;
	dewarper.setMap(feat_map);
	dewarper.setGrid(grid2);
	dewarper.dewarp();
	
	DFeatureMap<2> dewarped = dewarper.getMap();
	DFeatureMap<2>::const_iterator map_iter = dewarped.begin();
	TEST_REAL_EQUAL(map_iter->getPosition()[0],5.0);
	TEST_REAL_EQUAL(map_iter->getPosition()[1],2.0);
	
	map_iter++;
	TEST_REAL_EQUAL(map_iter->getPosition()[0],9.0);
	TEST_REAL_EQUAL(map_iter->getPosition()[1],4.0);
	
	map_iter++;
	TEST_REAL_EQUAL(map_iter->getPosition()[0],13.0);
	TEST_REAL_EQUAL(map_iter->getPosition()[1],6.0);	
		
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
	DFeatureMap<2> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);	
		
	TEST_EQUAL( map==dewarper.getMap(), true )
	
RESULT

CHECK((void setMap(MapType& elem)))
	DFeatureMap<2> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);	
		
	TEST_EQUAL( map==dewarper.getMap(), true )
	
RESULT

CHECK((const MapType& getMap() const))
	DFeatureMap<2> map;
	DMapDewarper<> dewarper;
	dewarper.setMap(map);		
	const  DFeatureMap<2> map2 = dewarper.getMap();
	
	TEST_EQUAL( map==map2, true )
	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
