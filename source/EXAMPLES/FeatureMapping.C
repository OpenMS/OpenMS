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

#include<OpenMS/KERNEL/DFeature.h>
#include<OpenMS/KERNEL/DFeatureMap.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include<OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>

#include<iostream.h>  

using namespace OpenMS;

/// The grid is simply a vector of cells.
typedef DGrid<2> Grid;
/// The feature pairs are computed by the feature matching class
typedef DFeaturePairVector<2> FeaturePairVector;

enum DimensionId
		{
			RT = DimensionDescription < LCMS_Tag >::RT,
			MZ = DimensionDescription < LCMS_Tag >::MZ
		};	

/** @brief An example of how the feature mapping classes can be
           obtained to generate a feature mapping application.

*/
int main()
{
	// In the real world, some class would provide us with
	// a list of feature pairs. In this simple example, we
	// generate our own pairs.
	DFeatureMap<2> feat_map;
	
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
	
	feat_map.push_back(feat1);
	feat_map.push_back(feat3);
	feat_map.push_back(feat5);
	
	// the grid consists of one cell only
	// including all features
	DGridCell<2> cell1(0,0,20,20);
		
	Grid the_grid;
	the_grid.push_back(cell1);
	
	// estimate the transform using linear regression
	DMapMatcherRegression<2> mmatcher;
	mmatcher.setGrid(the_grid);
	mmatcher.setFeaturePairs(pairs);
	
	mmatcher.estimateTransform();
	Grid grid2 = mmatcher.getGrid();
	DMapDewarper<2> dewarper;
	
	dewarper.setFeatures(feat_map);
	dewarper.setGrid(grid2);
	dewarper.dewarp();

	// show output
	std::cout << "Vor dewarping: " << std::endl;
	
	for (DFeatureMap<2>::const_iterator map_iter = feat_map.begin();
	     map_iter != feat_map.end();
			 ++map_iter) 
	{
		std::cout << map_iter->getPosition()[0] << " ";
		std::cout << map_iter->getPosition()[1] << " ";
		std::cout << map_iter->getIntensity() << std::endl;
	}
						 
	std::cout << "Nach dewarping: " << std::endl;	
	DFeatureMap<2> dewarped = dewarper.getFeatures();
	for (DFeatureMap<2>::const_iterator map_iter = dewarped.begin();
	     map_iter != dewarped.end();
			 ++map_iter)
	{
		std::cout << map_iter->getPosition()[0] << " ";
		std::cout << map_iter->getPosition()[1] << " ";
		std::cout << map_iter->getIntensity() << std::endl;
	}
	
	return 0;
}
