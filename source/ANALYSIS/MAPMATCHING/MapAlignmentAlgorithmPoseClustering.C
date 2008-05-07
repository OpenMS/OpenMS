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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>

namespace OpenMS
{

	MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmPoseClustering");
		
		PoseClusteringPairwiseMapMatcher<std::vector< ConsensusFeature< FeatureMap<> > > > pairwise_matcher;
		defaults_.insert("pair_matcher:",pairwise_matcher.getParameters());
		
		defaultsToParam_();
	}

	MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
	{
	}

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(std::vector< MSExperiment<> >& maps)
	{		
		//define reference map (the one with most peaks)
		UInt reference_map_index = 0;
		UInt max_count = 0;		
		for (UInt m=0; m<maps.size(); ++m)
		{
			//init getSize() by calling updateRanges
			maps[m].updateRanges(1);
			if (maps[m].getSize()>max_count)
			{
				max_count = maps[m].getSize();
				reference_map_index = m;
			}
		}
		
    // build a consensus map of the elements of the reference map (tank the 400 highest peaks)
    std::vector< ConsensusFeature< DPeakArray<RawDataPoint2D> > > cons_ref_map;
    { //new scope to get rid of the tmp variables
			DPeakArray<RawDataPoint2D> tmp;
			tmp.sortByIntensity(true);
			if (tmp.size()>400) tmp.resize(400);
			maps[reference_map_index].get2DData(tmp);
	    for (UInt i=0; i < tmp.size(); ++i)
	    {
	      ConsensusFeature< DPeakArray<RawDataPoint2D> > c(reference_map_index,i,tmp[i]);
	      cons_ref_map.push_back(c);
	    }
	  }
		
		//init mapmatcher with the reference map
		PoseClusteringPairwiseMapMatcher< std::vector<ConsensusFeature< DPeakArray<RawDataPoint2D > > > > pairwise_matcher;
		pairwise_matcher.setParameters(param_.copy("pair_matcher:",true));
		pairwise_matcher.setElementMap(0,cons_ref_map); //define model

		for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				std::vector< ConsensusFeature< DPeakArray<RawDataPoint2D> > > map;
		    DPeakArray<RawDataPoint2D> tmp;
		    maps[i].get2DData(tmp);
		    for (UInt i2=0; i2 < tmp.size(); ++i2)
		    {
		      ConsensusFeature< DPeakArray<RawDataPoint2D> >  c(tmp[i2].getPosition(),tmp[i2].getIntensity());
		      map.push_back(c);
		    }
				
				pairwise_matcher.setElementMap(1, map); //define scene
				pairwise_matcher.initGridTransformation(map);
				pairwise_matcher.run();

				// calculate the transformation
				LinearMapping trafo;
				if (pairwise_matcher.getElementPairs().size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(pairwise_matcher.getElementPairs());
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo =  pairwise_matcher.getGrid();
				}
				
				// apply transformation
				for(UInt j=0; j< maps[i].size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					trafo.apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
	}


	void MapAlignmentAlgorithmPoseClustering::alignFeatureMaps(std::vector< FeatureMap<> >& maps)
	{
		//define reference map (the one with most peaks)
		UInt reference_map_index = 0;
		UInt max_count = 0;		
		for (UInt m=0; m<maps.size(); ++m)
		{
			if (maps[m].size()>max_count)
			{
				max_count = maps[m].size();
				reference_map_index = m;
			}
		}
		
    // build a consensus map of the elements of the reference map
    std::vector< ConsensusFeature< FeatureMap<> > > cons_ref_map;
    const FeatureMap<>& map = maps[reference_map_index];
    for (UInt i=0; i < map.size(); ++i)
    {
      ConsensusFeature< FeatureMap<> > c(reference_map_index,i,map[i]);
      cons_ref_map.push_back(c);
    }
   
    //initialize the pairwise mapmatcher
		PoseClusteringPairwiseMapMatcher<std::vector< ConsensusFeature< FeatureMap<> > > > pairwise_matcher;
		pairwise_matcher.setParameters(param_.copy("pair_matcher:",true));
    pairwise_matcher.setElementMap(0,cons_ref_map); //define scene

    for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				std::vector< ConsensusFeature< FeatureMap<> > > map;
		    const FeatureMap<>& map2 = maps[i];
		    for (UInt i2=0; i2 < map2.size(); ++i2)
		    {
		      ConsensusFeature< FeatureMap<> >  c(map2[i2].getPosition(),map2[i2].getIntensity());
		      map.push_back(c);
		    }

				// compute a transformation for each grid cell and find pairs in the reference_map_ and map_i
				pairwise_matcher.setElementMap(1, map); //define model
				pairwise_matcher.initGridTransformation(map);
				pairwise_matcher.run();

				// calculate the transformation
				LinearMapping trafo;
				if (pairwise_matcher.getElementPairs().size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(pairwise_matcher.getElementPairs());
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo =  pairwise_matcher.getGrid();
				}

				// apply transformation
				for (UInt j = 0; j < map.size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					trafo.apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
	}

} //namespace 
