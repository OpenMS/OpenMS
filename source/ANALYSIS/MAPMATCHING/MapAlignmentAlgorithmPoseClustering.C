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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>


using namespace std;

namespace OpenMS
{

	MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering()
		: MapAlignmentAlgorithm()
	{
		setName("MapAlignmentAlgorithmPoseClustering");
	
		defaults_.insert("superimposer:",PoseClusteringAffineSuperimposer<ConsensusMap>().getParameters());
		defaults_.insert("pairfinder:",DelaunayPairFinder().getParameters());
		
		defaultsToParam_();
	}

	MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering()
	{
	}

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(vector< MSExperiment<> >& maps)
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
		
    // build a consensus map of the elements of the reference map (take the 400 highest peaks)
    ConsensusMap cons_ref_map;
    { //new scope to get rid of the tmp variables
			DPeakArray<RawDataPoint2D> tmp;
			tmp.sortByIntensity(true);
			if (tmp.size()>400) tmp.resize(400); //TODO make this a parameter
			maps[reference_map_index].get2DData(tmp);
	    for (UInt i=0; i < tmp.size(); ++i)
	    {
	    	Feature tmp_feature;
	    	tmp_feature.RawDataPoint2D::operator=(tmp[i]);
	      cons_ref_map.push_back(ConsensusFeature(reference_map_index,i,tmp_feature));
	    }
	  }
		
		//init superimposer and pairfinder with model and parameters
		PoseClusteringAffineSuperimposer<ConsensusMap> superimposer;
    superimposer.setModelMap(cons_ref_map);
    superimposer.setParameters(param_.copy("superimposer:",true));
    
    DelaunayPairFinder pairfinder;
		pairfinder.setModelMap(0, cons_ref_map);
    pairfinder.setParameters(param_.copy("pairfinder:",true));
		
		for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				ConsensusMap map;
		    DPeakArray<RawDataPoint2D> tmp;
		    maps[i].get2DData(tmp);
		    for (UInt i2=0; i2 < tmp.size(); ++i2)
		    {
		      map.push_back(ConsensusFeature(tmp[i2]));
		    }
				
				//run superimposer    
	      superimposer.setSceneMap(map);
	      LinearMapping si_trafo;
	      superimposer.run(si_trafo);
	      //run pairfinder
	      pairfinder.setTransformation(0, si_trafo);        
	      pairfinder.setSceneMap(1,map);
	      vector< ElementPair < ConsensusFeature > > element_pairs;
	      pairfinder.setElementPairs(element_pairs);
	      pairfinder.findElementPairs();

				// calculate the transformation
				LinearMapping trafo;
				if (element_pairs.size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(element_pairs);
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo = si_trafo;
				}
				
				// apply transformation
				for (UInt j=0; j< maps[i].size(); ++j)
				{
					DoubleReal rt = maps[i][j].getRT();
					trafo.apply(rt);
					maps[i][j].setRT(rt);
				}
			}
		}
	}


	void MapAlignmentAlgorithmPoseClustering::alignFeatureMaps(vector< FeatureMap<> >& maps)
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
		
    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
    FeatureMapType cons_ref_map;
    const FeatureMap<>& ref_map = maps[reference_map_index];
    for (UInt i=0; i < ref_map.size(); ++i)
    {
      cons_ref_map.push_back(ConsensusFeature (reference_map_index,i,ref_map[i]));
    }
   
		//init superimposer and pairfinder with model and parameters
		PoseClusteringAffineSuperimposer<FeatureMapType> superimposer;
    superimposer.setModelMap(cons_ref_map);
    superimposer.setParameters(param_.copy("superimposer:",true));
    
    DelaunayPairFinder pairfinder;
		pairfinder.setModelMap(0, cons_ref_map);
    pairfinder.setParameters(param_.copy("pairfinder:",true));

    for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				FeatureMapType map;
		    const FeatureMap<>& map2 = maps[i];
		    map.reserve(map2.size());
		    for (UInt i2=0; i2 < map2.size(); ++i2)
		    {
		      map.push_back(ConsensusFeature(map2[i2]));
		    }

				//run superimposer    
	      superimposer.setSceneMap(map);
	      LinearMapping si_trafo;
	      superimposer.run(si_trafo);
	      //run pairfinder
	      pairfinder.setTransformation(0, si_trafo);        
	      pairfinder.setSceneMap(1, map);
	      vector< ElementPair < ConsensusFeature > > element_pairs;
	      pairfinder.setElementPairs(element_pairs);
	      pairfinder.findElementPairs();

				// calculate the transformation
				LinearMapping trafo;
				if (element_pairs.size() > 2) // estimate for each grid cell a better transformation using the element pairs
				{
					trafo = calculateRegression_(element_pairs);
				}
				else // otherwise take the estimated transformation of the superimposer
				{
					trafo =  si_trafo;
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
