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
#include <OpenMS/ANALYSIS/MAPMATCHING/MapMatcherRegression.h>
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

	void MapAlignmentAlgorithmPoseClustering::alignPeakMaps(std::vector< MSExperiment<> >&)
	{
    std::cout << "Add alignment here!" << std::endl;
    //TODO wie bei features, nur N largest peaks => parameter zum einstellen
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
		
    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
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

    MapMatcherRegression<ConsensusFeature< FeatureMap<> > > lin_regression;
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

				// use the linear regression only if there are more than 2 pairs
				LinearMapping trafo;
				if (pairwise_matcher.getElementPairs().size() > 2)  
				{
					// estimate for each grid cell a better transformation using the element pairs
					lin_regression.setElementPairs(pairwise_matcher.getElementPairs());
					lin_regression.setGrid(pairwise_matcher.getGrid());
					lin_regression.setMinQuality(-1.0);
					lin_regression.estimateTransform();
					trafo = lin_regression.getGrid().begin()->getMappings()[RawDataPoint2D::RT];
				}
				// otherwise take the estimated transformation of the superimposer
				else
				{
					trafo =  pairwise_matcher.getGrid().begin()->getMappings()[RawDataPoint2D::RT];
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
