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

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>

namespace OpenMS
{

	FeatureGroupingAlgorithmUnlabeled::FeatureGroupingAlgorithmUnlabeled()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmUnlabeled");
		
		defaults_.insert("",DelaunayPairFinder< std::vector< ConsensusFeature< FeatureMap<> > > >().getParameters());
		
		defaultsToParam_();
	}

	FeatureGroupingAlgorithmUnlabeled::~FeatureGroupingAlgorithmUnlabeled()
	{
	}

	void FeatureGroupingAlgorithmUnlabeled::group(const std::vector< FeatureMap<> >& maps, ConsensusMap<>& out)
	{

		// define reference map (the one with most peaks)
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
		out.clear();
		const FeatureMap<>& ref_map = maps[reference_map_index];
    for (UInt i=0; i < ref_map.size(); ++i)
    {
      out.push_back(ConsensusFeature< FeatureMap<> > (reference_map_index,i,ref_map[i]));
    }
  
		// loop over all other maps, extend the groups
		for (UInt i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				//build a consensus map of map i
				std::vector< ConsensusFeature< FeatureMap<> > > map;
				const FeatureMap<>& map2 = maps[i];
				map.reserve(map2.size());
				for (UInt i2=0; i2 < map2.size(); ++i2)
				{
					map.push_back(ConsensusFeature< FeatureMap<> >(map2[i2].getPosition(),map2[i2].getIntensity()));
		    }
				
				// compute the consensus of the reference map and map i
				DelaunayPairFinder< std::vector< ConsensusFeature< FeatureMap<> > > > pair_finder;
				pair_finder.setParameters(param_.copy("",true));
				// std::cout << "pair_finder.getParameters()\n" << pair_finder.getParameters() << std::endl;
				pair_finder.computeConsensusMap(map,out);
			}
		}
	}

} // namespace OpenMS
