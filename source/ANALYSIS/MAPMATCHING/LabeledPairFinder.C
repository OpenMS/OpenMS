// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>

using namespace std;

namespace OpenMS
{
	const DoubleReal LabeledPairFinder::sqrt2_half_ = 0.5*sqrt(2);


	LabeledPairFinder::LabeledPairFinder()
		: BaseGroupFinder()
	{
		setName("LabeledPairFinder");
		
		defaults_.setValue("rt_pair_dist", 0.3, "optimal pair distance in RT [sec]");
		defaults_.setValue("rt_dev_low", 0.44, "maximum allowed deviation below optimal retention time distance");
		defaults_.setMinFloat("rt_dev_low",0.0);
		defaults_.setValue("rt_dev_high", 1.3, "maximum allowed deviation above optimal retention time distance");
		defaults_.setMinFloat("rt_dev_high",0.0);
		
		defaults_.setValue("mz_pair_dist", 4.0, "optimal pair distance in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)");
		defaults_.setValue("mz_dev", 0.05, "maximum allowed deviation from optimal m/z distance\n");
		defaults_.setMinFloat("mz_dev",0.0);

		defaultsToParam_();
	}

	void LabeledPairFinder::run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map) 
	{
		if (input_maps.size()!=1) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"exactly one input map required");
		if (result_map.getFileDescriptions().size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"two file descriptions required");
		if (result_map.getFileDescriptions().begin()->second.filename!=result_map.getFileDescriptions().rbegin()->second.filename) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"the two file descriptions have to contain the same file name");
		checkIds_(input_maps);
		
		//look up the light and heavy index
		UInt light_index = std::numeric_limits<UInt>::max();
		UInt heavy_index = std::numeric_limits<UInt>::max();	
		for (ConsensusMap::FileDescriptions::const_iterator it = result_map.getFileDescriptions().begin();
			 	 it!=result_map.getFileDescriptions().end();
			 	 ++it)
		{
			if (it->second.label=="heavy")
			{
				heavy_index = it->first;
			}
			else if (it->second.label=="light")
			{
				light_index = it->first;
			}
		}
		if (light_index == std::numeric_limits<UInt>::max() || heavy_index == std::numeric_limits<UInt>::max())
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"the input maps have to be labeled 'light' and 'heavy'");
		}
		
		result_map.clear();
		
		// sort consensus features by RT (and MZ) to speed up searching afterwards
		typedef ConstRefVector<ConsensusMap> RefMap;
		RefMap model_ref(input_maps[0].begin(),input_maps[0].end());
		model_ref.sortByPosition();
		
		//calculate matches
		ConsensusMap matches;
		//settings
		DoubleReal rt_pair_dist = param_.getValue("rt_pair_dist");
		DoubleReal rt_dev_low = param_.getValue("rt_dev_low");
		DoubleReal rt_dev_high = param_.getValue("rt_dev_high");
		DoubleReal mz_dev = param_.getValue("mz_dev");
		DoubleReal mz_pair_dist = param_.getValue("mz_pair_dist");
		// check each feature
		for (RefMap::const_iterator it=model_ref.begin(); it!=model_ref.end(); ++it)
		{
			RefMap::const_iterator range = lower_bound(model_ref.begin(),model_ref.end(),it->getRT()+rt_pair_dist - rt_dev_low, ConsensusFeature::NthPositionLess<0>());
			while (range!=model_ref.end() && range->getRT() <= it->getRT()+rt_pair_dist + rt_dev_high)
			{
				if (range->getCharge() == it->getCharge()
					&& range->getMZ() >=it->getMZ()+mz_pair_dist/it->getCharge()-mz_dev  
					&& range->getMZ() <=it->getMZ()+mz_pair_dist/it->getCharge()+mz_dev)
				{
					
					DoubleReal score = sqrt(
															  PValue_(range->getMZ() - it->getMZ(), mz_pair_dist/it->getCharge(), mz_dev, mz_dev)
															* PValue_(range->getRT() - it->getRT(), rt_pair_dist, rt_dev_low, rt_dev_high)
														);
					matches.push_back(ConsensusFeature(light_index,it->begin()->getElementIndex(),*it));
					matches.back().insert(heavy_index,range->begin()->getElementIndex(),*range);
					matches.back().setQuality(score);
					matches.back().computeConsensus();
				}
				++range;
			}
		}
		
		//compute best pairs
		// - sort matches by quality
		// - take highest-quality matches first (greedy) and mark them as used
		set<UInt> used_features;
		matches.sortByQuality(true);
		for (ConsensusMap::const_iterator match=matches.begin(); match!=matches.end(); ++match)
		{
			//check if features are not used yet
			if ( used_features.find(match->begin()->getElementIndex())==used_features.end() && 
					 used_features.find(match->rbegin()->getElementIndex())==used_features.end()
				 )
			{
				//if unused, add it to the final set of elements
				result_map.push_back(*match);
				used_features.insert(match->begin()->getElementIndex());
				used_features.insert(match->rbegin()->getElementIndex());
			}
		}
		
		// Very useful for checking the results, and the ids have no real meaning anyway
		result_map.sortByNthPosition(Peak2D::MZ);
	}
}
