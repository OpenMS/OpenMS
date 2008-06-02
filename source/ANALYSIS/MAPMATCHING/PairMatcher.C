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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>

using namespace std;

namespace OpenMS
{
	const DoubleReal PairMatcher::sqrt2_ = sqrt(2);


	PairMatcher::PairMatcher()
		: BasePairFinder(), 
			pairs_()
	{
		setName("PairMatcher");
		
		defaults_.setValue("rt_pair_dist", 0.3, "optimal pair distance in RT [sec]");
		defaults_.setValue("rt_stdev_low", 0.22, "standard deviation below optimal retention time distance");
		defaults_.setMinFloat("rt_stdev_low",0.0);
		defaults_.setValue("rt_stdev_high", 0.65, "standard deviation above optimal retention time distance");
		defaults_.setMinFloat("rt_stdev_high",0.0);
		
		defaults_.setValue("mz_pair_dist", 4.0, "optimal pair distance in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)");
		defaults_.setValue("mz_stdev", 0.025, "standard deviation from optimal m/z distance\n");
		defaults_.setMinFloat("mz_stdev",0.0);

		defaultsToParam_();
	}

	void PairMatcher::run(ConsensusMap& result_map) 
	{
		if (!model_map_) throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__,"model map not set");
		
		result_map.clear();
		
		// sort consensus features by RT (and MZ) to speed up searching afterwards
		typedef DPeakConstReferenceArray<ConsensusMap> RefMap;
		RefMap model_ref(model_map_->begin(),model_map_->end());
		model_ref.sortByPosition();
		
		//**********************************************************
		//calculate matches
		ConsensusMap matches;
		
		//settings
		DoubleReal rt_pair_dist = param_.getValue("rt_pair_dist");
		DoubleReal rt_stdev_low = param_.getValue("rt_stdev_low");
		DoubleReal rt_stdev_high = param_.getValue("rt_stdev_high");
		DoubleReal mz_stdev = param_.getValue("mz_stdev");
		DoubleReal mz_pair_dist = param_.getValue("mz_pair_dist");

		// check each feature
		for (RefMap::const_iterator it=model_ref.begin(); it!=model_ref.end(); ++it)
		{
			RefMap::const_iterator range = lower_bound(model_ref.begin(),model_ref.end(),it->getRT()+rt_pair_dist - 2.0*rt_stdev_low, ConsensusFeature::NthPositionLess<0>());
			while (range!=model_ref.end() && range->getRT() <= it->getRT()+rt_pair_dist + 2.0*rt_stdev_high)
			{
				if (range->getCharge() == it->getCharge()
					&& range->getMZ() >=it->getMZ()+mz_pair_dist/it->getCharge()-2.0*mz_stdev  
					&& range->getMZ() <=it->getMZ()+mz_pair_dist/it->getCharge()+2.0*mz_stdev) // TODO magic 2.0 ??
				{
					
					DoubleReal score =  PValue_(range->getMZ() - it->getMZ(), mz_pair_dist/it->getCharge(), mz_stdev, mz_stdev)
														* PValue_(range->getRT() - it->getRT(), rt_pair_dist, rt_stdev_low, rt_stdev_high);
					matches.push_back(ConsensusFeature(0,it->begin()->getElementIndex(),*it));
					matches.back().insert(1,range->begin()->getElementIndex(),*range);
					matches.back().setQuality(score);
					matches.back().computeConsensus();
				}
				++range;
			}
		}
		
		//**********************************************************
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
		// Very useful for checking the results, and the ids have no real meaning anyway :-) // TODO sort in algorithm?
		result_map.sortByNthPosition(RawDataPoint2D::MZ);
	}
}
