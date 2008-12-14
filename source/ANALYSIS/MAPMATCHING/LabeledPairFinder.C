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
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

using namespace std;
namespace OpenMS
{
	using namespace Math;

	LabeledPairFinder::LabeledPairFinder()
		: BaseGroupFinder()
	{
		setName("LabeledPairFinder");
		defaults_.setValue("rt_estimate", "true", "If 'true' the optimal RT pair distance and deviation are estimated by "
																							 "fitting a gaussian distribution to the histogram of pair distance. "
																							 "Note that this works only datasets with a significant amount of pairs! "
																							 "If 'false' the parameters 'rt_pair_dist', 'rt_dev_low' "
																							 "and 'rt_dev_high' define the optimal distance.");
		defaults_.setValidStrings("rt_estimate", StringList::create("true,false"));
		defaults_.setValue("rt_pair_dist", -20.0, "optimal pair distance in RT [sec] from light to heavy feature");
		defaults_.setValue("rt_dev_low", 15.0, "maximum allowed deviation below optimal retention time distance");
		defaults_.setMinFloat("rt_dev_low",0.0);
		defaults_.setValue("rt_dev_high", 15.0, "maximum allowed deviation above optimal retention time distance");
		defaults_.setMinFloat("rt_dev_high",0.0);
		
		defaults_.setValue("mz_pair_dist", 4.0, "optimal pair distance in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)");
		defaults_.setValue("mz_dev", 0.05, "maximum allowed deviation from optimal m/z distance\n");
		defaults_.setMinFloat("mz_dev",0.0);

		defaultsToParam_();
	}

	void LabeledPairFinder::run(const vector<ConsensusMap>& input_maps, ConsensusMap& result_map) 
	{
		if (input_maps.size()!=1) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"exactly one input map required");
		if (result_map.getFileDescriptions().size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"two file descriptions required");
		if (result_map.getFileDescriptions().begin()->second.filename!=result_map.getFileDescriptions().rbegin()->second.filename) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"the two file descriptions have to contain the same file name");
		checkIds_(input_maps);
		
		//look up the light and heavy index
		UInt light_index = numeric_limits<UInt>::max();
		UInt heavy_index = numeric_limits<UInt>::max();	
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
		if (light_index == numeric_limits<UInt>::max() || heavy_index == numeric_limits<UInt>::max())
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
		
		//estimate RT parameters
		if (param_.getValue("rt_estimate")=="true")
		{
			//find all possible RT distances of features with the same charge and a good m/z distance
			vector<DoubleReal> dists;
			dists.reserve(model_ref.size());
			for (RefMap::const_iterator it=model_ref.begin(); it!=model_ref.end(); ++it)
			{
				for (RefMap::const_iterator it2=model_ref.begin(); it2!=model_ref.end(); ++it2)
				{
					if (it2->getCharge() == it->getCharge()
						&& it2->getMZ() >=it->getMZ()+mz_pair_dist/it->getCharge()-mz_dev  
						&& it2->getMZ() <=it->getMZ()+mz_pair_dist/it->getCharge()+mz_dev)
					{
						dists.push_back(it2->getRT()-it->getRT());
					}
				}
			}
			if (dists.empty())
			{
				cout << "Warning: Could not find pairs for RT distance estimation. The manual settings are used!" << endl;
			}
			else
			{
				if (dists.size()<50)
				{
					cout << "Warning: Found only " << dists.size() << " pairs. The estimated shift and std deviation are probably not reliable!" << endl;
				}
				//--------------------------- estimate initial parameters of fit ---------------------------
				GaussFitter::GaussFitResult result;
				//first estimate of the optimal shift: median of the distances
				sort(dists.begin(),dists.end());
				UInt median_index = dists.size()/2; 
				result.x0 = dists[median_index]; 
				//create histogram of distances
				//consider only the maximum of pairs, centered around the optimal shift
				Int max_pairs = model_ref.size()/2;
				UInt start_index = max(0,(Int)(median_index) - max_pairs/2);
				UInt end_index = min((Int)(dists.size()-1),(Int)(median_index) + max_pairs/2);
				DoubleReal start_value = dists[start_index];
				DoubleReal end_value = dists[end_index];
				DoubleReal bin_step = fabs(end_value-start_value)/100;
				Math::Histogram<> hist(start_value,end_value,bin_step);
				//std::cout << "Histogram from " << start_value << " to " << end_value << " (bin size " << bin_step << ")" << endl;
				for (UInt i=start_index; i<=end_index; ++i)
				{
					//std::cout << "i: " << dists[i] << endl;		
					hist.inc(dists[i]);
				}
				//std::cout << hist << endl;
				dists.clear();
				//determine median of bins (uniform background distribution)
				vector<UInt> bins(hist.begin(),hist.end());
				sort(bins.begin(),bins.end());
				UInt bin_median = bins[bins.size()/2];
				bins.clear();
				//estimate scale A: maximum of the histogram
				UInt max_value = hist.maxValue();
				result.A = max_value-bin_median;
				//overwrite estimate of x0 with the position of the highest bin
				for (UInt i=0;i<hist.size();++i)
				{
					if (hist[i]==max_value)
					{
						result.x0 = hist.centerOfBin(i);
						break;
					}
				}
				//estimate sigma: first time the count is less or equal the median count in the histogram
				DoubleReal pos = result.x0;
				while (pos>start_value && hist.binValue(pos)>bin_median)
				{
					pos -= bin_step;
				}
				DoubleReal sigma_low =  result.x0 - pos;
				pos = result.x0;
				while (pos<end_value && hist.binValue(pos)>bin_median)
				{
					pos += bin_step;
				}
				DoubleReal sigma_high = pos - result.x0;		
				result.sigma = (sigma_high + sigma_low)/6.0;
				//cout << "estimated optimal RT distance (before fit): " << result.x0 << endl;
				//cout << "estimated allowed deviation (before fit): " << result.sigma*3.0 << endl;
				//--------------------------- do gauss fit ---------------------------
				vector<DPosition<2> > points(hist.size());
				for (UInt i=0;i<hist.size();++i)
				{
					points[i][0] = hist.centerOfBin(i);
					points[i][1] = max(0u,hist[i]);
				}
				GaussFitter fitter;
				fitter.setInitialParameters(result);
				result = fitter.fit(points);
				cout << "estimated optimal RT distance: " << result.x0 << endl;
				cout << "estimated allowed deviation: " << result.sigma*3.0 << endl;
				rt_pair_dist = result.x0;
				rt_dev_low = result.sigma*3.0;
				rt_dev_high = result.sigma*3.0;
			}
		}

		
		// check each feature
		for (RefMap::const_iterator it=model_ref.begin(); it!=model_ref.end(); ++it)
		{
			RefMap::const_iterator it2 = lower_bound(model_ref.begin(),model_ref.end(),it->getRT()+rt_pair_dist - rt_dev_low, ConsensusFeature::RTLess());
			while (it2!=model_ref.end() && it2->getRT() <= it->getRT()+rt_pair_dist + rt_dev_high)
			{
				if ( it2->getCharge() == it->getCharge() &&
						 it2->getMZ() >=it->getMZ()+mz_pair_dist/it->getCharge()-mz_dev &&
						 it2->getMZ() <=it->getMZ()+mz_pair_dist/it->getCharge()+mz_dev
					 )
				{
					
					DoubleReal score = sqrt(
																	PValue_(it2->getMZ() - it->getMZ(), mz_pair_dist/it->getCharge(), mz_dev, mz_dev) *
																	PValue_(it2->getRT() - it->getRT(), rt_pair_dist, rt_dev_low, rt_dev_high)
																 );
					matches.push_back(ConsensusFeature(light_index,it->begin()->getElementIndex(),*it));
					matches.back().clearMetaInfo();
					matches.back().insert(heavy_index,it2->begin()->getElementIndex(),*it2);
					matches.back().setQuality(score);
					matches.back().setCharge(it->getCharge());
					matches.back().computeConsensus();
				}
				++it2;
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
		
		//Add protein identifications to result map
		for (UInt i=0; i<input_maps.size(); ++i)
		{
			result_map.getProteinIdentifications().insert(result_map.getProteinIdentifications().end(),input_maps[i].getProteinIdentifications().begin(), input_maps[i].getProteinIdentifications().end());
		}

		//Add unassigned peptide identifications to result map
		for (UInt i=0; i<input_maps.size(); ++i)
		{
			result_map.getUnassignedPeptideIdentifications().insert(result_map.getUnassignedPeptideIdentifications().end(),input_maps[i].getUnassignedPeptideIdentifications().begin(), input_maps[i].getUnassignedPeptideIdentifications().end());
		}
		
		// Very useful for checking the results, and the ids have no real meaning anyway
		result_map.sortByMZ();
	}
}
