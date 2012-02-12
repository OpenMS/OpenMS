// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#ifdef Debug_StablePairFinder
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla) ;
#endif
#define VV_(bla) V_(""#bla": " << bla);

using namespace std;

namespace OpenMS
{

  StablePairFinder::StablePairFinder(): Base()
  {
    //set the name for DefaultParamHandler error messages
    Base::setName(getProductName());

    defaults_.setValue("second_nearest_gap", 2.0, "The distance to the second nearest neighbors must be larger by this factor than the distance to the matching element itself");
    defaults_.setMinFloat("second_nearest_gap", 1.0);

		defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account)");
		defaults_.setValidStrings("use_identifications", StringList::create("true,false"));

		defaults_.insert("", FeatureDistance().getDefaults());

    Base::defaultsToParam_();
  }

  void StablePairFinder::updateMembers_()
  {
    V_("@@@ StablePairFinder::updateMembers_()");

    second_nearest_gap_ = param_.getValue("second_nearest_gap");
		use_IDs_ = String(param_.getValue("use_identifications")) == "true";
  }

  void StablePairFinder::run(const std::vector<ConsensusMap>& input_maps,
														 ConsensusMap &result_map)
  {
    // empty output destination:
    result_map.clear(false);

		// sanity checks:
    if (input_maps.size() != 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "exactly two input maps required");
    }
    checkIds_(input_maps);

		// set up the distance functor:
		DoubleReal max_intensity = max(input_maps[0].getMaxInt(), 
																	 input_maps[1].getMaxInt());
		Param distance_params = param_.copy("");
		distance_params.remove("use_identifications");
		distance_params.remove("second_nearest_gap");
		FeatureDistance feature_distance(max_intensity, false);
		feature_distance.setParameters(distance_params);

		// keep track of pairing:
    std::vector<bool> is_singleton[2];
    is_singleton[0].resize(input_maps[0].size(), true);
    is_singleton[1].resize(input_maps[1].size(), true);

		typedef pair<DoubleReal, DoubleReal> DoublePair;
		DoublePair init = make_pair(FeatureDistance::infinity,
																FeatureDistance::infinity);

		// for every element in map 0:
		// - index of nearest neighbor in map 1:
    vector<UInt> nn_index_0(input_maps[0].size(), UInt(-1));
		// - distances to nearest and second-nearest neighbors in map 1:
    vector<DoublePair> nn_distance_0(input_maps[0].size(), init);

		// for every element in map 1:
		// - index of nearest neighbor in map 0:
    vector<UInt> nn_index_1(input_maps[1].size(), UInt(-1));
		// - distances to nearest and second-nearest neighbors in map 0:
    vector<DoublePair> nn_distance_1(input_maps[1].size(), init);

    // iterate over all feature pairs, find nearest neighbors:
    for (UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0)
    {
			const ConsensusFeature& feat0 = input_maps[0][fi0];

      for (UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1)
      {
				const ConsensusFeature& feat1 = input_maps[1][fi1];

				if (use_IDs_ && !compatibleIDs_(feat0, feat1)) // check peptide IDs
				{
					continue; // mismatch
				}

        pair<bool, DoubleReal> result = feature_distance(feat0, feat1);
				DoubleReal distance = result.second;
				// we only care if distance constraints are satisfied for "best
				// matches", not for second-best; this means that second-best distances
				// can become smaller than best distances!
				bool valid = result.first;

				// update entries for map 0:
				if (distance < nn_distance_0[fi0].second)
				{
					if (valid && (distance < nn_distance_0[fi0].first))
					{
						nn_distance_0[fi0].second = nn_distance_0[fi0].first;
						nn_distance_0[fi0].first = distance;
						nn_index_0[fi0] = fi1;
					}
					else nn_distance_0[fi0].second = distance;
				}
				// update entries for map 1:
				if (distance < nn_distance_1[fi1].second)
				{
					if (valid && (distance < nn_distance_1[fi1].first))
					{
						nn_distance_1[fi1].second = nn_distance_1[fi1].first;
						nn_distance_1[fi1].first = distance;
						nn_index_1[fi1] = fi0;
					}
					else nn_distance_1[fi1].second = distance;
				}
			}
    }

    // if features from the two maps are nearest neighbors of each other, they
    // can become a pair:
    for (UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0)
    {
			UInt fi1 = nn_index_0[fi0]; // nearest neighbor of "fi0" in map 1
			// cout << "index: " << fi0 << ", RT: " << input_maps[0][fi0].getRT()
			// 		 << ", MZ: " << input_maps[0][fi0].getMZ() << endl
			// 		 << "neighbor: " << fi1 << ", RT: " << input_maps[1][fi1].getRT()
			// 		 << ", MZ: " << input_maps[1][fi1].getMZ() << endl
			// 		 << "d(i,j): " << nn_distance_0[fi0].first << endl
			// 		 << "d2(i): " << nn_distance_0[fi0].second << endl
			// 		 << "d2(j): " << nn_distance_1[fi1].second << endl;		

			// criteria set by the parameters must be fulfilled:
			if ((nn_distance_0[fi0].first < FeatureDistance::infinity) && 
					(nn_distance_0[fi0].first * second_nearest_gap_ <= nn_distance_0[fi0].second))
			{
				// "fi0" satisfies constraints...
				if ((nn_index_1[fi1] == fi0) && 
						(nn_distance_1[fi1].first * second_nearest_gap_ <= nn_distance_1[fi1].second))
				{
					// ...nearest neighbor of "fi0" also satisfies constraints (yay!)
					// cout << "match!" << endl;
					result_map.push_back(ConsensusFeature());
					ConsensusFeature& f = result_map.back();
					
					f.insert(input_maps[0][fi0]);
					f.getPeptideIdentifications().insert(f.getPeptideIdentifications().end(), input_maps[0][fi0].getPeptideIdentifications().begin(), input_maps[0][fi0].getPeptideIdentifications().end());

					f.insert(input_maps[1][fi1]);
					f.getPeptideIdentifications().insert(f.getPeptideIdentifications().end(), input_maps[1][fi1].getPeptideIdentifications().begin(), input_maps[1][fi1].getPeptideIdentifications().end());

					f.computeConsensus();
					DoubleReal quality = 1.0 - nn_distance_0[fi0].first;
					DoubleReal quality0 = 1.0 - nn_distance_0[fi0].first * second_nearest_gap_ / nn_distance_0[fi0].second;
					DoubleReal quality1 = 1.0 - nn_distance_1[fi1].first * second_nearest_gap_ / nn_distance_1[fi1].second;
					quality = quality * quality0 * quality1; // TODO other formula?

					// incorporate existing quality values:
					Size size0 = max(input_maps[0][fi0].size(), size_t(1));
					Size size1 = max(input_maps[1][fi1].size(), size_t(1));
					// quality contribution from first map:
					quality0 = input_maps[0][fi0].getQuality() * (size0 - 1);
					// quality contribution from second map:
					quality1 = input_maps[1][fi1].getQuality() * (size1 - 1);
					f.setQuality((quality + quality0 + quality1) / (size0 + size1 - 1));
				
					is_singleton[0][fi0] = false;
					is_singleton[1][fi1] = false;
				}
			}
    }

    // write out unmatched consensus features
    for ( UInt input = 0; input <= 1; ++ input )
    {
      for ( UInt index = 0; index < input_maps[input].size(); ++index )
      {
        if ( is_singleton[input][index] )
        {
          result_map.push_back(input_maps[input][index]);
					if (result_map.back().size() < 2) // singleton consensus feature
					{
						result_map.back().setQuality(0.0);
					}
        }
      }
    }

    // canonical ordering for checking the results, and the ids have no real meaning anyway
    result_map.sortByMZ();

		// protein IDs and unassigned peptide IDs are added to the result by the
		// FeatureGroupingAlgorithm!
  }

	
	bool StablePairFinder::compatibleIDs_(const ConsensusFeature& feat1, const ConsensusFeature& feat2) const
	{
		vector<PeptideIdentification> pep1 = feat1.getPeptideIdentifications(),pep2 = feat2.getPeptideIdentifications();
		// a feature without identifications always matches:
		if (pep1.empty() || pep2.empty()) return true;
		set<AASequence> best1, best2;
		for (vector<PeptideIdentification>::iterator pep_it = pep1.begin(); pep_it != pep1.end(); ++pep_it)
		{
			if (pep_it->getHits().empty()) continue; // shouldn't be the case
			pep_it->sort();
			best1.insert(pep_it->getHits()[0].getSequence());
		}
		for (vector<PeptideIdentification>::iterator pep_it = pep2.begin(); pep_it != pep2.end(); ++pep_it)
		{
			if (pep_it->getHits().empty()) continue; // shouldn't be the case
			pep_it->sort();
			best2.insert(pep_it->getHits()[0].getSequence());
		}
		return (best1 == best2);
	}

}
