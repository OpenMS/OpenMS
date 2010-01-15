// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>


namespace OpenMS
{

	FeatureGroupingAlgorithmUnlabeled::FeatureGroupingAlgorithmUnlabeled()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmUnlabeled");
		defaults_.insert("",StablePairFinder().getParameters());
		defaultsToParam_();
	}

	FeatureGroupingAlgorithmUnlabeled::~FeatureGroupingAlgorithmUnlabeled()
	{
	}

	void FeatureGroupingAlgorithmUnlabeled::group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out)
	{
		// check that the number of maps is ok
		if (maps.size()<2)
		{
		  throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"At least two maps must be given!");
		}

		// define reference map (the one with most peaks)
		Size reference_map_index = 0;
		Size max_count = 0;
		for (Size m=0; m<maps.size(); ++m)
		{
			if (maps[m].size()>max_count)
			{
				max_count = maps[m].size();
				reference_map_index = m;
			}
		}

		std::vector<ConsensusMap> input(2);

    // build a consensus map of the elements of the reference map (contains only singleton consensus elements)
		ConsensusMap::convert(reference_map_index, maps[reference_map_index],
													input[0]);

		// loop over all other maps, extend the groups
		StablePairFinder pair_finder;
		pair_finder.setParameters(param_.copy("", true));

		for (Size i = 0; i < maps.size(); ++i)
		{
			if (i != reference_map_index)
			{
				ConsensusMap::convert( i, maps[i], input[1] );
				// compute the consensus of the reference map and map i
				ConsensusMap result;
				pair_finder.run(input, result);
				input[0].swap(result);
			}
		}

		// replace result with temporary map
		out.swap(input[0]);
		// copy back the input maps (they have been deleted while swapping)
		out.getFileDescriptions() = input[0].getFileDescriptions();

		// add protein IDs and unassigned peptide IDs to the result map here,
		// to keep the same order as the input maps (useful for output later)
		for (std::vector<FeatureMap<> >::const_iterator map_it = maps.begin();
				 map_it != maps.end(); ++map_it)
		{
			// add protein identifications to result map
			out.getProteinIdentifications().insert(
				out.getProteinIdentifications().end(),
				map_it->getProteinIdentifications().begin(),
				map_it->getProteinIdentifications().end());

			// add unassigned peptide identifications to result map
			out.getUnassignedPeptideIdentifications().insert(
				out.getUnassignedPeptideIdentifications().end(),
				map_it->getUnassignedPeptideIdentifications().begin(),
				map_it->getUnassignedPeptideIdentifications().end());
		}

    // canonical ordering for checking the results, and the ids have no real meaning anyway
#if 1 // the way this was done in DelaunayPairFinder and StablePairFinder
    out.sortByMZ();
#else
    out.sortByQuality();
    out.sortByMaps();
    out.sortBySize();
#endif

		return;
	}

} // namespace OpenMS
