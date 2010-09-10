// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTPairFinder.h>


namespace OpenMS
{

FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmQT");
		defaults_.insert("",QTPairFinder().getParameters());
		defaultsToParam_();
	}

FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT()
	{
	}

	void FeatureGroupingAlgorithmQT::group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out)
	{
		// check that the number of maps is ok
		if (maps.size()<2)
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"At least two maps must be given!");
		}
		// build a consensus map of the elements of the reference map (contains only singleton consensus elements)
		//		ConsensusMap::convert(reference_map_index, maps[reference_map_index],input[0]);

		// loop over all other maps, extend the groups
		QTPairFinder pair_finder;
		pair_finder.setParameters(param_.copy("", true));
		ConsensusMap result;
		pair_finder.run(maps, out);


		// add protein IDs and unassigned peptide IDs to the result map here,
		// to keep the same order as the input maps (useful for output later)
		for (std::vector<FeatureMap<> >::const_iterator map_it = maps.begin();map_it != maps.end(); ++map_it)
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

	void FeatureGroupingAlgorithmQT::group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)
	{
		// check that the number of maps is ok
		if (maps.size()<2)
		{
			throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"At least two maps must be given!");
		}

		// build a consensus map of the elements of the reference map (contains only singleton consensus elements)
		//		ConsensusMap::convert(reference_map_index, maps[reference_map_index],input[0]);

		// loop over all other maps, extend the groups
		QTPairFinder pair_finder;
		pair_finder.setParameters(param_.copy("", true));
		pair_finder.run(maps, out);

		out.getFileDescriptions().clear();

		// add protein IDs and unassigned peptide IDs to the result map here,
		// to keep the same order as the input maps (useful for output later)
		Size file_description_offset=0;
		bool keep_subelements = String(param_.getValue("keep_subelements")) == "true";
		for (std::vector<ConsensusMap>::const_iterator map_it = maps.begin();map_it != maps.end(); ++map_it)
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

			if (keep_subelements)
			{
				for (ConsensusMap::FileDescriptions::const_iterator description_it=map_it->getFileDescriptions().begin();description_it!=map_it->getFileDescriptions().end();++description_it)
				{
					out.getFileDescriptions()[description_it->first+file_description_offset]=description_it->second;
				}
				file_description_offset+=map_it->getFileDescriptions().size();
			}
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
