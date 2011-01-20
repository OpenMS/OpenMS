// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>

using namespace std;

namespace OpenMS
{

	FeatureGroupingAlgorithmQT::FeatureGroupingAlgorithmQT()
		: FeatureGroupingAlgorithm()
	{
		setName("FeatureGroupingAlgorithmQT");
		defaults_.insert("", QTClusterFinder().getParameters());
		defaultsToParam_();
	}

	FeatureGroupingAlgorithmQT::~FeatureGroupingAlgorithmQT()
	{
	}

	template <typename MapType>
	void FeatureGroupingAlgorithmQT::group_(const vector<MapType>& maps,
																					ConsensusMap& out)
	{
		// check that the number of maps is ok:
		if (maps.size() < 2)
		{
			throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
																			 "At least two maps must be given!");
		}

		QTClusterFinder cluster_finder;
		cluster_finder.setParameters(param_.copy("", true));
		ConsensusMap result;
		cluster_finder.run(maps, out);

		// add protein IDs and unassigned peptide IDs to the result map here,
		// to keep the same order as the input maps (useful for output later):
		for (typename vector<MapType>::const_iterator map_it = maps.begin();
				 map_it != maps.end(); ++map_it)
		{
			// add protein identifications to result map:
			out.getProteinIdentifications().insert(
				out.getProteinIdentifications().end(),
				map_it->getProteinIdentifications().begin(),
				map_it->getProteinIdentifications().end());

			// add unassigned peptide identifications to result map:
			out.getUnassignedPeptideIdentifications().insert(
				out.getUnassignedPeptideIdentifications().end(),
				map_it->getUnassignedPeptideIdentifications().begin(),
				map_it->getUnassignedPeptideIdentifications().end());
		}

		// canonical ordering for checking the results:
		out.sortByQuality();
		out.sortByMaps();
		out.sortBySize();
		return;
	}


	void FeatureGroupingAlgorithmQT::group(const std::vector<FeatureMap<> >& maps,
																				 ConsensusMap& out)
	{
		group_(maps, out);
	}


	void FeatureGroupingAlgorithmQT::group(const std::vector<ConsensusMap>& maps,
																				 ConsensusMap& out)
	{
		group_(maps, out);

		// maybe TODO:
		// move the code below to FeatureLinker or FeatureGroupingAlgorithm, and/or
		// move the handling of file descriptions from FeatureLinker
		bool keep_subelements = String(param_.getValue("keep_subelements")) == 
			"true";
		if (keep_subelements)
		{
			// components of the output map are not the input maps themselves, but the
			// components of the input maps!

			// accumulate file descriptions from the input maps:
			// cout << "Updating file descriptions..." << endl;
			out.getFileDescriptions().clear();
			out.getFileDescriptions() = maps[0].getFileDescriptions();
			vector<Size> offsets(maps.size(), 0);
			for (Size i = 1; i < maps.size(); ++i)
			{
				const ConsensusMap& consensus = maps[i];
				Size offset = offsets[i - 1] + consensus.getFileDescriptions().size();
				for (ConsensusMap::FileDescriptions::const_iterator desc_it =
							 consensus.getFileDescriptions().begin(); desc_it !=
							 consensus.getFileDescriptions().end(); ++desc_it)
				{
					out.getFileDescriptions()[desc_it->first + offset] = 
						desc_it->second;
				}
				offsets[i] = offset;
			}

			// look-up table: input map -> unique ID -> consensus feature
			// cout << "Creating look-up table..." << endl;
			vector<map<UInt64, ConsensusMap::ConstIterator> > id_lookup(maps.size());
			for (Size i = 0; i < maps.size(); ++i)
			{
				const ConsensusMap& consensus = maps[i];
				for (ConsensusMap::ConstIterator feat_it = consensus.begin();
						 feat_it != consensus.end(); ++feat_it)
				{
          // do NOT use 'id_lookup[i][feat_it->getUniqueId()]=feat_it;' here as you will get
          // "attempt to copy- construct an iterator from a singular iterator." in STL debug mode
					id_lookup[i].insert(std::pair<UInt64, ConsensusMap::ConstIterator>(feat_it->getUniqueId(), feat_it));
				}			
			}
			// adjust the consensus features:
			// cout << "Adjusting consensus features..." << endl;
			for (ConsensusMap::iterator cons_it = out.begin(); cons_it != out.end(); 
					 ++cons_it)
			{
				ConsensusFeature adjusted = ConsensusFeature(
					static_cast<BaseFeature>(*cons_it)); // remove sub-features
				for (ConsensusFeature::HandleSetType::const_iterator sub_it = 
							 cons_it->getFeatures().begin(); sub_it !=
							 cons_it->getFeatures().end(); ++sub_it)
				{
					UInt64 id = sub_it->getUniqueId();
					Size map_index = sub_it->getMapIndex();
					ConsensusMap::ConstIterator origin = id_lookup[map_index][id];
					for (ConsensusFeature::HandleSetType::const_iterator handle_it = 
								 origin->getFeatures().begin(); handle_it != 
								 origin->getFeatures().end(); ++handle_it)
					{
						FeatureHandle handle = *handle_it;
						handle.setMapIndex(handle.getMapIndex() + offsets[map_index]);
						adjusted.insert(handle);
					}
				}
				*cons_it = adjusted;
			}
		}
	}

} // namespace OpenMS
