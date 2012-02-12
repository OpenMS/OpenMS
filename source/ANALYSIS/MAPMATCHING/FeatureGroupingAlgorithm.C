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
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmIdentification.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

using namespace std;

namespace OpenMS
{
	//register products here
	void FeatureGroupingAlgorithm::registerChildren()
	{
		// deprecated:
		// Factory<FeatureGroupingAlgorithm>::registerProduct ( FeatureGroupingAlgorithmIdentification::getProductName(), &FeatureGroupingAlgorithmIdentification::create );
		Factory<FeatureGroupingAlgorithm>::registerProduct ( FeatureGroupingAlgorithmLabeled::getProductName(), &FeatureGroupingAlgorithmLabeled::create );
		Factory<FeatureGroupingAlgorithm>::registerProduct ( FeatureGroupingAlgorithmUnlabeled::getProductName(), &FeatureGroupingAlgorithmUnlabeled::create );
		Factory<FeatureGroupingAlgorithm>::registerProduct ( FeatureGroupingAlgorithmQT::getProductName(), &FeatureGroupingAlgorithmQT::create );
	}

	FeatureGroupingAlgorithm::FeatureGroupingAlgorithm()
		: DefaultParamHandler("FeatureGroupingAlgorithm")
	{
	}

	void FeatureGroupingAlgorithm::group(const vector<ConsensusMap>& maps, ConsensusMap& out)
	{
    LOG_WARN << "FeatureGroupingAlgorithm::group() does not support ConsensusMaps directly. Converting to FeatureMaps." << endl;

		vector< FeatureMap<> > maps_f;
    for (Size i=0; i<maps.size(); ++i)
    {
      FeatureMap<> fm;
      ConsensusMap::convert(maps[i], true, fm);
      maps_f.push_back(fm);
    }
    // call FeatureMap version of group()
    group(maps_f, out);
	}

	void FeatureGroupingAlgorithm::transferSubelements(
		const vector<ConsensusMap>& maps, ConsensusMap& out) const
	{
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
				// do NOT use "id_lookup[i][feat_it->getUniqueId()] = feat_it;" here as
				// you will get "attempt to copy- construct an iterator from a singular
				// iterator" in STL debug mode:
				id_lookup[i].insert(make_pair(feat_it->getUniqueId(), feat_it));
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

	FeatureGroupingAlgorithm::~FeatureGroupingAlgorithm()
	{
	}

}
