// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SeedListGenerator.h>

using namespace std;

namespace OpenMS
{
	SeedListGenerator::SeedListGenerator()
	{
	}

	
 	void SeedListGenerator::getSeedList(const MSExperiment<>& experiment,
																			SeedList& seeds)
	{
		seeds.clear();
		for (MSExperiment<>::ConstIterator exp_it = experiment.begin();
				 exp_it != experiment.end(); ++exp_it)
		{
			if (exp_it->getMSLevel() == 2) // MS2 spectrum -> look for precursor
			{
				MSExperiment<>::ConstIterator prec_it =
					experiment.getPrecursorSpectrum(exp_it);
				const vector<Precursor>& precursors = exp_it->getPrecursors();
				DPosition<2> point(prec_it->getRT(), precursors[0].getMZ());
				seeds.push_back(point);
			}
		}
	}


	void SeedListGenerator::getSeedList(const vector<PeptideIdentification>&
																			peptides, SeedList& seeds)
	{
		seeds.clear();
		for (vector<PeptideIdentification>::const_iterator pep_it =
					 peptides.begin(); pep_it != peptides.end(); ++pep_it)
		{
			DPosition<2> point(pep_it->getMetaValue("RT"),
												 pep_it->getMetaValue("MZ"));
			seeds.push_back(point);
		}
	}

	
	void SeedListGenerator::getSeedLists(const ConsensusMap& consensus,
																			 Map<UInt64, SeedList>& seed_lists)
	{
		seed_lists.clear();
		// iterate over all consensus features...
		for (ConsensusMap::ConstIterator cons_it = consensus.begin();
				 cons_it != consensus.end(); ++cons_it)
		{
			DPosition<2> point(cons_it->getRT(), cons_it->getMZ());
			// for each sub-map in the consensus map, add a seed at the position of
			// this consensus feature:
			for (ConsensusMap::FileDescriptions::ConstIterator file_it =
						 consensus.getFileDescriptions().begin(); file_it !=
						 consensus.getFileDescriptions().end(); ++file_it)
				seed_lists[file_it->first].push_back(point);
			// for each feature contained in the consensus feature, remove the seed of
			// the corresponding map:
			for (ConsensusFeature::HandleSetType::const_iterator feat_it =
						 cons_it->getFeatures().begin(); feat_it !=
						 cons_it->getFeatures().end(); ++feat_it)
			{
				seed_lists[feat_it->getMapIndex()].pop_back();
			}
			// this leaves seeds for maps where no feature was found near the
			// consensus position
		}
	}


	void SeedListGenerator::convert(const SeedList& seeds, FeatureMap<>& features)
	{
		features.clear(true); // "true" should really be a default value here...
		for (SeedList::const_iterator seed_it = seeds.begin();
				 seed_it != seeds.end(); ++seed_it)
		{
			Feature feature;
			feature.setRT(seed_it->getX());
			feature.setMZ(seed_it->getY());
			features.push_back(feature);
		}
	}


	void SeedListGenerator::convert(const FeatureMap<>& features, SeedList& seeds)
	{
		seeds.clear();
		for (FeatureMap<>::ConstIterator feat_it = features.begin();
				 feat_it != features.end(); ++feat_it)
		{
			DPosition<2> point(feat_it->getRT(), feat_it->getMZ());
			seeds.push_back(point);
		}
	}
	
} // namespace OpenMS
