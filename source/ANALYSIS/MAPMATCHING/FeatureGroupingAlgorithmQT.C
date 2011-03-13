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
// $Maintainer: Hendrik Weisser $
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
	}

} // namespace OpenMS
