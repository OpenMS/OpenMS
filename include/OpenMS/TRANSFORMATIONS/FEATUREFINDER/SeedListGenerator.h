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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SEEDLISTGENERATOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SEEDLISTGENERATOR_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
	/**
    @brief Generate seed lists for feature detection

		Seed lists specify locations in an MS experiment where features are expected. Currently, only the "centroided" FeatureFinder algorithm (class @p FeatureFinderAlgorithmPicked) supports custom seed lists (in featureXML format).

		@experimental This is a new class that still needs to be tested.
  */
	class OPENMS_DLLAPI SeedListGenerator
	{
	public:

		/// List of seed positions
		typedef std::vector<DPosition<2> > SeedList;

		
		/// Default constructor
		SeedListGenerator();

		
		/**
			 @brief Generate a seed list based on an MS experiment

			 This uses the locations of MS2 precurors as seed positions.
		*/
		void generateSeedList(const MSExperiment<>& experiment, SeedList& seeds);

		
		/**
			 @brief Generate a seed list based on a list of peptide identifications

			 This uses the locations of the corresponding MS2 precursors ("RT" and "MZ" MetaInfo entries) as seed positions.
		*/
		void generateSeedList(const std::vector<PeptideIdentification>& peptides,
													SeedList& seeds);

		
		/**
			 @brief Generate seed lists based on a consensus map

			 This creates one seed list per constituent map of the consensus map. For each constituent map, a seed is generated at the position of each consensus feature that does not contain a sub-feature from this map. (The idea is to fill "holes" in the consensus map by looking for features explicitly at the positions of the "holes".)
			 The seed lists are indexed by the IDs of the constituent maps in the consensus map.

			 Note that the resulting seed lists use the retention time scale of the consensus map, which might be different from the original time scales of the experiments if e.g. the MapAligner tool was used to perform retention time correction as part of the alignment process.
		*/		
		void generateSeedLists(const ConsensusMap& consensus,
													 Map<UInt64, SeedList>& seed_lists);

		
		/// Convert a list of seed positions to a feature map (expected format for FeatureFinder)
		void convertSeedList(const SeedList& seeds, FeatureMap<>& features);


		/// Convert a feature map with seed positions back to a simple list
		void convertSeedList(const FeatureMap<>& features, SeedList& seeds);
		
	};

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SEEDLISTGENERATOR_H
