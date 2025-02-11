// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
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

         This uses the locations of MS2 precursors as seed positions.
    */
    void generateSeedList(const PeakMap & experiment, SeedList & seeds);


    /**
         @brief Generate a seed list based on a list of peptide identifications

         This uses the retention time of the MS2 precursor ("RT" MetaInfo entry) together with either the precursor m/z value ("MZ" MetaInfo entry) or - if @p use_peptide_mass is true - the monoisotopic mass and charge of the best peptide hit to define the seed position for a peptide identification.

         The peptide hits in @p peptides will be sorted if @p use_peptide_mass is true.
    */
    void generateSeedList(std::vector<PeptideIdentification> & peptides,
                          SeedList & seeds, bool use_peptide_mass = false);


    /**
         @brief Generate seed lists based on a consensus map

         This creates one seed list per constituent map of the consensus map. For each constituent map, a seed is generated at the position of each consensus feature that does not contain a sub-feature from this map. (The idea is to fill "holes" in the consensus map by looking for features explicitly at the positions of the "holes".)
         The seed lists are indexed by the IDs of the constituent maps in the consensus map.

         Note that the resulting seed lists use the retention time scale of the consensus map, which might be different from the original time scales of the experiments if e.g. the MapAligner tool was used to perform retention time correction as part of the alignment process.
    */
    void generateSeedLists(const ConsensusMap & consensus,
                           std::map<UInt64, SeedList> & seed_lists);


    /// Convert a list of seed positions to a feature map (expected format for FeatureFinder)
    void convertSeedList(const SeedList & seeds, FeatureMap & features);


    /// Convert a feature map with seed positions back to a simple list
    void convertSeedList(const FeatureMap & features, SeedList & seeds);

  };

}

