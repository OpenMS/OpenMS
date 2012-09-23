// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    void generateSeedList(const MSExperiment<> & experiment, SeedList & seeds);


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
                           Map<UInt64, SeedList> & seed_lists);


    /// Convert a list of seed positions to a feature map (expected format for FeatureFinder)
    void convertSeedList(const SeedList & seeds, FeatureMap<> & features);


    /// Convert a feature map with seed positions back to a simple list
    void convertSeedList(const FeatureMap<> & features, SeedList & seeds);

  };

}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SEEDLISTGENERATOR_H
