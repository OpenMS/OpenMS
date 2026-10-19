// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <unordered_map>
#include <vector>

namespace OpenMS
{
  class AASequence;
  class DigestionEnzymeProtein;

  /**
     @brief Methods to generate isobaric decoy sequences for DDA target-decoy searches.
  */
  class OPENMS_DLLAPI DecoyGenerator
  {
    public:
      // initializes random generator
      DecoyGenerator();

      // destructor
      ~DecoyGenerator() = default;

      // random seed for shuffling
      void setSeed(UInt64 seed);

      /**
        @brief Generate a repeat-preserving de Bruijn decoy database.

        The mapping is built across the complete input database, so repeated
        (k + 1)-mers receive the same replacement in every protein. No enzyme
        cleavage rule is required, which also supports nonspecific searches.

        @param[in] proteins Target protein sequences
        @param[in] k Vertex k-mer length (must be positive; default: 2)
        @return One decoy sequence for each input protein, in input order
      */
      std::vector<AASequence> generateDeBruijnDecoys(const std::vector<AASequence>& proteins, Size k = 2);

      /* 
         @brief reverses the protein sequence. 
         note: modifications are discarded
      */
      AASequence reverseProtein(const AASequence& protein) const;
    
      /* 
          @brief reverses the protein's peptide sequences between enzymatic cutting positions. 
          note: modifications are discarded
      */
      AASequence reversePeptides(const AASequence& protein, const std::string& protease) const;

      /**
        @brief Generate decoy protein sequences using shuffle algorithm
        
        Digests the protein using the specified protease and shuffles each resulting peptide
        to minimize sequence identity with the target. For top-down proteomics, use "no cleavage"
        as the protease to shuffle the entire protein as a single sequence.

        @param[in] protein The protein sequence to generate decoys from
        @param[in] protease The enzyme name (e.g., "Trypsin", "Trypsin/P", "no cleavage")
        @param[in] decoy_factor Number of decoy variants to generate per target peptide (default: 1)
        @return Vector of shuffled decoy sequences (one entry per decoy variant)
        
        @note
          - Peptides <= 2 amino acids are kept unchanged
          - Each peptide uses a fresh DecoyGenerator with fixed seed to ensure
            identical peptides produce identical decoys across different proteins
          - Modifications are discarded
          - Returns decoy_factor number of complete protein sequences
      */
      std::vector<AASequence> shuffle(const AASequence& protein, const std::string& protease, int decoy_factor = 1);

      /*
          @brief shuffle the protein's peptide sequences between enzymatic cutting positions.
          each peptide is shuffled @p max_attempts times to minimize sequence identity.

          Note: 
            - Generated decoys are retrieved from a cache to prevent that same peptide (in different proteins) 
              leads to different decoys.
            - modifications are discarded 
      */
      AASequence shufflePeptides(
            const AASequence& aas,
            const std::string& protease,
            const int max_attempts = 100
            );
    
    private:
      // sequence identity by matching AAs
      static double SequenceIdentity_(const std::string& decoy, const std::string& target);

      // portable shuffle
      Math::RandomShuffler shuffler_;

      /// Seed used to make database-wide decoy generation reproducible
      UInt64 seed_;

      // ensures that shuffling same peptide (in different proteins) leads to same decoy
      std::unordered_map<std::string, std::string> td_cache_;
  };
}

