// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, kg290 $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <string>
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
        @brief Prepare repeat-preserving decoy generation for a protein database.

        Call this before adding target proteins. The same generator must receive every target
        sequence in the database before finalizeDeBruijn() is called, so repeated (k+1)-mers
        across proteins receive the same decoy residue.

        @param[in] k Length of the k-mer vertices used by the de Bruijn construction
        @param[in] seed Seed for deterministic edge relabeling
        @param[in] keep_residues Residues that should remain unchanged in decoys (e.g. "KR")
      */
      void startDeBruijn(Size k, UInt64 seed = 4711, const std::string& keep_residues = "");

      /**
        @brief Add one unmodified target protein to the database being prepared.

        @param[in] protein Target protein to add
      */
      void addProteinToDeBruijn(const AASequence& protein);

      /**
        @brief Assign replacement residues to all distinct (k+1)-mers collected so far.

        Replacement residues are sampled using the target database's amino-acid frequencies.
      */
      void finalizeDeBruijn();

      /**
        @brief Generate a repeat-preserving decoy for a protein added before finalization.

        @param[in] protein Target protein included before finalizeDeBruijn()
        @return A decoy sequence with the same length as the target
      */
      AASequence deBruijn(const AASequence& protein) const;

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

      // ensures that shuffling same peptide (in different proteins) leads to same decoy
      std::unordered_map<std::string, std::string> td_cache_;

      Size debruijn_k_ = 0;
      UInt64 debruijn_seed_ = 4711;
      bool debruijn_ready_ = false;
      std::string debruijn_keep_residues_;
      std::unordered_map<std::string, Size> debruijn_edge_counts_;
      std::unordered_map<char, Size> debruijn_residue_counts_;
      std::unordered_map<std::string, char> debruijn_edge_labels_;
  };
}

