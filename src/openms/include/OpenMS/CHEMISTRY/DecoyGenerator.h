// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Sam Wein $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <unordered_map>

namespace OpenMS
{
  class AASequence;
  class NASequence;
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

      /**
        @brief Reverse RNA products while preserving composition, terminal chemistry and RNase cleavage positions.
        @param[in] rna RNA parent
        @param[in] rnase RNase name (unspecific cleavage is unsupported)
        @param[in] keep_five_prime Additionally keep the parent's first nucleotide
        @param[in] keep_three_prime Additionally keep the parent's last nucleotide
        @return Cleavage-preserving RNA decoy; an untransformable product is retained unchanged
        @throws Exception::IllegalArgument RNase products do not partition the parent
      */
      NASequence reversePeptides(const NASequence& rna, const std::string& rnase,
                                bool keep_five_prime = false, bool keep_three_prime = false) const;

      /**
        @brief Portably shuffle RNA products without changing composition or modifications.

        RNase cut contexts and terminal chemistry are retained. The identity threshold is a
        best-effort target: return the lowest-identity valid permutation within the budget,
        and report an unattained target. No nucleotide or modification is ever mutated.
        A zero budget returns the original sequence.
        Each RNA product starts from the configured seed so identical products
        with the same fixed positions have identical decoys across parents.
        @param[in] rna RNA parent
        @param[in] rnase RNase name (unspecific cleavage is unsupported)
        @param[in] max_attempts Shuffle budget per RNase product
        @param[in] identity_threshold Best-effort maximum sequence identity in [0,1]
        @param[in] keep_five_prime Additionally keep the parent's first nucleotide
        @param[in] keep_three_prime Additionally keep the parent's last nucleotide
        @return Best valid, composition-preserving RNA decoy found
        @throws Exception::InvalidParameter Invalid budget or identity target
        @throws Exception::IllegalArgument RNase products do not partition the parent
      */
      NASequence shufflePeptides(const NASequence& rna, const std::string& rnase,
                                int max_attempts = 100, double identity_threshold = 0.5,
                                bool keep_five_prime = false, bool keep_three_prime = false);
    
    private:
      // sequence identity by matching AAs
      static double SequenceIdentity_(const std::string& decoy, const std::string& target);

      // portable shuffle
      Math::RandomShuffler shuffler_;
      UInt64 seed_ = 0;

      // ensures that shuffling same peptide (in different proteins) leads to same decoy
      std::unordered_map<std::string, std::string> td_cache_;
  };
}
