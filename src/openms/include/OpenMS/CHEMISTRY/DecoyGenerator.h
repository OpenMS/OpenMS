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

      /* 
         @brief reverses the protein sequence. 
         note: modifications are discarded
      */
      AASequence reverseProtein(const AASequence& protein) const;
    
      /* 
          @brief reverses the protein's peptide sequences between enzymatic cutting positions. 
          note: modifications are discarded
      */
      AASequence reversePeptides(const AASequence& protein, const String& protease) const;

      /* 
          @brief shuffle the protein's peptide sequences between enzymatic cutting positions.
          each peptide is shuffled @param max_attempts times to minimize sequence identity.
          note: modifications are discarded 
      */
      AASequence shufflePeptides(
            const AASequence& aas,
            const String& protease,
            const int max_attempts = 100
            );
    
    private:
      // sequence identity by matching AAs
      static double SequenceIdentity_(const String& decoy, const String& target);

      // portable shuffle
      Math::RandomShuffler shuffler_;
  };
}

