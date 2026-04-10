// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#pragma once
 
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <functional>
#include <sstream>
 
namespace OpenMS
{

    /*
      @ingroup Chemistry

      @brief This class is used for Hydrophobicity Profiling of Peptides
  */

    class OPENMS_DLLAPI HydrophobicityProfile {
 public:
    
    /// @brief calculates Hydrophobicity profile per residue
    /// @param seq: the Aminoacid Sequence for which to calculate the profile
    /// @param scale: the scale to use for the calculation 
    /// possible scales are:                                                                       
    /// KYTEDOOLITTLE,
    /// EISENBERG,
    /// HOPPWOODS,
    /// BULLBREESE,
    /// BLACKMOULD,
    /// GUY
    /// @return Hydrophobicity Profile
    std::vector<double> computeProfile(
        const AASequence& seq, 
        const HydrophobicityScaleNumber scale = KYTEDOOLITTLE
    );
    
    /// @brief calculates windowed Hydrophobicity profile
    /// @param seq the Aminoacid Sequence for which to calculate the profile
    /// @param window_size the size of the sliding window for the calculation
    /// @param scale the scale to use for the calculation
    /// @return windowed hydrophobicity profile
    std::vector<double> computeWindowedProfile(
        const AASequence& seq,
        Size window_size = 7,
        const HydrophobicityScaleNumber scale = KYTEDOOLITTLE
    );
    
    // GRAVY score
    /// @brief calculates GRAVY score
    /// @param seq the Aminoacid Sequence for which to calculate the GRAVY score
    /// @return the GRAVY score
    double computeGRAVY(const AASequence& seq);
    
    // Hydrophobic moment (Eisenberg et al., 1982)
    /// @brief calculates Hydrophobic moments
    /// @param seq the Amino acid sequence
    /// @param window_size the size of the sliding window for the calculation
    /// @return ?????????
    std::vector<double> computeHydrophobicMoment(
        const AASequence& seq,
        Size window_size = 11,
        double angle = 100.0  // degrees; 100=alpha-helix, 160=beta-sheet
    );  

    //just for testing will be removed later
    int testfunktion();
    
}; 
   
} // namespace OpenMS