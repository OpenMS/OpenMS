// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#pragma once
 
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/CommonEnums.h>
#include <OpenMS/CONCEPT/Types.h>
#include <functional>
#include <sstream>
 
namespace OpenMS
{
    /** 
    @brief This class is used for hydrophobicity profiling of peptides.
    Currently there are 6 scales available.
    The scales are implemented in:
    @ref HydrophobicityScaleMethod 
    */
  
    class OPENMS_DLLAPI HydrophobicityProfile {

 public:
    
    /// @brief Calculates hydrophobicity profile per residue
    /// @param seq The amino acid sequence for which to calculate the profile
    /// @param scale The scale to use for the calculation    
    /// @return Hydrophobicity profile
    std::vector<double> computeProfile(
        const AASequence& seq, 
        const HydrophobicityScaleMethod scale = HydrophobicityScaleMethod::KYTE_DOOLITTLE
    );
    
    /// @brief Calculates windowed hydrophobicity profile
    /// @param seqTthe Amino acid Sequence for which to calculate the profile
    /// @param window_size The size of the sliding window for the calculation
    /// @param scale The scale to use for the calculation
    /// @return Windowed hydrophobicity profile
    std::vector<double> computeWindowedProfile(
        const AASequence& seq,
        Size window_size = 7,
        const HydrophobicityScaleMethod scale = HydrophobicityScaleMethod::KYTE_DOOLITTLE
    );
    
    /// @brief Calculates the GRAVY score
    /// @param seq The amino acid sequence for which to calculate the GRAVY score
    /// @return The GRAVY score
    double computeGRAVY(const AASequence& seq);
    
    /// @brief Calculates hydrophobic moments
    /// @param seq The Amino acid sequence
    /// @param window_size The size of the sliding window for the calculation
    /// @return Hydrophobic moments of the amino acid sequence
    std::vector<double> computeHydrophobicMoment(
        const AASequence& seq,
        Size window_size = 11,
        double angle = 100.0  // degrees; 100=alpha-helix, 160=beta-sheet
    );  

    //just for testing will be removed later
    int testfunktion();
    
}; 
   
} // namespace OpenMS