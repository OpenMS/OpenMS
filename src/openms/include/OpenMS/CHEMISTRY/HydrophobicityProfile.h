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
    /* class HydrophobicityProfile {
public:
    // Per-residue and windowed profiles
    static std::vector<double> computeProfile(
        const AASequence& seq, 
        const HydrophobicityScaleNumber scale = KYTEDOOLITTLE
    );
    
    static std::vector<double> computeWindowedProfile(
        const AASequence& seq,
        Size window_size = 7,
        const HydrophobicityScaleNumber scale = KYTEDOOLITTLE
    );
    
    // GRAVY score
    static double computeGRAVY(const AASequence& seq);
    
    // Hydrophobic moment (Eisenberg et al., 1982)
    static std::vector<double> computeHydrophobicMoment(
        const AASequence& seq,
        Size window_size = 11,
        double angle = 100.0  // degrees; 100=alpha-helix, 160=beta-sheet
    );
}; */
   
} // namespace OpenMS