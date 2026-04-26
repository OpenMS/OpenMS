// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Markus Apel, Nora Heese $
// --------------------------------------------------------------------------
 
#pragma once
 
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/Types.h>
#include <functional>
#include <sstream>
 
namespace OpenMS
{
  /** 
  @ingroup Chemistry

  @brief This class is used for hydrophobicity profiling of peptides.

  The possible scales are listed here: @ref HydrophobicityScaleMethod

  References for the hydrophobicity scales:

  - Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 1982 May 5;157(1):105-32. doi: 10.1016/0022-2836(82)90515-0. PMID: 7108955.
  - Eisenberg D, Schwarz E, Komaromy M, Wall R. Analysis of membrane and surface protein sequences with the hydrophobic moment plot. J Mol Biol. 1984 Oct 15;179(1):125-42. doi: 10.1016/0022-2836(84)90309-7. PMID: 6502707.
  - T.P. Hopp, & K.R. Woods, Prediction of protein antigenic determinants from amino acid sequences., Proc. Natl. Acad. Sci. U.S.A. 78 (6) 3824-3828, https://doi.org/10.1073/pnas.78.6.3824 (1981).
  - Henry B. Bull, Keith Breese, Surface tension of amino acid solutions: A hydrophobicity scale of the amino acid residues, Archives of Biochemistry and Biophysics, Volume 161, Issue 2, 1974, Pages 665-670, ISSN 0003-9861, https://doi.org/10.1016/0003-9861(74)90352-X. (https://www.sciencedirect.com/science/article/pii/000398617490352X)
  - Black SD, Mould DR. Development of hydrophobicity parameters to analyze proteins which bear post- or cotranslational modifications. Anal Biochem. 1991 Feb 15;193(1):72-82. doi: 10.1016/0003-2697(91)90045-u. PMID: 2042744.
  - Guy, H. R. (1985). Amino acid side-chain partition energies and distribution of residues in soluble proteins. Biophysical journal, 47(1), 61-70.

  Paper for calculating the hydrophobic moments (Eisenberg et al. 1984):
  
  - Eisenberg D, Schwarz E, Komaromy M, Wall R. Analysis of membrane and surface protein sequences with the hydrophobic moment plot. J Mol Biol. 1984 Oct 15;179(1):125-42. doi: 10.1016/0022-2836(84)90309-7. PMID: 6502707.
  - Eisenberg D, Weiss RM, Terwilliger TC. The hydrophobic moment detects periodicity in protein hydrophobicity. Proc Natl Acad Sci U S A. 1984 Jan;81(1):140-4. doi: 10.1073/pnas.81.1.140. PMID: 6582470; PMCID: PMC344626.

  Referenced code for calculating the hydrophobic moment of a windowed peptide \n
  R package R.Peptides: https://github.com/dosorio/Peptides/blob/master/R/hmoment.R
  */
  
  class OPENMS_DLLAPI HydrophobicityProfile 
  {

public:

    /// @brief Calculates the GRAVY score
    /// @param seq The amino acid sequence for which to calculate the GRAVY score
    /// @return The GRAVY score of a given peptide
    /// @throw Exception::InvalidValue Throws an exception if the sequence is empty or contains invalid symbols (sequence can only contain the standard 20 Amino acids)
    double computeGRAVY(const AASequence& seq) const;
    
    /// @brief Calculates hydrophobicity profile per residue
    /// @param seq The amino acid sequence for which to calculate the profile
    /// @param scale The scale to use for the calculation    
    /// @return A vector which contains the hydrophobicity value for each amino acid
    /// @throw Exception::InvalidValue Throws an exception if the sequence is empty or contains invalid symbols (sequence can only contain the standard 20 Amino acids) or when an unknown scale is used
    std::vector<double> computeProfile(
      const AASequence& seq, 
      const HydrophobicityScaleMethod scale = HydrophobicityScaleMethod::KYTE_DOOLITTLE
    ) const;
    
    /// @brief Calculates windowed hydrophobicity profile
    /// @param seq The Amino acid Sequence for which to calculate the profile
    /// @param window_size The size of the sliding window for the calculation
    /// @param scale The scale to use for the calculation
    /// @return A vector which contains the average hydrophobicity for each window
    /// @throw Exception::InvalidValue Throws an exception if the sequence is empty or contains invalid symbols (sequence can only contain the standard 20 Amino acids) or when an unknown scale is used
    /// @throw Exception::InvalidSize Throws an exception when the window size is 0
    /// @warning When the window size is larger than the sequence the window size will be clamped to the sequence length
    std::vector<double> computeWindowedProfile(
      const AASequence& seq,
      const Size window_size = 7,
      const HydrophobicityScaleMethod scale = HydrophobicityScaleMethod::KYTE_DOOLITTLE
    ) const;
      
    /// @brief Calculates hydrophobic moments
    /// @param seq The Amino acid sequence
    /// @param window_size The size of the sliding window for the calculation
    /// @param angle The angle to use: 100 for alpha-helices, 160 for beta-sheets, other angles can be used
    /// @return A vector which contains the hydrophobic moment for each window
    /// @throw Exception::InvalidValue Throws an exception if the sequence is empty or contains invalid symbols (sequence can only contain the standard 20 Amino acids)
    /// @throw Exception::InvalidSize Throws an exception when the window size is 0
    /// @warning When the window size is larger than the sequence the window size will be clamped to the sequence length
    std::vector<double> computeHydrophobicMoment(
      const AASequence& seq,
      const Size window_size = 11,
      const double angle = 100.0  // degrees; 100 = alpha-helix, 160 = beta-sheet
    ) const;  
  }; 
} // namespace OpenMS