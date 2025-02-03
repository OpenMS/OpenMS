// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>


#include <string>
#include <utility> // for pair
#include <vector>
#include <map>

// #define DEBUG_MRMDECOY

namespace OpenMS
{
  /**
  @brief This class generates a TargetedExperiment object with decoys based on a TargetedExperiment object

  There are multiple methods to create the decoy transitions, the simplest ones
  are reverse and pseudo-reverse which reverse the sequence either completely or
  leaving the last (tryptic) AA untouched respectively.

  Another decoy generation method is "shuffle" which uses an algorithm similar
  to the one described in Lam, Henry, et al. (2010). "Artificial decoy spectral
  libraries for false discovery rate estimation in spectral library searching in
  proteomics".  Journal of Proteome Research 9, 605-610. It shuffles the amino
  acid sequence and shuffles the fragment ion intensities accordingly, however
  for this to work the fragment ions need to be matched to annotated before.


  First, the algorithm goes through all peptides and applies the decoy method to
  the target peptide sequence (pseudo-reverse, reverse or shuffle) in order to
  produce the decoy sequence. Then, for each peptide, the fragment ions in the
  target library are matched to their most likely origin (e.g. the ions are
  annotated with their ion series (a,b,y) and the fragment number and optionally
  a neutral loss (10 different neutral losses are currently implemented)). For
  each fragment ion from the target peptide, an equivalent ion is created for the
  decoy peptide with the same intensity (e.g. if the target peptide sequence has
  a b5 ion with a normalized intensity of 200, an equivalent b5 ion for the
  decoy sequence is created and assigned the intensity 200).
  Optionally, the m/z values are corrected to reflect the theoretical value rather
  than the experimental value in the library.

 */
  class OPENMS_DLLAPI MRMDecoy :
    public DefaultParamHandler,
    public ProgressLogger
  {

public:

    typedef std::vector<size_t> IndexType;

    MRMDecoy();

    /**
      @brief Generate decoys from a TargetedExperiment

      Will generate decoy peptides for each target peptide provided in exp and
      write them into the decoy experiment.

      Valid methods: shuffle, reverse, pseudo-reverse

      If theoretical is true, the target transitions will be returned but their
      masses will be adjusted to match the theoretical value of the fragment ion
      that is the most likely explanation for the product.

      mz_threshold is used for the matching of theoretical ion series to the observed one

      To generate decoys with different precursor mass, use the "switchKR" flag
      which switches terminal K/R (switches K to R and R to K). This generates
      different precursor m/z and ensures that the y ion series has a different
      mass. For a description of the procedure, see (supplemental material)

      Bruderer et al. Mol Cell Proteomics. 2017. 10.1074/mcp.RA117.000314.

    */
    void generateDecoys(const OpenMS::TargetedExperiment& exp,
                        OpenMS::TargetedExperiment& dec,
                        const String& method,
                        const double aim_decoy_fraction,
                        const bool switchKR,
                        const String& decoy_tag,
                        const int max_attempts,
                        const double identity_threshold,
                        const double precursor_mz_shift,
                        const double product_mz_shift,
                        const double product_mz_threshold,
                        const std::vector<String>& fragment_types,
                        const std::vector<size_t>& fragment_charges,
                        const bool enable_specific_losses,
                        const bool enable_unspecific_losses,
                        const int round_decPow = -4) const;

    /**
       @brief Switch the final Amino Acid of a tryptic peptide.
       E.g. If the last Amino Acid is "K" switch to "R" (and vice versa).

       @note If the last Amino Acid is neither "K" or "R", the last Amino Acid is changed to a random Amino Acid .
    */
    void switchKR(OpenMS::TargetedExperiment::Peptide& peptide) const;

    typedef std::vector<OpenMS::TargetedExperiment::Protein> ProteinVectorType;
    typedef std::vector<OpenMS::TargetedExperiment::Peptide> PeptideVectorType;
    typedef std::vector<OpenMS::ReactionMonitoringTransition> TransitionVectorType;

    typedef std::map<String, std::vector<const ReactionMonitoringTransition*> > PeptideTransitionMapType;

    /**
      @brief Compute relative identity (relative number of matches of amino
      acids at the same position) between two sequences.
    */
    float AASequenceIdentity(const String& sequence, const String& decoy) const;

    /**
      @brief Shuffle a peptide (with its modifications) sequence

      This function will shuffle the given peptide sequences and its
      modifications such that the resulting relative sequence identity is below
      identity_threshold.
    */
    OpenMS::TargetedExperiment::Peptide shufflePeptide(
                OpenMS::TargetedExperiment::Peptide peptide,
                const double identity_threshold,
                int seed = -1,
                const int max_attempts = 100) const;

    /**
      @brief Reverse a peptide sequence (with its modifications)

      @param peptide The peptide sequence and modifications
      @param keepN Whether to keep N terminus in place
      @param keepC Whether to keep C terminus in place
      @param const_pattern A list of AA to leave in place
    */
    static OpenMS::TargetedExperiment::Peptide reversePeptide(
                const OpenMS::TargetedExperiment::Peptide& peptide,
                const bool keepN,
                const bool keepC, 
                const String& const_pattern = String());

    /**
      @brief Find all residues in a sequence that should not be reversed / shuffled
      
      @param sequence The amino acid sequence
      @param keepN Whether to keep N terminus constant
      @param keepC Whether to keep C terminus constant
      @param keep_const_pattern A string containing the AA to not change (e.g. 'KRP')
    */
    static IndexType findFixedResidues(const std::string& sequence,
        bool keepN, bool keepC, const OpenMS::String& keep_const_pattern);

protected:

    /**
      @brief Check if a peptide has C or N terminal modifications
    */
    bool hasCNterminalMods_(const OpenMS::TargetedExperiment::Peptide& peptide, bool checkCterminalAA) const;

    /**
      @brief Find all K, R, P sites in a sequence to be set as fixed

      This method was adapted from the SpectraST decoy generator
    */
    IndexType findFixedResidues_(const std::string& sequence) const;

    /**
      @brief Find all K, R, P and C-/N-terminal sites in a sequence to be set as fixed

      This method was adapted from the SpectraST decoy generator
    */
    IndexType findFixedAndTermResidues_(const std::string& sequence) const;

    /**
      @brief Pseudo-reverse a peptide sequence (with its modifications)

      @note Pseudo reverses a peptide sequence, leaving the C terminus (the
            last AA constant)
    */
    OpenMS::TargetedExperiment::Peptide pseudoreversePeptide_(
      const OpenMS::TargetedExperiment::Peptide& peptide) const;

    /**
      @brief Reverse a peptide sequence (with its modifications)

      @note Does not keep N / C terminus in place.
    */
    OpenMS::TargetedExperiment::Peptide reversePeptide_(
      const OpenMS::TargetedExperiment::Peptide& peptide) const;

    /**
     @brief Convert a peptide to a string which contains the peptide sequence and modifications
    **/
    String getModifiedPeptideSequence_(const OpenMS::TargetedExperiment::Peptide& pep) const;

    /// Synchronize members with param class
    void updateMembers_() override;

    String keep_const_pattern_;
    bool keepN_;
    bool keepC_;
  };
}
