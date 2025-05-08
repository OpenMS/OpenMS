// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
//#include <numeric>

namespace OpenMS
{
  class OPENMS_DLLAPI OPXLDataStructs
  {

    public:

      /**
       * @brief The ProteinProteinCrossLinkType enum enumerates possible types of Protein-Protein cross-linking reaction results. Cross-link, Mono-link or Loop-link.
       */
      enum ProteinProteinCrossLinkType
      {
        CROSS = 0,
        MONO = 1,
        LOOP = 2,
        NUMBER_OF_CROSS_LINK_TYPES
      };

      /**
       * @brief The ProteinProteinCrossLink struct represents a cross-link between two peptides in OpenPepXL.

          This struct completely defines a cross-link linking two peptides. It contains the two peptides alpha and beta as AASequences,
          the positions that are linked in each peptide, the mass and name of the linker and terminal specificities to distinguish linkers attached to side chains
          of the terminal residues from linkers attached to the termini themselves.
          The beta peptide can be an empty AASequence object in case of a mono- or loop-link.
          Used to represent a theoretical candidate for generation of theoretical spectra and matching to experimental spectra.
       */
      struct ProteinProteinCrossLink
      {
        const AASequence *alpha = nullptr; ///< longer peptide
        const AASequence *beta = nullptr; ///< shorter peptide (empty for mono-link), tie breaker: mass then lexicographical
        std::pair<SignedSize, SignedSize> cross_link_position; ///< index in alpha, beta or between alpha, alpha in loop-links
        double cross_linker_mass = 0;
        String cross_linker_name;
        ResidueModification::TermSpecificity term_spec_alpha = ResidueModification::TermSpecificity::ANYWHERE;
        ResidueModification::TermSpecificity term_spec_beta = ResidueModification::TermSpecificity::ANYWHERE;
        int precursor_correction = 0;

        ProteinProteinCrossLinkType getType() const
        {
          if (beta && !beta->empty()) return CROSS;

          if (cross_link_position.second == -1) return MONO;

          return LOOP;
        }

        bool operator==(const ProteinProteinCrossLink & other) const
        {
          return alpha == other.alpha &&
                      beta == other.beta &&
                      cross_link_position == other.cross_link_position &&
                      cross_linker_mass == other.cross_linker_mass &&
                      cross_linker_name == other.cross_linker_name &&
                      term_spec_alpha == other.term_spec_alpha &&
                      term_spec_beta == other.term_spec_beta &&
                      precursor_correction == other.precursor_correction;
        }
      };

      /**
       * @brief The CrossLinkSpectrumMatch struct represents a PSM between a ProteinProteinCrossLink and a spectrum in OpenPepXL.

          This struct contains a ProteinProteinCrossLink and indices to one or two spectra.
          It also contains the results of a match between the ProteinProteinCrossLink and these spectra as scores and peak annotations.
          Used as a temporary container to collect results efficiently, since only a few top matches will be kept for each experimental spectrum for output.
       */
      struct CrossLinkSpectrumMatch
      {
        /// structure of the cross-link
        ProteinProteinCrossLink cross_link;

        /// reference to pair of spectra
        Size scan_index_light = 0;
        Size scan_index_heavy = 0;

        /// final score
        double score = 0;

        /// rank among the matches to the same spectrum
        Size rank = 0;

        /// counts, scores and other data for xQuest-like output
        double xquest_score = 0;
        double pre_score = 0;
        double percTIC = 0;
        double wTIC = 0;
        double wTICold = 0;
        double int_sum = 0;
        double intsum_alpha = 0;
        double intsum_beta = 0;
        double total_current = 0;
        double precursor_error_ppm = 0;

        double match_odds = 0;
        double match_odds_alpha = 0;
        double match_odds_beta = 0;
        double log_occupancy = 0;
        double log_occupancy_alpha = 0;
        double log_occupancy_beta = 0;
        double xcorrx_max = 0;
        double xcorrc_max = 0;
        Size matched_linear_alpha = 0;
        Size matched_linear_beta = 0;
        Size matched_xlink_alpha = 0;
        Size matched_xlink_beta = 0;

        double num_iso_peaks_mean = 0;
        double num_iso_peaks_mean_linear_alpha = 0;
        double num_iso_peaks_mean_linear_beta = 0;
        double num_iso_peaks_mean_xlinks_alpha = 0;
        double num_iso_peaks_mean_xlinks_beta = 0;

        double ppm_error_abs_sum_linear_alpha = 0;
        double ppm_error_abs_sum_linear_beta = 0;
        double ppm_error_abs_sum_xlinks_alpha = 0;
        double ppm_error_abs_sum_xlinks_beta = 0;
        double ppm_error_abs_sum_linear = 0;
        double ppm_error_abs_sum_xlinks = 0;
        double ppm_error_abs_sum_alpha = 0;
        double ppm_error_abs_sum_beta = 0;
        double ppm_error_abs_sum = 0;

        int precursor_correction = 0;

        double precursor_total_intensity = 0;
        double precursor_target_intensity = 0;
        double precursor_signal_proportion = 0;
        Size precursor_target_peak_count = 0;
        Size precursor_residual_peak_count = 0;

        std::vector<PeptideHit::PeakAnnotation> frag_annotations;

        Size peptide_id_index = 0;
      };

      /**
        * @brief Comparator to sort CrossLinkSpectrumMatches by the main score

       */
      struct CLSMScoreComparator
      {
        bool operator() (const CrossLinkSpectrumMatch& a, const CrossLinkSpectrumMatch& b)
        {
          if (a.score == b.score)
          {
            // in rare cases when the sequences are the same, multiple candidates with different cross-linked positions can have the same score
            // that leads to ambiguous sorting and may cause differences between compilers
            // in those cases we prefer higher positions (just like the score),
            // because the lower position might be an N-term link, which is usually less likely and all other positions are equal (because the score is equal)
            if (a.cross_link.cross_link_position.first == b.cross_link.cross_link_position.first)
            {
              return a.cross_link.cross_link_position.second < b.cross_link.cross_link_position.second;
            }
            return a.cross_link.cross_link_position.first < b.cross_link.cross_link_position.first;
          }
          return a.score < b.score;
        }
      };

      /**
       * @brief The XLPrecursor struct represents a cross-link candidate in the process of filtering candidates by precursor masses in OpenPepXL.

          Since the precursor mass does not change depending on the exact linked residues, one XLPrecursor can represent several ProteinProteinCrossLinks
          in the process of filtering by precursor mass. The precursor mass is the sum of the masses of the two peptides and the cross-linker.
          This struct also contains the indices of the two peptides in a vector, so that the two peptides can be identified again.
          This precursor mass is enumerated for all peptide pairs in the protein database given as input to OpenPepXL
          and is one of the major contributors to the memory usage of the tool because of the squared complexity of this enumeration.
          Therefore this should be kept as compact as possible.
       */
      struct XLPrecursor
      {
        float precursor_mass{};
        unsigned int alpha_index = 0;
        unsigned int beta_index = 0;
        String alpha_seq;
        String beta_seq;
      };

      // comparator for sorting XLPrecursor vectors and using upper_bound and lower_bound using only a precursor mass
      /**
       * @brief The XLPrecursorComparator is a comparator for XLPrecursors, that allows direct comparison of the XLPrecursor precursor mass with double numbers.

          This comparator can be used to sort a vector of XLPrecursors by precursor mass and to search for XLPrecursors with precursor masses within a double type mass range.
       */
      struct XLPrecursorComparator
      {
        bool operator() (const XLPrecursor& a, const XLPrecursor& b) const
        {
          return a.precursor_mass < b.precursor_mass;
        }
        bool operator() (const XLPrecursor& a, const double& b) const
        {
          return a.precursor_mass < b;
        }
        bool operator() (const double& a, const XLPrecursor& b) const
        {
          return a < b.precursor_mass;
        }
      };

      /**
       * @brief The PeptidePosition enum

          Used to record the positions of peptides in proteins after in silico digestion determine whether protein terminal modifications are possible on a peptide.
       */
      enum PeptidePosition
      {
        INTERNAL = 0,
        C_TERM = 1,
        N_TERM = 2
      };

      /**
       * @brief The AASeqWithMass struct represents a normal peptide with its precomputed mass.

          This struct stores information about a peptide as an AASequence and a PeptidePosition.
          It is used to enumerate pairs of peptides in OpenPepXL.
          Since the mass of every peptide is used many times, it is precomputed once and also stored in this struct.
          A vector of these structs is used to represent the digested protein database in OpenPepXL.
          An instance of this struct is created only once for each peptide in the digested database, so it does not contribute to memory usage
          as much as XLPrecursor does.
       */
      struct AASeqWithMass
      {
        double peptide_mass = 0;
        AASequence peptide_seq;
        PeptidePosition position = PeptidePosition::INTERNAL;
        String unmodified_seq;
      };

      /**
       * @brief The AASeqWithMassComparator is a comparator for AASeqWithMass objects.

          This comparator allows to sort AASeqWithMass objects by the precomputed peptide mass and search for AASeqWithMass objects within a double type mass range.
       */
      struct AASeqWithMassComparator
      {
        bool operator() (const AASeqWithMass& a, const AASeqWithMass& b) const
        {
          return a.peptide_mass < b.peptide_mass;
        }
        bool operator() (const AASeqWithMass& a, const double& b) const
        {
          return a.peptide_mass < b;
        }
        bool operator() (const double& a, const AASeqWithMass& b) const
        {
          return a < b.peptide_mass;
        }
      };

      /**
       * @brief The PreprocessedPairSpectra struct represents the result of comparing a light and a heavy labeled spectra to each other.

          OpenPepXL can use labeled cross-linkers to denoise MS2 spectra. The PeakMaps contained in this struct represent the result of this
          denoising process for a whole mzML input file.
       */
      struct PreprocessedPairSpectra
      {

        MSExperiment spectra_linear_peaks; ///< merge spectrum of linear peaks (present in both spectra)
        MSExperiment spectra_xlink_peaks; ///< Xlink peaks in the light spectrum (linear peaks between spectra_light_different and spectra heavy_to_light)
        MSExperiment spectra_all_peaks;

        // pre-initialize so we can simply std::swap the spectra (no synchronization in multi-threading context needed as we get no reallocation of the PeakMaps).
        PreprocessedPairSpectra(Size size)
        {
          for (Size i = 0; i != size; ++i)
          {
            spectra_linear_peaks.addSpectrum(PeakSpectrum());
            spectra_xlink_peaks.addSpectrum(PeakSpectrum());
            spectra_all_peaks.addSpectrum(PeakSpectrum());
          }
        }
      };

  }; // class
} // namespace OpenMS
