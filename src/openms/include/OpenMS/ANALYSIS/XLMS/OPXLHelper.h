// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <numeric>

namespace OpenMS
{
  /**
   * @brief The OPXLHelper class contains functions needed by OpenPepXL and OpenPepXLLF to reduce duplicated code
   */
  class OPENMS_DLLAPI OPXLHelper
  {
    public:

      /**
       * @brief A comparator for PeptideIdentifications that compares the scores in the first PeptideHit
       *
       */
      struct PeptideIDScoreComparator
      {
        bool operator() (const PeptideIdentification& a, const PeptideIdentification& b) const
        {
          if (!a.getHits().empty() && !b.getHits().empty())
          {
            return a.getHits()[0].getScore() < b.getHits()[0].getScore();
          }
          else
          {
            return false;
          }
        }
        bool operator() (const PeptideIdentification& a, const double& b) const
        {
          if (!a.getHits().empty())
          {
            return a.getHits()[0].getScore() < b;
          }
          else
          {
            return false;
          }
        }
        bool operator() (const double& a, const PeptideIdentification& b) const
        {
          if (!b.getHits().empty())
          {
            return a < b.getHits()[0].getScore();
          }
          else
          {
            return false;
          }
        }
      };

      /**
       * @brief Enumerates precursor masses for all candidates in an XL-MS search

          Assumes the list of peptides and the list of spectrum precursor masses are sorted by mass in ascending order,
          and the list of mono-link masses is sorted in descending order.

       * @param peptides The peptides with precomputed masses from the digestDatabase function
       * @param cross_link_mass_light Mass of the cross-linker, only the light one if a labeled linker is used
       * @param cross_link_mass_mono_link A list of possible masses for the cross-link, if it is attached to a peptide on one side
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param spectrum_precursors A vector of all MS2 precursor masses of the searched spectra. Used to filter out candidates.
       * @param precursor_correction_positions A vector of the position of the used precursor correction
       * @param precursor_mass_tolerance The precursor mass tolerance
       * @param precursor_mass_tolerance_unit_ppm The unit of the precursor mass tolerance ("Da" or "ppm")
       * @return A vector of XLPrecursors containing all possible candidate cross-links
       */
      static std::vector<OPXLDataStructs::XLPrecursor> enumerateCrossLinksAndMasses(const std::vector<OPXLDataStructs::AASeqWithMass>&  peptides, double cross_link_mass_light, const DoubleList& cross_link_mass_mono_link, const StringList& cross_link_residue1, const StringList& cross_link_residue2, const std::vector< double >& spectrum_precursors, std::vector< int >& precursor_correction_positions, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm);

      /**
       * @brief Digests a database with the given EnzymaticDigestion settings and precomputes masses for all peptides

          Also keeps track of the peptides at protein terminals and builds peptide candidates with all possible modification patterns
          according to the parameters.

       * @param fasta_db The protein database
       * @param digestor The object containing the digestion settings, e.g. the enzyme
       * @param min_peptide_length The minimal peptide length for the digestion
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param fixed_modifications The list of fixed modifications
       * @param variable_modifications The list of variable modifications
       * @param max_variable_mods_per_peptide The maximal number of variable modifications per peptide
       * @return A vector of AASeqWithMass containing the peptides, their masses and information about terminal peptides
       */
      static std::vector<OPXLDataStructs::AASeqWithMass> digestDatabase(std::vector<FASTAFile::FASTAEntry> fasta_db,
        const EnzymaticDigestion& digestor, Size min_peptide_length, const StringList& cross_link_residue1, const StringList& cross_link_residue2,
        const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
        const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
        Size max_variable_mods_per_peptide);

      /**
       * @brief Builds specific cross-link candidates with all possible combinations of linked positions from peptide pairs. Used to build candidates for the precursor mass window of a single MS2 spectrum.
       * @param candidates The XLPrecursors containing indices of two peptides
       * @param precursor_corrections Same size as @p candidates
       * @param precursor_correction_positions Same size as @p candidates; A vector of the position of the used precursor correction
       * @param peptide_masses The digested peptide database, that the indices in the XLPrecursors refer to
       * @param cross_link_residue1 A list of residues, to which the first side of the linker can react
       * @param cross_link_residue2 A list of residues, to which the second side of the linker can react
       * @param cross_link_mass mass of the cross-linker, only the light one if a labeled linker is used
       * @param cross_link_mass_mono_link A list of possible masses for the cross-link, if it is attached to a peptide on one side
       * @param spectrum_precursor_vector  The precursor masses of the experimental spectrum (used to filter out certain candidates, e.g. mono- and loop-links have a different mass)
       * @param allowed_error_vector Allowed mass error in Da for the respective precursor in @p spectrum_precursor_vector.
       * @param cross_link_name The name of the cross-linker
       * @return A vector of ProteinProteinCrossLink candidates containing all necessary information to generate theoretical spectra
       */
      static std::vector <OPXLDataStructs::ProteinProteinCrossLink> buildCandidates(const std::vector< OPXLDataStructs::XLPrecursor > & candidates,
                                                                                    const std::vector< int > & precursor_corrections,
                                                                                    const std::vector< int > & precursor_correction_positions,
                                                                                    const std::vector<OPXLDataStructs::AASeqWithMass> & peptide_masses,
                                                                                    const StringList & cross_link_residue1,
                                                                                    const StringList & cross_link_residue2,
                                                                                    double cross_link_mass,
                                                                                    const DoubleList & cross_link_mass_mono_link,
                                                                                    const std::vector< double >& spectrum_precursor_vector,
                                                                                    const std::vector< double >& allowed_error_vector,
                                                                                    const String& cross_link_name);

      /**
       * @brief Fills up the given FragmentAnnotation vector with annotations from a theoretical spectrum

          This function takes an alignment of a theoretical spectrum with meta information and an experimental spectrum
          and builds annotations taking the MZ and intensity values from the experimental spectrum and the ion names and charges from the theoretical spectrum
          to annotate matched experimental peaks.

       * @param frag_annotations The vector to fill. Does not have to be empty, as annotations from several alignments can just be added on to the same vector.
       * @param matching The alignment between the two spectra
       * @param theoretical_spectrum The theoretical spectrum with meta information
       * @param experiment_spectrum The experimental spectrum
       */
      static void buildFragmentAnnotations(std::vector<PeptideHit::PeakAnnotation> & frag_annotations, const std::vector< std::pair< Size, Size > > & matching, const PeakSpectrum & theoretical_spectrum, const PeakSpectrum & experiment_spectrum);

      /**
       * @brief Builds PeptideIdentifications and PeptideHits
       * @param peptide_ids The vector of PeptideIdentifications for the whole experiment. The created PepIds will be pushed on this one.
       * @param top_csms_spectrum All CrossLinkSpectrumMatches from the current spectrum to be written out
       * @param all_top_csms A vector of all CrossLinkSpectrumMatches of the experiment, that is also extended in this function
       * @param all_top_csms_current_index The index of the current spectrum in all_top_csms (some spectra have no matches, so this is not equal to the spectrum index)
       * @param spectra The searched spectra as a PeakMap
       * @param scan_index The index of the current spectrum
       * @param scan_index_heavy The index of the heavy spectrum in a spectrum pair with labeled linkers
       */
      static void buildPeptideIDs(std::vector<PeptideIdentification> & peptide_ids, const std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > & top_csms_spectrum, std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > > & all_top_csms, Size all_top_csms_current_index, const PeakMap & spectra, Size scan_index, Size scan_index_heavy);

      /**
       * @brief adds MetaValues for cross-link positions to PeptideHits
       * @param peptide_ids The vector of peptide_ids containing XL-MS search results with alpha and beta PeptideHits, after mapping of peptides to proteins
       */
      static void addProteinPositionMetaValues(std::vector< PeptideIdentification > & peptide_ids);

      /**
       * @brief adds xl_target_decoy MetaValue that combines alpha and beta target_decoy info
       * @param peptide_ids The vector of peptide_ids containing XL-MS search results with alpha and beta PeptideHits, after mapping of peptides to proteins
       */
      static void addXLTargetDecoyMV(std::vector< PeptideIdentification > & peptide_ids);

      /**
       * @brief adds accessions_beta MetaValue to alpha peptides for TOPPView visualization and CSV table output
       * @param peptide_ids The vector of peptide_ids containing XL-MS search results with alpha and beta PeptideHits, after mapping of peptides to proteins
       */
      static void addBetaAccessions(std::vector< PeptideIdentification > & peptide_ids);

      /**
       * @brief removes beta peptides from cross-link IDs, since all info is already contained in the alpha peptide hits
       * @param peptide_ids The vector of peptide_ids containing XL-MS search results with alpha and beta PeptideHits
       */
      static void removeBetaPeptideHits(std::vector< PeptideIdentification > & peptide_ids);

      /**
       * @brief Adds the list of features that percolator should use for OpenPepXL
       * @param prot_id Receives updated search parameters ('feature_extractor' and 'extra_features') as meta value for OpenPepXL
       */
      static void addPercolatorFeatureList(ProteinIdentification& prot_id);

      /**
       * @brief sorts PeptideHits for each PeptideIdentification by score and adds the delta score as a MetaValue
       * @param peptide_ids The vector of peptide_ids containing XL-MS search results without beta PeptideHits
       */
      static void computeDeltaScores(std::vector< PeptideIdentification >& peptide_ids);

      /**
       * @brief combines all hits to spectrum pairs with the same light spectrum into one ranked list
       *
       * This function is a post-processing step for OpenPepXL with labeled linkers.
       * This function collects PeptideIdentifications from all spectrum pairs with the same light spectrum,
       * then resorts them by the score, makes them unique in case of equal candidates and reduces their number down to the chosen number of reported top hits.
       *
       * @param peptide_ids PeptideIdentifications from a Cross-Linking MS search with labeled linkers
       * @param number_top_hits The chosen number of reported top hits
       */
      static std::vector< PeptideIdentification > combineTopRanksFromPairs(std::vector< PeptideIdentification > & peptide_ids, Size number_top_hits);

      /**
       * @brief Searches for cross-link candidates for a MS/MS spectrum

          This function uses enumerateCrossLinksAndMasses and buildCandidates to search for peptide pairs fitting to the given precursor mass_light
          and all considered precursor corrections.

       * @param precursor_correction_steps An IntList of integers as indices of isotopic peaks around the experimental precursor
       * @param precursor_mass The decharged precursor mass
       * @param precursor_mass_tolerance The precursor tolerance
       * @param precursor_mass_tolerance_unit_ppm The unit of the precursor tolerance. "ppm" if true, "Da" if false
       * @param filtered_peptide_masses A vector of AASeqWithMass containing the sorted (ascending) peptide database with precomputed peptide masses
       * @param cross_link_mass The mass of the cross-linker (light mass, if labeled)
       * @param cross_link_mass_mono_link A list of possible mono-link masses
       * @param cross_link_residue1 A list of one-letter-code residues, that the first side of the cross-linker can attach to
       * @param cross_link_residue2 A list of one-letter-code residues, that the second side of the cross-linker can attach to
       * @param cross_link_name The name of the cross-linker, e.g. "DSS" or "BS3"
       * @param use_sequence_tags Whether to use sequence tags to filter out candidates
       * @param tags The list of sequence tags that are used to filter candidate sequences. Only applied if use_sequence_tags = true
       */
      static std::vector <OPXLDataStructs::ProteinProteinCrossLink> collectPrecursorCandidates(const IntList& precursor_correction_steps,
                                                                                                double precursor_mass,
                                                                                                double precursor_mass_tolerance,
                                                                                                bool precursor_mass_tolerance_unit_ppm,
                                                                                                const std::vector<OPXLDataStructs::AASeqWithMass>& filtered_peptide_masses,
                                                                                                double cross_link_mass,
                                                                                                const DoubleList& cross_link_mass_mono_link,
                                                                                                const StringList& cross_link_residue1,
                                                                                                const StringList& cross_link_residue2,
                                                                                                String cross_link_name,
                                                                                                bool use_sequence_tags = false,
                                                                                                const std::vector<std::string>& tags = std::vector<std::string>());

      /**
       * @brief Computes the mass error of a precursor mass to a hit

       * @param csm The cross-link spectrum match containing the hit
       * @param precursor_mz The precursor mz of the MS/MS spectrum
       * @param precursor_charge The charge of the precursor
       */
      static double computePrecursorError(const OPXLDataStructs::CrossLinkSpectrumMatch& csm, double precursor_mz, int precursor_charge);

      /**
       * @brief Computes the mean of alpha, beta, xlinks-alpha and xlinks-beta respectively and store means in @p csm
       */
      static void isoPeakMeans(OPXLDataStructs::CrossLinkSpectrumMatch& csm, 
                               const DataArrays::IntegerDataArray& num_iso_peaks_array,
                               const std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha,
                               const std::vector< std::pair< Size, Size > >& matched_spec_linear_beta,
                               const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha,
                               const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta);

      /**
       * @brief Filters the list of candidates for cases that include at least one of the tags in at least one of the two sequences

       * @param candidates The list of XLPrecursors as enumerated by e.g. enumerateCrossLinksAndMasses
       * @param precursor_correction_positions 
       * @param tags The list of tags for the current spectrum produced by the Tagger
       */
      static void filterPrecursorsByTags(std::vector <OPXLDataStructs::XLPrecursor>& candidates, std::vector< int >& precursor_correction_positions, const std::vector<std::string>& tags);
  };
}
