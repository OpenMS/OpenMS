// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors: Raphael Förster $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>

namespace OpenMS
{

/**
  @brief Fragment-index-based peptide database search algorithm (experimental).

  Provides a self-contained search engine that matches MS/MS spectra against a protein
  database using an FI (Fragment Index). Typical usage:
  - Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
  - Call search() with an input mzML file and a FASTA database to populate identification
    outputs (ProteinIdentification and PeptideIdentificationList)
  - Intended for educational/prototyping use and to demonstrate FI-backed searching

  Notes:
  - Used by the PeptideDataBaseSearchFI TOPP tool
  - Experimental; interfaces and behavior may change
*/
/**
  @brief Fragment-index-based peptide database search algorithm (experimental).

  Provides a self-contained search engine that matches MS/MS spectra against a protein
  database using an FI (Fragment Index). Typical usage:
  - Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
  - Call search() with an input mzML file and a FASTA database to populate identification
    outputs (ProteinIdentification and PeptideIdentificationList)
  - Intended for educational/prototyping use and to demonstrate FI-backed searching

  Notes:
  - Used by the PeptideDataBaseSearchFI TOPP tool
  - Experimental; interfaces and behavior may change
*/
class OPENMS_DLLAPI PeptideSearchEngineFIAlgorithm :
  public DefaultParamHandler,
  public ProgressLogger
{
  public:
    PeptideSearchEngineFIAlgorithm(); 

    /// Exit codes
    enum class ExitCodes
    {
      EXECUTION_OK,
      INPUT_FILE_EMPTY,
      UNEXPECTED_RESULT,
      UNKNOWN_ERROR,
      ILLEGAL_PARAMETERS
    };

    /**
     * @brief Search spectra in an mzML file against a protein database using an FI-backed workflow.
     *
     * Populates protein and peptide identifications, including search meta data, PSM hits,
     * and search engine annotations. Parameters are taken from this instance (DefaultParamHandler).
     *
     * @param in_mzML Input path to the mzML file containing MS/MS spectra to search.
     * @param in_db   Input path to the protein sequence database in FASTA format.
     * @param prot_ids Output container receiving search meta data and protein-level information.
     * @param pep_ids  Output container receiving spectrum-level peptide identifications (PSMs).
     *
     * @return ExitCodes indicating success (EXECUTION_OK) or the encountered error condition.
     *
     * Side effects:
     *  - Updates ProgressLogger during processing.
     *  - Assigns identifiers and sets search engine name/version in prot_ids/pep_ids.
     *
     * Errors:
     *  - May signal invalid parameters via ILLEGAL_PARAMETERS exit code.
     *  - May propagate OpenMS exceptions (e.g., I/O or parse errors) from underlying components.
     */
    ExitCodes search(const String& in_mzML,
      const String& in_db,
      std::vector<ProteinIdentification>& prot_ids,
      PeptideIdentificationList& pep_ids) const;
  protected:
    void updateMembers_() override;

    /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
    struct AnnotatedHit_
    {
      AASequence sequence;
      /*
      StringView sequence;
      SignedSize peptide_mod_index; ///< enumeration index of the non-RNA peptide modification
      */
      double score = 0; ///< main score
      std::vector<PeptideHit::PeakAnnotation> fragment_annotations;      
      double prefix_fraction = 0; ///< fraction of annotated b-ions
      double suffix_fraction = 0; ///< fraction of annotated y-ions
      double mean_error = 0.0; ///< mean absolute fragment mass error

      static bool hasBetterScore(const AnnotatedHit_& a, const AnnotatedHit_& b)
      {
        if (a.score != b.score) return a.score > b.score;
        return a.sequence < b.sequence;
      }
    };

    /// @brief filter, deisotope, decharge spectra
    static void preprocessSpectra_(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm);

    /**
     * @brief Filter and annotate search results.
     *
     * Trims per-spectrum candidate hits to the top N and converts them into
     * PeptideIdentification objects, adding requested PSM annotations and
     * populating protein-level search metadata.
     *
     * @param exp Input MS experiment providing spectra/metadata for annotation.
     * @param annotated_hits Per-spectrum candidate hits (trimmed to @p top_hits in-place).
     * @param protein_ids Output container for protein-level identification and search metadata.
     * @param peptide_ids Output container for spectrum-level peptide identifications (PSMs).
     * @param top_hits Number of top-scoring hits to retain per spectrum (report_top_hits_).
     * @param modifications_fixed Fixed modifications (by name) used during the search.
     * @param modifications_variable Variable modifications (by name) used during the search.
     * @param peptide_missed_cleavages Allowed missed cleavages in digestion.
     * @param precursor_mass_tolerance Precursor mass tolerance value.
     * @param fragment_mass_tolerance Fragment mass tolerance value.
     * @param precursor_mass_tolerance_unit_ppm Precursor tolerance unit ("true"->ppm, "false"->Da).
     * @param fragment_mass_tolerance_unit_ppm Fragment tolerance unit ("true"->ppm, "false"->Da).
     * @param precursor_min_charge Minimum precursor charge considered.
     * @param precursor_max_charge Maximum precursor charge considered.
     * @param enzyme Digestion enzyme name.
     * @param database_name Database file name used for the search (stored in protein_ids).
     */
    void postProcessHits_(const PeakMap& exp,
      std::vector<std::vector<PeptideSearchEngineFIAlgorithm::AnnotatedHit_> >& annotated_hits,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      /* TODO check if ok if unused
      Size top_hits,
      const ModifiedPeptideGenerator::MapToResidueType& fixed_modifications,
      const ModifiedPeptideGenerator::MapToResidueType& variable_modifications,
      */
      Size max_variable_mods_per_peptide,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const String& precursor_mass_tolerance_unit_ppm,
      const String& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const String& enzyme,
      const String& database_name) const;

    double precursor_mass_tolerance_;
    String precursor_mass_tolerance_unit_;

    Size precursor_min_charge_;
    Size precursor_max_charge_;

    IntList precursor_isotopes_;

    double fragment_mass_tolerance_;

    String fragment_mass_tolerance_unit_;

    StringList modifications_fixed_;

    StringList modifications_variable_;

    Size modifications_max_variable_mods_per_peptide_;

    String enzyme_;

    bool decoys_;

    StringList annotate_psm_;

    Size peptide_min_size_;
    Size peptide_max_size_;
    Size peptide_missed_cleavages_;

    String peptide_motif_;

    Size report_top_hits_;
};

} // namespace