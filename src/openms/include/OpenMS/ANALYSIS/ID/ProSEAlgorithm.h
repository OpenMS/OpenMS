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

#include <OpenMS/ANALYSIS/ID/FragmentIndex.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <algorithm>   // std::min (used by inline computeModMatchTolerance_)
#include <vector>

namespace OpenMS
{

/**
  @brief Fragment-index-based peptide database search algorithm (experimental).

  Provides a self-contained search engine that matches MS/MS spectra against a protein
  database using an FI (Fragment Index). Typical usage:
  - Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
  - Call search() with an input spectrum file (mzML or Bruker .d) and a FASTA database to populate identification
    outputs (ProteinIdentification and PeptideIdentificationList)
  - Intended for educational/prototyping use and to demonstrate FI-backed searching

  Notes:
  - Used by the ProSE TOPP tool
  - Experimental; interfaces and behavior may change
*/
class OPENMS_DLLAPI ProSEAlgorithm :
  public DefaultParamHandler,
  public ProgressLogger
{
  public:
    ProSEAlgorithm(); 

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
     * @brief Comprehensive search result including modification analysis
     *
     * This structure contains all outputs from an open search including:
     * - Standard protein and peptide identifications
     * - Delta mass statistics table (histogram of mass shifts)
     * - PTM statistics table (mapped modifications with residue analysis)
     */
    struct SearchResult
    {
      ExitCodes exit_code = ExitCodes::EXECUTION_OK;
      std::vector<ProteinIdentification> protein_ids;
      PeptideIdentificationList peptide_ids;
      OpenSearchModificationAnalysis::OpenSearchAnalysisResult modification_analysis;
      bool is_open_search = false;
    };

    /**
     * @brief Multi-file search result bundle.
     *
     * Returned by the file-list searchWithModificationAnalysis() overloads. Holds
     * one SearchResult per input file (in @p per_file, in input order) and a
     * single @p aggregate result whose peptide_ids are the concatenation of all
     * per-file PSMs and whose modification_analysis is computed once on the
     * pooled set of PSMs.
     *
     * Special cases for @p aggregate:
     *  - When the input list contains exactly one file, @p aggregate is left
     *    almost-empty (only @c is_open_search and @c exit_code are set) — the
     *    single-file pooled aggregate would just duplicate @c per_file[0] and
     *    re-run modification analysis on the same PSMs. Callers should use
     *    @c per_file[0] for the result in this case.
     *  - When every per-file run failed, @p aggregate.exit_code is set to the
     *    first non-OK per-file exit code (so callers can inspect it without
     *    walking the @p per_file vector).
     *
     * The aggregate's @c protein_ids template is taken from the first
     * successful per-file result (search parameters are identical across files
     * by construction), with the primary MS run path overwritten to list every
     * input file.
     */
    struct MultiFileSearchResult
    {
      std::vector<SearchResult> per_file;
      SearchResult aggregate;
    };

    /**
     * @brief Prepared per-database state shared across multiple spectrum files.
     *
     * Holds the (decoy-augmented) protein database and the built FragmentIndex
     * so that searching N spectrum files against the same FASTA pays the index
     * build cost only once. Construct via prepareContext() and pass to the
     * context-taking search() overload.
     */
    struct SearchContext
    {
      std::vector<FASTAFile::FASTAEntry> db;
      FragmentIndex fragment_index;
    };

    /**
     * @brief Search spectra in a spectrum file (mzML or Bruker .d) against a protein database using an FI-backed workflow.
     *
     * Populates protein and peptide identifications, including search meta data, PSM hits,
     * and search engine annotations. Parameters are taken from this instance (DefaultParamHandler).
     *
     * @param[in] in_spectra Input path to the spectra file (mzML or Bruker .d) containing MS/MS spectra to search.
     * @param[in] in_db   Input path to the protein sequence database in FASTA format.
     * @param[out] prot_ids Output container receiving search meta data and protein-level information.
     * @param[out] pep_ids  Output container receiving spectrum-level peptide identifications (PSMs).
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
    ExitCodes search(const String& in_spectra,
      const String& in_db,
      std::vector<ProteinIdentification>& prot_ids,
      PeptideIdentificationList& pep_ids) const;

    /**
     * @brief Search with comprehensive results including modification analysis tables
     *
     * This method performs a peptide database search and additionally returns
     * structured modification analysis results for open search mode. This is the
     * recommended method for modification discovery workflows.
     *
     * The method automatically:
     * - Detects open search mode based on precursor tolerance
     * - Computes delta mass statistics
     * - Maps delta masses to known modifications
     * - Generates PTM statistics with residue localization
     * - Writes TSV output files if output_base_name is provided
     *
     * @param in_spectra Input path to the spectra file (mzML or Bruker .d) containing MS/MS spectra
     * @param in_db Input path to the protein sequence database in FASTA format
     * @param output_base_name Optional base name for output files (TSV tables)
     * @return SearchResult containing identifications and modification analysis
     *
     * Example usage:
     * @code
     * ProSEAlgorithm algo;
     * Param p = algo.getParameters();
     * p.setValue("precursor:mass_tolerance_lower", 500.0);  // Open search
     * p.setValue("precursor:mass_tolerance_upper", 500.0);
     * p.setValue("precursor:mass_tolerance_unit", "Da");
     * algo.setParameters(p);
     *
     * auto result = algo.searchWithModificationAnalysis("spectra.mzML", "database.fasta", "output");
     * if (result.exit_code == ExitCodes::EXECUTION_OK && result.is_open_search)
     * {
     *   // Access PTM statistics
     *   for (const auto& ptm : result.modification_analysis.ptm_stats.entries)
     *   {
     *     std::cout << ptm.name << ": " << ptm.count << " PSMs" << std::endl;
     *   }
     *
     *   // Access delta mass statistics
     *   for (const auto& dm : result.modification_analysis.delta_mass_stats.entries)
     *   {
     *     std::cout << dm.delta_mass << " Da: " << dm.count << " PSMs" << std::endl;
     *   }
     * }
     * @endcode
     */
    SearchResult searchWithModificationAnalysis(const String& in_spectra,
                                                const String& in_db,
                                                const String& output_base_name = "") const;

    /**
     * @brief In-memory search: search spectra against a protein database without file I/O.
     *
     * Same as the file-based search() but takes pre-loaded spectra and FASTA entries directly.
     * Spectra are preprocessed in-place (filtered, deisotoped, normalized).
     *
     * @param[in,out] spectra  MS/MS spectra to search (preprocessed in-place).
     * @param[in] fasta_db  Protein sequence database as FASTA entries.
     * @param[out] prot_ids  Output protein-level identifications.
     * @param[out] pep_ids   Output spectrum-level peptide identifications (PSMs).
     * @return ExitCodes indicating success or error.
     *
     * Internally this is a thin wrapper around prepareContext() + the
     * context-taking search() overload, so the FragmentIndex is rebuilt on every
     * call. For repeated searches against the same database, prefer calling
     * prepareContext() once and reusing the resulting SearchContext.
     */
    ExitCodes search(PeakMap& spectra,
                     const std::vector<FASTAFile::FASTAEntry>& fasta_db,
                     std::vector<ProteinIdentification>& prot_ids,
                     PeptideIdentificationList& pep_ids) const;

    /**
     * @brief Build a SearchContext (decoy-augmented database + FragmentIndex) for reuse.
     *
     * Performs the database preparation and FragmentIndex construction steps so
     * that subsequent calls to search(spectra, ctx, ...) can reuse the same
     * index across many spectrum files. If decoy generation is enabled
     * (parameter "decoys"), decoys are generated and shuffled into the
     * returned context's db member exactly once here.
     *
     * @param[in] fasta_db Protein sequence database as FASTA entries.
     * @return Prepared SearchContext containing the (possibly decoy-augmented)
     *         database and the built FragmentIndex.
     *
     * Thread-safety: the returned context's FragmentIndex is read-only during
     * subsequent search() calls; concurrent search() calls reading the same
     * SearchContext are safe (per FragmentIndex query thread-safety contract).
     * Do not call prepareContext() concurrently on the same algorithm instance.
     */
    SearchContext prepareContext(const std::vector<FASTAFile::FASTAEntry>& fasta_db) const;

    /**
     * @brief In-memory search using a pre-built SearchContext.
     *
     * Searches @p spectra against the database and FragmentIndex held in @p ctx.
     * The fragment index build cost (decoy generation, peptide/fragment
     * generation, sorting, bucketing) is paid by prepareContext() and is not
     * repeated here, making this overload the right choice when searching many
     * spectrum files against the same database.
     *
     * @param[in,out] spectra MS/MS spectra to search (preprocessed in-place).
     * @param[in,out] ctx Pre-built SearchContext from prepareContext(). Taken
     *            by non-const reference because the underlying FragmentIndex
     *            query API is non-const, even though the index content is not
     *            modified during the search; the @c db member is also handed
     *            non-const to the downstream PeptideIndexing step (which
     *            requires a non-const reference).
     * @param[out] prot_ids Output protein-level identifications.
     * @param[out] pep_ids  Output spectrum-level peptide identifications (PSMs).
     * @return ExitCodes indicating success or error.
     */
    ExitCodes search(PeakMap& spectra,
                     SearchContext& ctx,
                     std::vector<ProteinIdentification>& prot_ids,
                     PeptideIdentificationList& pep_ids) const;

    /**
     * @brief In-memory search with modification analysis: no file I/O required.
     *
     * Same as the file-based searchWithModificationAnalysis() but takes pre-loaded data.
     *
     * @param[in,out] spectra  MS/MS spectra (preprocessed in-place).
     * @param[in] fasta_db  Protein sequence database as FASTA entries.
     * @param[in] output_base_name  Optional base name for TSV output files.
     * @return SearchResult containing identifications and modification analysis.
     */
    SearchResult searchWithModificationAnalysis(PeakMap& spectra,
                                                const std::vector<FASTAFile::FASTAEntry>& fasta_db,
                                                const String& output_base_name = "") const;

    /**
     * @brief Multi-file search with modification analysis (in-memory FASTA).
     *
     * Builds a single SearchContext (decoy generation + FragmentIndex) from
     * @p fasta_db and reuses it across all input spectrum files. Each input
     * file produces its own SearchResult including a per-file modification
     * analysis (TSV written if a non-empty per-file base name is provided). An
     * additional aggregate SearchResult is computed by pooling all per-file
     * peptide identifications and running modification analysis once on the
     * pooled set.
     *
     * @param[in] in_spectra_files Spectrum file paths (mzML or Bruker .d).
     * @param[in] fasta_db Protein sequence database as FASTA entries.
     * @param[in] output_base_names Optional per-file base names for
     *            modification-analysis TSV outputs. Must be empty or have the
     *            same length as @p in_spectra_files. Empty entries skip TSV
     *            writing for that file.
     * @param[in] aggregate_base_name Optional base name for the aggregate
     *            modification-analysis TSV output. Empty disables aggregate
     *            TSV writing (the aggregate analysis is still computed).
     * @return MultiFileSearchResult bundling per-file results and aggregate.
     *
     * Errors:
     *  - Throws Exception::InvalidParameter if @p output_base_names is
     *    non-empty and its size differs from @p in_spectra_files.
     */
    MultiFileSearchResult searchWithModificationAnalysis(
        const std::vector<String>& in_spectra_files,
        const std::vector<FASTAFile::FASTAEntry>& fasta_db,
        const std::vector<String>& output_base_names = {},
        const String& aggregate_base_name = "") const;

    /**
     * @brief Multi-file search with modification analysis (FASTA file path).
     *
     * Convenience overload that loads the FASTA database from @p in_db and
     * delegates to the in-memory multi-file overload. The database file path
     * is recorded in each per-file ProteinIdentification's SearchParameters
     * (and on the aggregate result).
     *
     * @see searchWithModificationAnalysis(const std::vector<String>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<String>&, const String&) const
     */
    MultiFileSearchResult searchWithModificationAnalysis(
        const std::vector<String>& in_spectra_files,
        const String& in_db,
        const std::vector<String>& output_base_names = {},
        const String& aggregate_base_name = "") const;

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
      // Layout: doubles first, then floats, then int, then uint16_t — minimizes padding (40 bytes excluding AASequence)
      double score = 0; ///< main score
      double delta_mass = 0.0; ///< mass difference for open search (Da)
      float prefix_fraction = 0; ///< fraction of annotated prefix ions (a/b/c)
      float suffix_fraction = 0; ///< fraction of annotated suffix ions (x/y/z)
      float mean_error = 0.0f; ///< mean absolute fragment mass error
      int isotope_error = 0; ///< isotope offset used for this PSM
      uint16_t applied_charge = 0; ///< precursor charge used for this PSM
      uint16_t matched_prefix_ions = 0; ///< number of matched prefix ions (a/b/c)
      uint16_t matched_suffix_ions = 0; ///< number of matched suffix ions (x/y/z)

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
     * @param[in] exp Input MS experiment providing spectra/metadata for annotation.
     * @param[in,out] annotated_hits Per-spectrum candidate hits (trimmed to @p top_hits in-place).
     * @param[out] protein_ids Output container for protein-level identification and search metadata.
     * @param[out] peptide_ids Output container for spectrum-level peptide identifications (PSMs).
     * @param[in] top_hits Number of top-scoring hits to retain per spectrum (report_top_hits_).
     * @param[in] modifications_fixed Fixed modifications (by name) used during the search.
     * @param[in] modifications_variable Variable modifications (by name) used during the search.
     * @param[in] peptide_missed_cleavages Allowed missed cleavages in digestion.
     * @param[in] precursor_mass_tolerance Precursor mass tolerance value.
     * @param[in] fragment_mass_tolerance Fragment mass tolerance value.
     * @param[in] precursor_mass_tolerance_unit_ppm Precursor tolerance unit ("true"->ppm, "false"->Da).
     * @param[in] fragment_mass_tolerance_unit_ppm Fragment tolerance unit ("true"->ppm, "false"->Da).
     * @param[in] precursor_min_charge Minimum precursor charge considered.
     * @param[in] precursor_max_charge Maximum precursor charge considered.
     * @param[in] enzyme Digestion enzyme name.
     * @param[out] database_name Database file name used for the search (stored in protein_ids).
     */
    void postProcessHits_(const PeakMap& exp,
      std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      Size top_hits,
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

    /// Calibration overwrites these with the calibrated magnitudes for the duration of
    /// search(); pure runtime-state mutation that does not affect the logical const-ness
    /// of search(), matching the `mutable` pattern used by last_calibration_result_.
    mutable double precursor_mass_tolerance_lower_{20.0};   ///< positive magnitude
    mutable double precursor_mass_tolerance_upper_{20.0};   ///< positive magnitude
    String precursor_mass_tolerance_unit_{"ppm"};

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
    String decoy_prefix_;

    double fdr_psm_{0.0};
    double fdr_protein_{0.0};

    StringList annotate_psm_;

    Size peptide_min_size_;
    Size peptide_max_size_;
    Size peptide_missed_cleavages_;
    EnzymaticDigestion::Specificity peptide_enzyme_specificity_{EnzymaticDigestion::SPEC_FULL};

    String peptide_motif_;

    Size report_top_hits_;

    bool add_a_ions_{false};
    bool add_b_ions_{true};
    bool add_c_ions_{false};
    bool add_x_ions_{false};
    bool add_y_ions_{true};
    bool add_z_ions_{false};

    bool calibration_enabled_{false};
    double calibration_subset_ratio_{0.1};
    Size calibration_min_psms_{50};

    /**
     * @brief Result of a calibration pass.
     *
     * Holds the estimated precursor and fragment tolerances computed from
     * confident PSMs during the calibration pass. When @c success is false,
     * the tolerance values are undefined and should not be used.
     */
    struct CalibrationResult_
    {
      double precursor_shift{0};     ///< signed median of precursor errors (calibration bias)
      double precursor_spread{0};    ///< median(|e - shift|) + 3 * MAD(|e - shift|)
      double cal_lower{0};           ///< calibrated lower magnitude (valid iff !extreme_bias && success)
      double cal_upper{0};           ///< calibrated upper magnitude (valid iff !extreme_bias && success)
      double fragment_tolerance{0};  ///< estimated fragment tolerance (same unit as configured)
      double fragment_shift{0};      ///< reserved for future fragment m/z shift correction
      bool extreme_bias{false};      ///< |shift| >= spread — writeback skipped (test observability)
      bool success{false};           ///< true if enough PSMs were found for reliable estimation
    };

    /// Most recent calibration result (valid after any search that invoked runCalibrationPass_).
    /// Stored for test observability and diagnostics. Marked `mutable` because it is pure
    /// diagnostic/telemetry state that doesn't affect the logical const-ness of search().
    mutable CalibrationResult_ last_calibration_result_;

    /// Scalar tolerance passed to OpenSearchModificationAnalysis on the most recent
    /// search() call. Stored for test observability: because the calibration writeback
    /// restores the tolerance members on exit (to avoid per-file state leaks in the
    /// multi-file wrapper), tests that want to verify "the mod analyzer received the
    /// calibrated value, not the user-configured one" can't just read the members
    /// post-search — they need to see what was actually passed to the analyzer.
    /// Default -1.0 (sentinel: no search has run yet).
    mutable double last_mod_match_tolerance_used_{-1.0};

    /// Scalar tolerance passed to OpenSearchModificationAnalysis under asymmetric bounds.
    /// Uses the tighter of the two positive magnitudes — semantically correct for
    /// UniMod Δmass matching precision. OpenSearchModificationAnalysis internally clamps
    /// this at MAX_MOD_MAPPING_TOL_ = 0.02 Da; see spec §7 for rationale.
    ///
    /// Zero on one side is a legal one-sided window (e.g., [0, 500] Da = "search only
    /// positive mass shifts"). In that case std::min() would collapse to 0, passing
    /// a useless zero tolerance into the mod analyzer — masked in ppm mode by the
    /// internal clamp, but genuinely broken in Da mode. Fall back to the non-zero
    /// side so the mod-matching precision reflects the configured tolerance.
    double computeModMatchTolerance_() const
    {
      if (precursor_mass_tolerance_lower_ <= 0.0) return precursor_mass_tolerance_upper_;
      if (precursor_mass_tolerance_upper_ <= 0.0) return precursor_mass_tolerance_lower_;
      return std::min(precursor_mass_tolerance_lower_, precursor_mass_tolerance_upper_);
    }

    /**
     * @brief Run a fast calibration pass on a subset of spectra to estimate mass accuracy.
     *
     * Scores a TIC-ranked subset of spectra against the fragment index,
     * collects precursor and fragment mass errors from high-confidence PSMs,
     * and returns calibrated tolerances using median + 3*MAD estimation.
     *
     * @param[in] spectra  Preprocessed MS/MS spectra (subset is selected internally by TIC).
     * @param[in,out] fragment_index  Pre-built fragment index for candidate lookup.
     * @param[in] db  Protein database (for sequence reconstruction of candidates).
     * @return CalibrationResult_ with estimated tolerances, or success=false if insufficient PSMs.
     */
    CalibrationResult_ runCalibrationPass_(PeakMap& spectra,
                                           FragmentIndex& fragment_index,
                                           const std::vector<FASTAFile::FASTAEntry>& db) const;

    /// Helper: log the modification analysis summary (shared by in-memory and file-based paths)
    void logModificationAnalysisSummary_(const SearchResult& result,
                                         const String& output_base_name) const;

    /// Helper: log search summary statistics and per-run tolerance estimation
    void logSearchDiagnostics_(const PeakMap& spectra,
                               const std::vector<ProteinIdentification>& protein_ids,
                               const PeptideIdentificationList& peptide_ids) const;

    /// Helper function to determine if open search should be used based on tolerance
    bool isOpenSearchMode_() const
    {
      return FragmentIndex::isOpenSearchMode(precursor_mass_tolerance_lower_,
                                             precursor_mass_tolerance_upper_,
                                             precursor_mass_tolerance_unit_ == "ppm");
    }
};

} // namespace
