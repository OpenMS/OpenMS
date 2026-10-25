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
#include <OpenMS/ANALYSIS/ID/FragmentIonLikelihoodModel.h>
#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <algorithm>   // std::min (used by inline computeModMatchTolerance_)
#include <functional>  // std::hash (ProteinMapping_)
#include <iosfwd>      // std::ostream (renderRunSummary / renderModificationSummary)
#include <map>
#include <memory>      // std::shared_ptr (released_hits_)
#include <string>      // std::string (renderRunSummaryJson return / manifest)
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>     // std::pair (renderRunSummaryJson manifest)
#include <vector>

namespace OpenMS
{

class TheoreticalSpectrumGenerator;

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
     * @brief Per-run identification statistics for the end-of-search report.
     *
     * Populated by collectRunStatistics_() once a single spectrum file has been
     * searched (post-FDR), plus a few fields captured at well-defined points
     * during search() (target/decoy counts pre-FDR, achieved q-value, timing).
     * All counts refer to one input file. Cross-file/shared facts (database,
     * fragment index, configuration) live in SharedSearchStats instead.
     */
    struct RunStatistics
    {
      std::string input_file;                  ///< spectrum file this run searched (basename or path)
      Size ms2_spectra = 0;                    ///< number of MS2 spectra in the input
      Size matched_spectra = 0;                ///< spectra with >=1 retained PSM in the final IDs (after FDR, if applied)
      Size target_psms = 0;                    ///< target PSMs in the final IDs (after FDR, if applied)
      Size decoy_psms = 0;                     ///< decoy PSMs in the final IDs (after FDR, if applied)
      bool fdr_applied = false;                ///< true if PSM-level FDR filtering ran
      double achieved_psm_fdr = -1.0;          ///< max retained q-value after FDR (<0 = n/a)
      Size unique_peptides = 0;                ///< distinct peptide sequences among top hits
      Size unique_proteins = 0;                ///< distinct protein accessions among top hits
      std::map<Int, Size> charge_histogram;    ///< precursor charge -> PSM count
      std::map<Size, Size> missed_cleavage_histogram; ///< missed cleavages -> PSM count
      bool score_stats_valid = false;          ///< true if hyperscore_* below are meaningful
      double hyperscore_min = 0.0;
      double hyperscore_median = 0.0;
      double hyperscore_max = 0.0;
      bool prec_tol_valid = false;             ///< true if precursor-error estimate present
      double prec_err_median = 0.0, prec_err_mad = 0.0, prec_err_recommended = 0.0;
      bool frag_tol_valid = false;             ///< true if fragment-error estimate present
      double frag_err_median = 0.0, frag_err_mad = 0.0, frag_err_recommended = 0.0;
      double seconds_search = 0.0;             ///< scoring + post-processing wall time
      double seconds_calibration = 0.0;        ///< calibration pass wall time (0 if disabled)
      double seconds_fdr = 0.0;                ///< FDR filtering wall time (0 if not applied)
    };

    /**
     * @brief Configuration, database and fragment-index facts shared across all
     * input files of one ProSE invocation.
     *
     * Computed once (the fragment index is built once and reused), so these
     * costs/counts must NOT be summed per file. Populated by the multi-file
     * searchWithModificationAnalysis() overloads.
     */
    struct SharedSearchStats
    {
      std::string database_file;               ///< FASTA path (empty for in-memory db)
      std::string enzyme;
      double precursor_tol_lower = 0.0, precursor_tol_upper = 0.0;
      std::string precursor_tol_unit;
      double fragment_tol = 0.0;
      std::string fragment_tol_unit;
      Int min_charge = 0, max_charge = 0;
      Size missed_cleavages = 0;
      std::vector<std::string> fixed_mods, variable_mods, ion_series;
      bool open_search = false;
      bool calibration_enabled = false;
      bool peptide_mass_mode = false;
      bool chunked = false;
      std::string decoy_mode;                  ///< "generated" | "external" | "none (target-only)"
      double psm_fdr_threshold = 0.0, protein_fdr_threshold = 0.0;
      Size db_target_proteins = 0;             ///< target entries in the searched (augmented) db
      Size db_decoy_proteins = 0;              ///< decoy entries in the searched (augmented) db
      Size indexed_peptides = 0;               ///< peptides in the fragment index (summed over chunks)
      Size indexed_fragments = 0;              ///< theoretical fragments in the index (summed over chunks)
      double seconds_index_build = 0.0;        ///< decoy generation + fragment index build wall time
      double seconds_total = 0.0;              ///< whole-search wall time (set by the caller)
    };

    /**
     * @brief Outcome of the Percolator rescoring of one search run (see rescoreWithPercolator() and the parameter 'rescore').
     *
     * With 'rescore' off nothing is attempted: @c applied is false and @c skipped_reason empty.
     */
    struct PercolatorRescoring
    {
      /// True if the PSMs were rescored
      bool applied = false;
      /// Why the PSMs were not rescored, as a sentence to show to the user (empty if @c applied or not attempted)
      std::string skipped_reason;
      /// A caveat of an applied rescoring, as a sentence to show to the user (empty if none): Percolator's
      /// sanity check replaced the learned model by the best single feature (Percolator::usedDefaultDirection())
      std::string warning;
      /// With 'rescore_reuse:mode' other than per_file: the models that scored the PSMs, as a phrase to show to the user (empty: trained
      /// on the run itself)
      std::string model;
    };

    /**
     * @brief Options of Percolator's training for the rescoring (the parameters 'rescoring:...'; the defaults are
     *        ProSE's training: all PSMs, 10 iterations)
     */
    struct PercolatorTraining
    {
      int subset_max_train = 0;  ///< train on at most this many PSMs (drawn by spectrum, then all are scored; 0 = all)
      std::string subset_mode = "out_of_fold"; ///< with a cap: Percolator's 'subset_mode' ("out_of_fold" | "out_of_fold_avg" | "in_sample" | "learning_curve")
      int curve_start = 100000;                ///< learning_curve: rows of the first sample
      double curve_min_gain = 0.002;           ///< learning_curve: relative gain below which a fold stops
      int iterations = 10;       ///< iterations of the cross-validation training
      double tron_epsilon = 1e-7;  ///< stopping tolerance of the SVM solver (TRON)
      bool warm_start = false;     ///< SVM trainings start from the previous iteration's weights
      bool single_cneg = false;    ///< one Cneg candidate instead of three
      std::string solver = "tron"; ///< Percolator's SVM solver (parameter 'rescore_solver')
    };

    /**
     * @brief rescore_reuse: whether spectrum @p spectrum (its index in the run) of file @p file is drawn for the
     *        training of the models at the sampling fraction @p fraction (rescore_reuse:max_train)
     *
     * A seeded hash-uniform draw per spectrum (all PSMs of a spectrum together), independent of its position: a capped
     * training sample is spread at random over the whole run, not its first spectra. The 5 % of the spectra whose draw
     * is highest are calibration spectra and never drawn for training.
     */
    static bool reuseTrainingDraw(UInt64 seed, Size file, Size spectrum, double fraction);

    /// rescore_reuse: whether spectrum @p spectrum of file @p file is drawn as a calibration spectrum (whose decoys set
    /// the scale of a model; never trained on) at the calibration fraction @p fraction: from the 5 % of the spectra
    /// whose draw is highest, so never a training spectrum (reuseTrainingDraw()) at any fraction
    static bool reuseCalibrationDraw(UInt64 seed, Size file, Size spectrum, double fraction);

    /**
     * @brief Comprehensive search result including modification analysis
     *
     * This structure contains all outputs from an open search including:
     * - Standard protein and peptide identifications
     * - Delta mass statistics table (histogram of mass shifts)
     * - PTM statistics table (mapped modifications with residue analysis)
     * - Per-run identification statistics for the end-of-search report
     */
    struct SearchResult
    {
      ExitCodes exit_code = ExitCodes::EXECUTION_OK;
      std::vector<ProteinIdentification> protein_ids;
      PeptideIdentificationList peptide_ids;
      OpenSearchModificationAnalysis::OpenSearchAnalysisResult modification_analysis;
      bool is_open_search = false;
      RunStatistics stats;
      /// Outcome of the Percolator rescoring of this run (parameter 'rescore')
      PercolatorRescoring percolator;
      PeptideIdentificationList pre_fdr_peptide_ids; ///< the PSMs before rescoring and FDR:PSM; only filled with setKeepPreFdrPsms(true)
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
     *  - When the caller passed @c build_pooled_aggregate = false, the pooled
     *    PSM copy is not retained: @p aggregate is left almost-empty (only
     *    @c is_open_search and @c exit_code), exactly as in the single-file
     *    case. In OPEN search the pooling still happens internally because the
     *    aggregate modification analysis needs it, so
     *    @c aggregate.modification_analysis is populated as usual and only
     *    @c aggregate.protein_ids / @c aggregate.peptide_ids are released
     *    afterwards. In closed search nothing is pooled at all.
     *  - When Percolator (parameter 'rescore') rescored only some of the files,
     *    their runs report different search engines (Percolator / ProSE) and
     *    are not merged: nothing is pooled (a warning is logged), so
     *    @p aggregate is left almost-empty as above, also in open search.
     *
     * The aggregate's @c protein_ids template is taken from the first
     * successful per-file result (search parameters are identical across files
     * by construction), with the primary MS run path overwritten to list every
     * input file. Its PSMs carry the meta value 'id_merge_index', the input
     * position of their file (for inputs with distinct paths), also when other
     * files failed.
     */
    struct MultiFileSearchResult
    {
      std::vector<SearchResult> per_file;
      SearchResult aggregate;

      /// Effective decoy marker resolved from the shared database, for a
      /// caller-side merged-PSM protein-FDR step (e.g. the ProSE TOPP tool's
      /// -out_merged path). Empty when the search was target-only (decoys=ignore).
      std::string decoy_string;
      /// Position of @c decoy_string (true = prefix, false = suffix).
      bool decoy_is_prefix = true;
      /// True when the searched databases contained decoys (FDR possible).
      bool have_decoys = false;

      /// Configuration / database / fragment-index facts shared across all input
      /// files (the index is built once and reused), for the end-of-search report.
      SharedSearchStats shared;
    };

    /// The pooled evaluation of a multi-file search (the ProSE tool's -out_merged).
    struct MergedResult
    {
      std::vector<ProteinIdentification> protein_ids; ///< one pooled run; empty if no file succeeded
      /// PSMs of all files; meta value 'id_merge_index' = input position of the file (inputs with distinct paths)
      PeptideIdentificationList peptide_ids;
      /// picked-protein FDR applied and its filter at FDR:protein (and, with remove_decoys, the decoys removed)
      bool protein_fdr_applied = false;
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

      /// When true, the context-taking search() overload will release
      /// `fragment_index` (via clear()) after scoreSpectraAgainstIndex_
      /// returns, reclaiming its heap footprint before the PeptideIndexing
      /// Aho-Corasick pass — which is otherwise the RSS high-water mark
      /// on large databases. Set by single-use callers (e.g. the
      /// internal search(spectra, fasta_db, ...) wrapper that creates
      /// a throw-away ctx). Default false preserves the reuse contract
      /// for callers that built ctx via prepareContext() and want to
      /// run multiple search() calls against it.
      bool release_fragment_index_after_scoring = false;

      /// Effective decoy marker carried by `db` (auto-detected for external
      /// decoys, or the configured prefix for internally generated ones).
      /// Empty when the search is target-only. Feeds PeptideIndexing and the
      /// protein-level FDR so they recognise the same decoys that were searched.
      std::string decoy_string;
      /// Position of `decoy_string` (true = prefix, false = suffix).
      bool decoy_is_prefix = true;
      /// True when `db` contains decoy entries (generated or external), i.e.
      /// target-decoy FDR is possible.
      bool have_decoys = false;
      /// True when `fragment_index` also holds c and z+1 ions for electron-activated
      /// spectra (ions:by_activation); only those spectra are matched against them. See
      /// prepareContext(). search() does not add them to a context: for such spectra and
      /// a context without these ions, it builds a temporary index for the call.
      bool electron_ions = false;
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
    ExitCodes search(const std::string& in_spectra,
      const std::string& in_db,
      std::vector<ProteinIdentification>& prot_ids,
      PeptideIdentificationList& pep_ids) const;

    /**
     * @brief Finalize protein-level FDR on a COMPLETE protein set (a single input file,
     * or a merged cross-file aggregate).
     *
     * Runs protein inference (BasicProteinInferenceAlgorithm), picked-protein FDR
     * (Savitski et al. 2015), threshold filtering at @p protein_fdr, then removes decoys
     * and repairs indistinguishable-protein / protein-group / peptide-evidence references
     * so the result stores as schema-valid, decoy-free idXML.
     *
     * Shared by the file-based single-file search() above and the ProSE TOPP tool's
     * single-file finalization, so the exact IDFilter sequence and ordering live in one
     * place (they previously existed as two copies that could drift).
     *
     * Precondition: @p protein_ids is non-empty and the set is statistically complete —
     * picked-protein FDR does not compose across runs. The caller gates on "decoys present"
     * and "FDR requested"; this helper additionally skips (with a warning) when no decoy
     * proteins survive inference, to avoid target-only q-values.
     *
     * @param[in,out] protein_ids Protein identifications (operates on protein_ids[0]).
     * @param[in,out] peptide_ids PSMs feeding inference; decoy PSMs are removed.
     * @param[in] decoy_string Accession marker identifying decoy proteins.
     * @param[in] decoy_is_prefix Whether @p decoy_string is a prefix (true) or suffix (false).
     * @param[in] protein_fdr Picked-protein q-value threshold (expected > 0).
     */
    static void applyCompleteSetProteinFDR(std::vector<ProteinIdentification>& protein_ids,
                                           PeptideIdentificationList& peptide_ids,
                                           const std::string& decoy_string,
                                           bool decoy_is_prefix,
                                           double protein_fdr);

    /**
     * @brief Annotate PSM q-values by target-decoy competition of the best hit of each spectrum, as FDR:PSM does.
     *
     * Uses FalseDiscoveryRate (decoy hits are kept and annotated) and replaces hit scores by q-values; no filtering.
     * With FDR:PSM_groups = 'scored_charges' and multiply charged fragments scored (scoring:fragment_charges), PSMs
     * whose best hits were scored with different numbers of fragment charges compete separately, because HyperScore
     * grows with the number of theoretical ions. With one group, or a group without target or decoy hits, all PSMs
     * compete together.
     *
     * @param[in,out] peptide_ids PSMs with target_decoy annotations (PeptideIndexing).
     */
    void annotatePsmQValues(PeptideIdentificationList& peptide_ids) const;

    /**
     * @brief Rescore the PSMs of one search run in-process with Percolator, as the parameter 'rescore' does in every search.
     *
     * Trains Percolator (Percolator::rescorePSMs) on the standard PIN feature set
     * (PercolatorInfile::getStandardFeatureSet) plus the search engine features of the run (the
     * search parameters' 'extra_features', or 'score' without them), with train_best_positive,
     * target-decoy competition (post_processing_tdc), all PSMs for training (subset_max_train 0)
     * and nonparametric PEPs, as the percolator executable. Percolator q-values become the main
     * score (score type "q-value"; the SVM score and PEP are kept as meta values), only the best
     * hit per spectrum is kept, and the runs are marked as post-processed by Percolator (search
     * engine "Percolator", version "3.08-vendored", meta value 'percolator' and the
     * "Percolator:*" search parameter meta values PercolatorAdapter writes).
     *
     * The searches of this class call it for every run when the parameter 'rescore' is on (the
     * default), in place of their PSM-level FDR, which then follows on the rescored PSMs
     * (applyDeferredPsmFDR()). Call it directly only to rescore PSMs of a search with 'rescore'
     * off, whose FDR filtering must then be off too (FDR:PSM and FDR:protein 0), so that
     * Percolator gets all target and decoy PSMs.
     *
     * The PSMs are not rescored, and stay unchanged, with fewer than 100 PSMs, without decoy PSMs
     * or if Percolator fails; @c skipped_reason of the result then says why. Nothing is logged
     * about the skip, so that the caller can report it with its context (e.g. the input file).
     *
     * @param[in,out] protein_ids The runs of the search; the first one gives the enzyme and the
     *                            search engine features. Marked as post-processed by Percolator if
     *                            the PSMs were rescored.
     * @param[in,out] peptide_ids Target and decoy PSMs of the run (with 'target_decoy' annotations
     *                            and peptide evidences, as from PeptideIndexing).
     * @param[in] threads Threads Percolator may use (values above 9 use 9, values <= 0 use 3); never
     *                    more than OpenMP allows. Cross-validation trains 9 models at a time (3 folds x
     *                    3 Cneg values), so more than 9 do not help. The results do not depend on it.
     *
     * @return Whether the PSMs were rescored, and if not, why.
     */
    static PercolatorRescoring rescoreWithPercolator(std::vector<ProteinIdentification>& protein_ids,
                                                     PeptideIdentificationList& peptide_ids,
                                                     int threads);

    /**
     * @brief Apply the PSM-level FDR that a search with Percolator rescoring deferred.
     *
     * For PSMs that rescoreWithPercolator() rescored (@p rescored), the scores are Percolator
     * q-values already. Otherwise they get q-values by annotatePsmQValues() first (FDR:PSM_groups
     * of this instance), provided the run has decoy proteins (protein hits with 'target_decoy'
     * "decoy"). The PSMs are then filtered at @p psm_fdr (IDFilter::filterHitsByScore); target and
     * decoy PSMs are filtered alike, decoys are not removed (that is the protein FDR's
     * finalization, see applyCompleteSetProteinFDR()).
     *
     * The searches of this class apply it after rescoreWithPercolator() when the parameter
     * 'rescore' is on.
     *
     * @param[in] protein_ids The runs of the search; the first one tells whether there are decoys.
     * @param[in,out] peptide_ids PSMs of the run, rescored by rescoreWithPercolator() or with their
     *                            search engine scores.
     * @param[in] rescored Whether rescoreWithPercolator() rescored @p peptide_ids (its @c applied).
     * @param[in] psm_fdr PSM q-value threshold.
     *
     * @return False if nothing was filtered: @p psm_fdr <= 0, no runs or no PSMs, or PSMs that were
     *         not rescored without decoy proteins to estimate the FDR with.
     */
    bool applyDeferredPsmFDR(const std::vector<ProteinIdentification>& protein_ids,
                             PeptideIdentificationList& peptide_ids,
                             bool rescored,
                             double psm_fdr) const;

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
     * - Applies FDR:protein as search() does, after the modification analysis
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
    SearchResult searchWithModificationAnalysis(const std::string& in_spectra,
                                                const std::string& in_db,
                                                const std::string& output_base_name = "") const;

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
     * Thread-safety: concurrent search() calls reading the same SearchContext are
     * safe (per FragmentIndex query thread-safety contract) as long as
     * calibration:enabled is off. A successful calibration sets the calibrated
     * tolerances as query parameters of the context's FragmentIndex until the call
     * returns, so calibrated searches must not share a context concurrently.
     * Do not call prepareContext() concurrently on the same algorithm instance.
     */
    SearchContext prepareContext(const std::vector<FASTAFile::FASTAEntry>& fasta_db) const;

    /**
     * @brief Build a SearchContext whose FragmentIndex also holds c and z+1 ions.
     *
     * As prepareContext(fasta_db). With @p electron_ions, the index also holds c and z+1
     * ions, which ions:by_activation scores for electron-activated spectra (ETD, ECD, EThcD,
     * ETciD). Only those spectra are matched against them, so other spectra get the same
     * candidates as from prepareContext(fasta_db). Use it when the spectra to search contain
     * such spectra: search() otherwise builds a temporary index with these ions for each call.
     *
     * @param[in] fasta_db Protein sequence database as FASTA entries.
     * @param[in] electron_ions Also index c and z+1 ions.
     * @return Prepared SearchContext.
     */
    SearchContext prepareContext(const std::vector<FASTAFile::FASTAEntry>& fasta_db, bool electron_ions) const;

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
                                                const std::string& output_base_name = "") const;

    /**
     * @brief Multi-file search with modification analysis (in-memory FASTA).
     *
     * Builds a single SearchContext (decoy generation + FragmentIndex) from
     * @p fasta_db and reuses it across all input spectrum files. Each input
     * file produces its own SearchResult including a per-file modification
     * analysis (TSV written if a non-empty per-file base name is provided). An
     * additional aggregate SearchResult is computed by pooling all per-file
     * peptide identifications and running modification analysis once on the
     * pooled set (see @p build_pooled_aggregate to skip retaining that second
     * copy of every PSM).
     *
     * Progress: the search of each input file is a progress step of its own,
     * labelled "Searching &lt;file name&gt; (i of n)..." (i from 1, in input order).
     * Consumers identify the per-file step by this label. Every start of a
     * per-file step is balanced by its end, also if the search of the file
     * fails or throws; the steps nested in it are not exception-safe (one that
     * an exception leaves stays open). The steps of the search of the file
     * (index build, filtering, scoring, ...) nest one level below it (depth 0
     * for the per-file step, depth 1 inside), but this nesting holds only with
     * one OpenMP thread: with more, the background reading of the next mzML
     * file shares the static step depth of the ProgressLogger, and the depths
     * reported are not reliable. Chunk-major searches (database:chunk_size > 0
     * with a database of more proteins, decoys included, than the chunk size)
     * emit no per-file step.
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
     * @param[in] build_pooled_aggregate Whether the pooled PSM set is retained
     *            in @p aggregate. Pooling copies every per-file PSM a second
     *            time, so callers that never read @c aggregate.peptide_ids
     *            should pass false. With false, @c aggregate.protein_ids and
     *            @c aggregate.peptide_ids come back empty (only
     *            @c is_open_search and @c exit_code are set); the aggregate
     *            modification analysis is unaffected — in open search the
     *            pooled set is still built internally, analyzed, and then
     *            released. Default true (unchanged behavior).
     * @return MultiFileSearchResult bundling per-file results and aggregate.
     *
     * Errors:
     *  - Throws Exception::InvalidParameter if @p output_base_names is
     *    non-empty and its size differs from @p in_spectra_files.
     */
    MultiFileSearchResult searchWithModificationAnalysis(
        const std::vector<std::string>& in_spectra_files,
        const std::vector<FASTAFile::FASTAEntry>& fasta_db,
        const std::vector<std::string>& output_base_names = {},
        const std::string& aggregate_base_name = "",
        bool build_pooled_aggregate = true) const;

    /**
     * @brief Multi-file search with modification analysis (FASTA file path).
     *
     * Convenience overload that loads the FASTA database from @p in_db and
     * delegates to the in-memory multi-file overload. The database file path
     * is recorded in each per-file ProteinIdentification's SearchParameters
     * (and on the aggregate result).
     *
     * @p build_pooled_aggregate is forwarded unchanged; see the in-memory
     * overload for its exact effect on @c aggregate.
     *
     * @see searchWithModificationAnalysis(const std::vector<std::string>&, const std::vector<FASTAFile::FASTAEntry>&, const std::vector<std::string>&, const std::string&, bool) const
     */
    MultiFileSearchResult searchWithModificationAnalysis(
        const std::vector<std::string>& in_spectra_files,
        const std::string& in_db,
        const std::vector<std::string>& output_base_names = {},
        const std::string& aggregate_base_name = "",
        bool build_pooled_aggregate = true) const;

    /**
     * @brief Pools result (moving its PSMs out), infers proteins across files, applies FDR:protein of this instance
     * with picked-protein FDR if decoy proteins were identified and then, by default, removes the decoys. Throws
     * Exception::InvalidValue if Percolator rescored only some files.
     *
     * The evaluation of the ProSE tool's -out_merged, for the result of a multi-file
     * searchWithModificationAnalysis() with the parameters of this instance:
     *  - Pooling: the pooled aggregate of the search (@c build_pooled_aggregate) is taken as it is. Without one (e.g.
     *    @c build_pooled_aggregate false, or a search of a single file, which builds none), the runs of the files that
     *    were searched (exit code EXECUTION_OK) are merged with IDMergerAlgorithm: the PSMs of the files are moved,
     *    their runs are copied and stay with the files. 'id_merge_index' is the input position of a PSM's file (its
     *    index in @p in_spectra_files, for inputs with distinct paths), also when other files failed. The pooled run
     *    lists @p in_spectra_files as its primary MS run paths and gets the current date.
     *  - Cross-file protein inference (BasicProteinInferenceAlgorithm).
     *  - With FDR:protein > 0 and decoys in the searched database: picked-protein FDR (Savitski et al. 2015) and the
     *    filter at FDR:protein, provided the inference identified decoy proteins (otherwise a warning, and nothing is
     *    filtered). The finalization then removes the decoys and the references to removed proteins, as
     *    applyCompleteSetProteinFDR() does. With @p remove_decoys false the picked-protein FDR keeps the decoy proteins
     *    in the run, with q-values (FalseDiscoveryRate's add_decoy_proteins; the q-values of the targets are the
     *    same), and the filter applies to them as to the targets; the decoy PSMs and the decoy proteins that pass are
     *    kept, and only the references to the proteins the filter removed go (peptide evidences, indistinguishable
     *    proteins and protein groups), so that the result can still be stored, e.g. as idXML. Without the protein
     *    FDR the result keeps the targets and decoys.
     *
     * @param[in,out] result The search result. Its pooled aggregate's identifications are moved into the returned
     *                       result; without an aggregate, the PSMs of the files that were searched are, while the runs
     *                       of the files (@c per_file[i].protein_ids: search parameters, meta values such as
     *                       PeptDeepRescoring::NCE_META_VALUE, primary MS run paths) stay as they are. The statistics,
     *                       exit codes, Percolator outcomes and modification analyses stay, too. Unchanged if
     *                       Exception::IllegalArgument or Exception::InvalidValue is thrown.
     * @param[in] in_spectra_files The input files of the search, in its order: one per entry of @c result.per_file.
     * @param[in] remove_decoys Whether the finalization of the protein FDR removes the decoy PSMs and proteins (default,
     *                          as in the ProSE tool's -out_merged). False keeps those that pass, e.g. to show them;
     *                          @c protein_fdr_applied of the result still tells whether the FDR was applied. No effect
     *                          without the protein FDR.
     * @return The pooled run and PSMs; empty if no file was searched.
     *
     * @throws Exception::IllegalArgument if @p in_spectra_files has not one entry per file of @p result.
     * @throws Exception::InvalidValue if Percolator ('rescore') rescored only some of the files that were searched: their
     *         runs report different search engines and scores and cannot be pooled.
     */
    MergedResult mergeFileResults(MultiFileSearchResult& result, const StringList& in_spectra_files,
                                  bool remove_decoys = true) const;

    /// With @p keep = true, the multi-file searchWithModificationAnalysis() overloads also return the PSMs of each
    /// file before the rescoring ('rescore') and FDR:PSM, in SearchResult::pre_fdr_peptide_ids (e.g. for a
    /// Percolator .pin): the search engine scores, every reported hit and the decoys. Default: false (no copy).
    void setKeepPreFdrPsms(bool keep)
    {
      keep_pre_fdr_psms_ = keep;
      if (!keep) { pre_fdr_peptide_ids_ = PeptideIdentificationList(); }
    }

  protected:
    void updateMembers_() override;

    /// Slimmer structure as storing all scored candidates in PeptideHit objects takes too much space
    struct AnnotatedHit_
    {
      AASequence sequence;
      /*
      std::string_view sequence;
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
      UInt32 index_peptide = 0; ///< the candidate: its index in FragmentIndex::getPeptides() of the searched index (protein mapping)

      static bool hasBetterScore(const AnnotatedHit_& a, const AnnotatedHit_& b)
      {
        if (a.score != b.score) return a.score > b.score;
        return a.sequence < b.sequence;
      }
    };

    /**
     * @brief Running summary of the *complete* candidate pool of one spectrum.
     *
     * scoreSpectraAgainstIndex_() prunes each spectrum to max(report_top_hits_, 2)
     * candidates as soon as it has scored them, so by the time postProcessHits_()
     * runs the surviving hits are no longer a sample of the search space: with the
     * default report:top_hits=1 only two candidates remain, and derived features
     * such as `ln_num_candidates` or `hyperscore_zscore` would degenerate into a
     * binary flag and a rescaled delta score respectively.
     *
     * add() is therefore called for every candidate the moment it is scored --
     * before pruning and before zero-scoring candidates are dropped -- so the
     * summary reflects the full pool. Instances accumulate across chunks in the
     * chunked search paths, where each chunk contributes its own candidates for
     * the same spectrum. With peptide:deduplicate enabled for a fragment index,
     * repeated peptidoform/charge/isotope hypotheses in later chunks are skipped
     * before add(). The per-chunk candidate cap still determines which distinct
     * hypotheses enter the pool; chunked and unchunked pools can therefore differ.
     */
    struct CandidatePoolStats_
    {
      double sum = 0.0;         ///< sum of all candidate scores
      double sumsq = 0.0;       ///< sum of squared candidate scores
      double best = 0.0;        ///< best candidate score (HyperScore is non-negative, so 0 doubles as "none seen")
      double second_best = 0.0; ///< runner-up candidate score
      Size count = 0;           ///< number of candidates scored
      /// Exact peptide/charge/isotope keys already scored in earlier database chunks.
      /// Empty in non-chunked searches; full keys resolve hash collisions by equality.
      std::unordered_set<std::string> seen_candidates;

      /// Fold one freshly scored candidate into the summary.
      void add(double score)
      {
        if (score > best) { second_best = best; best = score; }
        else if (score > second_best) { second_best = score; }
        sum += score;
        sumsq += score * score;
        ++count;
      }

      /// Best minus runner-up over the full pool.
      double deltaScore() const { return best - second_best; }

      /**
       * @brief How much of an outlier the best score is within its own candidate pool.
       *
       * Standard score of `best` against the mean and (population) SD of all
       * `count` candidates -- a lightweight significance proxy, since HyperScore
       * (unlike e.g. MS-GF+'s SpecEValue) has no closed-form e-value.
       *
       * Standardising against the whole pool rather than against the `count - 1`
       * non-best candidates keeps the statistic bounded: by Samuelson's inequality
       * a pool member cannot deviate from its own mean by more than sqrt(count - 1)
       * SDs, so the feature tops out at sqrt(scoring:max_candidates_per_spectrum - 1)
       * (7 at the default cap of 50). The leave-one-out form has no such bound --
       * with two candidates the SD of the single remaining score is 0 by
       * construction, and the ratio diverges. That mattered in practice: on a
       * 70k-PSM Orbitrap run the leave-one-out form produced values up to 6.5e+07,
       * which is enough to dominate Percolator's feature standardisation.
       *
       * Returns 0 when fewer than two candidates were scored, and when every
       * candidate scored the same -- in both cases the best candidate does not
       * stand out from the pool.
       */
      double zScore() const;
    };

    /**
      @brief MS2 deisotoping rule of preprocessSpectra_().

      Filled from the parameters fragment:deisotope_* by updateMembers_(), whose defaults are the Sage-like rule
      (min_peaks = 2, charge_cap_precursor: charges up to the precursor charge, at most 3, and sum_intensity). The member
      defaults are the behaviour before these parameters existed (envelopes of at least three peaks, charges 1-3 in every
      spectrum, the monoisotopic peak keeps its own intensity); the preprocessSpectra_() overload without settings uses them.
    */
    struct DeisotopingSettings_
    {
      unsigned int min_peaks{3};          ///< fragment:deisotope_min_peaks
      bool charge_cap_precursor{false};   ///< fragment:deisotope_charge_cap = precursor
      bool sum_intensity{false};          ///< fragment:deisotope_sum_intensity
    };

    /**
      @brief filter, deisotope, decharge spectra

      @p deisotoping sets the envelope rule of the deisotoper (see DeisotopingSettings_).

      With @p ion_evidence, a compact copy of every spectrum's peaks is kept for the self-trained ion priors
      (annotate:self_trained_ion_priors): the peaks after deisotoping, before the window and top-N filters, or, with
      @p ion_evidence_scored_peaks, the peaks that are scored. The lists are indexed like @p exp after its sort by RT.

      With @p query_spectra, a copy of each spectrum after zero-intensity removal and m/z sorting, before
      deisotoping and local filtering, is kept for candidate retrieval (fragment:query_spectrum=raw); with
      @p evidence_spectra, a copy after deisotoping, before window/top-N filtering (annotate:local_fragment_evidence).
      Both are aligned with @p exp after its sort by RT.
    */
    static void preprocessSpectra_(PeakMap& exp,
                                   double fragment_mass_tolerance,
                                   bool fragment_mass_tolerance_unit_ppm,
                                   bool deisotope_requested,
                                   Size peaks_keep_n,
                                   Int peaks_window_top,
                                   const std::string& window_type,
                                   const DeisotopingSettings_& deisotoping,
                                   FragmentIonLikelihoodModel::PeakLists* ion_evidence = nullptr,
                                   bool ion_evidence_scored_peaks = false,
                                   PeakMap* evidence_spectra = nullptr,
                                   PeakMap* query_spectra = nullptr);

    /// preprocessSpectra_() with the DeisotopingSettings_ defaults (the behaviour before its parameters existed)
    static void preprocessSpectra_(PeakMap& exp,
                                   double fragment_mass_tolerance,
                                   bool fragment_mass_tolerance_unit_ppm,
                                   bool deisotope_requested,
                                   Size peaks_keep_n,
                                   Int peaks_window_top,
                                   const std::string& window_type = "auto");

    /**
     * @brief Learn per-run fragment ion likelihoods from confident PSMs and annotate every hit with them.
     *
     * Implements annotate:self_trained_ion_priors. The spectra are split into two folds by the parity of their
     * scan_index. Each fold trains its own FragmentIonLikelihoodModel on its rank-one target hits at target-decoy
     * competition q <= annotate:ion_prior_train_fdr of the native score, computed within the fold, with their
     * reversed sequences matched against the same spectra as noise. Every hit is then scored by the model of the
     * other fold (cross-fitting), so neither a PSM's own spectrum nor its label, nor any other label of its fold,
     * enters the model that scores it; targets and decoys are scored alike. Both folds need
     * annotate:ion_prior_min_psms training PSMs (a fold without a decoy hit has the estimate (0 + 1) / T); otherwise
     * all features of the file are 0 and a warning is logged. This run-level gate is the one place where the labels
     * of a fold reach its own features: at the boundary, a single label (in either fold) decides whether both folds
     * get their cross-fitted features or all zeros. The gate treats targets and decoys alike, but it is not formally
     * independent of the labels. Adds the feature names to the search parameters' extra_features and records the
     * settings that determine the features (annotate:ion_prior_*) and the training set (ion_prior:*). Native scores
     * and the retained candidates are unchanged. Call after PeptideIndexing (target/decoy labels) and before FDR
     * (which overwrites the native scores). Does nothing when the ion priors are off, which includes target-only
     * searches (decoys=ignore).
     *
     * @param[in] spectra Preprocessed spectra, indexed by the scan_index meta value of the PSMs (precursor
     *            activation selects the ion series).
     * @param[in] evidence Peak lists of the same spectra from preprocessSpectra_.
     * @param[in,out] protein_ids Search parameters of protein_ids[0] receive the feature names and training facts.
     * @param[in,out] peptide_ids Every hit receives the meta values of annotate:ion_prior_features (ion_prior_llr,
     *                ion_prior_explained, ion_prior_topk_observed).
     */
    void annotateIonPriors_(const PeakMap& spectra,
                            const FragmentIonLikelihoodModel::PeakLists& evidence,
                            std::vector<ProteinIdentification>& protein_ids,
                            PeptideIdentificationList& peptide_ids) const;

    /// Reversed sequence keeping the C-terminal residue and every residue's modification: the noise hypothesis of a confident PSM
    static AASequence reversedNoiseSequence_(const AASequence& sequence);

    /// Keep the strongest peaks in each non-overlapping 100 Da window, including a short final window.
    static void filterLocalPeaks_(MSSpectrum& spectrum, Size peaks_per_window);

    /// How decoys are obtained/recognised for a search (parameter "decoys").
    enum class DecoyMode_
    {
      AUTO,     ///< ensure decoys: use those already present, otherwise generate
      GENERATE, ///< always (re)generate decoys from the targets (strip existing first)
      IGNORE    ///< search targets only; remove any decoys present in the database
    };

    /**
     * @brief Resolved decoy handling for one concrete input database.
     *
     * Produced by resolveDecoyStrategy_() and consumed by
     * buildDecoyAugmentedDB_() and the downstream PeptideIndexing / FDR steps,
     * so the same decoys that are searched are also the ones scored.
     */
    struct DecoyStrategy_
    {
      bool generate{false};       ///< reverse target proteins to synthesise decoys
      bool strip_existing{false}; ///< drop pre-existing decoy entries before searching
      bool have_decoys{false};    ///< searched DB will contain decoys (FDR possible)
      std::string decoy_string;   ///< effective marker for PeptideIndexing + protein FDR
      bool is_prefix{true};       ///< position of decoy_string
      std::string strip_string;   ///< marker of pre-existing decoys to strip
      bool strip_is_prefix{true}; ///< position of strip_string
    };

    /**
     * @brief Decide how to obtain/recognise decoys for @p db.
     *
     * Detects pre-existing decoys via the common-marker heuristic
     * (OpenMS::DecoyHelper), falling back to the configured decoy_prefix for
     * custom markers, then maps the "decoys" mode (auto/generate/ignore) onto a
     * concrete DecoyStrategy_.
     */
    DecoyStrategy_ resolveDecoyStrategy_(
        const std::vector<FASTAFile::FASTAEntry>& db) const;

    /**
     * @brief ProSE parameters made safe to hand to a FragmentIndex.
     *
     * FragmentIndex declares its own (unused) boolean "decoys" flag with
     * restrictions {true,false}. ProSE's "decoys" parameter is the richer
     * auto/generate/ignore enum and manages decoys at the database level, so
     * forwarding it verbatim would trip FragmentIndex's validation. This returns
     * getParameters() with "decoys" overridden to "false".
     *
     * @param[in] electron_ions Also index c and z+1 ions, for electron-activated spectra
     *            (ions:by_activation).
     */
    Param fragmentIndexParameters_(bool electron_ions = false) const;

    /// Theoretical spectrum generators for the configured ion series and for electron-activated
    /// spectra (ions:by_activation); defined in the source file
    struct SpectrumGenerators_;

    /// Generators for scoring, annotation and calibration, configured from the ions:* parameters
    SpectrumGenerators_ spectrumGenerators_() const;

    /// Local fragment evidence of one candidate (annotate:local_fragment_evidence).
    struct LocalFragmentEvidence_
    {
      double chance_match_surprise = 0.0;
      double mass_competition_evidence = 0.0;
    };

    /// Peak counts per Da in a fixed +/-50 m/z window around each peak of a sorted spectrum.
    static std::vector<double> localPeakDensities_(const MSSpectrum& spectrum);

    /**
     * @brief Compute local fragment-match evidence on the spectrum before window/top-N peak filtering.
     *
     * Matches each intact theoretical ion to its nearest experimental peak. Surprise is
     * sum(max(0, -log(2 * tolerance_Da * local_density))). Competition evidence is
     * sum(1 / (1 + alternative_assignments + local_density)). Alternative hypotheses
     * include all configured intact ions and generic H2O/NH3 losses, divided by the
     * fragment charge. Losses only discount ambiguous matches; they are not scored as
     * additional matches and do not require sequence-specific loss eligibility.
     *
     * @param[in] spectrum Sorted experimental spectrum after zero-intensity removal and
     *                    optional deisotoping, before window/top-N peak filtering.
     * @param[in] theoretical Sorted intact ions, with the TSG's aligned integer charge array.
     * @param[in] densities Result of localPeakDensities_(spectrum), reused for all candidates.
     * @param[in] tolerance Positive, finite matching tolerance.
     * @param[in] ppm Whether tolerance is given in ppm rather than Da.
     * @return Two additive annotations; zero without matches. Neither is a calibrated PSM p-value.
     */
    static LocalFragmentEvidence_ localFragmentEvidence_(const MSSpectrum& spectrum,
                                                        const MSSpectrum& theoretical,
                                                        const std::vector<double>& densities,
                                                        double tolerance, bool ppm);

    /// localFragmentEvidence_() restricted to the ions of @p theoretical with charge <= @p max_charge, with a
    /// caller-owned buffer for the sorted alternative assignments (reused across candidates).
    static LocalFragmentEvidence_ localFragmentEvidence_(const MSSpectrum& spectrum,
                                                        const MSSpectrum& theoretical,
                                                        int max_charge,
                                                        const std::vector<double>& densities,
                                                        double tolerance, bool ppm,
                                                        std::vector<double>& alternatives);

    /// localFragmentEvidence_() with the charges of the peaks of @p theoretical given (@p charges, one per peak)
    static LocalFragmentEvidence_ localFragmentEvidence_(const MSSpectrum& spectrum,
                                                        const MSSpectrum& theoretical,
                                                        const std::vector<Int>& charges,
                                                        int max_charge,
                                                        const std::vector<double>& densities,
                                                        double tolerance, bool ppm,
                                                        std::vector<double>& alternatives);

    /// True if the precursor of @p spectrum was activated by electrons (ETD, ECD, EThcD or ETciD)
    static bool isElectronActivated_(const MSSpectrum& spectrum);

    /// Number of spectra that ions:by_activation also scores with c and z+1 ions (0 if it is off)
    Size countElectronActivated_(const PeakMap& spectra) const;

    /**
     * @brief Build the searched database according to @p strategy.
     *
     * Optionally strips pre-existing decoys and/or generates fresh decoys by
     * reversing the target proteins. Produces FASTA entries only — does not
     * construct a FragmentIndex.
     */
    std::vector<FASTAFile::FASTAEntry> buildDecoyAugmentedDB_(
        const std::vector<FASTAFile::FASTAEntry>& fasta_db,
        const DecoyStrategy_& strategy) const;

    /// As above, but takes the entries of @p fasta_db over instead of copying them
    std::vector<FASTAFile::FASTAEntry> buildDecoyAugmentedDB_(
        std::vector<FASTAFile::FASTAEntry>&& fasta_db,
        const DecoyStrategy_& strategy) const;

    /// The multi-file searchWithModificationAnalysis() overloads: searches @p fasta_db_given, or, with @p in_db set,
    /// the database read from *in_db (while the first spectrum file is read in the background), whose entries the
    /// search then takes over instead of copying them.
    MultiFileSearchResult searchWithModificationAnalysis_(
        const std::vector<std::string>& in_spectra_files,
        const std::vector<FASTAFile::FASTAEntry>& fasta_db_given,
        const std::string* in_db,
        const std::vector<std::string>& output_base_names,
        const std::string& aggregate_base_name,
        bool build_pooled_aggregate) const;

    /// prepareContext(fasta_db, electron_ions) that takes the entries of @p fasta_db over instead of copying them.
    /// With @p searched_spectra, the index holds only the peptides these spectra can reach (see FragmentIndex::build):
    /// the context then serves only one search of these spectra, without calibration.
    SearchContext prepareContext_(std::vector<FASTAFile::FASTAEntry>&& fasta_db, bool electron_ions,
                                  const std::function<const PeakMap*(Size)>& searched_spectra = {}) const;

    /// Whether a single-use context may index only the peptides the searched spectra can reach: not with
    /// calibration, which may change the precursor windows after the index is built.
    bool restrictIndexToSpectra_() const { return !calibration_enabled_; }

    /**
     * @brief Build a strided protein sample for chunked calibration.
     *
     * Sample size is tied to database_chunk_size_ so the calibration FI never
     * exceeds the user's declared memory budget. Crucial for immunopeptidomics
     * (non-specific digestion) where a fixed-size 5000-protein sample would
     * generate tens of GB of fragment index. Strided rather than first-N so
     * small chunk_size values don't starve the calibration pool (#9182).
     *
     * Note: very small chunk_size may still produce a calibration pool below
     * calibration_min_psms_; calibration then no-ops silently (follow-up
     * issue: pool PSMs across first-K chunks of the main search loop).
     */
    std::vector<FASTAFile::FASTAEntry> buildCalibrationSample_(
        const std::vector<FASTAFile::FASTAEntry>& full_db) const;

    /**
     * @brief Chunked database search implementation.
     *
     * Splits full_db into chunks of database_chunk_size_ proteins, builds a
     * FragmentIndex per chunk, scores all spectra against each chunk, and
     * accumulates hits before a single post-processing pass. full_db is the
     * already-resolved database (decoys generated/stripped per @p strategy);
     * @p strategy carries the effective decoy marker for PeptideIndexing and
     * FDR. Taken by non-const reference because PeptideIndexing::run() mutates it.
     *
     * @return EXECUTION_OK on success; INPUT_FILE_EMPTY / UNEXPECTED_RESULT /
     *         UNKNOWN_ERROR on PeptideIndexing failure.
     */
    ExitCodes searchChunked_(
        PeakMap& spectra,
        std::vector<FASTAFile::FASTAEntry>& full_db,
        const DecoyStrategy_& strategy,
        std::vector<ProteinIdentification>& protein_ids,
        PeptideIdentificationList& peptide_ids) const;

    /**
     * @brief Score all spectra against one FragmentIndex.
     *
     * Shared by the non-chunked and chunked search paths. Appends per-scan
     * AnnotatedHit_ entries to annotated_hits and prunes each spectrum to
     * max(report_top_hits_, 2) candidates to bound memory. Expects a pre-built
     * FragmentIndex whose parameters already reflect any calibrated tolerances
     * the caller wants to apply.
     *
     * @param[in] spectra Preprocessed MS2 spectra to score.
     * @param[in] fi Pre-built FragmentIndex to query for candidates.
     * @param[in] db Database the FragmentIndex was built from, used to reconstruct candidate sequences.
     * @param[in] generators Generators for the theoretical spectrum of each candidate, chosen per spectrum.
     * @param[in] effective_fragment_tol Fragment mass tolerance to score with (calibrated, if calibration ran).
     * @param[in] fragment_mass_tolerance_unit_ppm Whether @p effective_fragment_tol is in ppm rather than Da.
     * @param[in] open_search_mode Whether to record the precursor delta mass on each hit.
     * @param[in,out] annotated_hits Per-spectrum candidate hits, one entry per spectrum.
     * @param[in,out] pool_stats Per-spectrum summary of the full (unpruned) candidate
     *                pool, one entry per spectrum. Accumulated rather than overwritten,
     *                so chunked callers can pass the same vector for every chunk.
     * @param[in] progress_label Label shown by the progress logger for this scoring pass.
     * @param[in] query_spectra Optional spectra for candidate retrieval, aligned with @p spectra
     *                          (fragment:query_spectrum=raw); scoring always uses @p spectra.
     */
    void scoreSpectraAgainstIndex_(
        const PeakMap& spectra,
        FragmentIndex& fi,
        const std::vector<FASTAFile::FASTAEntry>& db,
        const SpectrumGenerators_& generators,
        double effective_fragment_tol,
        bool fragment_mass_tolerance_unit_ppm,
        bool open_search_mode,
        std::vector<std::vector<AnnotatedHit_>>& annotated_hits,
        std::vector<CandidatePoolStats_>& pool_stats,
        const std::string& progress_label,
        const PeakMap* query_spectra = nullptr) const;

    /**
     * @brief Filter and annotate search results.
     *
     * Trims per-spectrum candidate hits to the top N and converts them into
     * PeptideIdentification objects, adding requested PSM annotations and
     * populating protein-level search metadata.
     *
     * @param[in] exp Input MS experiment providing spectra/metadata for annotation.
     * @param[in,out] annotated_hits Per-spectrum candidate hits (trimmed to @p top_hits in-place).
     * @param[in] pool_stats Per-spectrum summary of the full candidate pool as collected
     *            by scoreSpectraAgainstIndex_(), used for the pool-derived PSM features
     *            (delta score, z-score, candidate count). Must match @p annotated_hits in size.
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
     * @param[in] evidence_spectra Peak lists retained before local/top-N filtering, aligned with @p exp;
     *                             required when annotate:local_fragment_evidence is enabled, otherwise unused.
     * @param[in] defer_peak_annotations Leave out the fragment annotations (annotate:PSM fragment_annotation): the caller
     *                                   adds them with annotateFragmentPeaks_() to the hits it keeps.
     */
    void postProcessHits_(const PeakMap& exp,
      std::vector<std::vector<ProSEAlgorithm::AnnotatedHit_> >& annotated_hits,
      const std::vector<CandidatePoolStats_>& pool_stats,
      std::vector<ProteinIdentification>& protein_ids,
      PeptideIdentificationList& peptide_ids,
      Size top_hits,
      const StringList& modifications_fixed,
      const StringList& modifications_variable,
      Int peptide_missed_cleavages,
      double precursor_mass_tolerance,
      double fragment_mass_tolerance,
      const std::string& precursor_mass_tolerance_unit_ppm,
      const std::string& fragment_mass_tolerance_unit_ppm,
      const Int precursor_min_charge,
      const Int precursor_max_charge,
      const std::string& enzyme,
      const std::string& database_name,
      const PeakMap* evidence_spectra = nullptr,
      bool defer_peak_annotations = false) const;

    /**
     * @brief The fragment annotations that postProcessHits_() leaves out with defer_peak_annotations, for the hits of
     * @p peptide_ids (by their "scan_index" in @p exp), computed as postProcessHits_() computes them.
     */
    void annotateFragmentPeaks_(const PeakMap& exp, PeptideIdentificationList& peptide_ids,
                                double fragment_mass_tolerance, const std::string& fragment_mass_tolerance_unit_ppm) const;

    /**
     * @brief Adds PeptDeep-prediction-derived rescoring features to the PSMs.
     *
     * No-op unless both @p peptdeep:ms2_model and @p peptdeep:rt_model are set. Warns
     * and does nothing when OpenMS was built without ONNX support, so an ini file that
     * requests the models still runs (without those features) on such a build.
     */
    void annotatePeptDeepFeatures_(const PeakMap& spectra,
                                   std::vector<ProteinIdentification>& protein_ids,
                                   PeptideIdentificationList& peptide_ids) const;

    bool peptdeep_enable_{false};
    std::string peptdeep_ms2_model_;
    std::string peptdeep_rt_model_;
    std::string peptdeep_instrument_{"QE"};
    double peptdeep_nce_{-1.0};
    std::string peptdeep_rt_model_type_{"b_spline"};
    int peptdeep_batch_size_{500};

    /// Calibration overwrites these with the calibrated magnitudes for the duration of
    /// search(); pure runtime-state mutation that does not affect the logical const-ness
    /// of search(), matching the `mutable` pattern used by last_calibration_result_.
    mutable double precursor_mass_tolerance_lower_{10.0};   ///< positive magnitude (default matches the param)
    mutable double precursor_mass_tolerance_upper_{10.0};   ///< positive magnitude (default matches the param)
    std::string precursor_mass_tolerance_unit_{"ppm"};

    Size precursor_min_charge_;
    Size precursor_max_charge_;

    double fragment_mass_tolerance_;

    std::string fragment_mass_tolerance_unit_;

    /// Resolved MS2 deisotoping request (param fragment:deisotope != "false").
    /// preprocessSpectra_ still gates on Deisotoper::isToleranceSupported() so the
    /// deisotoper is never called out of range (it would throw -> terminate in the
    /// OpenMP region). See OpenMS#9619.
    bool deisotope_requested_{true};
    bool scoring_multiple_charges_{false}; ///< Resolved scoring:fragment_charges: score charges above 1
    int scoring_max_charge_{2};            ///< Upper fragment charge bound (fragment:max_charge)
    bool mass_accuracy_score_{false};      ///< Resolved scoring:method: mass-accuracy weighted HyperScore
    double mass_error_sd_ppm_{7.0};        ///< Width of the mass-accuracy kernel in ppm (scoring:mass_error_sd)

    /// Fragment charges that score a candidate of the given precursor charge; also used by the calibration pass.
    int scoringMaxCharge_(int precursor_charge) const
    {
      return scoring_multiple_charges_ ? std::max(1, std::min(precursor_charge - 1, scoring_max_charge_)) : 1;
    }
    Size peaks_keep_n_{0};     ///< NLargest cap on MS2 peaks before scoring; 0 = resolution-aware auto (peaks:keep_n)
    Int peaks_window_top_{20}; ///< WindowMower peaks-per-100Da before scoring (peaks:window_top)
    std::string peaks_window_type_ {"auto"}; ///< Resolution-aware treatment of the final peak window
    DeisotopingSettings_ deisotoping_;        ///< MS2 deisotoping rule (fragment:deisotope_*)
    bool query_raw_spectrum_{false}; ///< Retrieve candidates with the spectrum before deisotoping and local filtering (fragment:query_spectrum)

    StringList modifications_fixed_;

    StringList modifications_variable_;

    Size modifications_max_variable_mods_per_peptide_;

    std::string enzyme_;

    DecoyMode_ decoy_mode_{DecoyMode_::AUTO};
    std::string decoy_prefix_;

    bool rescore_{true}; ///< rescore: Percolator rescoring of every run before its PSM-level FDR
    PercolatorTraining percolator_training_; ///< rescoring:...: options of Percolator's training
    Param percolator_extra_; ///< rescore_method 'reset': Percolator's RESET parameters (empty: Percolator)
    mutable std::string reset_run_name_; ///< spectra file of the run being rescored (names its RESET stem table)
    std::string reuse_mode_{"per_file"}; ///< rescore_reuse:mode
    Size reuse_max_train_{250000};       ///< rescore_reuse:max_train
    Size reuse_static_files_{1};         ///< rescore_reuse:static_files
    double reuse_max_shift_{0.5};        ///< rescore_reuse:max_shift
    double fdr_psm_{0.0};
    bool fdr_psm_by_scored_charges_{true}; ///< FDR:PSM_groups: separate competitions per number of scored fragment charges
    double fdr_protein_{0.0};

    StringList annotate_psm_;

    bool self_trained_ion_priors_{false};       ///< annotate:self_trained_ion_priors
    FragmentIonLikelihoodModel::ContextSet ion_prior_model_{FragmentIonLikelihoodModel::ContextSet::RICH}; ///< annotate:ion_prior_model
    bool ion_prior_scored_peaks_{false};        ///< annotate:ion_prior_peaks == scored
    int ion_prior_max_fragment_charge_{2};      ///< annotate:ion_prior_max_fragment_charge (capped at precursor charge - 1)
    double ion_prior_train_fdr_{0.01};          ///< annotate:ion_prior_train_fdr: TDC q-value threshold of the training PSMs
    Size ion_prior_min_psms_{100};              ///< annotate:ion_prior_min_psms: training PSMs required per fold
    bool ion_prior_feature_llr_{true};          ///< annotate:ion_prior_features holds ion_prior_llr
    bool ion_prior_feature_explained_{true};    ///< annotate:ion_prior_features holds ion_prior_explained
    bool ion_prior_feature_topk_{false};        ///< annotate:ion_prior_features holds ion_prior_topk_observed

    Size peptide_min_size_;
    Size peptide_max_size_;
    Size peptide_missed_cleavages_;
    EnzymaticDigestion::Specificity peptide_enzyme_specificity_{EnzymaticDigestion::SPEC_FULL};

    std::string peptide_motif_;

    Size report_top_hits_;
    /// report:isotope_error_convention: the PSMs report isotope_error as observed minus theoretical (true) or with the
    /// sign of the search offset, theoretical minus observed (false)
    bool isotope_error_observed_minus_theoretical_{true};

    bool add_a_ions_{false};
    bool add_b_ions_{true};
    bool add_c_ions_{false};
    bool add_x_ions_{false};
    bool add_y_ions_{true};
    bool add_z_ions_{false};
    bool add_zp1_ions_{false};
    bool ions_by_activation_{true}; ///< add c and z+1 ions for electron-activated spectra
    bool protein_mapping_from_index_{true}; ///< peptide:protein_mapping = index

    Size database_chunk_size_{0};  ///< 0 = disabled; >0 = chunk DB into groups of this many proteins

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

    /// Per-run statistics of the most recent search(spectra, ctx, ...) call.
    /// Bridges the const search() (which returns ExitCodes, not a SearchResult)
    /// to its callers, which copy this into SearchResult::stats. Reset at the
    /// start of each search() call. `mutable` for the same reason as above:
    /// pure diagnostic state, orthogonal to logical const-ness.
    mutable RunStatistics last_run_stats_;

    /// Percolator rescoring outcome of the most recent search(spectra, ctx, ...) call, bridged to
    /// SearchResult::percolator as last_run_stats_ is to SearchResult::stats.
    mutable PercolatorRescoring last_percolator_;
    bool keep_pre_fdr_psms_ = false; ///< see setKeepPreFdrPsms()
    /// The PSMs of the most recent search(spectra, ctx, ...) before FDR:PSM (only with keep_pre_fdr_psms_), bridged to
    /// the caller like last_run_stats_.
    mutable PeptideIdentificationList pre_fdr_peptide_ids_;
    /// The hits that the most recent rescoring dropped (all but the best per spectrum), destroyed by a helper thread
    /// meanwhile (see applyPsmFDR_()); the holder joins it when it goes (at the next rescoring of this object, or with
    /// the last copy of this object that shares it).
    mutable std::shared_ptr<void> released_hits_;

    /// rescore_reuse: how the PSMs of a run are scored instead of training Percolator on them (defined in the .cpp file)
    struct ReuseScoring_;
    /// rescore_reuse: the two cross-fitted models (defined in the .cpp file)
    struct ReuseModel_;
    /// rescore_reuse: a run that trains models: its index in the search (the random halves depend on it) and its PSMs
    struct ReuseRun_
    {
      Size file_index;
      const std::vector<ProteinIdentification>* protein_ids;
      PeptideIdentificationList* peptide_ids;
    };

    /// rescoreWithPercolator(); with @p dropped, the hits it drops are moved there (as memory allows) instead of
    /// destroyed; @p training: the options of Percolator's training (parameters 'rescoring:...' and 'rescore_solver');
    /// @p extra: further Percolator parameters (rescore_method 'reset': the RESET parameters). With @p reuse, the PSMs
    /// are scored as it says instead of by a model trained on them.
    static PercolatorRescoring rescoreWithPercolator_(std::vector<ProteinIdentification>& protein_ids,
                                                      PeptideIdentificationList& peptide_ids,
                                                      int threads,
                                                      std::vector<std::vector<PeptideHit>>* dropped,
                                                      const PercolatorTraining& training,
                                                      const Param& extra,
                                                      const ReuseScoring_* reuse = nullptr);

    /// rescore_reuse: trains the two cross-fitted models on PSMs of @p runs (see the parameter); throws if Percolator cannot
    /// train them
    std::shared_ptr<const ReuseModel_> trainReuseModel_(const std::vector<ReuseRun_>& runs, int threads) const;

    /// rescore_reuse: rescores the PSMs of a run (file @p file_index of the search) with @p model, as rescoreWithPercolator_();
    /// with @p check_shift, a run whose PSMs (targets and decoys, labels unused) differ from the training PSMs by more than reuse_max_shift_ is rescored by
    /// Percolator trained on it instead
    PercolatorRescoring rescoreWithReuseModel_(std::vector<ProteinIdentification>& protein_ids,
                                               PeptideIdentificationList& peptide_ids,
                                               Size file_index,
                                               const ReuseModel_& model,
                                               bool check_shift,
                                               int threads,
                                               std::vector<std::vector<PeptideHit>>* dropped) const;

    /// FDR:PSM on the PSMs of a run after its rescoring (see applyPsmFDR_())
    void applyRescoredPsmFDR_(const std::vector<ProteinIdentification>& protein_ids,
                              PeptideIdentificationList& peptide_ids,
                              bool has_decoys,
                              RunStatistics& stats,
                              bool rescored) const;

    /// rescore_reuse joint/static: trains the models on the files of @p deferred (indices into @p mfres.per_file, searched
    /// with reuse_defer_), then rescores them and applies FDR:PSM
    void finishDeferredRescoring_(MultiFileSearchResult& mfres, const std::vector<Size>& deferred) const;

    /// rescore_reuse: the models of a multi-file search ('static': trained on its first files; null otherwise)
    mutable std::shared_ptr<const ReuseModel_> reuse_model_;
    /// rescore_reuse: the current search defers its rescoring and FDR:PSM to finishDeferredRescoring_()
    mutable bool reuse_defer_ = false;
    /// rescore_reuse: the index of the file of the current search in a multi-file search (0 otherwise)
    mutable Size reuse_file_index_ = 0;

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
     * @param[in] query_spectra Optional spectra for candidate retrieval, aligned with @p spectra.
     * @return CalibrationResult_ with estimated tolerances, or success=false if insufficient PSMs.
     */
    CalibrationResult_ runCalibrationPass_(PeakMap& spectra,
                                           FragmentIndex& fragment_index,
                                           const std::vector<FASTAFile::FASTAEntry>& db,
                                           const PeakMap* query_spectra = nullptr) const;

    /**
     * @brief The protein occurrences of the sequences of a search's hits (peptide:protein_mapping = index).
     *
     * buildProteinMapping_() collects them while the fragment index exists, applyProteinMapping_() writes what
     * PeptideIndexing::run() writes. PeptideIndexing (without mismatches and I/L equivalence; missed cleavages
     * ignored; the initial methionine and the residue after it may be cleaved off) maps an unmodified hit
     * sequence S to every (protein, position) at which S occurs, where up to aaa_max ambiguous residues (B, J, Z, X)
     * of the protein may stand for residues of S, and isValidProduct() accepts the span. These are:
     *  1. the spans of the digest of the index with the residues S (FragmentIndex::getProteinOccurrences(), which
     *     includes the occurrences removed by peptide:deduplicate): they lie at enzymatic cleavage sites, have the
     *     missed cleavages, length and mass of S and no ambiguous residue;
     *  2. spans at position 1 or 2 of a protein that starts with M, which isValidProduct() takes as N-terminal;
     *  3. spans over 1 to aaa_max ambiguous residues, each one standing for a residue that AhoCorasickAmbiguous
     *     matches it with.
     * The spans of 2. and 3. are checked with isValidProduct(). Configurations, databases and hits for which this
     * does not reproduce PeptideIndexing exactly leave a reason (fallback_reason); then PeptideIndexing runs.
     */
    struct ProteinMapping_
    {
      /// Number of hash tables that sequences are spread over (by their hash), so that they are filled in parallel
      static constexpr Size SHARDS = 64;
      /// The unmodified hit sequences (views into the database), one entry each
      std::vector<std::string_view> sequences;
      /// Per entry: its occurrences {protein index, position}, ascending (PeptideIndexing's evidence order)
      std::vector<std::vector<std::pair<UInt32, UInt32>>> occurrences;
      /// Sequence -> entry, in SHARDS hash tables
      std::vector<std::unordered_map<std::string_view, Size>> entries;
      /// Why PeptideIndexing has to run instead; empty if the mapping applies
      std::string fallback_reason;

      /// The entry of @p sequence, or sequences.size() if there is none
      Size find(std::string_view sequence) const
      {
        if (entries.empty()) return sequences.size();
        const auto& table = entries[std::hash<std::string_view>{}(sequence) % SHARDS];
        const auto entry = table.find(sequence);
        return entry == table.end() ? sequences.size() : entry->second;
      }
    };

    /**
     * @brief Collects the protein occurrences of the sequences of the index peptides @p candidates.
     *
     * @param[in] index The searched fragment index (built from @p db); must still hold its peptides
     * @param[in] db The searched database
     * @param[in] candidates Indices into FragmentIndex::getPeptides(): the candidates of the hits (in any order,
     *            repeats allowed)
     * @param[in] indexer_parameters The parameters PeptideIndexing would run with
     */
    static ProteinMapping_ buildProteinMapping_(const FragmentIndex& index,
                                                const std::vector<FASTAFile::FASTAEntry>& db,
                                                std::vector<UInt32> candidates,
                                                const Param& indexer_parameters);

    /**
     * @brief Writes what PeptideIndexing::run(@p db, @p protein_ids, @p peptide_ids) with @p indexer_parameters
     * writes: the evidences, target_decoy and protein_references of every hit, the protein hits and the
     * PeptideIndexer:* search parameters.
     *
     * @return false, with @p fallback_reason set and nothing changed, if @p mapping cannot reproduce PeptideIndexing
     *         for these identifications (e.g. a hit whose sequence it does not hold); then PeptideIndexing has to run.
     */
    static bool applyProteinMapping_(const ProteinMapping_& mapping,
                                     const std::vector<FASTAFile::FASTAEntry>& db,
                                     const Param& indexer_parameters,
                                     std::vector<ProteinIdentification>& protein_ids,
                                     PeptideIdentificationList& peptide_ids,
                                     std::string& fallback_reason);

    /// Helper: does @p accession carry the decoy @p marker at the given position?
    /// Empty marker → false. Pure std::string (no String dependency).
    static bool accessionHasDecoyMarker_(const std::string& accession,
                                         const std::string& marker, bool is_prefix);

    /// Helper: pool the runs of the files of @p per_file that were searched (exit code EXECUTION_OK) with
    /// IDMergerAlgorithm into one run (@p protein_ids, one entry) that lists @p in_spectra_files as primary MS run
    /// paths, and their PSMs (@p peptide_ids). Copies the identifications, or with @p move_psms moves the PSMs out of
    /// @p per_file (leaving them empty) and copies the runs, which stay with the files. 'id_merge_index' of a PSM is
    /// the input position of its file (its index in @p per_file and @p in_spectra_files, for inputs with distinct
    /// paths), also when other files failed. Shared by the multi-file search's pooled aggregate and
    /// mergeFileResults().
    static void poolSearchedFiles_(std::vector<SearchResult>& per_file, bool move_psms, const StringList& in_spectra_files,
                                   std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids);

    /// Helper: capture statistics that must be read BEFORE FDR filtering, namely
    /// target/decoy PSM counts (via the "target_decoy" meta set by PeptideIndexing;
    /// "target+decoy" counts as target, matching OpenMS FDR semantics) and the
    /// HyperScore distribution (FDR overwrites each hit's score with its q-value).
    /// Fills stats.target_psms, stats.decoy_psms, stats.hyperscore_* and
    /// stats.score_stats_valid.
    static void capturePreFdrStats_(const PeptideIdentificationList& peptide_ids,
                                    RunStatistics& stats);

    /// Helper: maximum top-hit score among retained PSMs. After PSM-FDR the score
    /// is the q-value, so this is the achieved FDR of the filtered set
    /// (-1.0 if empty).
    static double maxRetainedScore_(const PeptideIdentificationList& peptide_ids);

    /// Helper: the PSM-level FDR of one search run, the last step of every search path. With
    /// 'rescore' on, rescores the PSMs with Percolator first (rescoreWithPercolator()) and then
    /// applies FDR:PSM with applyDeferredPsmFDR(); otherwise applies FDR:PSM to the search engine
    /// scores (if @p has_decoys). Fills the FDR fields of @p stats. If @p pre_fdr_psms is given, it
    /// receives the PSMs before the rescoring and FDR:PSM (see setKeepPreFdrPsms()).
    /// @return The outcome of the rescoring (not applied, no reason with 'rescore' off).
    PercolatorRescoring applyPsmFDR_(std::vector<ProteinIdentification>& protein_ids,
                                     PeptideIdentificationList& peptide_ids,
                                     bool has_decoys,
                                     RunStatistics& stats,
                                     PeptideIdentificationList* pre_fdr_psms = nullptr) const;

    /// Helper: fill a RunStatistics with identification counts, distributions and
    /// per-run tolerance estimation for one searched file (post-FDR). Silent;
    /// rendering is done separately (ProSE TOPP tool or renderRunSummary()).
    /// Fields filled at other points (target/decoy counts, achieved FDR, timing)
    /// are left untouched.
    void collectRunStatistics_(const PeakMap& spectra,
                               const std::vector<ProteinIdentification>& protein_ids,
                               const PeptideIdentificationList& peptide_ids,
                               RunStatistics& stats) const;

  public:
    /// Recompute result-level statistics (matched count, unique peptide/protein
    /// counts, charge & missed-cleavage histograms, target/decoy counts and
    /// achieved PSM FDR) from a FINAL, post-processed PSM list. Use this when the
    /// identifications are mutated AFTER search() returns — e.g. the ProSE TOPP
    /// tool's Percolator rescoring + deferred FDR — so the report reflects what
    /// the user actually receives. Does NOT touch @c ms2_spectra, the HyperScore
    /// distribution (captured pre-FDR during the search) or any timing field.
    /// @p fdr_applied records whether a PSM-level FDR filter was applied; when
    /// true @c achieved_psm_fdr is set to the maximum retained q-value.
    static void updateFinalStats(RunStatistics& stats,
                                 const PeptideIdentificationList& peptide_ids,
                                 const std::string& enzyme,
                                 bool fdr_applied);

    /// Free the protein and peptide identifications of @p result with the OpenMP threads (destroyed one after the
    /// other, the hits of a file take tens of ms at the end of a run). Leaves them empty; the statistics are kept.
    static void releaseInParallel(SearchResult& result);

    /// Render a human-readable single-run summary block (configuration recap,
    /// database/index stats, identification results, tolerance estimate, timing
    /// and — for open searches — modification discovery) to @p os.
    /// @p shared carries the configuration/db/index context (built once).
    static void renderRunSummary(const RunStatistics& stats,
                                 const SharedSearchStats& shared,
                                 const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
                                 bool is_open_search,
                                 std::ostream& os);

    /// Render the modification-discovery section (top PTMs, unknown delta masses)
    /// to @p os. Shared by the per-file and aggregate report blocks.
    static void renderModificationSummary(const OpenSearchModificationAnalysis::OpenSearchAnalysisResult& mod_analysis,
                                          std::ostream& os);

    /// Serialize the complete end-of-search report to a machine-readable YAML string:
    /// shared configuration/database/index facts, per-file identification statistics,
    /// the output-file @p manifest (label -> written paths) and failed/total file counts.
    /// Hand-rolled (no YAML library dependency); every string scalar is double-quoted so
    /// values containing ':' (e.g. Windows paths) or other YAML metacharacters round-trip
    /// safely, and non-finite numbers are emitted as null.
    static std::string renderRunSummaryYaml(
        const MultiFileSearchResult& mfres,
        const std::vector<std::pair<std::string, std::vector<std::string>>>& manifest,
        Size files_failed,
        Size files_total);

  private:

    /// Helper function to determine if open search should be used based on tolerance
    bool isOpenSearchMode_() const
    {
      return FragmentIndex::isOpenSearchMode(precursor_mass_tolerance_lower_,
                                             precursor_mass_tolerance_upper_,
                                             precursor_mass_tolerance_unit_ == "ppm");
    }
};

} // namespace
