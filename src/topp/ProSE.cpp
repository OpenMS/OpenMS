// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/openms_data_path.h>
#ifdef WITH_ONNX
#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepUtils.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <new>
#include <numeric>
#include <optional>
#include <sstream>

#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ProSE ProSE

@brief Identifies peptides in MS/MS spectra.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; ProSE &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML or Bruker .d format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

@em This search engine is mainly for educational/benchmarking/prototyping use cases.
It lacks behind in speed and/or quality of results when compared to state-of-the-art search engines.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.
@note Open-search mode is automatically determined by the precursor mass tolerance: enabled when tolerance exceeds 1 Da or 1000 ppm. No explicit open-search parameter is needed. This is logged at runtime and recorded in the output search parameters as UserParam 'open_search'.
@note Decoy handling is controlled by '-Search:decoys'. The default 'auto' ensures decoys are available for target-decoy FDR: it reuses decoys already present in the FASTA (the marker is auto-detected, prefix or suffix, e.g. from DecoyDatabase) or generates them internally (prefixing accessions with '-Search:decoy_prefix', default "DECOY_") if none are found. Use 'generate' to always (re)build decoys from the targets, or 'ignore' to search the targets only.
@note Decoy reporting is tied to protein-level FDR, not to PSM-level FDR. Setting '-Search:FDR:protein' > 0 signals "finalize this result": picked-protein FDR is applied and decoys are removed. PSM-level FDR ('-Search:FDR:PSM') filters target and decoy PSMs alike by the q-value threshold; it does no decoy-specific stripping, so decoys that pass the threshold are kept. With '-Search:FDR:protein' = 0 (the default), decoys are retained in every output (target+decoy evidence with scores). The '-out_pin' files are an exception to both: they are written before any FDR filter and keep the decoys. To obtain a clean, decoy-free result without protein FDR, run @ref TOPP_IDFilter with '-remove_decoys' downstream.
@note Protein-level FDR scope: picked-protein FDR does not compose across runs, so it is applied (and decoys removed) only on a @em complete protein set — a single input file, or the pooled '-out_merged' set of a multi-file run. For a multi-file run with '-Search:FDR:protein' > 0 but no '-out_merged', protein FDR is NOT applied (no output represents a complete experiment); ProSE warns and leaves the per-file outputs as intermediates (decoys retained).
@note Percolator rescoring: by default ('-Search:rescore true'), the PSMs of each input file are rescored separately with OpenMS' built-in (in-process) Percolator implementation, the same backend PercolatorAdapter uses by default. No percolator executable is needed (the former '-percolator_executable' parameter and '-rescore' flag were removed). Training uses the standard Percolator feature set plus ProSE's search engine features, '-train_best_positive' and target-decoy competition ('-post_processing_tdc'); Percolator q-values become the main score (the SVM score and PEP are kept as meta values MS:1001492/MS:1001493, the HyperScore as meta value 'ln(hyperscore)'), only the best hit per spectrum is kept, and the runs report 'Percolator' as search engine. The PSM and protein FDR thresholds are then applied on the rescored PSMs. A file with fewer than 100 PSMs or without decoys, or for which Percolator fails, is @em not rescored: ProSE prints a warning, keeps its HyperScores (FDR is then computed from them) and continues; the end-of-search report states how many files were rescored. Set '-Search:rescore false' to report ProSE's own scores. '-out_pin' always holds ProSE's own scores (see there).
@note `Search:peptide:clip_nterm_methionine=true` also searches the mature N-terminal
peptide after loss of exactly the initial methionine of an M-leading protein.
Both forms are searched with the configured length, mass and missed-cleavage limits.
Protein N-terminal variable modifications also apply after clipping. Set false to
restore the previous search space; internal methionines are not clipped. Generated
decoys retain the initial M when enabled; supplied decoys are searched as provided.
@note Local peak filtering: `Search:peaks:window_type=auto` retains up to `Search:peaks:window_top`
peaks in every 100 Da window, including the short final window, for high-resolution fragments
(tolerance <= 0.1 Da or <= 100 ppm). Low-resolution auto retains the legacy width-scaled final
quota. Use `jump` for legacy filtering at every resolution, or `jump_full` to always retain
the full quota. This changes high-resolution preprocessing, not the scoring formula.
@note Self-trained fragment ion priors: with `Search:annotate:self_trained_ion_priors` (on by default) each file learns
how its fragment ions appear: presence, intensity rank among all peaks after deisotoping (before the
`Search:peaks:window_top` and `Search:peaks:keep_n` filters) and mass error, per ion series, precursor and fragment
charge, relative cleavage position, residues at the cleavage site (before P, after D or E) and presence of the
complementary ion (`Search:annotate:ion_prior_model` 'rich'; 'basic' keeps series, charges and position, and the rank
only). The model learns from the confident target PSMs (target-decoy competition q <= `Search:annotate:ion_prior_train_fdr`
of the native score) against their reversed sequences on the same spectra. The spectra are split into two halves by scan
parity, and every PSM is scored by the model of the other half (cross-fitting), so a PSM's own spectrum and label never
enter the model that scores it. Every PSM gets the Percolator features of `Search:annotate:ion_prior_features` (ion_prior_llr and
ion_prior_explained); native scores and the reported candidates are unchanged. Nothing is pre-trained: a file with fewer
than `Search:annotate:ion_prior_min_psms` confident PSMs in either half gets zeros (and a warning); at that boundary a
single label decides between the features and zeros for the whole file. Target-only searches
(`Search:decoys` ignore) learn and write nothing. The peak lists are kept in a compact form (m/z and rank, 5 bytes per
peak) until the PSMs are annotated.

@note Ion mobility: the drift time of an MS2 spectrum is written to its identification. An ion mobility (e.g. 1/K0 of
Bruker timsTOF DDA-PASEF data) is written as meta value 'IM', with its unit as 'IM' on the run; a FAIMS compensation voltage
(Thermo FAIMS data, from .raw or mzML) is written as 'FAIMS_CV', which the FAIMS-aware tools (e.g. @ref TOPP_ProteomicsLFQ,
@ref TOPP_MS1LabeledWorkflow, @ref TOPP_FeatureFinderIdentification) read.
@note Extended QPX output: '-qpx_extended' (with '-out_qpx') adds optional columns to the QPX psm view and writes the
QPX mz, run, sample and provenance views, so that a viewer can show a search from the QPX files alone. QPX readers ignore
columns and keys they do not know; the names are those proposed in bigbio/qpx#355 and #259. Each row of the per-input
psm.parquet files and of quantms.psm.parquet gets: 'additional_scores' with every numeric score and meta value of the
PSM in snake_case instead of the known scores only ('global_qvalue': the PSM-level q-value, from target-decoy
competition or Percolator, computed over the PSMs of the input file; 'ln_hyperscore'; Percolator's 'percolator_score',
'percolator_qvalue', 'percolator_pep'; 'precursor_mz_error_ppm' (corrected for the isotope error), 'isotope_error' (as
the column), 'delta_mass' (open search) and the Percolator features), with 'higher_better' where it is known;
'spectrum_reference' (string, the native ID of the spectrum; null if the mz view has no row of it, see below); 'isotope_error'
(int16: the 13C peaks from the monoisotopic to the selected precursor peak, i.e. observed minus theoretical, whatever
'Search:report:isotope_error_convention' says); 'protein_mappings' (list of struct accession, start, end, aa_before,
aa_after: each protein of the peptide with the peptide's 1-based first and last position and the residues before and
after it, '-' at a protein terminus). The per-input psm files also get 'predicted_ion_type_array', 'predicted_mz_array',
'predicted_intensity_array': the b and y ions of charge 1 and (for precursors of charge 2 or more) 2, named as in
'ion_type_array' ("b3+", "y5++"), with their theoretical m/z and the intensity PeptDeep predicts, relative to the most
intense of them; ions of 0.01 or less are left out. The predictions use the settings of 'Search:peptdeep:*' (models,
instrument) at the normalised collision energy of the file: the one PeptDeep rescoring used ('Search:peptdeep:enable';
predicted again at that NCE), else 'Search:peptdeep:nce' if set, else the one PeptDeep rescoring selects, run on a
copy of the written PSMs; PSMs of precursor charge 1 to 6 are predicted, others get null. Each input is read once more
after its search for the spectra; a failed re-read is a warning and leaves out its mz view (and the
'spectrum_reference' of its PSMs null), a failed prediction a warning and null prediction columns; without ONNX
Runtime the prediction columns are null. Each input gets the QPX mz view (runs/<basename>.mz.parquet) of the MS2 spectra of
its written PSMs only (one row per spectrum, in the order of its first PSM; of several spectra of one native ID only
the first, as 'id' is the view's key, so the PSMs of the others have no spectrum there and a null 'spectrum_reference':
they are unresolved, as a join on their native ID would give the first spectrum, and a reader must not look them up by
'scan' either; 'id' the native ID, 'run_file_name',
'scan', 'ms_level', 'centroid', 'scan_start_time' (minutes), 'inverse_ion_mobility', 'ion_injection_time' (null if
unknown), 'total_ion_current' (of all its peaks), 'precursors', 'cv_params' (null)), whose 'mz' and 'intensity' are
filtered: the '-qpx_peak_top' most intense peaks of each '-qpx_peak_window' Th window (jumping windows: each starts at
its first peak; the last one too keeps '-qpx_peak_top' peaks) and the peaks of the fragment annotations 'mz_array' of
its PSMs ('-qpx_peak_top 0' keeps all); of deisotoped spectra, whose annotation has every fragment at its 1+ m/z, also
the peak of the fragment at its own charge 2 or 3, so that every annotated fragment keeps its peak (the charges the
search's deisotoping tried: up to 3, and with 'Search:fragment:deisotope_charge_cap precursor' at most the spectrum's
recorded precursor charge, whatever the charge of the PSM; up to 3 if that is unknown). Its footer records the filter:
'prose_peak_window', 'prose_peak_top', 'prose_peak_annotated' (kept). A PSM's spectrum is the mz row whose 'id' is its
'spectrum_reference'; a null one means no stored spectrum (the re-read failed, the spectrum was not found again, or
another spectrum of its native ID has the row). A per-input psm file and its mz view are written as soon as the input is
done (neither is kept for the rest of the run), in row groups of 256 rows, so that a viewer can read the spectrum and
the prediction of one PSM alone; the psm file's footer has the keys of its prediction columns:
'prose_prediction_ms2_model', 'prose_prediction_instrument'. quantms.psm.parquet has the rows of the inputs whose
files were written, in input order, without the prediction columns (and no rows if there are none). The output
directory is one QPX dataset (prefix 'quantms'): the quantms.* files of the whole run on top, the per-input files
(<basename>.psm/mz/pg.parquet) in its subfolder runs/, so an input may have any name, 'quantms' too.
At the end of the run three QPX metadata views describe it:
quantms.run.parquet (one row per input, in input order): 'run_accession' and 'run_file_name' (the input's stem),
'file_name', 'samples' (one sample per run: the stem, label 'LFQ'; ProSE knows no design), 'fraction', 'instrument',
'dissociation_method', 'acquisition_method', 'additional_terms' (null: unknown to ProSE), 'enzymes',
'modification_parameters' (the search's fixed and variable modifications: UNIMOD accession, name, fixed, position,
target residue); and the optional columns 'search_status' (EXECUTION_OK, INPUT_FILE_EMPTY, UNEXPECTED_RESULT,
UNKNOWN_ERROR, ILLEGAL_PARAMETERS), 'ms2_spectra_count', 'identified_spectra_count' (target PSMs after the FDR
filters), 'decoy_psm_count', 'psm_rows' (its rows in quantms.psm.parquet), int64; 'rescoring_applied' (bool; null:
the search did not rescore), 'rescoring_skipped_reason', 'rescoring_warning' (string); 'prediction_nce' (double; null:
no prediction); 'precursor_error_median_ppm', 'precursor_error_mad_ppm', 'precursor_tolerance_recommended_ppm',
'fragment_error_median_ppm', 'fragment_error_mad_ppm', 'fragment_tolerance_recommended_ppm' (double; null: no
estimate); 'search_seconds' (double).
quantms.sample.parquet: the sample of each run ('sample_accession' the stem; 'organism' and 'organism_part' "not
available", as SDRF has an unknown value).
quantms.provenance.parquet: one row per step, in order: 'database_search' 'peptide_search' (tool ProSE; parameters
database_file, enzyme, enzyme_specificity, missed_cleavages, peptide_min_length, peptide_max_length,
precursor_mass_tolerance_lower, precursor_mass_tolerance_upper, precursor_mass_tolerance_unit,
fragment_mass_tolerance, fragment_mass_tolerance_unit, min_charge, max_charge, fixed_modifications,
variable_modifications, clip_nterm_methionine, decoy_mode, decoy_string, decoy_position, open_search, calibration,
chunked, isotope_error_convention (the sign of ProSE's own isotope_error meta value, which Percolator saw:
observed_minus_theoretical or theoretical_minus_observed), extra_features (the Percolator features beyond the standard
ones, by their names in 'additional_scores'), db_target_proteins, db_decoy_proteins, indexed_peptides,
indexed_fragments, seconds_index_build, seconds_total, and in open search with -out_mod_analysis_dir
mod_analysis_delta_masses and mod_analysis_ptms, the file names of the run's tables in that directory; 'config': every
parameter of this tool as one JSON object); 'rescoring' 'spectrum_prediction' (tool OpenMS, PeptDeep's models; only if
some input was predicted; parameters method, ms2_model, instrument, rescoring_features); 'rescoring' 'psm_rescoring'
(tool OpenMS, its Percolator; with 'Search:rescore'; parameter method); 'fdr_filtering' 'target_decoy_fdr' (tool
ProSE; psm_fdr, protein_fdr, psm_fdr_method, protein_fdr_method, protein_fdr_applied: the protein FDR was applied to a
complete protein set, one input or the pooled set of -out_merged); 'protein_inference' 'protein_grouping' (tool
ProSE; algorithm, pooled). 'tool_version' is the version of OpenMS.
Protein groups: with '-qpx_extended' a run without protein inference gets one group per protein in its pg.parquet
(without a global_qvalue), and with '-out_merged' and several inputs quantms.pg.parquet holds the groups of the pooled
inference (with its protein FDR) instead of the concatenated groups of the inputs.
@note Decoys: '-keep_decoys' keeps the decoy PSMs and decoy proteins that pass the FDR thresholds where the protein FDR
would remove them (a single input, and -out_merged), e.g. for showing target and decoy score distributions. With
'-out_merged' and '-out_parquet', the pooled proteins and protein groups (cross-file inference, protein FDR) are also
written as merged.proteins.parquet and merged.pg.parquet in the -out_parquet directory.
@note Several input files: give all files that share the database and the search parameters to one call ('-in a.mzML b.mzML ...'). The fragment index is then built once instead of once per file, and with more than one thread the next mzML file is read while the current one is searched. Each file gets the same identifications as when it is searched alone, and the whole run usually takes less time than one call per file. Compare search engines in the same mode.
@note Threads: '-threads' sets the OpenMP threads of the search; mzML files are read with a share of them. With many threads on a shared machine, set the environment variable OMP_WAIT_POLICY=passive: idle OpenMP threads then sleep instead of spinning, which saves CPU time at about the same wall time. Do not bind threads (OMP_PROC_BIND), which made ProSE slower.
@note Memory in chunked multi-file runs: '-Search:database:chunk_size' bounds the fragment-index memory only. With multiple '-in' files and chunking active, the chunk-major schedule keeps every input file's preprocessed MS2 spectra in memory for the whole search (each chunk's index is built once and scored against all files). Budget roughly the sum of all files' MS2 peak data on top of one chunk's index, or split very large cohorts across separate invocations (see the sharded-FDR workflow below).
@note Deferred / distributed (sharded) FDR: to search shards on separate nodes and control FDR globally afterwards, run each shard with '-Search:FDR:protein' = 0 (the default), optionally with '-Search:FDR:PSM' > 0 for per-run PSM filtering. Per-file outputs retain the full target+decoy set, so you can pool them and apply FDR once downstream — e.g. @ref TOPP_IDMerger &rarr; @ref TOPP_ProteinInference / @ref TOPP_Epifany &rarr; @ref TOPP_FalseDiscoveryRate / @ref TOPP_IDFilter (idXML route) — or run a single ProSE process over all shards with '-out_merged'.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ProSE.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ProSE.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class ProSE :
    public TOPPBase
{
  public:
    ProSE() :
      TOPPBase("ProSE",
        "Annotates bottom-up MS/MS spectra using ProSE.")
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFileList_("in", "<files>", StringList(), "Input spectrum file(s). Multiple files are searched against the same database; the fragment index is built once and reused.");
      setValidFormats_("in", { "mzML",
#ifdef WITH_OPENTIMS
        "d",
#endif
#ifdef WITH_THERMO_RAW
        "raw",
#endif
      });

      registerInputFile_("database", "<file>", "", "Input protein sequence database in FASTA format.");
      setValidFormats_("database", ListUtils::create<std::string>("fasta"));

      registerOutputFileList_("out_idxml", "<files>", StringList(), "Output idXML identification file(s). Must have the same number of entries as -in.", false);
      setValidFormats_("out_idxml", ListUtils::create<std::string>("idXML"));

      registerOutputDir_("out_qpx", "<dir>", "", "Output directory for QPX exchange format Parquet files. Writes per-input runs/<basename>.psm.parquet and runs/<basename>.pg.parquet, and quantms.psm.parquet and quantms.pg.parquet of the whole run (see also -qpx_extended): the folder is one QPX dataset (prefix quantms). The pg files are written only when protein groups were actually inferred, which needs 'FDR:protein' > 0 together with decoys, or '-out_merged' with more than one input; an identification-only run produces no pg file rather than an empty one.", false, true);

      registerOutputDir_("out_parquet", "<dir>", "", "Output directory for OpenMS internal format Parquet files. Writes per-input runs/<basename>.psm/proteins/pg/search_params.parquet, and openms.psm/proteins/pg/search_params.parquet of the whole run: the folder is one dataset (prefix openms). With -out_merged and several inputs also merged.proteins.parquet and merged.pg.parquet (the pooled proteins).", false, true);

      registerOutputFile_("out_merged", "<file>", "", "Optional merged output file containing all PSMs pooled across input files with cross-file protein inference (BasicProteinInferenceAlgorithm) and optional picked-protein FDR (FDR:protein). Per-file outputs (-out_idxml / -out_qpx / -out_parquet) retain run-level information. Only useful with multiple -in files.", false);
      setValidFormats_("out_merged", ListUtils::create<std::string>("idXML"));

      registerOutputFileList_("out_pin", "<files>", StringList(), "Output Percolator input (.pin/.tsv) file(s) for external rescoring. Must have the same number of entries as -in. It holds the PSMs before the rescoring and before any FDR filter, with ProSE's own scores and all reported hits: 'Search:rescore' and 'Search:FDR:*' apply to the other outputs only.", false, true);
      setValidFormats_("out_pin", ListUtils::create<std::string>("tsv"));

      registerOutputDir_("out_mod_analysis_dir", "<dir>", "", "Optional directory to write modification-analysis tables (delta-mass, PTM stats) when running in open-search mode. When set, per-file tables are written using each input file's basename and an additional aggregate table is written across all input files. Has no effect in closed-search mode.", false, true);

      registerOutputFile_("summary_out", "<file>", "", "Optional YAML file capturing the end-of-search report (per-file identification statistics, shared configuration/database/index facts, timing and output manifest) for pipeline ingestion. The same information is printed to the command line unless -no_summary is set.", false, true);
      setValidFormats_("summary_out", ListUtils::create<std::string>("yaml"));

      registerFlag_("no_summary", "Suppress the end-of-search report printed to the command line. The -summary_out YAML file (if requested) is still written.", true);

      registerFlag_("qpx_extended", "With -out_qpx: add optional columns to the QPX psm view (all scores in 'additional_scores', 'spectrum_reference', 'isotope_error', 'protein_mappings'; in the per-input files also PeptDeep's predicted b/y intensities) and write the QPX mz view of the spectra of each input's PSMs (filtered peaks, see -qpx_peak_top) and the QPX run, sample and provenance views (quantms.run/sample/provenance.parquet: the inputs with their statistics, the steps of the search with their parameters); a run without protein inference gets one protein group per protein. See the tool documentation.", true);
      registerDoubleOption_("qpx_peak_window", "<Th>", 100.0, "-qpx_extended: the m/z window of the peak filter of the mz view.", false, true);
      setMinFloat_("qpx_peak_window", 1.0);
      registerIntOption_("qpx_peak_top", "<n>", 100, "-qpx_extended: the most intense peaks kept per -qpx_peak_window (with every annotated peak); 0 keeps all peaks.", false, true);
      setMinInt_("qpx_peak_top", 0);
      registerFlag_("keep_decoys", "Keep the decoy PSMs and proteins that pass the FDR thresholds in the outputs; the protein FDR (one input, or -out_merged) otherwise removes them.", true);


      // put search algorithm parameters at Search: subtree of parameters
      Param search_algo_params_with_subsection;
      search_algo_params_with_subsection.insert("Search:", ProSEAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    /// The spectrum of one PSM row of -qpx_extended: its prediction columns (null when @p has_prediction is false), and
    /// whether the mz view has the row of its spectrum (else its 'spectrum_reference' is null)
    struct PSMSpectrum_
    {
      bool has_prediction = false;
      std::vector<std::string> ion_type;
      std::vector<float> mz, intensity;
      /// the mz view has a row of its spectrum, under its native ID; if not (the re-read failed, the spectrum was not
      /// found again, or another spectrum of its native ID has the row), its 'spectrum_reference' is null: no stored
      /// spectrum, as a join on its native ID would give nothing or another spectrum's peaks
      bool has_mz_row = false;
    };

    /// The rows per row group of the per-input psm and mz files of -qpx_extended: a viewer reads the peaks and the
    /// prediction of one PSM by decoding its row group (with ~600 peaks per spectrum about 1.2 MB of peaks)
    static constexpr Int64 QPX_SPECTRA_ROW_GROUP = 256;

    /// The peak filter of -qpx_extended: the @p top most intense peaks per @p window Th (0: all), and the annotated ones
    struct PeakFilter_
    {
      double window = 100.0;
      Size top = 100;
      /// The fragment annotation is of deisotoped spectra ('Search:fragment:deisotope', a high-resolution fragment
      /// tolerance): the deisotoper moved each fragment of charge z (up to @p max_charge, or the precursor's charge with
      /// @p cap_at_precursor) to its 1+ m/z, which is no peak of the spectrum as the file has it.
      bool deisotoped = false;
      int max_charge = 3;
      bool cap_at_precursor = true;
    };

    /// The deisotoping of ProSE's search for PeakFilter_ (as ProSEAlgorithm decides it, from its parameters @p search_params)
    static void setDeisotoping_(const Param& search_params, PeakFilter_& filter)
    {
      const auto value = [&](const std::string& key, const std::string& fallback)
      {
        return search_params.exists(key) ? search_params.getValue(key).toString() : fallback;
      };
      const double tolerance = search_params.exists("fragment:mass_tolerance")
                                 ? static_cast<double>(search_params.getValue("fragment:mass_tolerance")) : 20.0;
      filter.deisotoped = value("fragment:deisotope", "auto") != "false"
                          && Deisotoper::isToleranceSupported(tolerance, value("fragment:mass_tolerance_unit", "ppm") == "ppm");
      filter.cap_at_precursor = value("fragment:deisotope_charge_cap", "precursor") == "precursor";
    }

    /**
      @brief The indices of the peaks of @p spectrum (sorted by m/z) that the filter keeps: the most intense ones per
      window (jumping windows, as WindowMower's but with the full quota in the last one) and those of the fragment
      annotation @p annotations. All peaks if @p filter.top is 0.

      An annotated m/z keeps the peak nearest to it (within 0.01 Th). Of deisotoped spectra (@p filter.deisotoped) it
      also keeps, for each charge z from 2 to the deisotoper's maximum (3, or @p precursor_charge if lower and
      @p filter.cap_at_precursor; @p precursor_charge is the spectrum's recorded one, as in the search's preprocessing,
      0 if unknown), the peak nearest to (m/z + (z - 1) * proton) / z (within 0.01 Th): the monoisotopic
      peak of a fragment of charge z, which the deisotoper moved to the annotated 1+ m/z. (That conversion is exact, so
      the raw peak lies at that m/z.) Every annotated fragment thus keeps its peak in the spectrum as the file has it; a
      1+ fragment may keep an unrelated peak at its 2+ or 3+ position too.
    */
    static std::vector<Size> keptPeaks_(const MSSpectrum& spectrum, const std::vector<PeptideHit::PeakAnnotation>& annotations,
                                        int precursor_charge, const PeakFilter_& filter)
    {
      std::vector<bool> keep(spectrum.size(), filter.top == 0);
      auto indexOf = [&spectrum](double mz) -> Size // the peak nearest to mz
      {
        const auto above = std::lower_bound(spectrum.begin(), spectrum.end(), mz,
                                            [](const Peak1D& p, double m) { return p.getMZ() < m; });
        Size i = static_cast<Size>(above - spectrum.begin());
        if (i == spectrum.size() || (i > 0 && mz - spectrum[i - 1].getMZ() < above->getMZ() - mz)) { --i; }
        return i;
      };
      if (filter.top > 0 && !spectrum.empty())
      {
        // Jumping windows as WindowMower's (a window starts at its first peak; the next at the first peak past it),
        // but the last window keeps the full quota too: WindowMower scales it by the window's width, which would
        // drop the highest-m/z peaks of nearly every spectrum.
        std::vector<Size> order;
        for (Size begin = 0, end = 0; begin < spectrum.size(); begin = end)
        {
          while (end < spectrum.size() && spectrum[end].getMZ() - spectrum[begin].getMZ() < filter.window) { ++end; }
          order.resize(end - begin);
          std::iota(order.begin(), order.end(), begin);
          const Size n = std::min(filter.top, order.size());
          std::partial_sort(order.begin(), order.begin() + n, order.end(), [&spectrum](Size a, Size b)
          {
            return spectrum[a].getIntensity() > spectrum[b].getIntensity() || (spectrum[a].getIntensity() == spectrum[b].getIntensity() && a < b);
          });
          for (Size k = 0; k < n; ++k) { keep[order[k]] = true; }
        }
        int max_charge = filter.deisotoped ? filter.max_charge : 1;
        if (filter.deisotoped && filter.cap_at_precursor && precursor_charge > 0) { max_charge = std::min(max_charge, precursor_charge); }
        for (const PeptideHit::PeakAnnotation& a : annotations)
        {
          for (int z = 1; z <= std::max(1, max_charge); ++z)
          {
            const double mz = (a.mz + (z - 1) * Constants::PROTON_MASS_U) / z;
            const Size i = indexOf(mz);
            if (std::fabs(spectrum[i].getMZ() - mz) <= 0.01) { keep[i] = true; }
          }
        }
      }
      std::vector<Size> kept;
      for (Size i = 0; i < keep.size(); ++i) { if (keep[i]) { kept.push_back(i); } }
      return kept;
    }

    /// The spectra and predictions of the PSMs of one searched file (collectPSMSpectra_())
    struct SpectraOfFile_
    {
      /// the QPX mz view of the file's MS2 spectra of its PSMs (null if they could not be read)
      std::shared_ptr<arrow::Table> mz;
      /// the predictions: one entry per element of the file's peptide identifications
      std::vector<PSMSpectrum_> psms;
    };

    /**
      @brief The QPX mz view of the MS2 spectra of the PSMs of one searched file, with filtered peaks, and the PeptDeep
      predictions of its PSMs (-qpx_extended)

      Never throws because of the spectra or the predictions: a failed re-read of @p in_file or a failed prediction is
      logged as a warning and leaves no mz view or null prediction columns; the search result of the file is not
      affected. Only std::bad_alloc propagates: out of memory is not taken for an unreadable file.

      The spectra are read again from @p in_file (MS2 spectra, by native ID; of several spectra with a PSM's native ID
      the one at its RT). The mz view has one row per spectrum, in the order of the first PSM of each (in the order of
      @p peptide_ids), its peaks filtered with @p filter (keeping the annotated peaks of every PSM of the spectrum); a
      PSM whose spectrum is not found has none, nor has one whose spectrum shares its native ID with an earlier spectrum
      (the view's key is the native ID). Only a PSM whose spectrum has the row of its native ID is marked
      PSMSpectrum_::has_mz_row; the others (also all after a failed re-read) get a null spectrum_reference. The
      predictions are PeptDeep's at the NCE of the file: the one PeptDeep rescoring stored on the run ('peptdeep_nce'),
      else 'peptdeep:nce' if set, else the one PeptDeepRescoring::annotate() selects on a copy of @p peptide_ids; the
      NCE is then stored on the run.
    */
    SpectraOfFile_ collectPSMSpectra_(const std::string& in_file, const Param& search_params,
                                      std::vector<ProteinIdentification>& protein_ids,
                                      const PeptideIdentificationList& peptide_ids, const PeakFilter_& filter)
    {
      SpectraOfFile_ result;
      std::vector<PSMSpectrum_>& out = result.psms;
      out.resize(peptide_ids.size());

      // -- peaks: the MS2 spectra of the PSMs, as the file has them --
      std::unordered_map<std::string, std::vector<Size>> psms_of_spectrum; // native ID -> PSMs
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (!peptide_ids[i].getHits().empty()) { psms_of_spectrum[peptide_ids[i].getSpectrumReference()].push_back(i); }
      }
      PeakMap spectra; // the spectra of the PSMs (several, if a native ID is not unique in the file)
      bool read = true;
      try
      {
        MSDataTransformingConsumer keep;
        keep.setSpectraProcessingFunc([&](MSSpectrum& spectrum)
        {
          if (spectrum.getMSLevel() == 2 && psms_of_spectrum.count(spectrum.getNativeID()) > 0)
          {
            if (!spectrum.isSorted()) { spectrum.sortByPosition(); }
            spectra.addSpectrum(spectrum);
          }
        });
        FileHandler fh;
        fh.getOptions().setMSLevels({2});
        fh.transformExperiment(in_file, &keep, {}, log_type_);
      }
      catch (const std::bad_alloc&)
      {
        throw; // out of memory is no unreadable file: the caller fails the input
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_WARN << "-qpx_extended: cannot read the spectra of " << in_file << " again (" << e.what()
                        << "); no mz view is written for it, and the spectrum_reference of its PSMs is null." << endl;
        spectra.clear(true);
        read = false;
      }
      std::unordered_map<std::string, std::vector<Size>> spectra_of_id;
      for (Size s = 0; s < spectra.size(); ++s) { spectra_of_id[spectra[s].getNativeID()].push_back(s); }
      // the spectrum of PSM i: by its native ID, of several the one at its RT; null if none
      auto spectrumOf = [&](Size i) -> const MSSpectrum*
      {
        const auto found = spectra_of_id.find(peptide_ids[i].getSpectrumReference());
        if (found == spectra_of_id.end()) { return nullptr; }
        if (found->second.size() == 1) { return &spectra[found->second.front()]; }
        for (Size s : found->second)
        {
          if (std::fabs(spectra[s].getRT() - peptide_ids[i].getRT()) <= 1e-6) { return &spectra[s]; }
        }
        return nullptr;
      };

      // -- the mz view: a row per spectrum, in the order of its first PSM --
      using Column = ArrowIOHelpers::ParquetColumn;
      using Value = Column::StructValue;
      auto column = [](const std::string& name, Column::Type type, bool nullable = true)
      {
        Column c; c.name = name; c.type = type; c.nullable = nullable; return c;
      };
      Column id = column("id", Column::Type::UTF8, false);
      Column run_file = column("run_file_name", Column::Type::UTF8);
      Column scan = column("scan", Column::Type::INT32_LIST);
      Column ms_level = column("ms_level", Column::Type::INT32, false);
      Column centroid = column("centroid", Column::Type::BOOL, false);
      Column start_time = column("scan_start_time", Column::Type::FLOAT, false);
      Column mobility = column("inverse_ion_mobility", Column::Type::FLOAT);
      Column injection = column("ion_injection_time", Column::Type::FLOAT);
      Column tic = column("total_ion_current", Column::Type::FLOAT, false);
      Column precursors = column("precursors", Column::Type::STRUCT_LIST);
      precursors.struct_fields = {{"selected_ion_mz", Column::Type::FLOAT, false}, {"selected_ion_charge", Column::Type::INT32},
                                  {"selected_ion_intensity", Column::Type::FLOAT}, {"isolation_window_target", Column::Type::FLOAT},
                                  {"isolation_window_lower", Column::Type::FLOAT}, {"isolation_window_upper", Column::Type::FLOAT},
                                  {"spectrum_ref", Column::Type::UTF8}};
      Column mz = column("mz", Column::Type::FLOAT_LIST, false);
      Column intensity = column("intensity", Column::Type::FLOAT_LIST, false);
      Column cv_params = column("cv_params", Column::Type::STRUCT_LIST);
      cv_params.struct_fields = {{"cv_name", Column::Type::UTF8, false}, {"cv_value", Column::Type::UTF8, false}};
      const std::string run_file_name = ArrowIOHelpers::qpxRunFileName(in_file);
      auto optional_float = [](Column& c, double v, bool valid)
      {
        c.double_values.push_back(valid ? v : 0.0);
        c.null_rows.push_back(!valid);
      };
      // the spectrum with a row of each native ID: one row per native ID ('id' is the mz view's primary key), of
      // several spectra of one native ID (at other RTs) the first's; the PSMs of the others have none
      std::map<std::string, const MSSpectrum*> written;
      Size peaks_read = 0, peaks_kept = 0;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (peptide_ids[i].getHits().empty()) { continue; }
        const MSSpectrum* spectrum = spectrumOf(i);
        if (spectrum == nullptr) { continue; }
        if (!written.emplace(spectrum->getNativeID(), spectrum).second) { continue; } // its native ID has a row
        // its peaks: those the filter keeps for any of its PSMs (each keeps its own annotated peaks), with the
        // deisotoper's charge cap of the search's preprocessing: the spectrum's recorded precursor charge (0: unknown,
        // then up to 3), not the charge of a PSM (a spectrum of unknown charge was deisotoped up to 3+, whatever the
        // charge of its PSMs)
        const int precursor_charge = spectrum->getPrecursors().empty() ? 0 : spectrum->getPrecursors()[0].getCharge();
        std::vector<bool> keep(spectrum->size(), false);
        for (Size j : psms_of_spectrum[peptide_ids[i].getSpectrumReference()])
        {
          if (spectrumOf(j) != spectrum) { continue; }
          const PeptideHit& hit = peptide_ids[j].getHits()[0];
          for (Size k : keptPeaks_(*spectrum, hit.getPeakAnnotations(), precursor_charge, filter)) { keep[k] = true; }
        }
        std::vector<float> peak_mz, peak_intensity;
        for (Size k = 0; k < keep.size(); ++k)
        {
          if (!keep[k]) { continue; }
          peak_mz.push_back(static_cast<float>((*spectrum)[k].getMZ()));
          peak_intensity.push_back((*spectrum)[k].getIntensity());
        }
        peaks_read += spectrum->size();
        peaks_kept += peak_mz.size();

        id.utf8_values.push_back(spectrum->getNativeID());
        run_file.utf8_values.push_back(run_file_name);
        scan.int32_lists.push_back(ArrowIOHelpers::qpxScanComponents(spectrum->getNativeID()));
        ms_level.int64_values.push_back(static_cast<Int64>(spectrum->getMSLevel()));
        SpectrumSettings::SpectrumType type = spectrum->getType(false);
        if (type == SpectrumSettings::SpectrumType::UNKNOWN) { type = spectrum->getType(true); }
        centroid.int64_values.push_back(type == SpectrumSettings::SpectrumType::CENTROID ? 1 : 0);
        start_time.double_values.push_back(spectrum->getRT() / 60.0);
        optional_float(mobility, spectrum->getDriftTime(),
                       spectrum->getDriftTimeUnit() == DriftTimeUnit::VSSC && spectrum->getDriftTime() >= 0.0);
        const bool has_injection = !spectrum->getAcquisitionInfo().empty()
                                   && spectrum->getAcquisitionInfo()[0].metaValueExists("MS:1000927");
        optional_float(injection, has_injection ? static_cast<double>(spectrum->getAcquisitionInfo()[0].getMetaValue("MS:1000927")) : 0.0,
                       has_injection);
        tic.double_values.push_back(spectrum->calculateTIC());
        std::vector<Column::StructEntry> ions;
        for (const Precursor& p : spectrum->getPrecursors())
        {
          ions.push_back({p.getMZ(), p.getCharge() > 0 ? Value(static_cast<Int32>(p.getCharge())) : Value(),
                          p.getIntensity() > 0 ? Value(static_cast<double>(p.getIntensity())) : Value(), Value(),
                          p.getIsolationWindowLowerOffset() > 0 ? Value(p.getIsolationWindowLowerOffset()) : Value(),
                          p.getIsolationWindowUpperOffset() > 0 ? Value(p.getIsolationWindowUpperOffset()) : Value(), Value()});
        }
        precursors.struct_lists.push_back(std::move(ions));
        mz.float_lists.push_back(std::move(peak_mz));
        intensity.float_lists.push_back(std::move(peak_intensity));
        cv_params.struct_lists.emplace_back();
        cv_params.null_rows.push_back(true);
      }
      const Size mz_rows = id.utf8_values.size();
      if (read)
      {
        result.mz = ArrowIOHelpers::tableOfColumns({id, run_file, scan, ms_level, centroid, start_time, mobility, injection, tic,
                                                    precursors, mz, intensity, cv_params}, mz_rows);
      }
      // Each PSM against the rows, once all are made (a later PSM may make the row of a native ID): a PSM has its
      // spectrum's row only if the row of its native ID is that spectrum's. Otherwise (no spectrum of its native ID and RT,
      // no mz view, or another spectrum of its native ID has the row) its spectrum_reference becomes null: a join on its
      // native ID would give no peaks or another spectrum's
      Size without_spectrum = 0, unresolved = 0;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (peptide_ids[i].getHits().empty()) { continue; }
        const auto row = written.find(peptide_ids[i].getSpectrumReference());
        if (row == written.end()) { ++without_spectrum; } // no spectrum of its native ID (and RT), and no row to join
        else if (row->second != spectrumOf(i)) { ++unresolved; }
        else { out[i].has_mz_row = result.mz != nullptr; }
      }
      if (without_spectrum > 0 && read)
      {
        OPENMS_LOG_WARN << "-qpx_extended: " << without_spectrum << " PSM(s) of " << in_file
                        << " have no MS2 spectrum of their native ID (and RT) in the file; the mz view has none for them,"
                        << " and their spectrum_reference is null." << endl;
      }
      if (unresolved > 0)
      {
        OPENMS_LOG_WARN << "-qpx_extended: " << unresolved << " PSM(s) of " << in_file << " are of a spectrum that shares its"
                        << " native ID with another spectrum (at another RT), which has the mz row of that native ID: the mz"
                        << " view has none for them, and their spectrum_reference is null (unresolved)." << endl;
      }
      if (read && !result.mz)
      {
        OPENMS_LOG_WARN << "-qpx_extended: the mz view of " << in_file << " could not be built (see the error above); no mz"
                        << " view is written for it, and the spectrum_reference of its PSMs is null." << endl;
      }

      // -- predictions --
#ifdef WITH_ONNX
      try
      {
        predictPSMSpectra_(in_file, search_params, protein_ids, peptide_ids, spectra, out);
      }
      catch (const std::bad_alloc&)
      {
        throw; // out of memory is no failed prediction: the caller fails the input
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_WARN << "-qpx_extended: PeptDeep could not predict the spectra of " << in_file << " (" << e.what()
                        << "); the predicted_* columns of its PSMs are null." << endl;
        for (PSMSpectrum_& entry : out)
        {
          entry.has_prediction = false;
          entry.ion_type.clear(); entry.mz.clear(); entry.intensity.clear();
        }
      }
#else
      (void)search_params;
      (void)protein_ids;
      static bool told = false;
      if (!told)
      {
        OPENMS_LOG_INFO << "-qpx_extended: this OpenMS was built without ONNX Runtime, so PeptDeep cannot predict; "
                        << "the spectra are written and the predicted_* columns are null." << endl;
        told = true;
      }
#endif
      OPENMS_LOG_INFO << "-qpx_extended: " << in_file << ": " << mz_rows << " spectra of its PSMs, " << peaks_kept
                      << " of " << peaks_read << " peaks kept";
      if (filter.top > 0)
      {
        OPENMS_LOG_INFO << " (the " << filter.top << " most intense per " << filter.window << " Th and the annotated ones)";
      }
      OPENMS_LOG_INFO << endl;
      return result;
    }

#ifdef WITH_ONNX
    /// The predictions of collectPSMSpectra_(), into @p out. Throws on failure (the caller logs it).
    void predictPSMSpectra_(const std::string& in_file, const Param& search_params,
                            std::vector<ProteinIdentification>& protein_ids,
                            const PeptideIdentificationList& peptide_ids, const PeakMap& spectra,
                            std::vector<PSMSpectrum_>& out)
    {
      if (protein_ids.empty()) { return; }
      // As ProSEAlgorithm resolves the models of 'peptdeep:enable' (relative to share/OpenMS of the build or install).
      const StringList model_dirs = {std::string(OPENMS_BINARY_PATH) + "/share/OpenMS"};
      const std::string ms2_model = File::find(search_params.getValue("peptdeep:ms2_model").toString(), model_dirs);
#ifdef _OPENMP
      const int threads = std::max(1, omp_get_max_threads());
#else
      const int threads = 1;
#endif
      const int batch_size = search_params.getValue("peptdeep:batch_size");
      const std::string instrument = search_params.getValue("peptdeep:instrument").toString();

      // The NCE: the one PeptDeep rescoring used, else 'peptdeep:nce' if set, else the one PeptDeepRescoring::annotate()
      // selects, run on a copy of the PSMs only for that (its features and RT calibration are discarded).
      ProteinIdentification& run = protein_ids[0];
      std::string nce_source = "from PeptDeep rescoring";
      const double fixed_nce = search_params.getValue("peptdeep:nce");
      if (!run.metaValueExists(PeptDeepRescoring::NCE_META_VALUE) && fixed_nce > 0.0)
      {
        run.setMetaValue(PeptDeepRescoring::NCE_META_VALUE, fixed_nce);
        nce_source = "set";
      }
      else if (!run.metaValueExists(PeptDeepRescoring::NCE_META_VALUE))
      {
        PeptDeepRescoring rescoring;
        Param p = rescoring.getParameters();
        p.setValue("ms2_model", ms2_model);
        p.setValue("rt_model", File::find(search_params.getValue("peptdeep:rt_model").toString(), model_dirs));
        p.setValue("instrument", instrument);
        p.setValue("rt_model_type", search_params.getValue("peptdeep:rt_model_type"));
        p.setValue("batch_size", batch_size);
        p.setValue("threads", threads);
        rescoring.setParameters(p);
        rescoring.setLogType(log_type_);
        std::vector<ProteinIdentification> runs = protein_ids;
        PeptideIdentificationList psms = peptide_ids;
        OPENMS_LOG_INFO << "-qpx_extended: selecting the NCE of " << in_file
                        << " with PeptDeep rescoring on a copy of its PSMs" << endl;
        rescoring.annotate(spectra, runs, psms);
        if (!runs[0].metaValueExists(PeptDeepRescoring::NCE_META_VALUE))
        {
          OPENMS_LOG_WARN << "-qpx_extended: no NCE selected for " << in_file << "; the predicted_* columns are null." << endl;
          return;
        }
        run.setMetaValue(PeptDeepRescoring::NCE_META_VALUE, runs[0].getMetaValue(PeptDeepRescoring::NCE_META_VALUE));
        nce_source = "selected";
      }
      const double nce = run.getMetaValue(PeptDeepRescoring::NCE_META_VALUE);

      // The distinct (peptide, charge) pairs of the rows, of the charges the model predicts (1 to 6).
      std::map<std::pair<std::string, int>, Size> key_of;
      std::vector<std::string> peptides;
      std::vector<float> charges;
      std::vector<Size> row_key(peptide_ids.size(), std::numeric_limits<Size>::max());
      Size outside_charges = 0;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (peptide_ids[i].getHits().empty()) { continue; }
        const PeptideHit& hit = peptide_ids[i].getHits()[0];
        const int z = hit.getCharge();
        if (z < 1 || z > 6) { ++outside_charges; continue; }
        if (hit.getSequence().size() < 2) { continue; }
        const auto [it, inserted] = key_of.emplace(std::make_pair(hit.getSequence().toString(), z), peptides.size());
        if (inserted)
        {
          peptides.push_back(it->first.first);
          charges.push_back(static_cast<float>(z));
        }
        row_key[i] = it->second;
      }
      if (outside_charges > 0)
      {
        OPENMS_LOG_INFO << "-qpx_extended: " << outside_charges << " PSM(s) of " << in_file
                        << " have a precursor charge outside 1 to 6, which PeptDeep does not predict; their predicted_* columns are null." << endl;
      }
      if (peptides.empty()) { return; }

      PeptDeepMS2Inference ms2(ms2_model, threads, static_cast<size_t>(batch_size));
      const int64_t instrument_index = ML::instrumentIndex(instrument);
      std::vector<std::vector<float>> predicted;
      try
      {
        predicted = ms2.predictMS2(peptides, charges, std::vector<float>(peptides.size(), static_cast<float>(nce)),
                                   std::vector<int64_t>(peptides.size(), instrument_index));
      }
      catch (const Exception::BaseException&)
      {
        // A peptide the model cannot take fails the whole call: predict one by one, and leave out those that fail.
        predicted.assign(peptides.size(), {});
        Size failed = 0;
        for (Size k = 0; k < peptides.size(); ++k)
        {
          try
          {
            predicted[k] = ms2.predictMS2({peptides[k]}, {charges[k]}, {static_cast<float>(nce)}, {instrument_index}).at(0);
          }
          catch (const Exception::BaseException&) { ++failed; }
        }
        OPENMS_LOG_WARN << "-qpx_extended: PeptDeep could not predict " << failed << " of " << peptides.size()
                        << " peptides of " << in_file << "; their predicted_* columns are null." << endl;
      }

      // The b and y ions of charge 1 (and 2 for precursors of charge 2 or more) of each pair, relative to the most
      // intense of them, without those of 0.01 or less, with their m/z, by m/z.
      TheoreticalSpectrumGenerator generator;
      Param tsg = generator.getParameters();
      tsg.setValue("add_b_ions", "true");
      tsg.setValue("add_y_ions", "true");
      tsg.setValue("add_first_prefix_ion", "true");
      tsg.setValue("add_metainfo", "true");
      generator.setParameters(tsg);
      constexpr Size channels = 8; // [b_z1, b_z2, y_z1, y_z2, b_modloss_z1, b_modloss_z2, y_modloss_z1, y_modloss_z2]
      std::vector<PSMSpectrum_> ions_of_key(peptides.size());
      Size other_layout = 0;
      for (Size k = 0; k < peptides.size(); ++k)
      {
        const AASequence sequence = AASequence::fromString(peptides[k]);
        const Size length = sequence.size();
        const std::vector<float>& flat = predicted[k];
        if (flat.empty()) { continue; } // not predicted (counted above)
        if (flat.size() != channels * (length - 1)) { ++other_layout; continue; }
        const int z = static_cast<int>(charges[k]);
        // (type, ordinal, charge) -> intensity; cleavage `pos` carries b(pos + 1) and y(length - 1 - pos)
        std::map<std::tuple<char, int, int>, double> intensity;
        double base_peak = 0.0;
        for (Size pos = 0; pos + 1 < length; ++pos)
        {
          const int b = static_cast<int>(pos + 1), y = static_cast<int>(length - 1 - pos);
          intensity[{'b', b, 1}] = flat[pos * channels + 0];
          intensity[{'y', y, 1}] = flat[pos * channels + 2];
          if (z > 1)
          {
            intensity[{'b', b, 2}] = flat[pos * channels + 1];
            intensity[{'y', y, 2}] = flat[pos * channels + 3];
          }
        }
        for (const auto& [ion, value] : intensity) { base_peak = std::max(base_peak, value); }
        if (!(base_peak > 0.0)) { continue; }

        PeakSpectrum theoretical;
        generator.getSpectrum(theoretical, sequence, 1, z > 1 ? 2 : 1);
        const auto& names = theoretical.getStringDataArrays().at(0);
        const auto& ion_charges = theoretical.getIntegerDataArrays().at(0);
        PSMSpectrum_& ions = ions_of_key[k];
        ions.has_prediction = true;
        for (Size t = 0; t < theoretical.size(); ++t)
        {
          const std::string& name = names[t]; // e.g. "b3+", "y12++"
          Size end = 1;
          while (end < name.size() && name[end] >= '0' && name[end] <= '9') { ++end; }
          if (end == 1) { continue; }
          const auto found = intensity.find({name[0], std::stoi(name.substr(1, end - 1)), ion_charges[t]});
          if (found == intensity.end()) { continue; }
          const double relative = found->second / base_peak;
          if (!(relative > 0.01)) { continue; }
          ions.ion_type.push_back(name);
          ions.mz.push_back(static_cast<float>(theoretical[t].getMZ()));
          ions.intensity.push_back(static_cast<float>(relative));
        }
      }
      if (other_layout > 0)
      {
        OPENMS_LOG_WARN << "-qpx_extended: the MS2 model gave " << other_layout << " of " << peptides.size()
                        << " peptides of " << in_file << " another number of values than " << channels
                        << " per cleavage (an unexpected model layout); their predicted_* columns are null." << endl;
      }
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (row_key[i] == std::numeric_limits<Size>::max()) { continue; }
        const PSMSpectrum_& ions = ions_of_key[row_key[i]];
        out[i].has_prediction = ions.has_prediction;
        out[i].ion_type = ions.ion_type;
        out[i].mz = ions.mz;
        out[i].intensity = ions.intensity;
      }
      OPENMS_LOG_INFO << "-qpx_extended: " << in_file << ": PeptDeep predictions of " << peptides.size()
                      << " peptide(s) at NCE " << nce << " (" << nce_source << ")" << endl;
    }
#endif

    /// A score or meta value name in the snake_case of QPX 'additional_scores' ("ln(hyperscore)" -> "ln_hyperscore",
    /// "DeltaMass" -> "delta_mass"); the scores of Percolator by their names of bigbio/qpx#355; "" for the PEP, which
    /// has its own column
    static std::string qpxScoreName_(const std::string& name)
    {
      static const std::map<std::string, std::string> known = {
        {"q-value", "global_qvalue"},
        {"MS:1001491", "percolator_qvalue"}, {"percolator_q_value", "percolator_qvalue"},
        {"MS:1001492", "percolator_score"}, {"svm", "percolator_score"},
        {"MS:1001493", "percolator_pep"},
        {"ln(hyperscore)", "ln_hyperscore"}, {"ln(hyperscore)_score", "ln_hyperscore"},
        {"Posterior Error Probability", ""}, {"Posterior Error Probability_score", ""},
        {"IsotopeError", "isotope_error"}};
      const auto found = known.find(name);
      if (found != known.end()) { return found->second; }
      std::string out;
      for (Size i = 0; i < name.size(); ++i)
      {
        const unsigned char c = static_cast<unsigned char>(name[i]);
        const unsigned char before = i > 0 ? static_cast<unsigned char>(name[i - 1]) : 0;
        if (std::isupper(c) && (std::islower(before) || std::isdigit(before))) { out += '_'; }
        if (std::isalnum(c)) { out += static_cast<char>(std::tolower(c)); }
        else if (!out.empty() && out.back() != '_') { out += '_'; }
      }
      while (!out.empty() && out.back() == '_') { out.pop_back(); }
      return out;
    }

    /// Whether a higher value of the score @p name (a qpxScoreName_()) is better; none if not known
    static std::optional<bool> higherBetter_(const std::string& name)
    {
      if (name == "global_qvalue" || name == "percolator_qvalue" || name == "percolator_pep") { return false; }
      if (name == "percolator_score" || name == "ln_hyperscore") { return true; }
      return std::nullopt;
    }

    /// Whether ProSE's 'isotope_error' meta values of the run @p protein_ids are observed minus theoretical (the record
    /// 'isotope_error_convention' of its search parameters; without it they are theoretical minus observed)
    static bool isotopeObservedMinusTheoretical_(const std::vector<ProteinIdentification>& protein_ids)
    {
      if (protein_ids.empty()) { return true; } // no PSMs: no isotope errors
      const ProteinIdentification::SearchParameters& sp = protein_ids[0].getSearchParameters();
      return sp.metaValueExists("isotope_error_convention")
             && sp.getMetaValue("isotope_error_convention").toString() == "observed_minus_theoretical";
    }

    /// The columns -qpx_extended adds to the QPX psm view of a file (extendedColumns_()): those of every QPX file of
    /// the run, and the spectrum columns of the per-input files only
    struct ExtendedColumns_
    {
      std::vector<ArrowIOHelpers::ParquetColumn> rows;
      std::vector<ArrowIOHelpers::ParquetColumn> spectra;
    };

    /**
      @brief The columns -qpx_extended adds to the QPX psm view of a file, one entry per row of the view (the best hit
      of each element of @p peptide_ids with hits): 'additional_scores' (replaced: every numeric score and meta value
      of the PSM in snake_case), 'spectrum_reference' (null without PSMSpectrum_::has_mz_row), 'isotope_error',
      'protein_mappings' (ExtendedColumns_::rows);
      the peaks and the predictions of @p spectra (ExtendedColumns_::spectra), which are moved out of @p spectra.
    */
    static ExtendedColumns_ extendedColumns_(const std::vector<ProteinIdentification>& protein_ids,
                                             const PeptideIdentificationList& peptide_ids,
                                             std::vector<PSMSpectrum_>& spectra)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
      Column scores = column(QPXPSMSchema::ADDITIONAL_SCORES, Column::Type::SCORE_LIST);
      Column reference = column("spectrum_reference", Column::Type::UTF8);
      Column isotope = column("isotope_error", Column::Type::INT16);
      Column mappings = column("protein_mappings", Column::Type::STRUCT_LIST);
      mappings.struct_fields = {{"accession", Column::Type::UTF8, false}, {"start", Column::Type::INT32}, {"end", Column::Type::INT32},
                                {"aa_before", Column::Type::UTF8}, {"aa_after", Column::Type::UTF8}};
      Column ion_type = column("predicted_ion_type_array", Column::Type::UTF8_LIST);
      Column mz = column("predicted_mz_array", Column::Type::FLOAT_LIST);
      Column intensity = column("predicted_intensity_array", Column::Type::FLOAT_LIST);

      // The isotope error as QPX has it (bigbio/qpx#355: the 13C peaks from the monoisotopic to the selected peak, i.e.
      // observed minus theoretical), in the column and in 'additional_scores'; ProSE's meta value follows
      // 'report:isotope_error_convention' (the provenance parameter 'isotope_error_convention' says which).
      const int isotope_sign = isotopeObservedMinusTheoretical_(protein_ids) ? 1 : -1;
      // meta values that are no scores (columns of their own, or attributes of the spectrum)
      static const std::set<std::string> not_scores = {"target_decoy", "predicted_RT", "predicted_rt", "ion_mobility", "IM",
                                                       "FAIMS_CV", "scan", "reference_file_name", "id_merge_index"};
      auto residue = [](char aa) -> std::optional<std::string>
      {
        if (aa == PeptideEvidence::N_TERMINAL_AA || aa == PeptideEvidence::C_TERMINAL_AA) { return std::string("-"); }
        if (aa == PeptideEvidence::UNKNOWN_AA) { return std::nullopt; }
        return std::string(1, aa);
      };
      auto position = [](Int p) -> std::optional<Int32>
      {
        if (p == PeptideEvidence::UNKNOWN_POSITION || p < 0) { return std::nullopt; }
        return static_cast<Int32>(p + 1); // 1-based
      };

      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        const PeptideIdentification& id = peptide_ids[i];
        if (id.getHits().empty()) { continue; } // no QPX row
        const PeptideHit& hit = id.getHits()[0];

        std::vector<Column::Score> row_scores;
        std::set<std::string> seen;
        auto add = [&](const std::string& key, double value, std::optional<bool> higher_better)
        {
          const std::string name = qpxScoreName_(key);
          if (name.empty() || !seen.insert(name).second) { return; }
          const std::optional<bool> known = higherBetter_(name);
          if (name == "isotope_error" && value != 0.0) { value *= isotope_sign; } // (no -0)
          row_scores.push_back({name, value, known ? known : higher_better});
        };
        add(id.getScoreType(), hit.getScore(), id.isHigherScoreBetter());
        std::vector<std::string> keys;
        hit.getKeys(keys);
        for (const std::string& key : keys)
        {
          if (not_scores.count(key) > 0) { continue; }
          const DataValue& value = hit.getMetaValue(key);
          if (value.valueType() == DataValue::INT_VALUE || value.valueType() == DataValue::DOUBLE_VALUE)
          {
            add(key, static_cast<double>(value), std::nullopt);
          }
        }
        scores.score_lists.push_back(std::move(row_scores));

        // (null: the mz view has no row of its spectrum, no stored spectrum; see collectPSMSpectra_())
        PSMSpectrum_& s = spectra[i];
        reference.utf8_values.push_back(s.has_mz_row ? id.getSpectrumReference() : std::string());
        reference.null_rows.push_back(!s.has_mz_row || id.getSpectrumReference().empty());

        const std::string isotope_key = hit.metaValueExists("isotope_error") ? "isotope_error" : "IsotopeError";
        isotope.int64_values.push_back(hit.metaValueExists(isotope_key) ? isotope_sign * static_cast<int>(hit.getMetaValue(isotope_key)) : 0);
        isotope.null_rows.push_back(!hit.metaValueExists(isotope_key));

        std::vector<Column::StructEntry> row_mappings;
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          row_mappings.push_back({evidence.getProteinAccession(), position(evidence.getStart()), position(evidence.getEnd()),
                                  residue(evidence.getAABefore()), residue(evidence.getAAAfter())});
        }
        mappings.struct_lists.push_back(std::move(row_mappings));

        ion_type.utf8_lists.push_back(std::move(s.ion_type));
        ion_type.null_rows.push_back(!s.has_prediction);
        mz.float_lists.push_back(std::move(s.mz));
        mz.null_rows.push_back(!s.has_prediction);
        intensity.float_lists.push_back(std::move(s.intensity));
        intensity.null_rows.push_back(!s.has_prediction);
      }
      spectra.clear();
      spectra.shrink_to_fit();
      ExtendedColumns_ out;
      out.rows = {std::move(scores), std::move(reference), std::move(isotope), std::move(mappings)};
      out.spectra = {std::move(ion_type), std::move(mz), std::move(intensity)};
      return out;
    }

    /// @p value as a JSON string literal
    static std::string jsonString_(const std::string& value)
    {
      static const char* const hex = "0123456789abcdef";
      std::string out = "\"";
      for (char ch : value)
      {
        const unsigned char c = static_cast<unsigned char>(ch);
        if (c == '"' || c == '\\') { out += '\\'; out += ch; }
        else if (c == '\n') { out += "\\n"; }
        else if (c == '\t') { out += "\\t"; }
        else if (c == '\r') { out += "\\r"; }
        else if (c < 0x20) { out += "\\u00"; out += hex[c >> 4]; out += hex[c & 0x0F]; }
        else { out += ch; }
      }
      return out + "\"";
    }

    /// @p value as a locale-independent number, the shortest of 15 to 17 digits that reads back as @p value; non-finite
    /// as @p non_finite
    static std::string number_(double value, const std::string& non_finite = "null")
    {
      if (!std::isfinite(value)) { return non_finite; }
      std::string text;
      for (int digits = 15; digits <= 17; ++digits)
      {
        std::ostringstream o;
        o.imbue(std::locale::classic());
        o << std::setprecision(digits) << value;
        text = o.str();
        std::istringstream in(text);
        in.imbue(std::locale::classic());
        double back = 0.0;
        if ((in >> back) && back == value) { break; }
      }
      return text;
    }

    /// The footer keys of the mz.parquet files of -qpx_extended: the peak filter of their peaks
    static std::map<std::string, std::string> peakMetadata_(const PeakFilter_& filter)
    {
      return {{"prose_peak_window", number_(filter.window, "0")}, {"prose_peak_top", std::to_string(filter.top)},
              {"prose_peak_annotated", "kept"}};
    }

    /// The footer keys of the per-input psm.parquet files of -qpx_extended: the model of their 'predicted_*_array'
    /// columns; the rest is in the run and provenance views
    static std::map<std::string, std::string> predictionMetadata_(const Param& search_params)
    {
      std::map<std::string, std::string> meta;
      if (search_params.exists("peptdeep:ms2_model"))
      {
        meta["prose_prediction_ms2_model"] = File::basename(search_params.getValue("peptdeep:ms2_model").toString());
      }
      if (search_params.exists("peptdeep:instrument"))
      {
        meta["prose_prediction_instrument"] = search_params.getValue("peptdeep:instrument").toString();
      }
      return meta;
    }

    /// @p value (a Param value) as JSON
    static std::string jsonOfValue_(const ParamValue& value)
    {
      auto list = [](const auto& items, auto&& one)
      {
        std::string out = "[";
        for (Size k = 0; k < items.size(); ++k) { out += (k > 0 ? "," : "") + one(items[k]); }
        return out + "]";
      };
      switch (value.valueType())
      {
        case ParamValue::INT_VALUE: return std::to_string(static_cast<int>(value));
        case ParamValue::DOUBLE_VALUE: return number_(static_cast<double>(value));
        case ParamValue::STRING_LIST: return list(value.toStringVector(), [](const std::string& s) { return jsonString_(s); });
        case ParamValue::INT_LIST: return list(value.toIntVector(), [](int i) { return std::to_string(i); });
        case ParamValue::DOUBLE_LIST: return list(value.toDoubleVector(), [](double d) { return number_(d); });
        case ParamValue::EMPTY_VALUE: return "null";
        default: return jsonString_(value.toString());
      }
    }

    /// @p param as one JSON object of its full parameter names ("Search:precursor:mass_tolerance_lower") and values
    static std::string jsonOfParam_(const Param& param)
    {
      std::string out = "{";
      bool first = true;
      for (auto it = param.begin(); it != param.end(); ++it)
      {
        out += (first ? "" : ",") + jsonString_(it.getName()) + ":" + jsonOfValue_(it->value);
        first = false;
      }
      return out + "}";
    }

    /// The QPX MODIFICATION struct fields (run view: modification_parameters)
    static std::vector<ArrowIOHelpers::ParquetColumn::StructField> modificationFields_()
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      return {{"accession", Column::Type::UTF8, false}, {"name", Column::Type::UTF8, true}, {"fixed", Column::Type::BOOL, true},
              {"position", Column::Type::UTF8, true}, {"target_amino_acid", Column::Type::UTF8, true}};
    }

    /// The modifications of the search (OpenMS names, e.g. "Carbamidomethyl (C)") as QPX MODIFICATION structs
    static std::vector<ArrowIOHelpers::ParquetColumn::StructEntry> modificationParameters_(const std::vector<std::string>& fixed,
                                                                                            const std::vector<std::string>& variable)
    {
      static const std::map<ResidueModification::TermSpecificity, std::string> positions = {
        {ResidueModification::ANYWHERE, "Anywhere"}, {ResidueModification::N_TERM, "Any N-term"},
        {ResidueModification::C_TERM, "Any C-term"}, {ResidueModification::PROTEIN_N_TERM, "Protein N-term"},
        {ResidueModification::PROTEIN_C_TERM, "Protein C-term"}};
      using Value = ArrowIOHelpers::ParquetColumn::StructValue;
      std::vector<ArrowIOHelpers::ParquetColumn::StructEntry> out;
      for (const auto& [names, is_fixed] : {std::make_pair(&fixed, true), std::make_pair(&variable, false)})
      {
        for (const std::string& name : *names)
        {
          const ResidueModification* mod = nullptr;
          try { mod = ModificationsDB::getInstance()->getModification(name); }
          catch (const Exception::BaseException&) { mod = nullptr; }
          if (mod == nullptr)
          {
            out.push_back({name, name, is_fixed, Value(), Value()});
            continue;
          }
          std::string accession = mod->getUniModAccession(); // "UniMod:4"
          if (StringUtils::hasPrefix(accession, "UniMod:")) { accession = "UNIMOD:" + accession.substr(7); }
          if (accession.empty()) { accession = mod->getFullId(); }
          const auto position = positions.find(mod->getTermSpecificity());
          const char origin = mod->getOrigin();
          out.push_back({accession, mod->getId(), is_fixed,
                         position != positions.end() ? Value(position->second) : Value(),
                         std::isupper(static_cast<unsigned char>(origin)) && origin != 'X' ? Value(std::string(1, origin)) : Value()});
        }
      }
      return out;
    }

    /// The names of ProSE's exit codes (ProSEAlgorithm::ExitCodes), the run view's 'search_status'
    static std::string searchStatus_(ProSEAlgorithm::ExitCodes code)
    {
      switch (code)
      {
        case ProSEAlgorithm::ExitCodes::EXECUTION_OK: return "EXECUTION_OK";
        case ProSEAlgorithm::ExitCodes::INPUT_FILE_EMPTY: return "INPUT_FILE_EMPTY";
        case ProSEAlgorithm::ExitCodes::UNEXPECTED_RESULT: return "UNEXPECTED_RESULT";
        case ProSEAlgorithm::ExitCodes::ILLEGAL_PARAMETERS: return "ILLEGAL_PARAMETERS";
        default: return "UNKNOWN_ERROR";
      }
    }

    /**
      @brief The QPX run view of -qpx_extended (quantms.run.parquet): one row per input, in input order, with the
      columns of the spec (run_accession = run_file_name: the input's stem, file_name, samples: one sample of the
      stem with label LFQ, enzymes, modification_parameters) and ProSE's optional ones (see the tool documentation)
    */
    static std::shared_ptr<arrow::Table> runTable_(const ProSEAlgorithm::MultiFileSearchResult& mfres,
                                                   const StringList& in_list, const std::vector<Int64>& psm_rows, bool rescore)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type, bool nullable = true)
      {
        Column c; c.name = name; c.type = type; c.nullable = nullable; return c;
      };
      Column accession = column("run_accession", Column::Type::UTF8, false);
      Column run_file = column("run_file_name", Column::Type::UTF8, false);
      Column file_name = column("file_name", Column::Type::UTF8);
      Column samples = column("samples", Column::Type::STRUCT_LIST, false);
      samples.struct_fields = {{"sample_accession", Column::Type::UTF8, false}, {"label", Column::Type::UTF8, false},
                               {"biological_replicate", Column::Type::INT32, true}, {"technical_replicate", Column::Type::INT32, true}};
      // the spec's columns that ProSE cannot fill: present (QPX's validator requires them), every value null
      Column fraction = column("fraction", Column::Type::UTF8);
      Column instrument = column("instrument", Column::Type::UTF8);
      Column dissociation = column("dissociation_method", Column::Type::UTF8);
      Column acquisition = column("acquisition_method", Column::Type::UTF8);
      Column terms = column("additional_terms", Column::Type::STRUCT_LIST);
      terms.struct_fields = {{"accession", Column::Type::UTF8, false}, {"name", Column::Type::UTF8},
                             {"properties", Column::Type::STRUCT_LIST, true,
                              {{"key", Column::Type::UTF8, false}, {"value", Column::Type::UTF8, false}}}};
      Column enzymes = column("enzymes", Column::Type::UTF8_LIST);
      Column mods = column("modification_parameters", Column::Type::STRUCT_LIST);
      mods.struct_fields = modificationFields_();
      Column status = column("search_status", Column::Type::UTF8);
      Column ms2 = column("ms2_spectra_count", Column::Type::INT64);
      Column identified = column("identified_spectra_count", Column::Type::INT64);
      Column decoys = column("decoy_psm_count", Column::Type::INT64);
      Column rows = column("psm_rows", Column::Type::INT64);
      Column applied = column("rescoring_applied", Column::Type::BOOL);
      Column skipped = column("rescoring_skipped_reason", Column::Type::UTF8);
      Column warning = column("rescoring_warning", Column::Type::UTF8);
      Column nce = column("prediction_nce", Column::Type::DOUBLE);
      Column seconds = column("search_seconds", Column::Type::DOUBLE);
      std::vector<Column> errors;
      for (const char* name : {"precursor_error_median_ppm", "precursor_error_mad_ppm", "precursor_tolerance_recommended_ppm",
                                      "fragment_error_median_ppm", "fragment_error_mad_ppm", "fragment_tolerance_recommended_ppm"})
      {
        errors.push_back(column(name, Column::Type::DOUBLE));
      }
      const ProSEAlgorithm::SharedSearchStats& sh = mfres.shared;
      auto optional_text = [](Column& c, const std::string& s) { c.utf8_values.push_back(s); c.null_rows.push_back(s.empty()); };
      auto optional_number = [](Column& c, bool valid, double v)
      {
        c.double_values.push_back(valid && std::isfinite(v) ? v : 0.0);
        c.null_rows.push_back(!(valid && std::isfinite(v)));
      };
      for (Size i = 0; i < in_list.size(); ++i)
      {
        const std::string stem = File::stemName(in_list[i]);
        accession.utf8_values.push_back(stem);
        run_file.utf8_values.push_back(stem);
        file_name.utf8_values.push_back(File::basename(in_list[i]));
        samples.struct_lists.push_back({{stem, "LFQ", Column::StructValue(), Column::StructValue()}});
        for (Column* unknown : {&fraction, &instrument, &dissociation, &acquisition})
        {
          unknown->utf8_values.emplace_back();
          unknown->null_rows.push_back(true);
        }
        terms.struct_lists.emplace_back();
        terms.null_rows.push_back(true);
        enzymes.utf8_lists.push_back({sh.enzyme});
        enzymes.null_rows.push_back(sh.enzyme.empty());
        mods.struct_lists.push_back(modificationParameters_(sh.fixed_mods, sh.variable_mods));
        const bool searched = i < mfres.per_file.size();
        const ProSEAlgorithm::SearchResult* pf = searched ? &mfres.per_file[i] : nullptr;
        status.utf8_values.push_back(pf ? searchStatus_(pf->exit_code) : "UNKNOWN_ERROR");
        ms2.int64_values.push_back(pf ? static_cast<Int64>(pf->stats.ms2_spectra) : 0);
        identified.int64_values.push_back(pf ? static_cast<Int64>(pf->stats.target_psms) : 0);
        decoys.int64_values.push_back(pf ? static_cast<Int64>(pf->stats.decoy_psms) : 0);
        rows.int64_values.push_back(i < psm_rows.size() ? psm_rows[i] : 0);
        // the rescoring outcome: null without rescoring ('Search:rescore' off)
        applied.int64_values.push_back(pf && pf->percolator.applied ? 1 : 0);
        applied.null_rows.push_back(!rescore);
        optional_text(skipped, pf ? pf->percolator.skipped_reason : std::string());
        optional_text(warning, pf ? pf->percolator.warning : std::string());
        double run_nce = std::numeric_limits<double>::quiet_NaN();
        if (pf && !pf->protein_ids.empty() && pf->protein_ids[0].metaValueExists("peptdeep_nce"))
        {
          run_nce = static_cast<double>(pf->protein_ids[0].getMetaValue("peptdeep_nce"));
        }
        optional_number(nce, true, run_nce);
        seconds.double_values.push_back(pf ? pf->stats.seconds_search : 0.0);
        const ProSEAlgorithm::RunStatistics empty_stats;
        const ProSEAlgorithm::RunStatistics& st = pf ? pf->stats : empty_stats;
        optional_number(errors[0], st.prec_tol_valid, st.prec_err_median);
        optional_number(errors[1], st.prec_tol_valid, st.prec_err_mad);
        optional_number(errors[2], st.prec_tol_valid, st.prec_err_recommended);
        optional_number(errors[3], st.frag_tol_valid, st.frag_err_median);
        optional_number(errors[4], st.frag_tol_valid, st.frag_err_mad);
        optional_number(errors[5], st.frag_tol_valid, st.frag_err_recommended);
      }
      // (the spec's columns in its order, then ProSE's)
      std::vector<Column> columns = {accession, run_file, file_name, samples, fraction, instrument, enzymes, dissociation,
                                     acquisition, terms, mods, status, ms2, identified, decoys, rows, applied, skipped,
                                     warning, nce};
      columns.insert(columns.end(), errors.begin(), errors.end());
      columns.push_back(seconds);
      return ArrowIOHelpers::tableOfColumns(columns, in_list.size());
    }

    /// The QPX sample view of -qpx_extended (quantms.sample.parquet): the sample of each run (its stem), whose
    /// organism and organism part ProSE does not know ("not available", as SDRF has it)
    static std::shared_ptr<arrow::Table> sampleTable_(const StringList& in_list)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      std::vector<Column> columns(3);
      const char* names[] = {"sample_accession", "organism", "organism_part"};
      for (Size k = 0; k < 3; ++k) { columns[k].name = names[k]; columns[k].type = Column::Type::UTF8; columns[k].nullable = false; }
      for (const std::string& in : in_list)
      {
        columns[0].utf8_values.push_back(File::stemName(in));
        columns[1].utf8_values.push_back("not available");
        columns[2].utf8_values.push_back("not available");
      }
      return ArrowIOHelpers::tableOfColumns(columns, in_list.size());
    }

    /// A step of the QPX provenance view (provenanceTable_())
    struct ProvenanceStep_
    {
      std::string category, name, tool;
      std::vector<std::pair<std::string, std::string>> parameters;
      std::string config;
      std::vector<std::string> output_views;
    };

    /**
      @brief The QPX provenance view of -qpx_extended (quantms.provenance.parquet): one row per step of the run --
      database_search (ProSE, with the complete parameters of the tool as JSON in 'config'), PeptDeep's predictions
      (category rescoring, step 'spectrum_prediction'), rescoring (Percolator), fdr_filtering, protein_inference --
      with the parameters a viewer shows (see the tool documentation)
    */
    static std::shared_ptr<arrow::Table> provenanceTable_(const std::vector<ProvenanceStep_>& steps)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type, bool nullable = true)
      {
        Column c; c.name = name; c.type = type; c.nullable = nullable; return c;
      };
      Column order = column("step_order", Column::Type::INT32, false);
      Column category = column("step_category", Column::Type::UTF8, false);
      Column name = column("step_name", Column::Type::UTF8, false);
      Column tool = column("tool_name", Column::Type::UTF8, false);
      Column version = column("tool_version", Column::Type::UTF8);
      Column uri = column("tool_uri", Column::Type::UTF8);
      Column parameters = column("parameters", Column::Type::STRUCT_LIST);
      parameters.struct_fields = {{"key", Column::Type::UTF8, false}, {"value", Column::Type::UTF8, false}};
      Column config = column("config", Column::Type::UTF8);
      Column views = column("output_views", Column::Type::UTF8_LIST);
      for (Size k = 0; k < steps.size(); ++k)
      {
        const ProvenanceStep_& s = steps[k];
        order.int64_values.push_back(static_cast<Int64>(k + 1));
        category.utf8_values.push_back(s.category);
        name.utf8_values.push_back(s.name);
        tool.utf8_values.push_back(s.tool);
        version.utf8_values.push_back(VersionInfo::getVersion());
        uri.utf8_values.push_back("");
        uri.null_rows.push_back(true);
        std::vector<Column::StructEntry> entries;
        for (const auto& [key, value] : s.parameters) { entries.push_back({key, value}); }
        parameters.struct_lists.push_back(std::move(entries));
        config.utf8_values.push_back(s.config);
        config.null_rows.push_back(s.config.empty());
        views.utf8_lists.push_back(s.output_views);
        views.null_rows.push_back(s.output_views.empty());
      }
      return ArrowIOHelpers::tableOfColumns({order, category, name, tool, version, uri, parameters, config, views}, steps.size());
    }

    /// The steps of the run for provenanceTable_()
    std::vector<ProvenanceStep_> provenanceSteps_(const ProSEAlgorithm::MultiFileSearchResult& mfres, const Param& search_params,
                                                  bool rescore, bool protein_fdr_applied, bool pooled,
                                                  const std::string& mod_delta_masses, const std::string& mod_ptms) const
    {
      const ProSEAlgorithm::SharedSearchStats& sh = mfres.shared;
      auto flag = [](bool b) -> std::string { return b ? "true" : "false"; };
      auto value = [&](const std::string& key) -> std::string
      {
        return search_params.exists(key) ? search_params.getValue(key).toString() : std::string();
      };
      std::vector<ProvenanceStep_> steps;

      ProvenanceStep_ search{"database_search", "peptide_search", "ProSE", {}, jsonOfParam_(getParam_()), {"psm"}};
      auto& p = search.parameters;
      p.emplace_back("database_file", sh.database_file);
      p.emplace_back("enzyme", sh.enzyme);
      std::string specificity = value("peptide:enzyme_specificity");
      bool clip = false;
      std::vector<std::string> features;
      bool observed_minus_theoretical = true;
      for (const auto& pf : mfres.per_file)
      {
        if (pf.protein_ids.empty()) { continue; }
        const ProteinIdentification::SearchParameters& sp = pf.protein_ids[0].getSearchParameters();
        specificity = EnzymaticDigestion::NamesOfSpecificity[sp.enzyme_term_specificity];
        clip = sp.metaValueExists("peptide:clip_nterm_methionine") && sp.getMetaValue("peptide:clip_nterm_methionine").toString() == "true";
        observed_minus_theoretical = isotopeObservedMinusTheoretical_(pf.protein_ids);
        if (sp.metaValueExists("extra_features"))
        {
          for (const std::string& feature : ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString()))
          {
            const std::string name = qpxScoreName_(feature);
            if (!name.empty() && feature != "score") { features.push_back(name); }
          }
        }
        break;
      }
      p.emplace_back("enzyme_specificity", specificity);
      p.emplace_back("missed_cleavages", std::to_string(sh.missed_cleavages));
      p.emplace_back("peptide_min_length", value("peptide:min_size"));
      p.emplace_back("peptide_max_length", value("peptide:max_size"));
      p.emplace_back("precursor_mass_tolerance_lower", number_(sh.precursor_tol_lower, "0"));
      p.emplace_back("precursor_mass_tolerance_upper", number_(sh.precursor_tol_upper, "0"));
      p.emplace_back("precursor_mass_tolerance_unit", sh.precursor_tol_unit);
      p.emplace_back("fragment_mass_tolerance", number_(sh.fragment_tol, "0"));
      p.emplace_back("fragment_mass_tolerance_unit", sh.fragment_tol_unit);
      p.emplace_back("min_charge", std::to_string(sh.min_charge));
      p.emplace_back("max_charge", std::to_string(sh.max_charge));
      p.emplace_back("fixed_modifications", ListUtils::concatenate(sh.fixed_mods, ","));
      p.emplace_back("variable_modifications", ListUtils::concatenate(sh.variable_mods, ","));
      p.emplace_back("clip_nterm_methionine", flag(clip));
      p.emplace_back("decoy_mode", sh.decoy_mode);
      p.emplace_back("decoy_string", mfres.decoy_string);
      p.emplace_back("decoy_position", mfres.decoy_is_prefix ? "prefix" : "suffix");
      p.emplace_back("open_search", flag(sh.open_search));
      p.emplace_back("calibration", flag(sh.calibration_enabled));
      p.emplace_back("chunked", flag(sh.chunked));
      p.emplace_back("isotope_error_convention", observed_minus_theoretical ? "observed_minus_theoretical" : "theoretical_minus_observed");
      p.emplace_back("extra_features", ListUtils::concatenate(features, ","));
      p.emplace_back("db_target_proteins", std::to_string(sh.db_target_proteins));
      p.emplace_back("db_decoy_proteins", std::to_string(sh.db_decoy_proteins));
      p.emplace_back("indexed_peptides", std::to_string(sh.indexed_peptides));
      p.emplace_back("indexed_fragments", std::to_string(sh.indexed_fragments));
      p.emplace_back("seconds_index_build", number_(sh.seconds_index_build, "0"));
      p.emplace_back("seconds_total", number_(sh.seconds_total, "0"));
      if (!mod_delta_masses.empty()) { p.emplace_back("mod_analysis_delta_masses", mod_delta_masses); }
      if (!mod_ptms.empty()) { p.emplace_back("mod_analysis_ptms", mod_ptms); }
      steps.push_back(std::move(search));

      // PeptDeep: its features for the rescoring ('Search:peptdeep:enable') and the predicted spectra of the psm files
      // (-qpx_extended), at the NCE of each run (run view: prediction_nce); category 'rescoring', as the closest one
      bool predicted = false;
      for (const auto& pf : mfres.per_file)
      {
        if (!pf.protein_ids.empty() && pf.protein_ids[0].metaValueExists("peptdeep_nce")) { predicted = true; }
      }
      if (predicted)
      {
        // (tool OpenMS: PeptDeep's models run in OpenMS' ONNX implementation, of OpenMS' version)
        ProvenanceStep_ prediction{"rescoring", "spectrum_prediction", "OpenMS", {}, "", {"psm"}};
        prediction.parameters = {{"method", "peptdeep_ms2"},
                                 {"ms2_model", File::basename(value("peptdeep:ms2_model"))},
                                 {"instrument", value("peptdeep:instrument")},
                                 {"rescoring_features", value("peptdeep:enable")}};
        steps.push_back(std::move(prediction));
      }
      if (rescore)
      {
        // (tool OpenMS: its own Percolator implementation, of OpenMS' version)
        steps.push_back({"rescoring", "psm_rescoring", "OpenMS", {{"method", "percolator"}}, "", {"psm"}});
      }
      steps.push_back({"fdr_filtering", "target_decoy_fdr", "ProSE",
                       {{"psm_fdr", number_(sh.psm_fdr_threshold, "0")}, {"protein_fdr", number_(sh.protein_fdr_threshold, "0")},
                        {"psm_fdr_method", "target_decoy_competition"}, {"protein_fdr_method", "picked_protein"},
                        {"protein_fdr_applied", flag(protein_fdr_applied)}},
                       "", {"psm", "pg"}});
      steps.push_back({"protein_inference", "protein_grouping", "ProSE",
                       {{"algorithm", "basic_protein_inference"}, {"pooled", flag(pooled)}}, "", {"pg"}});
      return steps;
    }

    ExitCodes main_(int, const char**) override
    {
      const StringList in_list = getStringList_("in");
      const std::string database = getStringOption_("database");
      const StringList out_idxml_list = getStringList_("out_idxml");
      const StringList out_pin_list = getStringList_("out_pin");
      const std::string out_merged = getStringOption_("out_merged");
      const std::string out_qpx_dir = getOutputDirOption("out_qpx");
      const std::string out_parquet_dir = getOutputDirOption("out_parquet");
      // getOutputDirOption auto-creates the directory if it doesn't exist (and returns "" if unset).
      const std::string out_mod_analysis_dir = getOutputDirOption("out_mod_analysis_dir");
      const bool keep_decoys = getFlag_("keep_decoys");
      const bool qpx_extended = getFlag_("qpx_extended");
      PeakFilter_ peak_filter;
      peak_filter.window = getDoubleOption_("qpx_peak_window");
      peak_filter.top = static_cast<Size>(getIntOption_("qpx_peak_top"));

      // The per-input Parquet files go to runs/ below the output directory, the files of the whole run stay on top:
      // an input may be named like them (quantms, openms, merged) without overwriting them.
      const std::string qpx_runs_dir = out_qpx_dir + "/runs";
      const std::string parquet_runs_dir = out_parquet_dir + "/runs";

      if (in_list.empty())
      {
        OPENMS_LOG_ERROR << "No input files provided (-in)." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // At least one output must be specified
      if (out_idxml_list.empty() && out_pin_list.empty() && out_qpx_dir.empty() && out_parquet_dir.empty() && out_merged.empty())
      {
        OPENMS_LOG_ERROR << "No output specified. Provide at least one of -out_idxml, -out_pin, -out_qpx, -out_parquet, or -out_merged." << endl;
        return ILLEGAL_PARAMETERS;
      }
      if (qpx_extended && out_qpx_dir.empty())
      {
        OPENMS_LOG_ERROR << "-qpx_extended extends the QPX output: it needs -out_qpx." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // -out_idxml count must match -in count
      if (!out_idxml_list.empty() && in_list.size() != out_idxml_list.size())
      {
        OPENMS_LOG_ERROR << "Number of output files (-out_idxml, " << out_idxml_list.size()
                         << ") must match number of input files (-in, " << in_list.size() << ")." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // -out_pin count must match -in count
      if (!out_pin_list.empty() && in_list.size() != out_pin_list.size())
      {
        OPENMS_LOG_ERROR << "Number of output files (-out_pin, " << out_pin_list.size()
                         << ") must match number of input files (-in, " << in_list.size() << ")." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // The stem of the files of the whole run in the -out_qpx directory (quantms.psm.parquet, ...), and in the
      // -out_parquet directory (openms.psm.parquet, ...; with -out_merged and several inputs also the pooled proteins,
      // merged.proteins.parquet and merged.pg.parquet), next to the subfolder runs/ of the per-input files
      static const std::string QPX_AGGREGATE_STEM = "quantms";
      static const std::string PARQUET_AGGREGATE_STEM = "openms";
      static const std::string PARQUET_POOLED_STEM = "merged";
      // Validate basenames are unique (for parquet directory output)
      if (!out_qpx_dir.empty() || !out_parquet_dir.empty())
      {
        std::set<std::string> basenames;
        for (const auto& in_file : in_list)
        {
          std::string bn = File::stemName(in_file);
          if (!basenames.insert(bn).second)
          {
            OPENMS_LOG_ERROR << "Duplicate input basename '" << bn
                             << "'. Parquet directory output requires unique basenames." << endl;
            return ILLEGAL_PARAMETERS;
          }
        }
      }

      // Create runs/ (the output directories themselves exist already, see getOutputDirOption).
      for (const auto& [option, dir] : {std::pair<std::string, std::string>{"-out_qpx", out_qpx_dir},
                                        std::pair<std::string, std::string>{"-out_parquet", out_parquet_dir}})
      {
        if (!dir.empty() && !File::makeDir(dir + "/runs"))
        {
          OPENMS_LOG_ERROR << "Cannot create the output directory '" << dir << "/runs' (" << option << ")." << endl;
          return CANNOT_WRITE_OUTPUT_FILE;
        }
      }

      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      Param search_params = getParam_().copy("Search:", true);
      setDeisotoping_(search_params, peak_filter);
      // Percolator rescoring (Search:rescore) runs inside the search, per input file, before its
      // PSM-level FDR; the protein FDR below then works on the rescored PSMs as well.
      const bool rescore = search_params.getValue("rescore").toString() == "true";
      const double user_protein_fdr = static_cast<double>(search_params.getValue("FDR:protein"));

      ProSEAlgorithm sse;
      sse.setParameters(search_params);
      // -out_pin is the input of a later rescoring: the PSMs before the rescoring and before any FDR filter. Both
      // Search:rescore (Percolator q-values as scores, only the best hit per spectrum) and FDR:PSM (removes PSMs,
      // q-values as scores) change the PSMs in the search, so the search also hands back each file's PSMs before them.
      const double user_psm_fdr = static_cast<double>(search_params.getValue("FDR:PSM"));
      const bool pin_from_pre_fdr = !out_pin_list.empty() && (user_psm_fdr > 0.0 || rescore);
      sse.setKeepPreFdrPsms(pin_from_pre_fdr);

      // Determine open-search mode the same way the algorithm does, so we can give
      // the user actionable feedback if -out_mod_analysis_dir is set in closed-search.
      const double precursor_tol_lower = static_cast<double>(search_params.getValue("precursor:mass_tolerance_lower"));
      const double precursor_tol_upper = static_cast<double>(search_params.getValue("precursor:mass_tolerance_upper"));
      const std::string precursor_tol_unit = search_params.getValue("precursor:mass_tolerance_unit").toString();
      const bool open_search_mode = FragmentIndex::isOpenSearchMode(precursor_tol_lower,
                                                                    precursor_tol_upper,
                                                                    precursor_tol_unit == "ppm");

      if (!out_mod_analysis_dir.empty() && !open_search_mode)
      {
        OPENMS_LOG_WARN << "-out_mod_analysis_dir was set but the search is in CLOSED-search mode "
                        << "(precursor tolerance [-" << precursor_tol_lower
                        << ", +" << precursor_tol_upper << "] " << precursor_tol_unit
                        << "). Modification-analysis tables will NOT be written." << endl;
      }

      // Build per-file modification-analysis base names if -out_mod_analysis_dir is set.
      // Each per-file base name is "<dir>/<input_basename_without_ext>" so the algorithm
      // appends suffixes like _ModificationAnalysis_DeltaMassStats.tsv.
      // Use a long, unlikely-to-collide aggregate stem so an input named "aggregate.mzML" (or "aggregate.d")
      // doesn't silently overwrite the aggregate output (or vice versa).
      static const std::string AGGREGATE_STEM = "_aggregate_across_files";
      std::vector<std::string> mod_analysis_base_names;
      std::string aggregate_base_name;
      if (!out_mod_analysis_dir.empty() && open_search_mode)
      {
        // Detect collisions: two inputs with the same File::stemName, OR an input whose
        // stem matches the aggregate stem. We refuse to run rather than silently overwrite.
        std::map<std::string, Size> stem_to_first_index;
        for (Size i = 0; i < in_list.size(); ++i)
        {
          const std::string stem = File::stemName(in_list[i]);
          if (stem == AGGREGATE_STEM)
          {
            OPENMS_LOG_ERROR << "Input file '" << in_list[i] << "' has basename '" << stem
                             << "' which would collide with the reserved aggregate output name. "
                             << "Rename the file or omit -out_mod_analysis_dir." << endl;
            return ILLEGAL_PARAMETERS;
          }
          auto [it, inserted] = stem_to_first_index.emplace(stem, i);
          if (!inserted)
          {
            OPENMS_LOG_ERROR << "Two -in files share basename '" << stem
                             << "': '" << in_list[it->second] << "' and '" << in_list[i]
                             << "'. Per-file modification-analysis TSVs would overwrite each other. "
                             << "Rename one of them or omit -out_mod_analysis_dir." << endl;
            return ILLEGAL_PARAMETERS;
          }
        }

        mod_analysis_base_names.reserve(in_list.size());
        for (const std::string& in_file : in_list)
        {
          mod_analysis_base_names.push_back(out_mod_analysis_dir + "/" + File::stemName(in_file));
        }
        aggregate_base_name = out_mod_analysis_dir + "/" + AGGREGATE_STEM;
      }

      // Single-shot multi-file search: builds the fragment index once and iterates over -in.
      // Timer spans search (incl. Percolator and PSM FDR) + protein FDR + output writing (stopped at report time).
      StopWatch sw_total; sw_total.start();
      // Keep the algorithm's pooled aggregate only where the -out_merged block below consumes it.
      const bool build_pooled_aggregate = !out_merged.empty() && in_list.size() > 1;
      ProSEAlgorithm::MultiFileSearchResult mfres =
        sse.searchWithModificationAnalysis(in_list, database, mod_analysis_base_names, aggregate_base_name,
                                           build_pooled_aggregate);

      if (mfres.per_file.size() != in_list.size())
      {
        OPENMS_LOG_ERROR << "Internal error: per-file result count (" << mfres.per_file.size()
                         << ") does not match input file count (" << in_list.size() << ")." << endl;
        return INTERNAL_ERROR;
      }

      // .pin output feeds an EXTERNAL Percolator, which also needs decoys for target/decoy
      // competition. Warn once if the run is target-only (decoys=ignore, or no decoys in the
      // database) so the user knows the .pin will not be usable for FDR. (decoys=auto/generate
      // ensure decoys are present, so this only fires for an explicit target-only search.)
      if (!out_pin_list.empty() && !mfres.have_decoys)
      {
        OPENMS_LOG_WARN << "-out_pin was requested but the search ran target-only (decoys=ignore, "
                        << "or no decoys in the database); the .pin file will contain no decoys and "
                        << "external Percolator cannot estimate FDR from it. Use '-Search:decoys auto' "
                        << "/ 'generate' or supply a decoy FASTA." << endl;
      }

      // Percolator rescoring outcome per file (Search:rescore): the search rescored the PSMs of each
      // file before its PSM-level FDR; a file it could not rescore keeps its search engine scores
      // (and if only some were rescored, -out_merged cannot pool their runs, see below).
      if (rescore)
      {
        Size n_searched = 0, n_rescored = 0;
        for (Size i = 0; i < in_list.size(); ++i)
        {
          const auto& result = mfres.per_file[i];
          if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
          ++n_searched;
          if (result.percolator.applied)
          {
            ++n_rescored;
            if (!result.percolator.warning.empty())
            {
              OPENMS_LOG_WARN << "Percolator rescoring of " << in_list[i] << ": " << result.percolator.warning << endl;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "Skipping Percolator rescoring for " << in_list[i] << ": " << result.percolator.skipped_reason
                            << " Keeping the search engine scores." << endl;
          }
        }
        if (n_searched > 0 && n_rescored == 0)
        {
          OPENMS_LOG_WARN << "Percolator rescoring is on (Search:rescore) but no input file could be rescored; "
                          << "all outputs carry the search engine scores." << endl;
        }
      }

      // -out_pin, before the protein FDR below, which removes the decoys of a single input: the PSMs before the
      // rescoring and FDR:PSM, from the copy the search keeps with setKeepPreFdrPsms().
      std::vector<bool> pin_failed(in_list.size(), false);
      std::vector<std::string> written_pin; // paths actually written, for the output manifest in the report
      for (Size i = 0; i < out_pin_list.size(); ++i)
      {
        const std::string& in_file = in_list[i];
        ProSEAlgorithm::SearchResult& result = mfres.per_file[i];
        if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue; // reported with the other outputs
        try
        {
          if (result.protein_ids.empty())
          {
            OPENMS_LOG_ERROR << "Cannot write .pin for " << in_file
                             << ": no ProteinIdentification/search parameters available." << endl;
            pin_failed[i] = true;
          }
          else
          {
            const auto& sp = result.protein_ids.front().getSearchParameters();
            const std::string enz_str = sp.digestion_enzyme.getName();
            const auto colon = sp.charges.find(':');
            const int min_charge = (colon != std::string::npos)
                                     ? StringUtils::toInt32(sp.charges.substr(0, colon))
                                     : StringUtils::toInt32(sp.charges);
            const int max_charge = (colon != std::string::npos)
                                     ? StringUtils::toInt32(sp.charges.substr(colon + 1))
                                     : min_charge;

            // Standard columns (SpecId/Label/ScanNr + mass/charge/enzyme features)
            // come from the centralized helper so .pin output matches PercolatorAdapter
            // and is directly consumable by the percolator CLI.
            StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
            if (sp.metaValueExists("extra_features"))
            {
              StringList extra = ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString());
              feature_set.insert(feature_set.end(), extra.begin(), extra.end());
            }
            if (std::find(feature_set.begin(), feature_set.end(), "score") == feature_set.end())
            {
              feature_set.push_back("score");
            }
            feature_set.push_back("Peptide");
            feature_set.push_back("Proteins");
            PercolatorInfile::store(out_pin_list[i], pin_from_pre_fdr ? result.pre_fdr_peptide_ids : result.peptide_ids,
                                   feature_set, enz_str, min_charge, max_charge);
            written_pin.push_back(out_pin_list[i]);
          }
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "Failed to write .pin output for " << in_file
                           << " -> " << out_pin_list[i] << ": " << e.what() << endl;
          pin_failed[i] = true;
        }
        result.pre_fdr_peptide_ids = PeptideIdentificationList(); // release the copy
      }

      // Protein-level FDR. FDR does not compose across runs, so picked-protein FDR is
      // statistically valid only on a COMPLETE protein set: a single input file (handled
      // here) or the pooled -out_merged set (handled below). Filtering each file of a
      // multi-file run separately would not control the FDR of the combined protein list
      // (false positives accumulate across files), so that is deliberately never done.
      // Whether the protein FDR was applied to a complete protein set (one input here, the pooled set below).
      bool protein_fdr_applied = false;
      if (user_protein_fdr > 0.0)
      {
        if (in_list.size() == 1)
        {
          ProSEAlgorithm::SearchResult& result = mfres.per_file[0];
          if (result.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK
              && !result.protein_ids.empty() && !result.peptide_ids.empty())
          {
            if (mfres.have_decoys)
            {
              // The single file is the complete experiment, so picked-protein FDR is valid. Shared
              // finalization (inference + picked FDR + decoy removal + ref cleanup) lives in
              // ProSEAlgorithm so this and the library search() path can't drift; it uses the
              // resolved decoy marker (prefix or suffix) and skips internally if no decoy proteins
              // survived inference, which it reports (the protein score type cannot tell: the inference already
              // scores the proteins with the PSMs' q-values).
              protein_fdr_applied = ProSEAlgorithm::applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids,
                mfres.decoy_string, mfres.decoy_is_prefix, user_protein_fdr, /*remove_decoys=*/ !keep_decoys);
            }
            else
            {
              OPENMS_LOG_WARN << "FDR:protein is set but the search ran target-only (decoys=ignore, "
                              << "or no decoys in the database); skipping protein FDR. Use "
                              << "'-Search:decoys auto' / 'generate' or supply a decoy FASTA." << endl;
            }
          }
        }
        else if (out_merged.empty())
        {
          // Multiple inputs but no complete (aggregate) protein set to finalize. Per-file
          // protein FDR would not control the combined FDR, so it is NOT applied and no
          // output represents the requested protein FDR. Per-file outputs are left as
          // run-level intermediate evidence (target+decoy retained), so they can still be
          // pooled and FDR-controlled downstream. Warn loudly.
          OPENMS_LOG_WARN << "FDR:protein=" << user_protein_fdr << " was requested for "
                          << in_list.size() << " input files but -out_merged was not set. "
                          << "Protein-level FDR is NOT applied: filtering each run separately "
                          << "would not control the FDR of the pooled protein list, and no "
                          << "single output here represents a complete experiment. Per-file "
                          << "outputs (-out_idxml/-out_qpx/-out_parquet) retain target+decoy "
                          << "evidence as intermediates. Use -out_merged to obtain protein FDR "
                          << "over the aggregated set, or pool the per-file outputs and apply "
                          << "FDR downstream (e.g. IDMerger -> ProteinInference -> FalseDiscoveryRate)." << endl;
        }
      }

      // Protein-FDR finalization mutates the per-file identifications after the search produced
      // RunStatistics. Recompute the
      // result-level stats so the report reflects the final identifications the user
      // receives; the HyperScore distribution and per-phase timing captured during the
      // search are preserved by updateFinalStats().
      for (Size i = 0; i < in_list.size(); ++i)
      {
        auto& result = mfres.per_file[i];
        if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) { continue; }
        ProSEAlgorithm::updateFinalStats(result.stats, result.peptide_ids, mfres.shared.enzyme,
                                         /*fdr_applied=*/ result.stats.fdr_applied);
      }

      // Write per-file outputs and track failures.
      // Accumulate Arrow tables for merged parquet output.
      std::vector<std::shared_ptr<arrow::Table>> qpx_psm_tables, qpx_pg_tables;
      // Per-run QPX scan_format tokens, reconciled for the merged file below. A built table
      // no longer carries the native IDs, so this has to be collected as we go.
      std::vector<std::string> qpx_psm_scan_formats;
      std::vector<std::shared_ptr<arrow::Table>> oms_psm_tables, oms_prot_tables, oms_pg_tables, oms_sp_tables;

      Size failed_count = 0;
      // Paths actually written, for an accurate output manifest in the report.
      std::vector<std::string> written_idxml;
      // -qpx_extended: the rows of each input in quantms.psm.parquet (its footer says), the footer keys of the per-input
      // psm files, and how they are written: optional columns allowed, and small row groups, so that a viewer can read
      // the peaks of one PSM without decoding those of the whole file
      std::vector<Int64> qpx_psm_rows(in_list.size(), 0);
      const std::map<std::string, std::string> qpx_prediction_keys =
        qpx_extended ? predictionMetadata_(search_params) : std::map<std::string, std::string>{};
      ParquetWriteConfig qpx_extended_config;
      qpx_extended_config.qpx_extra_columns = true;
      ParquetWriteConfig qpx_spectra_config = qpx_extended_config;
      qpx_spectra_config.row_group_size = QPX_SPECTRA_ROW_GROUP;
      // m/z and intensities are mostly distinct: a dictionary per row group would make the files 40 % larger
      qpx_spectra_config.plain_columns = {"predicted_mz_array", "predicted_intensity_array"};
      ParquetWriteConfig qpx_mz_config;
      qpx_mz_config.row_group_size = QPX_SPECTRA_ROW_GROUP;
      qpx_mz_config.plain_columns = {"mz", "intensity"};
      for (Size i = 0; i < in_list.size(); ++i)
      {
        const std::string& in_file = in_list[i];
        ProSEAlgorithm::SearchResult& result = mfres.per_file[i];

        if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_ERROR << "Search failed for " << in_file
                           << " (algorithm exit code " << static_cast<int>(result.exit_code) << "). Skipping output." << endl;
          ++failed_count;
          continue;
        }

        const std::string basename = File::stemName(in_file);

        // In test mode, replace the absolute MS run path with file://<basename> for reproducible diffs.
        if (getFlag_("test") && !result.protein_ids.empty())
        {
          result.protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_file)});
        }

        // Output modes are independent: a failure in one mode must not skip the
        // sibling modes for the same input. Each mode contributes at most one
        // failure to input_failed; the input counts as failed if any mode failed.
        bool input_failed = pin_failed[i];

        // -- the peaks and predictions of -qpx_extended --
        // Collected before the other outputs: selecting the NCE stores it on the run, which they write. A failed re-read
        // or prediction leaves columns null (collectPSMSpectra_ logs it); running out of memory fails the input.
        SpectraOfFile_ spectra_of_file;
        if (qpx_extended)
        {
          try
          {
            spectra_of_file = collectPSMSpectra_(in_file, search_params, result.protein_ids, result.peptide_ids, peak_filter);
          }
          catch (const std::exception& e)
          {
            OPENMS_LOG_ERROR << "Failed to collect the spectra of the PSMs of " << in_file << ": " << e.what() << endl;
            spectra_of_file = SpectraOfFile_();
            spectra_of_file.psms.assign(result.peptide_ids.size(), PSMSpectrum_());
            input_failed = true;
          }
        }

        // -- idXML output --
        if (!out_idxml_list.empty())
        {
          try
          {
            FileHandler().storeIdentifications(out_idxml_list[i], result.protein_ids, result.peptide_ids, {FileTypes::IDXML});
            written_idxml.push_back(out_idxml_list[i]);
          }
          catch (const Exception::BaseException& e)
          {
            OPENMS_LOG_ERROR << "Failed to write idXML output for " << in_file
                             << " -> " << out_idxml_list[i] << ": " << e.what() << endl;
            input_failed = true;
          }
        }

        // -- QPX directory output --
        if (!out_qpx_dir.empty())
        {
          // Wrapped like the sibling modes above: the QPX exporters refuse input they cannot key
          // (a merged run whose PSMs lack 'id_merge_index', ambiguous feature identities), and an
          // escaping exception would abort the whole per-file loop rather than failing this mode.
          try
          {
          // PSM — build table once, write to file, accumulate same table for merge.
          const std::string qpx_psm_file = qpx_runs_dir + "/" + basename + ".psm.parquet";
          auto qpx_psm_table = QPXFile::exportPSMsToQPXArrow(result.protein_ids, result.peptide_ids, /*export_all_psms=*/false);
          std::vector<std::string> spec_refs;
          spec_refs.reserve(result.peptide_ids.size());
          for (const auto& pid : result.peptide_ids)
          {
            if (!pid.getSpectrumReference().empty()) { spec_refs.push_back(pid.getSpectrumReference()); }
          }
          const std::string psm_scan_format = ArrowIOHelpers::qpxScanFormat(spec_refs);
          if (qpx_psm_table && qpx_extended)
          {
            // The per-input file gets the prediction columns and is written now, with the mz view of its spectra (neither
            // is kept for the rest of the run); quantms.psm.parquet gets the rows without the predictions, and only of
            // the inputs whose files were written.
            ExtendedColumns_ columns = extendedColumns_(result.protein_ids, result.peptide_ids, spectra_of_file.psms);
            auto rows_table = ArrowIOHelpers::withColumns(qpx_psm_table, columns.rows);
            columns.rows.clear();
            auto file_table = rows_table ? ArrowIOHelpers::withColumns(rows_table, columns.spectra) : nullptr;
            columns.spectra.clear();
            std::map<std::string, std::string> keys = qpx_prediction_keys;
            if (!psm_scan_format.empty()) { keys["scan_format"] = psm_scan_format; }
            const auto meta = ArrowIOHelpers::qpxFileMetadata("psm_file", qpx_spectra_config, keys);
            if (!file_table)
            {
              OPENMS_LOG_ERROR << "Failed to add the -qpx_extended columns for " << in_file << endl;
              input_failed = true;
            }
            else if (!meta || !ArrowIOHelpers::concatenateAndWriteToParquet({file_table}, qpx_psm_file, qpx_spectra_config, meta))
            {
              OPENMS_LOG_ERROR << "Failed to write QPX PSM parquet for " << in_file << " -> " << qpx_psm_file << endl;
              input_failed = true;
            }
            else
            {
              qpx_psm_tables.push_back(rows_table);
              qpx_psm_scan_formats.push_back(psm_scan_format);
              qpx_psm_rows[i] = ArrowIOHelpers::tableRowCount(rows_table);
              // the mz view of the spectra of its PSMs (none if they could not be read again)
              if (spectra_of_file.mz)
              {
                const std::string qpx_mz_file = qpx_runs_dir + "/" + basename + ".mz.parquet";
                std::map<std::string, std::string> mz_keys = peakMetadata_(peak_filter);
                if (!psm_scan_format.empty()) { mz_keys["scan_format"] = psm_scan_format; }
                const auto mz_meta = ArrowIOHelpers::qpxFileMetadata("mz_file", qpx_mz_config, mz_keys);
                if (!mz_meta || !ArrowIOHelpers::concatenateAndWriteToParquet({spectra_of_file.mz}, qpx_mz_file, qpx_mz_config, mz_meta))
                {
                  OPENMS_LOG_ERROR << "Failed to write the QPX mz view for " << in_file << " -> " << qpx_mz_file << endl;
                  input_failed = true;
                }
                spectra_of_file.mz.reset();
              }
            }
          }
          else if (qpx_psm_table)
          {
            qpx_psm_tables.push_back(qpx_psm_table);
            qpx_psm_scan_formats.push_back(psm_scan_format);
            if (!QPXFile::exportToParquet(qpx_psm_table, qpx_psm_file, ParquetWriteConfig{}, psm_scan_format))
            {
              OPENMS_LOG_ERROR << "Failed to write QPX PSM parquet for " << in_file << " -> " << qpx_psm_file << endl;
              input_failed = true;
            }
          }
          else if (!input_failed)
          {
            OPENMS_LOG_WARN << "QPX PSM table build returned null for " << in_file << " — skipping " << qpx_psm_file << endl;
          }

          // Protein groups — independent of PSM result.
          const std::string qpx_pg_file = qpx_runs_dir + "/" + basename + ".pg.parquet";
          std::shared_ptr<arrow::Table> qpx_pg_table;
          if (qpx_extended && !result.protein_ids.empty() && result.protein_ids[0].getIndistinguishableProteins().empty())
          {
            // -qpx_extended: without protein inference, one group per protein (without a q-value)
            std::vector<ProteinIdentification> runs = result.protein_ids;
            for (const ProteinHit& protein : runs[0].getHits())
            {
              ProteinIdentification::ProteinGroup group;
              group.probability = std::numeric_limits<double>::quiet_NaN();
              group.accessions = {protein.getAccession()};
              runs[0].getIndistinguishableProteins().push_back(group);
            }
            qpx_pg_table = ProteinGroupArrowExport::exportToArrow(runs, result.peptide_ids);
          }
          else
          {
            qpx_pg_table = ProteinGroupArrowExport::exportToArrow(result.protein_ids, result.peptide_ids);
          }
          // Row count via the library helper, not arrow::Table::num_rows(): Arrow is confined to
          // libOpenMS' implementation, and on Windows its symbols are dllimport, so calling the
          // member here links on Linux but not on MSVC.
          if (qpx_pg_table && ArrowIOHelpers::tableRowCount(qpx_pg_table) > 0)
          {
            qpx_pg_tables.push_back(qpx_pg_table);
            if (!ProteinGroupArrowExport::exportToParquet(qpx_pg_table, qpx_pg_file))
            {
              OPENMS_LOG_ERROR << "Failed to write QPX PG parquet for " << in_file << " -> " << qpx_pg_file << endl;
              input_failed = true;
            }
          }
          else if (qpx_pg_table)
          {
            // No protein groups, so no pg file. Rows come only from getIndistinguishableProteins(),
            // and ProSE infers proteins in just two places: single-file finalization, which needs
            // 'FDR:protein' > 0 AND decoys, and the merged path, which needs '-out_merged' with more
            // than one input. Identification-only runs open neither, so an empty table is the normal
            // outcome rather than a failure. Writing it anyway produced a schema-valid file saying
            // "no protein groups", which a consumer cannot tell apart from "inference ran and found
            // none" - so the file is omitted and the reason logged instead.
            OPENMS_LOG_INFO << "No protein groups inferred for " << in_file << " — not writing "
                            << qpx_pg_file << ". Protein inference requires 'FDR:protein' > 0 with "
                            << "decoys, or '-out_merged' with several inputs." << endl;
          }
          else
          {
            OPENMS_LOG_WARN << "QPX PG table build returned null for " << in_file << " — skipping " << qpx_pg_file << endl;
          }
          }
          catch (const Exception::BaseException& e)
          {
            OPENMS_LOG_ERROR << "Failed to write QPX output for " << in_file
                             << " -> " << out_qpx_dir << ": " << e.what() << endl;
            input_failed = true;
          }
          catch (const std::exception& e)
          {
            // Unlike the sibling modes, this one calls into Arrow/Parquet, whose exceptions do
            // not derive from Exception::BaseException. The write helpers check Status rather
            // than throwing, but allocation and internal parquet failures still surface here,
            // and letting one escape would skip the remaining outputs for this input.
            OPENMS_LOG_ERROR << "Failed to write QPX output for " << in_file
                             << " -> " << out_qpx_dir << " (Arrow/Parquet): " << e.what() << endl;
            input_failed = true;
          }
        }

        // -- OpenMS internal directory output --
        if (!out_parquet_dir.empty())
        {
          // PSM (internal PSMSchema) — build once, write, accumulate for merge.
          const std::string oms_psm_file = parquet_runs_dir + "/" + basename + ".psm.parquet";
          auto oms_psm_table = QPXFile::exportToArrow(result.protein_ids, result.peptide_ids, /*export_all_psms=*/false);
          if (oms_psm_table)
          {
            oms_psm_tables.push_back(oms_psm_table);
            if (!ArrowIOHelpers::writeTableToParquet(oms_psm_table, oms_psm_file))
            {
              OPENMS_LOG_ERROR << "Failed to write internal PSM parquet for " << in_file << " -> " << oms_psm_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "Internal PSM table build returned null for " << in_file << " — skipping " << oms_psm_file << endl;
          }

          // Proteins
          const std::string oms_prot_file = parquet_runs_dir + "/" + basename + ".proteins.parquet";
          auto prot_table = ProteinIdentificationArrowIO::exportProteinsToArrow(result.protein_ids);
          if (prot_table)
          {
            oms_prot_tables.push_back(prot_table);
            if (!ProteinIdentificationArrowIO::exportProteinsToParquet(result.protein_ids, oms_prot_file))
            {
              OPENMS_LOG_ERROR << "Failed to write proteins parquet for " << in_file << " -> " << oms_prot_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "Proteins table build returned null for " << in_file << " — skipping " << oms_prot_file << endl;
          }

          // Protein groups
          const std::string oms_pg_file = parquet_runs_dir + "/" + basename + ".pg.parquet";
          auto pg_table = ProteinIdentificationArrowIO::exportProteinGroupsToArrow(result.protein_ids);
          if (pg_table)
          {
            oms_pg_tables.push_back(pg_table);
            if (!ProteinIdentificationArrowIO::exportProteinGroupsToParquet(result.protein_ids, oms_pg_file))
            {
              OPENMS_LOG_ERROR << "Failed to write protein groups parquet for " << in_file << " -> " << oms_pg_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "Protein groups table build returned null for " << in_file << " — skipping " << oms_pg_file << endl;
          }

          // Search params
          const std::string oms_sp_file = parquet_runs_dir + "/" + basename + ".search_params.parquet";
          const auto sp_definitions = ModificationDefinitionIO::encodeByRun(
            result.protein_ids, ModificationDefinitionIO::collect(result.protein_ids, result.peptide_ids));
          auto sp_table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(result.protein_ids, sp_definitions);
          if (sp_table)
          {
            oms_sp_tables.push_back(sp_table);
            if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(result.protein_ids, oms_sp_file, ParquetWriteConfig{}, sp_definitions))
            {
              OPENMS_LOG_ERROR << "Failed to write search params parquet for " << in_file << " -> " << oms_sp_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "Search params table build returned null for " << in_file << " — skipping " << oms_sp_file << endl;
          }
        }

        if (input_failed) { ++failed_count; }
      }

      // -- Write merged Parquet files --
      // Always write merged files (even for single input — provides a canonical name). The QPX ones are written after
      // the merged idXML below: with -qpx_extended they need its protein inference and the statistics of the whole run.
      bool merged_ok = true;
      // -qpx_extended with several inputs and -out_merged: the protein groups of the pooled inference for quantms.pg.parquet
      std::shared_ptr<arrow::Table> pooled_pg_table;
      const bool pooled_pg = qpx_extended && !out_merged.empty() && in_list.size() > 1;
      auto write_merged_qpx = [&](const std::map<std::string, std::string>& footer)
      {
        // Concatenation inherits the first input's schema metadata, i.e. that per-run file's
        // uuid/creation_date. Stamp a fresh identity so the merged file is its own QPX file.
        // scan_format only describes the merged file if every contributing run agreed on a
        // convention; otherwise it is omitted rather than taken from an arbitrary run.
        std::map<std::string, std::string> psm_extra = footer;
        std::set<std::string> distinct_formats(qpx_psm_scan_formats.begin(), qpx_psm_scan_formats.end());
        if (distinct_formats.size() == 1 && !distinct_formats.begin()->empty())
        {
          psm_extra["scan_format"] = *distinct_formats.begin();
        }
        const ParquetWriteConfig psm_config = qpx_extended ? qpx_extended_config : ParquetWriteConfig{};
        auto psm_meta = ArrowIOHelpers::qpxFileMetadata("psm_file", psm_config, psm_extra);
        auto pg_meta  = ArrowIOHelpers::qpxFileMetadata("pg_file");
        merged_ok = psm_meta && pg_meta && merged_ok;
        if (psm_meta)
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_psm_tables, out_qpx_dir + "/" + QPX_AGGREGATE_STEM + ".psm.parquet", psm_config, psm_meta) && merged_ok;
        }
        // Same rule as the per-file pg above: with no group anywhere there is nothing to merge, and
        // an empty merged table would be indistinguishable from an inference that found nothing.
        if (pooled_pg)
        {
          if (pg_meta && pooled_pg_table && ArrowIOHelpers::tableRowCount(pooled_pg_table) > 0)
          {
            merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
              {pooled_pg_table}, out_qpx_dir + "/" + QPX_AGGREGATE_STEM + ".pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
          }
        }
        else if (pg_meta && !qpx_pg_tables.empty())
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_pg_tables, out_qpx_dir + "/" + QPX_AGGREGATE_STEM + ".pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
        }
      };

      if (!out_parquet_dir.empty())
      {
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_psm_tables,  out_parquet_dir + "/" + PARQUET_AGGREGATE_STEM + ".psm.parquet") && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_prot_tables, out_parquet_dir + "/" + PARQUET_AGGREGATE_STEM + ".proteins.parquet") && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_pg_tables,   out_parquet_dir + "/" + PARQUET_AGGREGATE_STEM + ".pg.parquet") && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_sp_tables,   out_parquet_dir + "/" + PARQUET_AGGREGATE_STEM + ".search_params.parquet") && merged_ok;
      }

      // Write merged idXML: ProSEAlgorithm::mergeFileResults() pools all per-file PSMs (possibly
      // Percolator-rescored) and runs cross-file protein inference + optional protein FDR. Wrapped in
      // try/catch so any failure here does NOT discard per-file outputs that already completed above.
      // The merge moves the identifications out of mfres; past this block only the per-file stats /
      // exit codes / Percolator outcomes / modification analyses are read (the search report and
      // renderRunSummaryYaml), never the identifications. A merged idXML that cannot be written fails
      // the run after the report, like the merged Parquet files above.
      ExitCodes merged_idxml_error = EXECUTION_OK;
      if (!out_merged.empty() && in_list.size() > 1)
      {
        try
        {
          ProSEAlgorithm::MergedResult merged;
          try
          {
            merged = sse.mergeFileResults(mfres, in_list, /*remove_decoys=*/ !keep_decoys);
            protein_fdr_applied = merged.protein_fdr_applied;
          }
          catch (const Exception::InvalidValue& e)
          {
            // The rescored runs report Percolator as search engine and carry its q-values, the others
            // ProSE's scores: the algorithm did not pool them, and mergeFileResults() rejects them.
            // (Only the merge is covered here: an InvalidValue of the store below is a write failure.)
            OPENMS_LOG_ERROR << "Cannot write the merged output -> " << out_merged << ": Percolator rescored only some "
                             << "of the input files (see the warnings above), and their PSMs cannot be pooled with the "
                             << "search engine scores of the others. Per-file outputs were written. Use "
                             << "'-Search:rescore false', or leave out the files that could not be rescored. ("
                             << e.what() << ")" << endl;
            merged_idxml_error = INCOMPATIBLE_INPUT_DATA;
          }

          if (merged_idxml_error == EXECUTION_OK)
          {
            if (getFlag_("test") && !merged.protein_ids.empty())
            {
              StringList basenames;
              for (const auto& f : in_list) basenames.push_back("file://" + File::basename(f));
              merged.protein_ids[0].setPrimaryMSRunPath(basenames);
            }

            FileHandler().storeIdentifications(out_merged, merged.protein_ids, merged.peptide_ids, {FileTypes::IDXML});
            if (pooled_pg && !out_qpx_dir.empty())
            {
              pooled_pg_table = ProteinGroupArrowExport::exportToArrow(merged.protein_ids, merged.peptide_ids);
            }

            // The pooled proteins and protein groups next to the per-file ones: the per-file runs have no protein
            // inference across files, so these are the only proteins with the protein FDR of a multi-file run.
            if (!out_parquet_dir.empty())
            {
              merged_ok = ProteinIdentificationArrowIO::exportProteinsToParquet(
                            merged.protein_ids, out_parquet_dir + "/" + PARQUET_POOLED_STEM + ".proteins.parquet") && merged_ok;
              merged_ok = ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
                            merged.protein_ids, out_parquet_dir + "/" + PARQUET_POOLED_STEM + ".pg.parquet") && merged_ok;
            }
          }
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "Failed to write merged idXML output -> " << out_merged
                           << ": " << e.what()
                           << ". Per-file outputs were written (check above for any per-file errors)." << endl;
          // Do NOT propagate yet; per-file outputs are already on disk and the report follows.
          merged_idxml_error = CANNOT_WRITE_OUTPUT_FILE;
        }
      }
      else if (!out_merged.empty() && in_list.size() == 1)
      {
        OPENMS_LOG_WARN << "-out_merged is only useful with multiple -in files. Skipping merged output." << endl;
      }

      // -- QPX: the merged files; with -qpx_extended also the run, sample and provenance views --
      if (!out_qpx_dir.empty())
      {
        std::map<std::string, std::string> footer;
        if (qpx_extended)
        {
          mfres.shared.seconds_total = sw_total.getClockTime();
          std::string mod_delta_masses, mod_ptms;
          if (!out_mod_analysis_dir.empty() && open_search_mode)
          {
            const std::string stem = in_list.size() > 1 ? AGGREGATE_STEM : File::stemName(in_list[0]);
            mod_delta_masses = stem + "_ModificationAnalysis_DeltaMassStats.tsv";
            mod_ptms = stem + "_ModificationAnalysis_PTMStats.tsv";
          }
          // (the footer keys of the spectrum columns go to the per-input files, which have the columns)
          if (qpx_psm_tables.empty())
          {
            // no input has rows: quantms.psm.parquet is written anyway, without rows
            std::vector<PSMSpectrum_> none;
            auto empty = QPXFile::exportPSMsToQPXArrow({}, {}, /*export_all_psms=*/false);
            if (empty) { empty = ArrowIOHelpers::withColumns(empty, extendedColumns_({}, {}, none).rows); }
            if (empty) { qpx_psm_tables.push_back(empty); }
          }
          // the run, sample and provenance views: the inputs and the steps of the run
          const std::vector<std::pair<std::string, std::shared_ptr<arrow::Table>>> views = {
            {"run", runTable_(mfres, in_list, qpx_psm_rows, rescore)},
            {"sample", sampleTable_(in_list)},
            {"provenance", provenanceTable_(provenanceSteps_(mfres, search_params, rescore, protein_fdr_applied,
                                                             pooled_pg, mod_delta_masses, mod_ptms))}};
          for (const auto& [view, table] : views)
          {
            const std::string file = out_qpx_dir + "/" + QPX_AGGREGATE_STEM + "." + view + ".parquet";
            const auto meta = ArrowIOHelpers::qpxFileMetadata(view + "_file");
            if (!table || !meta || !ArrowIOHelpers::concatenateAndWriteToParquet({table}, file, ParquetWriteConfig{}, meta))
            {
              OPENMS_LOG_ERROR << "Failed to write the QPX " << view << " view -> " << file << endl;
              merged_ok = false;
            }
          }
        }
        write_merged_qpx(footer);
      }

      // =====================================================================
      // End-of-search report (command line + optional YAML).
      // =====================================================================
      {
        sw_total.stop();
        mfres.shared.seconds_total = sw_total.getClockTime();
        const std::string summary_out = getStringOption_("summary_out");
        const bool no_summary = getFlag_("no_summary");
        const ProSEAlgorithm::SharedSearchStats& sh = mfres.shared;

        // -- Build the output manifest: every file actually written. --
        std::vector<std::pair<std::string, std::vector<std::string>>> manifest;
        auto add_manifest = [&manifest](const std::string& label, const std::vector<std::string>& paths) {
          std::vector<std::string> nonempty;
          for (const auto& p : paths) { if (!p.empty()) { nonempty.push_back(p); } }
          if (!nonempty.empty()) { manifest.emplace_back(label, std::move(nonempty)); }
        };
        add_manifest("idXML", written_idxml);
        add_manifest("pin", written_pin);
        if (!out_merged.empty() && in_list.size() > 1 && merged_idxml_error == EXECUTION_OK) { add_manifest("merged idXML", {out_merged}); }
        if (!out_qpx_dir.empty()) { add_manifest("QPX parquet", {out_qpx_dir}); }
        if (!out_parquet_dir.empty()) { add_manifest("OpenMS parquet", {out_parquet_dir}); }
        if (!out_mod_analysis_dir.empty() && open_search_mode) { add_manifest("mod-analysis tables", {out_mod_analysis_dir}); }

        // -- Command-line report. --
        if (!no_summary)
        {
          std::ostringstream r;
          r << "\n[ProSE] ===================== Search Report =====================\n";

          if (in_list.size() == 1)
          {
            // Single file: full block.
            const auto& pf = mfres.per_file[0];
            if (pf.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK)
            {
              r << "[ProSE] " << pf.stats.input_file << ": FAILED (exit code "
                << static_cast<int>(pf.exit_code) << ")\n";
            }
            else
            {
              ProSEAlgorithm::renderRunSummary(pf.stats, sh, pf.modification_analysis, pf.is_open_search, r);
            }
          }
          else
          {
            // Multi-file: shared config/db/index once, then a per-file table + totals.
            r << "[ProSE] Database     : " << (sh.database_file.empty() ? std::string("(in-memory)") : sh.database_file)
              << "  (" << sh.db_target_proteins << " target";
            if (sh.db_decoy_proteins > 0) { r << " + " << sh.db_decoy_proteins << " decoy"; }
            r << " proteins, decoys: " << (sh.decoy_mode.empty() ? std::string("n/a") : sh.decoy_mode) << ")\n";
            r << "[ProSE] Config       : " << sh.enzyme << ", " << sh.missed_cleavages << " MC | prec [-"
              << sh.precursor_tol_lower << ", +" << sh.precursor_tol_upper << "] " << sh.precursor_tol_unit
              << " | frag " << sh.fragment_tol << " " << sh.fragment_tol_unit
              << " | z " << sh.min_charge << "-" << sh.max_charge
              << " | mode: " << (sh.open_search ? "open" : "closed")
              << (sh.chunked ? " | chunked" : "") << "\n";
            r << "[ProSE] Fragment idx : " << sh.indexed_peptides << " peptides / " << sh.indexed_fragments
              << " fragments (build " << std::fixed << std::setprecision(1) << sh.seconds_index_build << " s)\n";
            r << "[ProSE] -------------------------------------------------------------\n";
            r << "[ProSE] " << std::left << std::setw(24) << "File" << std::right
              << std::setw(8) << "MS2" << std::setw(9) << "Matched" << std::setw(8) << "ID%"
              << std::setw(9) << "Peptide" << std::setw(8) << "Prot" << std::setw(8) << "Decoy" << "\n";

            Size tot_ms2 = 0, tot_matched = 0, tot_decoy = 0;
            for (const auto& pf : mfres.per_file)
            {
              std::string name = pf.stats.input_file;
              if (name.size() > 23) { name = "..." + name.substr(name.size() - 20); }
              if (pf.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK)
              {
                r << "[ProSE] " << std::left << std::setw(24) << name << std::right
                  << std::setw(8) << "FAILED" << "\n";
                continue;
              }
              const auto& st = pf.stats;
              const double idr = st.ms2_spectra > 0 ? (100.0 * st.matched_spectra / st.ms2_spectra) : 0.0;
              r << "[ProSE] " << std::left << std::setw(24) << name << std::right
                << std::setw(8) << st.ms2_spectra << std::setw(9) << st.matched_spectra
                << std::setw(7) << std::fixed << std::setprecision(1) << idr << "%"
                << std::setw(9) << st.unique_peptides << std::setw(8) << st.unique_proteins
                << std::setw(8) << st.decoy_psms << "\n";
              tot_ms2 += st.ms2_spectra; tot_matched += st.matched_spectra;
              tot_decoy += st.decoy_psms;
            }
            const double tot_idr = tot_ms2 > 0 ? (100.0 * tot_matched / tot_ms2) : 0.0;
            r << "[ProSE] " << std::left << std::setw(24) << "TOTAL" << std::right
              << std::setw(8) << tot_ms2 << std::setw(9) << tot_matched
              << std::setw(7) << std::fixed << std::setprecision(1) << tot_idr << "%"
              << std::setw(9) << "-" << std::setw(8) << "-" << std::setw(8) << tot_decoy << "\n";
            if (sh.decoy_mode == "none (target-only)")
            {
              r << "[ProSE] FDR          : n/a (target-only database)\n";
            }
            // Aggregate modification discovery (pooled PSMs) for open search.
            if (sh.open_search)
            {
              ProSEAlgorithm::renderModificationSummary(mfres.aggregate.modification_analysis, r);
            }
            r << "[ProSE] Total time   : " << std::fixed << std::setprecision(1) << sh.seconds_total << " s\n";
          }

          // -- Percolator rescoring status. --
          if (rescore)
          {
            Size n_ok = 0;
            for (const auto& pf : mfres.per_file) { if (pf.exit_code == ProSEAlgorithm::ExitCodes::EXECUTION_OK && pf.percolator.applied) { ++n_ok; } }
            r << "[ProSE] Percolator   : rescored " << n_ok << " / " << in_list.size() << " file(s)\n";
          }

          // -- Output manifest. --
          r << "[ProSE] ------------------------- Outputs ----------------------------\n";
          if (manifest.empty())
          {
            r << "[ProSE]   (no output files written)\n";
          }
          else
          {
            for (const auto& [label, paths] : manifest)
            {
              r << "[ProSE]   " << std::left << std::setw(16) << label << std::right << "-> ";
              for (size_t i = 0; i < paths.size(); ++i) { r << (i ? ", " : "") << paths[i]; }
              r << "\n";
            }
          }
          if (failed_count > 0)
          {
            r << "[ProSE] WARNING: " << failed_count << " of " << in_list.size() << " file(s) failed.\n";
          }
          r << "[ProSE] ==============================================================\n";
          OPENMS_LOG_INFO << r.str() << endl;
        }

        // -- Optional machine-readable YAML. Serialization lives in the OpenMS library
        //    (ProSEAlgorithm::renderRunSummaryYaml, hand-rolled, no YAML dependency); the
        //    tool just writes the returned string. --
        if (!summary_out.empty())
        {
          try
          {
            std::ofstream ofs(summary_out);
            ofs << ProSEAlgorithm::renderRunSummaryYaml(mfres, manifest, failed_count, in_list.size()) << std::endl;
            if (!ofs.good())
            {
              OPENMS_LOG_ERROR << "[ProSE] Failed to write summary YAML -> " << summary_out
                               << " (stream error)." << endl;
              return CANNOT_WRITE_OUTPUT_FILE;
            }
            OPENMS_LOG_INFO << "[ProSE] Wrote summary YAML -> " << summary_out << endl;
          }
          catch (const std::exception& e)
          {
            OPENMS_LOG_ERROR << "[ProSE] Failed to write summary YAML -> " << summary_out << ": " << e.what() << endl;
            return CANNOT_WRITE_OUTPUT_FILE;
          }
        }
      }

      if (failed_count > 0)
      {
        OPENMS_LOG_ERROR << "ProSE finished with " << failed_count
                         << " file(s) failing out of " << in_list.size() << "." << endl;
        return INTERNAL_ERROR;
      }
      if (!merged_ok)
      {
        OPENMS_LOG_ERROR << "ProSE failed to write one or more merged Parquet files." << endl;
        return INTERNAL_ERROR;
      }
      if (merged_idxml_error != EXECUTION_OK)
      {
        OPENMS_LOG_ERROR << "ProSE failed to write the merged idXML output -> " << out_merged << "." << endl;
        return merged_idxml_error;
      }
      return EXECUTION_OK;
    }
};

int main(int argc, const char** argv)
{
  ProSE tool;
  return tool.main(argc, argv);
}

///@endcond
