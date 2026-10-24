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
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
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
#include <OpenMS/FORMAT/QPXDatasetExport.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <algorithm>
#include <fstream>
#include <future>
#include <cmath>
#include <exception>
#include <iomanip>
#include <limits>
#include <new>
#include <sstream>

#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <unordered_set>
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
@note Protein inference of the protein FDR: as in ProteomicsLFQ, the inference and the picked-protein FDR use the PSMs up to a looser q-value than '-Search:FDR:PSM', '-Search:FDR:protein_inference_PSM' (default 0.1), which leaves more decoy proteins for the picked-protein FDR, and a protein scores with the best posterior error probability of its PSMs (Percolator's PEP; without rescoring the best HyperScore). The reported PSMs are still those that pass '-Search:FDR:PSM'; proteins and protein groups without one of them are removed after the protein FDR.
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
predicted again at that NCE), else 'Search:peptdeep:nce' if set, else the one PeptDeep rescoring would select for the
written PSMs (its selection alone, without its features); PSMs of precursor charge 1 to 6 are predicted, others get
null. Each input is read once more after its search for the spectra of its PSMs (as the search reads it: in one pass,
parsed in parallel); a failed re-read is a warning and leaves out its mz view (and the 'spectrum_reference' of its PSMs
null), a failed prediction a warning and null prediction columns; without ONNX Runtime the prediction columns are
null. The work of each input is logged in phase lines, each written at once when its phase starts, for a progress
display: "-qpx_extended: INPUT: reading spectra", "-qpx_extended: INPUT: selecting NCE" (only if the NCE is selected),
"-qpx_extended: INPUT: predicting N peptides at NCE X (SOURCE)" (N the distinct peptide and charge pairs, SOURCE
'from PeptDeep rescoring', 'set' or 'selected'; not without ONNX Runtime or without PSMs) and
"-qpx_extended: INPUT: writing" (the outputs of the input), INPUT the path of the input as given; other lines of that
prefix are summaries or warnings. The NCE selection logs no "[PeptDeepRescoring]" lines of its own (except warnings
about its input), so they stand for PeptDeep rescoring only. Each input gets the QPX mz view (runs/BASENAME.mz.parquet, BASENAME the input's basename) of the MS2 spectra of
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
another spectrum of its native ID has the row). A per-input mz view is written while the PSMs of the input are
predicted, its psm and pg files right after the input's other outputs, each file on a thread of its own while the run
goes on (an mz view whose psm file cannot be written is removed again); psm and mz files are in row groups of 256 rows,
so that a viewer can read the spectrum and the prediction of one PSM alone; the psm file's footer has the keys of its prediction columns:
'prose_prediction_ms2_model', 'prose_prediction_instrument'. quantms.psm.parquet has the rows of the inputs whose
psm file was written, in input order, without the prediction columns (and no rows if there are none). The output
directory is one QPX dataset (prefix 'quantms'): the quantms.* files of the whole run on top, the per-input files
(BASENAME.psm/mz/pg.parquet) in its subfolder runs/, so an input may be named 'quantms' too; two inputs of one
basename are still refused, as they would share runs/BASENAME.*.
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
estimate, fewer than 10 PSMs): the mass errors of the best hits of the run's target PSMs after the FDR filters (decoys
left out), signed, observed minus theoretical m/z in ppm: of the precursor, its m/z corrected by the hit's isotope error
(the 13C peak the search matched), in an open search of the PSMs without a delta mass (within 0.05 Da of zero) only;
of the fragments, every annotated b and y peak (without neutral losses) of these PSMs. The median is that of the signed errors (the bias of the run), the MAD their median absolute deviation from the
median (unscaled), and the recommended tolerance the symmetric one that covers the bias and three standard deviations
of a normal distribution estimated robustly: ceil(|median| + 3 * 1.4826 * MAD). The fragment errors are those of the
matched peaks, so within the searched fragment tolerance: a recommended fragment tolerance near it says that the
searched one may be too narrow rather than right (the report's "Tolerance estimate" is the same); 'search_seconds'
(double).
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

    /// The PeptDeep prediction of one PSM row of -qpx_extended: its 'predicted_*' columns (null when @p has_prediction
    /// is false)
    struct PSMPrediction_
    {
      bool has_prediction = false;
      std::vector<std::string> ion_type;
      std::vector<float> mz, intensity;
    };

    /// The rows per row group of the per-input psm and mz files of -qpx_extended: a viewer reads the peaks and the
    /// prediction of one PSM by decoding its row group (with ~600 peaks per spectrum about 1.2 MB of peaks)
    static constexpr Int64 QPX_SPECTRA_ROW_GROUP = 256;

    /// The deisotoping of ProSE's search for the peak filter of the mz view (as ProSEAlgorithm decides it, from its
    /// parameters @p search_params)
    static void setDeisotoping_(const Param& search_params, QPXDatasetExport::PeakFilter& filter)
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

    /// A phase line of the -qpx_extended work on the input @p in_file, "-qpx_extended: <in_file>: <what>", written at once
    /// (std::endl flushes the log): a progress parser follows the work by these lines (see the tool documentation)
    static void logPhase_(const std::string& in_file, const std::string& what)
    {
      OPENMS_LOG_INFO << "-qpx_extended: " << in_file << ": " << what << endl;
    }

    /// Where collectPSMSpectra_() writes the mz view of an input: the folder and prefix of its QPX dataset
    /// (runs/<basename>.mz.parquet), the writer configuration and the footer keys
    struct MzTarget_
    {
      std::string directory, prefix;
      ParquetWriteConfig config;
      std::map<std::string, std::string> keys;
      std::string file() const { return directory + "/" + prefix + ".mz.parquet"; }
    };

    /// The native IDs of @p peptide_ids that have one (for ArrowIOHelpers::qpxScanFormat())
    static std::vector<std::string> spectrumReferences_(const PeptideIdentificationList& peptide_ids)
    {
      std::vector<std::string> references;
      references.reserve(peptide_ids.size());
      for (const PeptideIdentification& id : peptide_ids)
      {
        std::string reference = id.getSpectrumReference();
        if (!reference.empty()) { references.push_back(std::move(reference)); }
      }
      return references;
    }

    /// The spectra and predictions of the PSMs of one searched file (collectPSMSpectra_())
    struct SpectraOfFile_
    {
      /// the QPX mz view of the file's MS2 spectra of its PSMs, without its table (written already; see @c mz_rows), and
      /// per element of the file's peptide identifications whether it has the row of its spectrum (for its
      /// 'spectrum_reference')
      QPXDatasetExport::MzView mz;
      /// the rows of the mz view; 0 if it could not be built (or the spectra could not be read)
      Size mz_rows = 0;
      /// the mz view was built, but could not be written (logged)
      bool mz_write_failed = false;
      /// the predictions: one entry per row of the file's psm view (per element of its peptide identifications with hits)
      std::vector<PSMPrediction_> predictions;
    };

    /**
      @brief The QPX mz view of the MS2 spectra of the PSMs of one searched file, with filtered peaks, and the PeptDeep
      predictions of its PSMs (-qpx_extended)

      Never throws because of the spectra or the predictions: a failed re-read of @p in_file or a failed prediction is
      logged as a warning and leaves no mz view or null prediction columns; the search result of the file is not
      affected. Only std::bad_alloc propagates: out of memory is not taken for an unreadable file.

      The spectra are read again from @p in_file, as the search reads them (ProSEAlgorithm::loadSpectra(): one pass,
      parsed in parallel), keeping only the MS2 spectra of the PSMs' native IDs; the mz view is QPXDatasetExport::mzView()
      of them, with @p filter, built and written to @p mz_target while the predictions are made (on a thread of its
      own): the caller removes its file if the psm file of the input cannot be written. After a failed re-read no PSM
      has its spectrum (a null spectrum_reference). The predictions
      are PeptDeep's at the NCE of the file: the one PeptDeep rescoring stored on the run ('peptdeep_nce'), else
      'peptdeep:nce' if set, else the one PeptDeepRescoring::selectNCE() selects for @p peptide_ids (as
      PeptDeepRescoring::annotate() would); the NCE is then stored on the run.
    */
    SpectraOfFile_ collectPSMSpectra_(const std::string& in_file, const Param& search_params,
                                      std::vector<ProteinIdentification>& protein_ids,
                                      const PeptideIdentificationList& peptide_ids, const QPXDatasetExport::PeakFilter& filter,
                                      const MzTarget_& mz_target)
    {
      SpectraOfFile_ result;
      std::vector<PSMPrediction_> out(peptide_ids.size()); // per element of peptide_ids

      // -- peaks: the MS2 spectra of the PSMs, as the file has them --
      std::unordered_set<std::string> native_ids; // of the PSMs
      for (const PeptideIdentification& id : peptide_ids)
      {
        if (!id.getHits().empty()) { native_ids.insert(id.getSpectrumReference()); }
      }
      PeakMap spectra; // the spectra of the PSMs (several, if a native ID is not unique in the file)
      bool read = true;
      logPhase_(in_file, "reading spectra");
      try
      {
        spectra = ProSEAlgorithm::loadSpectra(in_file, [&native_ids](const MSSpectrum& spectrum)
        {
          return spectrum.getMSLevel() == 2 && native_ids.count(spectrum.getNativeID()) > 0;
        });
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

      // -- the mz view: a row per spectrum, in the order of its first PSM; built and written on a thread of its own,
      //    while the predictions are made below --
      struct MzOfFile
      {
        QPXDatasetExport::MzView mz;
        Size rows = 0;
        bool table = false; // (built)
        bool write_failed = false;
      };
      std::future<MzOfFile> mz_of_file;
      if (read)
      {
#ifdef _OPENMP
        const int threads = omp_get_max_threads(); // (a new thread starts with the default of the process)
#endif
        mz_of_file = std::async(std::launch::async, [&]()
        {
#ifdef _OPENMP
          omp_set_num_threads(threads);
#endif
          MzOfFile built;
          built.mz = QPXDatasetExport::mzView(spectra, peptide_ids, ArrowIOHelpers::qpxRunFileName(in_file), filter);
          built.rows = ArrowIOHelpers::tableRowCount(built.mz.table);
          built.table = built.mz.table != nullptr;
          if (built.mz.table)
          {
            // (writeDataset() logs a failure)
            built.write_failed = !QPXDatasetExport::writeDataset(mz_target.directory, mz_target.prefix,
                                                                 {{"mz", {{built.mz.table}, mz_target.config, mz_target.keys}}});
            built.mz.table.reset();
          }
          return built;
        });
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
        out.assign(peptide_ids.size(), PSMPrediction_());
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
      // -- the mz view, once built and written --
      QPXDatasetExport::MzView& mz = result.mz;
      bool table = false;
      if (mz_of_file.valid())
      {
        MzOfFile built = mz_of_file.get();
        table = built.table;
        mz = std::move(built.mz);
        result.mz_rows = built.rows;
        result.mz_write_failed = built.write_failed;
      }
      else
      {
        mz.has_spectrum.assign(peptide_ids.size(), false);
      }
      if (mz.without_spectrum > 0)
      {
        OPENMS_LOG_WARN << "-qpx_extended: " << mz.without_spectrum << " PSM(s) of " << in_file
                        << " have no MS2 spectrum of their native ID (and RT) in the file; the mz view has none for them,"
                        << " and their spectrum_reference is null." << endl;
      }
      if (mz.unresolved > 0)
      {
        OPENMS_LOG_WARN << "-qpx_extended: " << mz.unresolved << " PSM(s) of " << in_file << " are of a spectrum that shares its"
                        << " native ID with another spectrum (at another RT), which has the mz row of that native ID: the mz"
                        << " view has none for them, and their spectrum_reference is null (unresolved)." << endl;
      }
      if (read && !table)
      {
        OPENMS_LOG_WARN << "-qpx_extended: the mz view of " << in_file << " could not be built (see the error above); no mz"
                        << " view is written for it, and the spectrum_reference of its PSMs is null." << endl;
      }
      OPENMS_LOG_INFO << "-qpx_extended: " << in_file << ": " << result.mz_rows
                      << " spectra of its PSMs, " << mz.peaks_kept << " of " << mz.peaks_read << " peaks kept";
      if (filter.top > 0)
      {
        OPENMS_LOG_INFO << " (the " << filter.top << " most intense per " << filter.window << " Th and the annotated ones)";
      }
      OPENMS_LOG_INFO << endl;
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (!peptide_ids[i].getHits().empty()) { result.predictions.push_back(std::move(out[i])); } // (a row)
      }
      return result;
    }
#ifdef WITH_ONNX
    /// The predicted ions of the PeptDeep prediction @p flat of @p peptide at precursor charge @p z: its b and y ions of
    /// charge 1 (and 2 for precursors of charge 2 or more), relative to the most intense of them, without those of 0.01
    /// or less, with their m/z (of @p generator), by m/z. None (no prediction) if @p flat is empty (not predicted), has
    /// another layout (@p layout false) or no intensity.
    static PSMPrediction_ predictedIons_(const TheoreticalSpectrumGenerator& generator, const std::string& peptide, int z,
                                         const std::vector<float>& flat, bool& layout)
    {
      PSMPrediction_ ions;
      if (flat.empty()) { return ions; } // not predicted (counted by the caller)
      const AASequence sequence = AASequence::fromString(peptide);
      const Size length = sequence.size();
      if (flat.size() != PeptDeepMS2Inference::ION_CHANNELS * (length - 1)) { layout = false; return ions; }
      const int max_charge = z > 1 ? 2 : 1;
      double base_peak = 0.0;
      for (Size ordinal = 1; ordinal < length; ++ordinal)
      {
        for (int charge = 1; charge <= max_charge; ++charge)
        {
          base_peak = std::max({base_peak, static_cast<double>(flat[PeptDeepMS2Inference::ionIndex(length, 'b', ordinal, charge)]),
                                static_cast<double>(flat[PeptDeepMS2Inference::ionIndex(length, 'y', ordinal, charge)])});
        }
      }
      if (!(base_peak > 0.0)) { return ions; }

      PeakSpectrum theoretical;
      generator.getSpectrum(theoretical, sequence, 1, max_charge);
      const auto& names = theoretical.getStringDataArrays().at(0);
      const auto& ion_charges = theoretical.getIntegerDataArrays().at(0);
      ions.has_prediction = true;
      for (Size t = 0; t < theoretical.size(); ++t)
      {
        const std::string& name = names[t]; // e.g. "b3+", "y12++"
        Size ordinal = 0;
        for (Size i = 1; i < name.size() && name[i] >= '0' && name[i] <= '9'; ++i) { ordinal = ordinal * 10 + static_cast<Size>(name[i] - '0'); }
        if (ion_charges[t] > max_charge) { continue; }
        const Size index = PeptDeepMS2Inference::ionIndex(length, name[0], ordinal, ion_charges[t]);
        if (index == PeptDeepMS2Inference::NO_ION) { continue; }
        const double relative = static_cast<double>(flat[index]) / base_peak;
        if (!(relative > 0.01)) { continue; }
        ions.ion_type.push_back(name);
        ions.mz.push_back(static_cast<float>(theoretical[t].getMZ()));
        ions.intensity.push_back(static_cast<float>(relative));
      }
      return ions;
    }

    /// The session of PeptDeep's MS2 model for the predictions of -qpx_extended: loaded once per run, at its first use
    /// (the model, threads and batch size are those of the search parameters, the same for every input)
    std::unique_ptr<PeptDeepMS2Inference> ms2_session_;

    PeptDeepMS2Inference& ms2Session_(const std::string& ms2_model, int threads, int batch_size)
    {
      if (!ms2_session_) { ms2_session_ = std::make_unique<PeptDeepMS2Inference>(ms2_model, threads, static_cast<size_t>(batch_size)); }
      return *ms2_session_;
    }

    /// The predictions of collectPSMSpectra_(), into @p out. Throws on failure (the caller logs it).
    void predictPSMSpectra_(const std::string& in_file, const Param& search_params,
                            std::vector<ProteinIdentification>& protein_ids,
                            const PeptideIdentificationList& peptide_ids, const PeakMap& spectra,
                            std::vector<PSMPrediction_>& out)
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

      // The NCE: the one PeptDeep rescoring used, else 'peptdeep:nce' if set, else the one PeptDeepRescoring::selectNCE()
      // selects (as annotate() would).
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
        // The NCE PeptDeep rescoring would have selected (its selection, without predicting all PSMs and computing its
        // features); it logs no "[PeptDeepRescoring]" lines of its own, so that a progress parser does not take this
        // for the rescoring
        logPhase_(in_file, "selecting NCE");
        PeptDeepRescoring rescoring;
        Param p = rescoring.getParameters();
        p.setValue("ms2_model", ms2_model);
        p.setValue("instrument", instrument);
        rescoring.setParameters(p);
        rescoring.setLogType(log_type_);
        rescoring.selectNCE(spectra, protein_ids, peptide_ids, ms2Session_(ms2_model, threads, batch_size));
        if (!run.metaValueExists(PeptDeepRescoring::NCE_META_VALUE))
        {
          OPENMS_LOG_WARN << "-qpx_extended: no NCE selected for " << in_file << "; the predicted_* columns are null." << endl;
          return;
        }
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

      std::ostringstream predicting;
      predicting << "predicting " << peptides.size() << " peptides at NCE " << nce << " (" << nce_source << ")";
      logPhase_(in_file, predicting.str());
      PeptDeepMS2Inference& ms2 = ms2Session_(ms2_model, threads, batch_size);
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

      // The b and y ions of each pair, independently of each other: in parallel, with a generator per thread
      // (exceptions must not leave the parallel region).
      Param tsg = TheoreticalSpectrumGenerator().getParameters();
      tsg.setValue("add_b_ions", "true");
      tsg.setValue("add_y_ions", "true");
      tsg.setValue("add_first_prefix_ion", "true");
      tsg.setValue("add_metainfo", "true");
      std::vector<PSMPrediction_> ions_of_key(peptides.size());
      Size other_layout = 0;
      std::exception_ptr failure;
#pragma omp parallel reduction(+ : other_layout)
      {
        TheoreticalSpectrumGenerator generator;
        generator.setParameters(tsg);
#pragma omp for schedule(dynamic, 64)
        for (SignedSize k = 0; k < static_cast<SignedSize>(peptides.size()); ++k)
        {
          try
          {
            bool layout = true;
            const Size key = static_cast<Size>(k);
            ions_of_key[key] = predictedIons_(generator, peptides[key], static_cast<int>(charges[key]), predicted[key], layout);
            if (!layout) { ++other_layout; }
          }
          catch (...)
          {
#pragma omp critical (ProSE_predictPSMSpectra)
            if (!failure) { failure = std::current_exception(); }
          }
        }
      }
      if (failure) { std::rethrow_exception(failure); }
      if (other_layout > 0)
      {
        OPENMS_LOG_WARN << "-qpx_extended: the MS2 model gave " << other_layout << " of " << peptides.size()
                        << " peptides of " << in_file << " another number of values than " << PeptDeepMS2Inference::ION_CHANNELS
                        << " per cleavage (an unexpected model layout); their predicted_* columns are null." << endl;
      }
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (row_key[i] == std::numeric_limits<Size>::max()) { continue; }
        const PSMPrediction_& ions = ions_of_key[row_key[i]];
        out[i].has_prediction = ions.has_prediction;
        out[i].ion_type = ions.ion_type;
        out[i].mz = ions.mz;
        out[i].intensity = ions.intensity;
      }
    }
#endif

    /// Whether ProSE's 'isotope_error' meta values of the run @p protein_ids are observed minus theoretical (the record
    /// 'isotope_error_convention' of its search parameters; without it they are theoretical minus observed)
    static bool isotopeObservedMinusTheoretical_(const std::vector<ProteinIdentification>& protein_ids)
    {
      if (protein_ids.empty()) { return true; } // no PSMs: no isotope errors
      const ProteinIdentification::SearchParameters& sp = protein_ids[0].getSearchParameters();
      return sp.metaValueExists("isotope_error_convention")
             && sp.getMetaValue("isotope_error_convention").toString() == "observed_minus_theoretical";
    }

    /// The prediction columns of the per-input psm files of -qpx_extended ('predicted_ion_type_array',
    /// 'predicted_mz_array', 'predicted_intensity_array'), one entry per row of @p predictions, which are moved out
    static std::vector<ArrowIOHelpers::ParquetColumn> predictionColumns_(std::vector<PSMPrediction_>& predictions)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
      Column ion_type = column("predicted_ion_type_array", Column::Type::UTF8_LIST);
      Column mz = column("predicted_mz_array", Column::Type::FLOAT_LIST);
      Column intensity = column("predicted_intensity_array", Column::Type::FLOAT_LIST);
      for (PSMPrediction_& p : predictions)
      {
        ion_type.utf8_lists.push_back(std::move(p.ion_type));
        ion_type.null_rows.push_back(!p.has_prediction);
        mz.float_lists.push_back(std::move(p.mz));
        mz.null_rows.push_back(!p.has_prediction);
        intensity.float_lists.push_back(std::move(p.intensity));
        intensity.null_rows.push_back(!p.has_prediction);
      }
      predictions.clear();
      predictions.shrink_to_fit();
      // (moved into the result: returning an initializer list would copy every column)
      std::vector<Column> columns;
      for (Column* c : {&ion_type, &mz, &intensity}) { columns.push_back(std::move(*c)); }
      return columns;
    }

    /// The footer keys of the mz.parquet files of -qpx_extended: the peak filter of their peaks
    static std::map<std::string, std::string> peakMetadata_(const QPXDatasetExport::PeakFilter& filter)
    {
      return {{"prose_peak_window", QPXDatasetExport::formatNumber(filter.window, "0")}, {"prose_peak_top", std::to_string(filter.top)},
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

    /// The runs of the QPX run view of -qpx_extended (quantms.run.parquet): one per input, in input order
    /// (run_accession = run_file_name: the input's stem; one sample of the stem with label LFQ; the search's enzyme and
    /// modifications)
    static std::vector<QPXDatasetExport::Run> runs_(const ProSEAlgorithm::SharedSearchStats& sh, const StringList& in_list)
    {
      std::vector<QPXDatasetExport::Run> runs;
      for (const std::string& in : in_list)
      {
        QPXDatasetExport::Run run;
        run.accession = File::stemName(in);
        run.run_file_name = run.accession;
        run.file_name = File::basename(in);
        QPXDatasetExport::RunSample sample;
        sample.accession = run.accession;
        run.samples.push_back(sample);
        if (!sh.enzyme.empty()) { run.enzymes.push_back(sh.enzyme); }
        run.fixed_modifications = sh.fixed_mods;
        run.variable_modifications = sh.variable_mods;
        runs.push_back(std::move(run));
      }
      return runs;
    }

    /// ProSE's optional columns of the QPX run view (runs_(); see the tool documentation): the search status, counts,
    /// rescoring, prediction NCE, error estimates and seconds of each input (with @p test_mode 0 seconds, so the view
    /// can be compared with a reference, as TOPPBase fixes the time of the data processing in test mode)
    static std::vector<ArrowIOHelpers::ParquetColumn> runExtraColumns_(const ProSEAlgorithm::MultiFileSearchResult& mfres,
                                                                       const StringList& in_list, const std::vector<Int64>& psm_rows,
                                                                       bool rescore, bool test_mode)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
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
      auto optional_text = [](Column& c, const std::string& s) { c.utf8_values.push_back(s); c.null_rows.push_back(s.empty()); };
      auto optional_number = [](Column& c, bool valid, double v)
      {
        c.double_values.push_back(valid && std::isfinite(v) ? v : 0.0);
        c.null_rows.push_back(!(valid && std::isfinite(v)));
      };
      for (Size i = 0; i < in_list.size(); ++i)
      {
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
        seconds.double_values.push_back(pf && !test_mode ? pf->stats.seconds_search : 0.0);
        const ProSEAlgorithm::RunStatistics empty_stats;
        const ProSEAlgorithm::RunStatistics& st = pf ? pf->stats : empty_stats;
        optional_number(errors[0], st.prec_tol_valid, st.prec_err_median);
        optional_number(errors[1], st.prec_tol_valid, st.prec_err_mad);
        optional_number(errors[2], st.prec_tol_valid, st.prec_err_recommended);
        optional_number(errors[3], st.frag_tol_valid, st.frag_err_median);
        optional_number(errors[4], st.frag_tol_valid, st.frag_err_mad);
        optional_number(errors[5], st.frag_tol_valid, st.frag_err_recommended);
      }
      std::vector<Column> columns = {status, ms2, identified, decoys, rows, applied, skipped, warning, nce};
      columns.insert(columns.end(), errors.begin(), errors.end());
      columns.push_back(seconds);
      return columns;
    }

    /// The steps of the run for the QPX provenance view of -qpx_extended (quantms.provenance.parquet), with the
    /// parameters a viewer shows (see the tool documentation)
    std::vector<QPXDatasetExport::ProvenanceStep> provenanceSteps_(const ProSEAlgorithm::MultiFileSearchResult& mfres,
                                                                   const Param& search_params, bool rescore,
                                                                   bool protein_fdr_applied, bool pooled,
                                                                   const std::string& mod_delta_masses,
                                                                   const std::string& mod_ptms) const
    {
      const ProSEAlgorithm::SharedSearchStats& sh = mfres.shared;
      auto flag = [](bool b) -> std::string { return b ? "true" : "false"; };
      auto value = [&](const std::string& key) -> std::string
      {
        return search_params.exists(key) ? search_params.getValue(key).toString() : std::string();
      };
      std::vector<QPXDatasetExport::ProvenanceStep> steps;

      QPXDatasetExport::ProvenanceStep search{"database_search", "peptide_search", "ProSE", {},
                                              QPXDatasetExport::paramToJSON(getParam_()), {"psm"}};
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
            const std::string name = QPXDatasetExport::scoreName(feature);
            if (!name.empty() && feature != "score") { features.push_back(name); }
          }
        }
        break;
      }
      p.emplace_back("enzyme_specificity", specificity);
      p.emplace_back("missed_cleavages", std::to_string(sh.missed_cleavages));
      p.emplace_back("peptide_min_length", value("peptide:min_size"));
      p.emplace_back("peptide_max_length", value("peptide:max_size"));
      p.emplace_back("precursor_mass_tolerance_lower", QPXDatasetExport::formatNumber(sh.precursor_tol_lower, "0"));
      p.emplace_back("precursor_mass_tolerance_upper", QPXDatasetExport::formatNumber(sh.precursor_tol_upper, "0"));
      p.emplace_back("precursor_mass_tolerance_unit", sh.precursor_tol_unit);
      p.emplace_back("fragment_mass_tolerance", QPXDatasetExport::formatNumber(sh.fragment_tol, "0"));
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
      p.emplace_back("seconds_index_build", QPXDatasetExport::formatNumber(sh.seconds_index_build, "0"));
      p.emplace_back("seconds_total", QPXDatasetExport::formatNumber(sh.seconds_total, "0"));
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
        QPXDatasetExport::ProvenanceStep prediction{"rescoring", "spectrum_prediction", "OpenMS", {}, "", {"psm"}};
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
                       {{"psm_fdr", QPXDatasetExport::formatNumber(sh.psm_fdr_threshold, "0")}, {"protein_fdr", QPXDatasetExport::formatNumber(sh.protein_fdr_threshold, "0")},
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
      QPXDatasetExport::PeakFilter peak_filter;
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
      // the outputs from here on, with the -qpx_extended work (the report's timing)
      const double seconds_searching = sw_total.getClockTime();
      StopWatch sw_qpx_extended;

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
              // survived inference, which it reports.
              protein_fdr_applied = ProSEAlgorithm::applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids,
                mfres.decoy_string, mfres.decoy_is_prefix, user_protein_fdr, /*remove_decoys=*/ !keep_decoys,
                result.inference_only_peptide_ids);
              result.inference_only_peptide_ids = PeptideIdentificationList(); // release the copy
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
          // the PSMs that only the protein inference of the protein FDR would use
          for (ProSEAlgorithm::SearchResult& result : mfres.per_file) { result.inference_only_peptide_ids = PeptideIdentificationList(); }
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
      std::vector<bool> input_failed_of(in_list.size(), false); // (counted in failed_count)
      // -qpx_extended: the writes of the inputs' psm and pg files, each on a thread of its own while the run goes on
      // (finished below, after quantms.*: the rows of an input are in quantms.psm.parquet unless its psm file fails)
      struct InputWrite
      {
        Size input = 0;
        std::string in_file, mz_file;          // (mz_file: its mz view, written already)
        std::shared_ptr<arrow::Table> rows;    // its rows of quantms.psm.parquet
        std::future<std::pair<bool, std::set<std::string>>> result; // written, and the views that failed
      };
      std::vector<InputWrite> input_writes;
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
        // The protein groups of the QPX output, built on a thread of their own while the input's other outputs are (from
        // a copy of its runs, on which -qpx_extended stores the NCE of its predictions; its PSMs are only read)
        std::future<std::shared_ptr<arrow::Table>> pg_build;
        if (!out_qpx_dir.empty())
        {
          pg_build = std::async(std::launch::async, [&result, qpx_extended, runs = result.protein_ids]() mutable
          {
            if (qpx_extended && !runs.empty() && runs[0].getIndistinguishableProteins().empty())
            {
              // -qpx_extended: without protein inference, one group per protein (without a q-value)
              for (const ProteinHit& protein : runs[0].getHits())
              {
                ProteinIdentification::ProteinGroup group;
                group.probability = std::numeric_limits<double>::quiet_NaN();
                group.accessions = {protein.getAccession()};
                runs[0].getIndistinguishableProteins().push_back(group);
              }
            }
            return ProteinGroupArrowExport::exportToArrow(runs, result.peptide_ids);
          });
        }

        SpectraOfFile_ spectra_of_file;
        std::string mz_file;     // -qpx_extended: the mz view of the input,
        bool mz_written = false; // written already (removed below if its psm file cannot be written)
        if (qpx_extended)
        {
          sw_qpx_extended.resume();
          try
          {
            MzTarget_ mz_target{qpx_runs_dir, basename, qpx_mz_config, peakMetadata_(peak_filter)};
            const std::string scan_format = ArrowIOHelpers::qpxScanFormat(spectrumReferences_(result.peptide_ids));
            if (!scan_format.empty()) { mz_target.keys["scan_format"] = scan_format; }
            mz_file = mz_target.file();
            spectra_of_file = collectPSMSpectra_(in_file, search_params, result.protein_ids, result.peptide_ids, peak_filter,
                                                 mz_target);
            mz_written = spectra_of_file.mz_rows > 0 && !spectra_of_file.mz_write_failed;
            if (spectra_of_file.mz_write_failed) { input_failed = true; } // (logged)
          }
          catch (const std::exception& e)
          {
            OPENMS_LOG_ERROR << "Failed to collect the spectra of the PSMs of " << in_file << ": " << e.what() << endl;
            spectra_of_file = SpectraOfFile_();
            spectra_of_file.mz.has_spectrum.assign(result.peptide_ids.size(), false);
            spectra_of_file.predictions.assign(std::count_if(result.peptide_ids.begin(), result.peptide_ids.end(),
                                                             [](const PeptideIdentification& id) { return !id.getHits().empty(); }),
                                               PSMPrediction_());
            input_failed = true;
          }
          sw_qpx_extended.stop();
          logPhase_(in_file, "writing");
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
        bool psm_file_written = false; // -qpx_extended: the input's psm file
        if (!out_qpx_dir.empty())
        {
          // Wrapped like the sibling modes above: the QPX exporters refuse input they cannot key
          // (a merged run whose PSMs lack 'id_merge_index', ambiguous feature identities), and an
          // escaping exception would abort the whole per-file loop rather than failing this mode.
          try
          {
          // PSM — build table once, write to file, accumulate same table for merge.
          const std::string qpx_psm_file = qpx_runs_dir + "/" + basename + ".psm.parquet";
          // -qpx_extended: the optional columns of the psm view, built on a thread of their own while the view is
          std::future<std::vector<ArrowIOHelpers::ParquetColumn>> extra_columns;
          if (qpx_extended)
          {
            extra_columns = std::async(std::launch::async, [&]()
            {
              return QPXDatasetExport::psmColumns(result.peptide_ids, spectra_of_file.mz.has_spectrum,
                                                  isotopeObservedMinusTheoretical_(result.protein_ids));
            });
          }
          auto qpx_psm_table = QPXFile::exportPSMsToQPXArrow(result.protein_ids, result.peptide_ids, /*export_all_psms=*/false);
          const std::string psm_scan_format = ArrowIOHelpers::qpxScanFormat(spectrumReferences_(result.peptide_ids));
          // -qpx_extended: the views of the input's own QPX dataset (runs/<basename>.psm/mz/pg.parquet), written below
          std::map<std::string, QPXDatasetExport::View> run_views;
          std::shared_ptr<arrow::Table> psm_rows_table; // (its rows in quantms.psm.parquet, once its psm file is written)
          if (qpx_psm_table && qpx_extended)
          {
            // The per-input file gets the prediction columns and is written below with its mz and pg views (it is not
            // kept for the rest of the run); quantms.psm.parquet gets the rows without the predictions, and only of the
            // inputs whose psm file was written.
            auto rows_table = ArrowIOHelpers::withColumns(qpx_psm_table, extra_columns.get());
            auto file_table = rows_table ? ArrowIOHelpers::withColumns(rows_table, predictionColumns_(spectra_of_file.predictions)) : nullptr;
            std::map<std::string, std::string> keys = qpx_prediction_keys;
            if (!psm_scan_format.empty()) { keys["scan_format"] = psm_scan_format; }
            if (!file_table)
            {
              OPENMS_LOG_ERROR << "Failed to add the -qpx_extended columns for " << in_file << endl;
              input_failed = true;
            }
            else
            {
              // written below with its pg view, at the same time (its mz view is written already)
              run_views["psm"] = {{file_table}, qpx_spectra_config, keys};
              psm_rows_table = rows_table;
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

          // Protein groups — independent of PSM result (built above).
          const std::string qpx_pg_file = qpx_runs_dir + "/" + basename + ".pg.parquet";
          std::shared_ptr<arrow::Table> qpx_pg_table = pg_build.get();
          // Row count via the library helper, not arrow::Table::num_rows(): Arrow is confined to
          // libOpenMS' implementation, and on Windows its symbols are dllimport, so calling the
          // member here links on Linux but not on MSVC.
          if (qpx_pg_table && ArrowIOHelpers::tableRowCount(qpx_pg_table) > 0)
          {
            qpx_pg_tables.push_back(qpx_pg_table);
            if (qpx_extended)
            {
              run_views["pg"] = {{qpx_pg_table}, ParquetWriteConfig{}, {}};
            }
            else if (!ProteinGroupArrowExport::exportToParquet(qpx_pg_table, qpx_pg_file))
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

          // (writeDataset() logs the view that failed)
          if (!run_views.empty())
          {
            // (writeDataset() logs the view that failed)
            auto views = std::make_shared<std::map<std::string, QPXDatasetExport::View>>(std::move(run_views));
            InputWrite input_write;
            input_write.input = i;
            input_write.in_file = in_file;
            input_write.mz_file = mz_written ? mz_file : std::string();
            input_write.rows = psm_rows_table;
            input_write.result = std::async(std::launch::async, [views, directory = qpx_runs_dir, basename]() mutable
            {
              std::pair<bool, std::set<std::string>> written;
              written.first = QPXDatasetExport::writeDataset(directory, basename, *views, &written.second);
              views.reset(); // (the tables, once written)
              return written;
            });
            input_writes.push_back(std::move(input_write));
          }
          if (psm_rows_table)
          {
            // (unless its psm file fails, see input_writes)
            qpx_psm_tables.push_back(psm_rows_table);
            qpx_psm_scan_formats.push_back(psm_scan_format);
            qpx_psm_rows[i] = ArrowIOHelpers::tableRowCount(psm_rows_table);
            psm_file_written = true;
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
        if (mz_written && !psm_file_written)
        {
          // (its rows are keyed by the spectra of the PSMs)
          OPENMS_LOG_ERROR << "Not keeping the QPX mz view, as the psm view of " << in_file << " was not written -> "
                           << mz_file << endl;
          File::remove(mz_file);
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
        input_failed_of[i] = input_failed;
      }

      // -- Write merged Parquet files --
      // Always write merged files (even for single input — provides a canonical name). With -qpx_extended the QPX ones
      // are written after the merged idXML below: they need its protein inference and the statistics of the whole run.
      bool merged_ok = true;
      // scan_format only describes the merged psm file if every contributing run agreed on a
      // convention; otherwise it is omitted rather than taken from an arbitrary run.
      std::map<std::string, std::string> qpx_psm_keys;
      if (!out_qpx_dir.empty())
      {
        std::set<std::string> distinct_formats(qpx_psm_scan_formats.begin(), qpx_psm_scan_formats.end());
        if (distinct_formats.size() == 1 && !distinct_formats.begin()->empty())
        {
          qpx_psm_keys["scan_format"] = *distinct_formats.begin();
        }
      }
      if (!out_qpx_dir.empty() && !qpx_extended)
      {
        // Concatenation inherits the first input's schema metadata, i.e. that per-run file's
        // uuid/creation_date. Stamp a fresh identity so the merged file is its own QPX file.
        auto psm_meta = ArrowIOHelpers::qpxFileMetadata("psm_file", ParquetWriteConfig{}, qpx_psm_keys);
        auto pg_meta  = ArrowIOHelpers::qpxFileMetadata("pg_file");
        merged_ok = psm_meta && pg_meta && merged_ok;
        if (psm_meta)
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_psm_tables, out_qpx_dir + "/" + QPX_AGGREGATE_STEM + ".psm.parquet", ParquetWriteConfig{}, psm_meta) && merged_ok;
        }
        // Same rule as the per-file pg above: with no group anywhere there is nothing to merge, and
        // an empty merged table would be indistinguishable from an inference that found nothing.
        if (pg_meta && !qpx_pg_tables.empty())
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_pg_tables, out_qpx_dir + "/" + QPX_AGGREGATE_STEM + ".pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
        }
      }
      // -qpx_extended with several inputs and -out_merged: the protein groups of the pooled inference for quantms.pg.parquet
      std::shared_ptr<arrow::Table> pooled_pg_table;
      const bool pooled_pg = qpx_extended && !out_merged.empty() && in_list.size() > 1;

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

      // -- -qpx_extended: the QPX dataset of the whole run, quantms.psm/pg/run/sample/provenance.parquet --
      auto write_aggregate = [&]() -> bool
      {
        mfres.shared.seconds_total = sw_total.getClockTime(); // (for the provenance view; the report's below)
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
          auto empty = QPXFile::exportPSMsToQPXArrow({}, {}, /*export_all_psms=*/false);
          if (empty) { empty = ArrowIOHelpers::withColumns(empty, QPXDatasetExport::psmColumns({}, {}, true)); }
          if (empty) { qpx_psm_tables.push_back(empty); }
        }
        // the protein groups: of the pooled inference (if it found any), else those of the inputs
        std::vector<std::shared_ptr<arrow::Table>> pg_tables = qpx_pg_tables;
        if (pooled_pg)
        {
          pg_tables.clear();
          if (pooled_pg_table && ArrowIOHelpers::tableRowCount(pooled_pg_table) > 0) { pg_tables.push_back(pooled_pg_table); }
        }
        // the run, sample and provenance views: the inputs and the steps of the run
        auto run_view = QPXDatasetExport::runView(runs_(mfres.shared, in_list));
        if (run_view) { run_view = ArrowIOHelpers::withColumns(run_view, runExtraColumns_(mfres, in_list, qpx_psm_rows, rescore, test_mode_)); }
        std::vector<QPXDatasetExport::Sample> samples(in_list.size());
        for (Size i = 0; i < in_list.size(); ++i) { samples[i].accession = File::stemName(in_list[i]); }
        const std::map<std::string, QPXDatasetExport::View> views = {
          {"psm", {qpx_psm_tables, qpx_extended_config, qpx_psm_keys}},
          {"pg", {pg_tables, ParquetWriteConfig{}, {}}},
          {"run", {{run_view}, ParquetWriteConfig{}, {}}},
          {"sample", {{QPXDatasetExport::sampleView(samples)}, ParquetWriteConfig{}, {}}},
          {"provenance", {{QPXDatasetExport::provenanceView(provenanceSteps_(mfres, search_params, rescore, protein_fdr_applied,
                                                                              pooled_pg, mod_delta_masses, mod_ptms))},
                          ParquetWriteConfig{}, {}}}};
        return QPXDatasetExport::writeDataset(out_qpx_dir, QPX_AGGREGATE_STEM, views);
      };
      // The writes of the inputs' files, once done: an input whose psm file failed fails, its mz view is not kept and its
      // rows are taken out of quantms.psm.parquet (which is then written again)
      auto finish_input_writes = [&]() -> bool
      {
        bool rows_removed = false;
        for (InputWrite& input_write : input_writes)
        {
          std::pair<bool, std::set<std::string>> written;
          try
          {
            written = input_write.result.get();
          }
          catch (const std::exception& e)
          {
            OPENMS_LOG_ERROR << "Failed to write QPX output for " << input_write.in_file << " -> " << out_qpx_dir
                             << " (Arrow/Parquet): " << e.what() << endl;
            written = {false, {"psm"}};
          }
          if (!written.first)
          {
            OPENMS_LOG_ERROR << "Failed to write the QPX files of " << input_write.in_file << " -> " << qpx_runs_dir << endl;
            if (!input_failed_of[input_write.input]) { input_failed_of[input_write.input] = true; ++failed_count; }
          }
          if (written.second.count("psm") == 0 || !input_write.rows) { continue; }
          const auto rows = std::find(qpx_psm_tables.begin(), qpx_psm_tables.end(), input_write.rows);
          if (rows != qpx_psm_tables.end())
          {
            qpx_psm_scan_formats.erase(qpx_psm_scan_formats.begin() + (rows - qpx_psm_tables.begin()));
            qpx_psm_tables.erase(rows);
          }
          qpx_psm_rows[input_write.input] = 0;
          if (!input_write.mz_file.empty())
          {
            // (its rows are keyed by the spectra of the PSMs)
            OPENMS_LOG_ERROR << "Not keeping the QPX mz view, as the psm view of " << input_write.in_file
                             << " could not be written -> " << input_write.mz_file << endl;
            File::remove(input_write.mz_file);
          }
          rows_removed = true;
        }
        input_writes.clear();
        return rows_removed;
      };
      if (qpx_extended)
      {
        // (written while the inputs' files are still written)
        const bool aggregate_written = write_aggregate();
        if (finish_input_writes())
        {
          qpx_psm_keys.clear();
          std::set<std::string> distinct_formats(qpx_psm_scan_formats.begin(), qpx_psm_scan_formats.end());
          if (distinct_formats.size() == 1 && !distinct_formats.begin()->empty()) { qpx_psm_keys["scan_format"] = *distinct_formats.begin(); }
          merged_ok = write_aggregate() && merged_ok;
        }
        else
        {
          merged_ok = aggregate_written && merged_ok;
        }
      }

      // =====================================================================
      // End-of-search report (command line + optional YAML).
      // =====================================================================
      {
        sw_total.stop();
        mfres.shared.seconds_total = sw_total.getClockTime();
        mfres.shared.seconds_output = mfres.shared.seconds_total - seconds_searching;
        mfres.shared.seconds_qpx_extended = sw_qpx_extended.getClockTime();
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
            r << "[ProSE] Total time   : " << ProSEAlgorithm::renderTotalTime(sh) << "\n";
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
