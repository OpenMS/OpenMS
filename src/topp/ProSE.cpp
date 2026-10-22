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
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>
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
#include <optional>
#include <sstream>

#include <map>
#include <set>
#include <tuple>
#include <unordered_map>

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
@note Extended QPX output: '-qpx_extended' (with '-out_qpx') adds optional columns to the QPX psm view and the
statistics of the run to the footer of its Parquet files, so that a viewer can show a search from the QPX files alone.
QPX readers ignore columns and keys they do not know; the names are those proposed in bigbio/qpx#355 and #259. Each
row of the per-input psm.parquet files and of quantms.psm.parquet gets: 'additional_scores' with every numeric
score and meta value of the PSM in snake_case instead of the known scores only ('global_qvalue': the PSM q-value of its
input file, from target-decoy competition or Percolator; 'ln_hyperscore'; Percolator's 'percolator_score',
'percolator_qvalue', 'percolator_pep'; 'precursor_mz_error_ppm' (corrected for the isotope error), 'delta_mass' (open
search) and the Percolator features), with 'higher_better' where it is known; 'spectrum_reference' (string, the native
ID of the spectrum); 'isotope_error' (int16: the 13C peaks from the monoisotopic to the selected precursor peak, i.e.
observed minus theoretical); 'protein_mappings' (list of struct accession, start, end, aa_before, aa_after: each protein
of the peptide with the peptide's 1-based first and last position and the residues before and after it, '-' at a
protein terminus); 'peak_mz_array' and 'peak_intensity_array' (float32, by m/z: the peaks of the spectrum as the input
file has them, filtered: the '-qpx_peak_top' most intense peaks of each '-qpx_peak_window' Th window (WindowMower,
jumping windows) and every peak of the fragment annotation 'mz_array'; '-qpx_peak_top 0' keeps all); and
'predicted_ion_type_array', 'predicted_mz_array', 'predicted_intensity_array': the b and y ions of charge 1 and (for
precursors of charge 2 or more) 2, named as in 'ion_type_array' ("b3+", "y5++"), with their theoretical m/z and the
intensity PeptDeep predicts, relative to the most intense of them; ions of 0.01 or less are left out. The predictions
use the settings of 'Search:peptdeep:*' (models, instrument) at the normalised collision energy of the file: the one
PeptDeep rescoring used ('Search:peptdeep:enable'; predicted again at that NCE), else 'Search:peptdeep:nce' if set,
else the one PeptDeep rescoring selects, run on a copy of the written PSMs; PSMs of precursor charge 1 to 6 are
predicted, others get null. Each input is read once more after its search for the peaks; a failed re-read or
prediction is a warning and leaves the columns null; without ONNX Runtime the prediction columns are null. The footer
of every psm.parquet of the run (they are written at its end) has the same keys about the whole run:
'prose_database_file', 'prose_enzyme', 'prose_enzyme_specificity' (none, semi, full), 'prose_clip_nterm_methionine',
'prose_decoy_mode', 'prose_decoy_string', 'prose_decoy_position' (prefix, suffix), 'prose_open_search',
'prose_calibration', 'prose_chunked' (true, false), 'prose_db_target_proteins', 'prose_db_decoy_proteins',
'prose_indexed_peptides', 'prose_indexed_fragments', 'prose_seconds_index_build', 'prose_seconds_total',
'prose_min_charge', 'prose_max_charge', 'prose_psm_fdr', 'prose_protein_fdr', 'prose_protein_fdr_applied' (the protein
FDR was applied to a complete protein set: one input, or the pooled set of -out_merged), 'prose_extra_features' (the
Percolator features beyond the standard ones, by their names in 'additional_scores', comma-separated),
'prose_peak_window', 'prose_peak_top' (the peak filter), 'prose_mod_analysis_delta_masses' and 'prose_mod_analysis_ptms'
(open search with -out_mod_analysis_dir: the file names of the tables of the run in that directory, of one input its
own, of several the aggregate), and 'prose_files', a JSON array with an object per input in input order: input_file,
exit_code, ms2_spectra, target_psms, decoy_psms, percolator (applied, skipped_reason, warning), peptdeep_nce (null:
none), precursor_error and fragment_error (median_ppm, mad_ppm, recommended_ppm; null: no estimate), seconds_search.
Protein groups: with '-qpx_extended' a run without protein inference gets one group per protein in its pg.parquet
(global_qvalue NaN), and with '-out_merged' and several inputs quantms.pg.parquet holds the groups of the pooled
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

      registerOutputDir_("out_qpx", "<dir>", "", "Output directory for QPX exchange format Parquet files. Writes per-input <basename>.psm.parquet and <basename>.pg.parquet, plus merged quantms.psm.parquet and quantms.pg.parquet (see also -qpx_extended). The pg files are written only when protein groups were actually inferred, which needs 'FDR:protein' > 0 together with decoys, or '-out_merged' with more than one input; an identification-only run produces no pg file rather than an empty one.", false, true);

      registerOutputDir_("out_parquet", "<dir>", "", "Output directory for OpenMS internal format Parquet files. Writes per-input <basename>.psm/proteins/pg/search_params.parquet, plus merged openms.* files.", false, true);

      registerOutputFile_("out_merged", "<file>", "", "Optional merged output file containing all PSMs pooled across input files with cross-file protein inference (BasicProteinInferenceAlgorithm) and optional picked-protein FDR (FDR:protein). Per-file outputs (-out_idxml / -out_qpx / -out_parquet) retain run-level information. Only useful with multiple -in files.", false);
      setValidFormats_("out_merged", ListUtils::create<std::string>("idXML"));

      registerOutputFileList_("out_pin", "<files>", StringList(), "Output Percolator input (.pin/.tsv) file(s) for external rescoring. Must have the same number of entries as -in. It holds the PSMs before the rescoring and before any FDR filter, with ProSE's own scores and all reported hits: 'Search:rescore' and 'Search:FDR:*' apply to the other outputs only.", false, true);
      setValidFormats_("out_pin", ListUtils::create<std::string>("tsv"));

      registerOutputDir_("out_mod_analysis_dir", "<dir>", "", "Optional directory to write modification-analysis tables (delta-mass, PTM stats) when running in open-search mode. When set, per-file tables are written using each input file's basename and an additional aggregate table is written across all input files. Has no effect in closed-search mode.", false, true);

      registerOutputFile_("summary_out", "<file>", "", "Optional YAML file capturing the end-of-search report (per-file identification statistics, shared configuration/database/index facts, timing and output manifest) for pipeline ingestion. The same information is printed to the command line unless -no_summary is set.", false, true);
      setValidFormats_("summary_out", ListUtils::create<std::string>("yaml"));

      registerFlag_("no_summary", "Suppress the end-of-search report printed to the command line. The -summary_out YAML file (if requested) is still written.", true);

      registerFlag_("qpx_extended", "With -out_qpx: add optional columns to the QPX psm view (all scores in 'additional_scores', 'spectrum_reference', 'isotope_error', 'protein_mappings', the filtered peaks of the spectrum and PeptDeep's predicted b/y intensities) and the statistics of the run to the footer of its psm.parquet files ('prose_' keys); a run without protein inference gets one protein group per protein. See the tool documentation.", true);
      registerDoubleOption_("qpx_peak_window", "<Th>", 100.0, "-qpx_extended: the m/z window of the peak filter.", false, true);
      setMinFloat_("qpx_peak_window", 1.0);
      registerIntOption_("qpx_peak_top", "<n>", 300, "-qpx_extended: the most intense peaks kept per -qpx_peak_window (with every annotated peak); 0 keeps all peaks.", false, true);
      setMinInt_("qpx_peak_top", 0);
      registerFlag_("keep_decoys", "Keep the decoy PSMs and proteins that pass the FDR thresholds in the outputs; the protein FDR (one input, or -out_merged) otherwise removes them.", true);


      // put search algorithm parameters at Search: subtree of parameters
      Param search_algo_params_with_subsection;
      search_algo_params_with_subsection.insert("Search:", ProSEAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    /// The spectrum columns of one PSM row of -qpx_extended; a column is null when its flag is false
    struct PSMSpectrum_
    {
      bool has_peaks = false;
      std::vector<float> peak_mz, peak_intensity;
      bool has_prediction = false;
      std::vector<std::string> ion_type;
      std::vector<float> mz, intensity;
    };

    /// The peak filter of -qpx_extended: the @p top most intense peaks per @p window Th (0: all), and the annotated ones
    struct PeakFilter_
    {
      double window = 100.0;
      Size top = 300;
    };

    /**
      @brief The indices of the peaks of @p spectrum (sorted by m/z) that the filter keeps: the most intense ones per
      window (WindowMower, jumping windows) and those nearest to the m/z of the fragment annotation @p annotations
      (within 0.01 Th). All peaks if @p filter.top is 0.
    */
    static std::vector<Size> keptPeaks_(const MSSpectrum& spectrum, const std::vector<PeptideHit::PeakAnnotation>& annotations,
                                        const PeakFilter_& filter)
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
        MSSpectrum mowed = spectrum;
        WindowMower mower;
        Param p = mower.getParameters();
        p.setValue("windowsize", filter.window);
        p.setValue("peakcount", static_cast<int>(filter.top));
        p.setValue("movetype", "jump");
        mower.setParameters(p);
        mower.filterPeakSpectrum(mowed);
        for (const Peak1D& peak : mowed) { keep[indexOf(peak.getMZ())] = true; }
        for (const PeptideHit::PeakAnnotation& a : annotations)
        {
          const Size i = indexOf(a.mz);
          if (std::fabs(spectrum[i].getMZ() - a.mz) <= 0.01) { keep[i] = true; }
        }
      }
      std::vector<Size> kept;
      for (Size i = 0; i < keep.size(); ++i) { if (keep[i]) { kept.push_back(i); } }
      return kept;
    }

    /**
      @brief The filtered peaks and PeptDeep predictions of the PSMs of one searched file (-qpx_extended): one entry
      per element of @p peptide_ids, for its best hit (the row of the QPX psm view).

      Never throws because of the spectra or the predictions: a failed re-read of @p in_file or a failed prediction is
      logged as a warning and leaves the columns null; the search result of the file is not affected. Only
      std::bad_alloc propagates: out of memory is not taken for an unreadable file.

      The peaks are read again from @p in_file (MS2 spectra, by native ID; of several spectra with the PSM's native
      ID the one at its RT) and filtered with @p filter. The predictions are PeptDeep's at the NCE of the file: the one
      PeptDeep rescoring stored on the run ('peptdeep_nce'), else 'peptdeep:nce' if set, else the one
      PeptDeepRescoring::annotate() selects on a copy of @p peptide_ids; the NCE is then stored on the run.
    */
    std::vector<PSMSpectrum_> collectPSMSpectra_(const std::string& in_file, const Param& search_params,
                                                 std::vector<ProteinIdentification>& protein_ids,
                                                 const PeptideIdentificationList& peptide_ids, const PeakFilter_& filter)
    {
      std::vector<PSMSpectrum_> out(peptide_ids.size());

      // -- peaks: the MS2 spectra of the PSMs, as the file has them --
      std::unordered_map<std::string, std::vector<Size>> psms_of_spectrum; // native ID -> PSMs
      for (Size i = 0; i < peptide_ids.size(); ++i)
      {
        if (!peptide_ids[i].getHits().empty()) { psms_of_spectrum[peptide_ids[i].getSpectrumReference()].push_back(i); }
      }
      PeakMap spectra; // the spectra of the PSMs (several, if a native ID is not unique in the file)
      Size with_peaks = 0, peaks_read = 0, peaks_kept = 0;
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
                        << "); the peak columns of its PSMs are null." << endl;
        spectra.clear(true);
      }
      std::unordered_map<std::string, std::vector<Size>> spectra_of_id;
      for (Size s = 0; s < spectra.size(); ++s) { spectra_of_id[spectra[s].getNativeID()].push_back(s); }
      Size without_spectrum = 0;
      for (const auto& [native_id, psms] : psms_of_spectrum)
      {
        const auto found = spectra_of_id.find(native_id);
        for (Size i : psms)
        {
          const MSSpectrum* spectrum = nullptr;
          if (found != spectra_of_id.end() && found->second.size() == 1) { spectrum = &spectra[found->second.front()]; }
          else if (found != spectra_of_id.end())
          {
            for (Size s : found->second)
            {
              if (std::fabs(spectra[s].getRT() - peptide_ids[i].getRT()) <= 1e-6) { spectrum = &spectra[s]; break; }
            }
          }
          if (spectrum == nullptr) { ++without_spectrum; continue; }
          PSMSpectrum_& entry = out[i];
          entry.has_peaks = true;
          ++with_peaks;
          const std::vector<Size> kept = keptPeaks_(*spectrum, peptide_ids[i].getHits()[0].getPeakAnnotations(), filter);
          peaks_read += spectrum->size();
          peaks_kept += kept.size();
          entry.peak_mz.reserve(kept.size());
          entry.peak_intensity.reserve(kept.size());
          for (Size k : kept)
          {
            entry.peak_mz.push_back(static_cast<float>((*spectrum)[k].getMZ()));
            entry.peak_intensity.push_back((*spectrum)[k].getIntensity());
          }
        }
      }
      if (without_spectrum > 0 && !spectra.empty())
      {
        OPENMS_LOG_WARN << "-qpx_extended: " << without_spectrum << " PSM(s) of " << in_file
                        << " have no MS2 spectrum of their native ID (and RT) in the file; their peaks are null." << endl;
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
                        << "the peaks are written and the predicted_* columns are null." << endl;
        told = true;
      }
#endif
      OPENMS_LOG_INFO << "-qpx_extended: " << in_file << ": peaks of " << with_peaks << " PSM(s), " << peaks_kept
                      << " of " << peaks_read << " peaks kept";
      if (filter.top > 0)
      {
        OPENMS_LOG_INFO << " (the " << filter.top << " most intense per " << filter.window << " Th and the annotated ones)";
      }
      OPENMS_LOG_INFO << endl;
      return out;
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

    /**
      @brief The columns -qpx_extended adds to the QPX psm view of a file, one entry per row of the view (the best hit
      of each element of @p peptide_ids with hits): 'additional_scores' (replaced: every numeric score and meta value
      of the PSM in snake_case), 'spectrum_reference', 'isotope_error', 'protein_mappings', the peaks and the
      predictions of @p spectra.
    */
    static std::vector<ArrowIOHelpers::ParquetColumn> extendedColumns_(const std::vector<ProteinIdentification>& protein_ids,
                                                                       const PeptideIdentificationList& peptide_ids,
                                                                       const std::vector<PSMSpectrum_>& spectra)
    {
      using Column = ArrowIOHelpers::ParquetColumn;
      auto column = [](const std::string& name, Column::Type type) { Column c; c.name = name; c.type = type; return c; };
      Column scores = column(QPXPSMSchema::ADDITIONAL_SCORES, Column::Type::SCORE_LIST);
      Column reference = column("spectrum_reference", Column::Type::UTF8);
      Column isotope = column("isotope_error", Column::Type::INT16);
      Column mappings = column("protein_mappings", Column::Type::PROTEIN_MAPPING_LIST);
      Column peak_mz = column("peak_mz_array", Column::Type::FLOAT_LIST);
      Column peak_intensity = column("peak_intensity_array", Column::Type::FLOAT_LIST);
      Column ion_type = column("predicted_ion_type_array", Column::Type::UTF8_LIST);
      Column mz = column("predicted_mz_array", Column::Type::FLOAT_LIST);
      Column intensity = column("predicted_intensity_array", Column::Type::FLOAT_LIST);

      // The isotope error as QPX has it (bigbio/qpx#355: the 13C peaks from the monoisotopic to the selected peak, i.e.
      // observed minus theoretical); ProSE's meta value follows 'report:isotope_error_convention'.
      int isotope_sign = 1;
      if (!protein_ids.empty() && protein_ids[0].getSearchParameters().metaValueExists("isotope_error_convention")
          && protein_ids[0].getSearchParameters().getMetaValue("isotope_error_convention").toString() == "theoretical_minus_observed")
      {
        isotope_sign = -1;
      }
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

        reference.utf8_values.push_back(id.getSpectrumReference());
        reference.null_rows.push_back(id.getSpectrumReference().empty());

        const std::string isotope_key = hit.metaValueExists("isotope_error") ? "isotope_error" : "IsotopeError";
        isotope.int64_values.push_back(hit.metaValueExists(isotope_key) ? isotope_sign * static_cast<int>(hit.getMetaValue(isotope_key)) : 0);
        isotope.null_rows.push_back(!hit.metaValueExists(isotope_key));

        std::vector<Column::ProteinMapping> row_mappings;
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          row_mappings.push_back({evidence.getProteinAccession(), position(evidence.getStart()), position(evidence.getEnd()),
                                  residue(evidence.getAABefore()), residue(evidence.getAAAfter())});
        }
        mappings.protein_mapping_lists.push_back(std::move(row_mappings));

        const PSMSpectrum_& s = spectra[i];
        peak_mz.float_lists.push_back(s.peak_mz);
        peak_mz.null_rows.push_back(!s.has_peaks);
        peak_intensity.float_lists.push_back(s.peak_intensity);
        peak_intensity.null_rows.push_back(!s.has_peaks);
        ion_type.utf8_lists.push_back(s.ion_type);
        ion_type.null_rows.push_back(!s.has_prediction);
        mz.float_lists.push_back(s.mz);
        mz.null_rows.push_back(!s.has_prediction);
        intensity.float_lists.push_back(s.intensity);
        intensity.null_rows.push_back(!s.has_prediction);
      }
      return {scores, reference, isotope, mappings, peak_mz, peak_intensity, ion_type, mz, intensity};
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

    /**
      @brief The footer metadata of the psm.parquet files of -qpx_extended: the 'prose_' keys about the whole run (see
      the tool documentation)
    */
    static std::map<std::string, std::string> runMetadata_(const ProSEAlgorithm::MultiFileSearchResult& mfres,
                                                           const PeakFilter_& filter, bool protein_fdr_applied,
                                                           const std::string& mod_delta_masses, const std::string& mod_ptms)
    {
      const ProSEAlgorithm::SharedSearchStats& sh = mfres.shared;
      auto flag = [](bool b) -> std::string { return b ? "true" : "false"; };
      std::map<std::string, std::string> meta;
      meta["prose_database_file"] = sh.database_file;
      meta["prose_enzyme"] = sh.enzyme;
      meta["prose_decoy_mode"] = sh.decoy_mode;
      meta["prose_decoy_string"] = mfres.decoy_string;
      meta["prose_decoy_position"] = mfres.decoy_is_prefix ? "prefix" : "suffix";
      meta["prose_open_search"] = flag(sh.open_search);
      meta["prose_calibration"] = flag(sh.calibration_enabled);
      meta["prose_chunked"] = flag(sh.chunked);
      meta["prose_db_target_proteins"] = std::to_string(sh.db_target_proteins);
      meta["prose_db_decoy_proteins"] = std::to_string(sh.db_decoy_proteins);
      meta["prose_indexed_peptides"] = std::to_string(sh.indexed_peptides);
      meta["prose_indexed_fragments"] = std::to_string(sh.indexed_fragments);
      meta["prose_seconds_index_build"] = number_(sh.seconds_index_build, "0");
      meta["prose_seconds_total"] = number_(sh.seconds_total, "0");
      meta["prose_min_charge"] = std::to_string(sh.min_charge);
      meta["prose_max_charge"] = std::to_string(sh.max_charge);
      meta["prose_psm_fdr"] = number_(sh.psm_fdr_threshold, "0");
      meta["prose_protein_fdr"] = number_(sh.protein_fdr_threshold, "0");
      meta["prose_protein_fdr_applied"] = flag(protein_fdr_applied);
      meta["prose_peak_window"] = number_(filter.window, "0");
      meta["prose_peak_top"] = std::to_string(filter.top);
      if (!mod_delta_masses.empty()) { meta["prose_mod_analysis_delta_masses"] = mod_delta_masses; }
      if (!mod_ptms.empty()) { meta["prose_mod_analysis_ptms"] = mod_ptms; }

      // the search settings of the first searched file's run (all files share them)
      for (const auto& pf : mfres.per_file)
      {
        if (pf.protein_ids.empty()) { continue; }
        const ProteinIdentification::SearchParameters& sp = pf.protein_ids[0].getSearchParameters();
        meta["prose_enzyme_specificity"] = EnzymaticDigestion::NamesOfSpecificity[sp.enzyme_term_specificity];
        meta["prose_clip_nterm_methionine"] =
          sp.metaValueExists("peptide:clip_nterm_methionine") ? sp.getMetaValue("peptide:clip_nterm_methionine").toString() : "false";
        std::vector<std::string> features;
        if (sp.metaValueExists("extra_features"))
        {
          for (const std::string& feature : ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString()))
          {
            const std::string name = qpxScoreName_(feature);
            if (!name.empty() && feature != "score") { features.push_back(name); }
          }
        }
        meta["prose_extra_features"] = ListUtils::concatenate(features, ",");
        break;
      }

      // per input file, in input order
      auto error = [&](bool valid, double median, double mad, double recommended) -> std::string
      {
        if (!valid) { return "null"; }
        return "{\"median_ppm\":" + number_(median) + ",\"mad_ppm\":" + number_(mad) + ",\"recommended_ppm\":" + number_(recommended) + "}";
      };
      std::string files = "[";
      for (Size i = 0; i < mfres.per_file.size(); ++i)
      {
        const auto& pf = mfres.per_file[i];
        const auto& st = pf.stats;
        double nce = std::numeric_limits<double>::quiet_NaN();
        if (!pf.protein_ids.empty() && pf.protein_ids[0].metaValueExists("peptdeep_nce"))
        {
          nce = static_cast<double>(pf.protein_ids[0].getMetaValue("peptdeep_nce"));
        }
        files += (i > 0 ? "," : "");
        files += "{\"input_file\":" + jsonString_(st.input_file)
               + ",\"exit_code\":" + std::to_string(static_cast<int>(pf.exit_code))
               + ",\"ms2_spectra\":" + std::to_string(st.ms2_spectra)
               + ",\"target_psms\":" + std::to_string(st.target_psms)
               + ",\"decoy_psms\":" + std::to_string(st.decoy_psms)
               + ",\"percolator\":{\"applied\":" + flag(pf.percolator.applied)
               + ",\"skipped_reason\":" + jsonString_(pf.percolator.skipped_reason)
               + ",\"warning\":" + jsonString_(pf.percolator.warning) + "}"
               + ",\"peptdeep_nce\":" + number_(nce)
               + ",\"precursor_error\":" + error(st.prec_tol_valid, st.prec_err_median, st.prec_err_mad, st.prec_err_recommended)
               + ",\"fragment_error\":" + error(st.frag_tol_valid, st.frag_err_median, st.frag_err_mad, st.frag_err_recommended)
               + ",\"seconds_search\":" + number_(st.seconds_search, "0") + "}";
      }
      meta["prose_files"] = files + "]";
      return meta;
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

      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      Param search_params = getParam_().copy("Search:", true);
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
              // survived inference.
              // The picked-protein FDR is applied unless no decoy protein is left for it, i.e. no PSM refers to one (the
              // inference keeps the proteins the PSMs refer to). The protein score type cannot tell: the inference
              // already scores the proteins with the PSMs' q-values.
              for (const PeptideIdentification& id : result.peptide_ids)
              {
                for (const PeptideHit& hit : id.getHits())
                {
                  for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
                  {
                    const std::string& accession = evidence.getProteinAccession();
                    protein_fdr_applied = protein_fdr_applied || (mfres.decoy_is_prefix
                      ? StringUtils::hasPrefix(accession, mfres.decoy_string) : StringUtils::hasSuffix(accession, mfres.decoy_string));
                  }
                }
              }
              ProSEAlgorithm::applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids,
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
      // -qpx_extended: the per-input QPX psm files, written at the end of the run with the run's statistics in their
      // footer (the merged file needs the tables anyway, so they are in memory either way)
      struct PendingPSMFile { Size input; std::string file; std::string scan_format; std::shared_ptr<arrow::Table> table; };
      std::vector<PendingPSMFile> pending_psm_files;
      std::vector<bool> input_counted_failed(in_list.size(), false);
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
        std::vector<PSMSpectrum_> spectra_of_file;
        if (qpx_extended)
        {
          try
          {
            spectra_of_file = collectPSMSpectra_(in_file, search_params, result.protein_ids, result.peptide_ids, peak_filter);
          }
          catch (const std::exception& e)
          {
            OPENMS_LOG_ERROR << "Failed to collect the spectra of the PSMs of " << in_file << ": " << e.what() << endl;
            spectra_of_file.assign(result.peptide_ids.size(), PSMSpectrum_());
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
          const std::string qpx_psm_file = out_qpx_dir + "/" + basename + ".psm.parquet";
          auto qpx_psm_table = QPXFile::exportPSMsToQPXArrow(result.protein_ids, result.peptide_ids, /*export_all_psms=*/false);
          if (qpx_psm_table && qpx_extended)
          {
            qpx_psm_table = ArrowIOHelpers::withColumns(qpx_psm_table,
                              extendedColumns_(result.protein_ids, result.peptide_ids, spectra_of_file));
            if (!qpx_psm_table)
            {
              OPENMS_LOG_ERROR << "Failed to add the -qpx_extended columns for " << in_file << endl;
              input_failed = true;
            }
          }
          if (qpx_psm_table)
          {
            qpx_psm_tables.push_back(qpx_psm_table);
            std::vector<std::string> spec_refs;
            spec_refs.reserve(result.peptide_ids.size());
            for (const auto& pid : result.peptide_ids)
            {
              if (!pid.getSpectrumReference().empty()) { spec_refs.push_back(pid.getSpectrumReference()); }
            }
            const std::string psm_scan_format = ArrowIOHelpers::qpxScanFormat(spec_refs);
            qpx_psm_scan_formats.push_back(psm_scan_format);
            if (qpx_extended)
            {
              pending_psm_files.push_back({i, qpx_psm_file, psm_scan_format, qpx_psm_table});
            }
            else if (!QPXFile::exportToParquet(qpx_psm_table, qpx_psm_file, ParquetWriteConfig{}, psm_scan_format))
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
          const std::string qpx_pg_file = out_qpx_dir + "/" + basename + ".pg.parquet";
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
          const std::string oms_psm_file = out_parquet_dir + "/" + basename + ".psm.parquet";
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
          const std::string oms_prot_file = out_parquet_dir + "/" + basename + ".proteins.parquet";
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
          const std::string oms_pg_file = out_parquet_dir + "/" + basename + ".pg.parquet";
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
          const std::string oms_sp_file = out_parquet_dir + "/" + basename + ".search_params.parquet";
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

        if (input_failed) { ++failed_count; input_counted_failed[i] = true; }
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
        auto psm_meta = ArrowIOHelpers::qpxFileMetadata("psm_file", ParquetWriteConfig{}, psm_extra);
        auto pg_meta  = ArrowIOHelpers::qpxFileMetadata("pg_file");
        merged_ok = psm_meta && pg_meta && merged_ok;
        if (psm_meta)
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_psm_tables, out_qpx_dir + "/quantms.psm.parquet", ParquetWriteConfig{}, psm_meta) && merged_ok;
        }
        // Same rule as the per-file pg above: with no group anywhere there is nothing to merge, and
        // an empty merged table would be indistinguishable from an inference that found nothing.
        if (pooled_pg)
        {
          if (pg_meta && pooled_pg_table && ArrowIOHelpers::tableRowCount(pooled_pg_table) > 0)
          {
            merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
              {pooled_pg_table}, out_qpx_dir + "/quantms.pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
          }
        }
        else if (pg_meta && !qpx_pg_tables.empty())
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_pg_tables, out_qpx_dir + "/quantms.pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
        }
      };

      if (!out_parquet_dir.empty())
      {
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_psm_tables,  out_parquet_dir + "/openms.psm.parquet")           && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_prot_tables, out_parquet_dir + "/openms.proteins.parquet")      && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_pg_tables,   out_parquet_dir + "/openms.pg.parquet")            && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_sp_tables,   out_parquet_dir + "/openms.search_params.parquet") && merged_ok;
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
                            merged.protein_ids, out_parquet_dir + "/merged.proteins.parquet") && merged_ok;
              merged_ok = ProteinIdentificationArrowIO::exportProteinGroupsToParquet(
                            merged.protein_ids, out_parquet_dir + "/merged.pg.parquet") && merged_ok;
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

      // -- QPX: the per-input psm files of -qpx_extended, with the statistics of the run in their footer, and the
      //    merged files --
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
          footer = runMetadata_(mfres, peak_filter, protein_fdr_applied, mod_delta_masses, mod_ptms);
        }
        for (const PendingPSMFile& pending : pending_psm_files)
        {
          std::map<std::string, std::string> keys = footer;
          if (!pending.scan_format.empty()) { keys["scan_format"] = pending.scan_format; }
          const auto meta = ArrowIOHelpers::qpxFileMetadata("psm_file", ParquetWriteConfig{}, keys);
          if (!meta || !ArrowIOHelpers::concatenateAndWriteToParquet({pending.table}, pending.file, ParquetWriteConfig{}, meta))
          {
            OPENMS_LOG_ERROR << "Failed to write QPX PSM parquet for " << in_list[pending.input] << " -> " << pending.file << endl;
            if (!input_counted_failed[pending.input])
            {
              ++failed_count;
              input_counted_failed[pending.input] = true;
            }
          }
        }
        pending_psm_files.clear();
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
            ofs << ProSEAlgorithm::renderRunSummaryYaml(mfres, manifest, failed_count, in_list.size());
            // Whether the protein FDR was applied: to the single input, or to the pooled set of -out_merged.
            ofs << "protein_fdr_applied: " << (protein_fdr_applied ? "true" : "false") << std::endl;
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
