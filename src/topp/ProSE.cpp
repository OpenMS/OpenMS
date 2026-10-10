// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/BasicProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <map>
#include <set>

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
@note Decoy reporting is tied to protein-level FDR, not to PSM-level FDR. Setting '-Search:FDR:protein' > 0 signals "finalize this result": picked-protein FDR is applied and decoys are removed. PSM-level FDR ('-Search:FDR:PSM') filters target and decoy PSMs alike by the q-value threshold; it does no decoy-specific stripping, so decoys that pass the threshold are kept. With '-Search:FDR:protein' = 0 (the default), decoys are retained in every output (target+decoy evidence with scores). To obtain a clean, decoy-free result without protein FDR, run @ref TOPP_IDFilter with '-remove_decoys' downstream.
@note Protein-level FDR scope: picked-protein FDR does not compose across runs, so it is applied (and decoys removed) only on a @em complete protein set — a single input file, or the pooled '-out_merged' set of a multi-file run. For a multi-file run with '-Search:FDR:protein' > 0 but no '-out_merged', protein FDR is NOT applied (no output represents a complete experiment); ProSE warns and leaves the per-file outputs as intermediates (decoys retained).
@note Percolator rescoring: with '-rescore', the PSMs of each input file are rescored separately with OpenMS' built-in (in-process) Percolator implementation, the same backend PercolatorAdapter uses by default. No percolator executable is needed (the former '-percolator_executable' parameter was removed). Training uses the standard Percolator feature set plus ProSE's search engine features, '-train_best_positive' and target-decoy competition ('-post_processing_tdc'); Percolator q-values become the main score (the SVM score and PEP are kept as meta values MS:1001492/MS:1001493, the HyperScore as meta value 'ln(hyperscore)'), only the best hit per spectrum is kept, and the runs report 'Percolator' as search engine. The PSM and protein FDR thresholds are then applied on the rescored PSMs. A file with fewer than 100 PSMs or without decoys, or for which Percolator fails, is @em not rescored: ProSE prints a warning, keeps its HyperScores (FDR is then computed from them) and continues; the end-of-search report states how many files were rescored.
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
    /// Identification output format: idparquet for a .idparquet name, otherwise idXML. Names
    /// without a known extension (e.g. TOPPAS '.unknown' outputs) are written as idXML; other
    /// recognised extensions are rejected by FileHandler::storeIdentifications.
    static FileTypes::Type identificationOutputType_(const std::string& path)
    {
      return FileHandler::getTypeByFileName(path) == FileTypes::IDPARQUET ? FileTypes::IDPARQUET : FileTypes::IDXML;
    }

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

      registerOutputFileList_("out", "<files>", StringList(), "Output identification file(s), one per -in: idXML or an idparquet directory bundle (format by extension). Must have the same number of entries as -in.", false);
      setValidFormats_("out", ListUtils::create<std::string>("idXML,idparquet"));

      registerOutputFileList_("out_idxml", "<files>", StringList(), "Deprecated, use -out (same behaviour). Output identification file(s): idXML or an idparquet directory bundle (format by extension). Must have the same number of entries as -in. Cannot be combined with -out.", false, true);
      setValidFormats_("out_idxml", ListUtils::create<std::string>("idXML,idparquet"));

      registerOutputDir_("out_qpx", "<dir>", "", "Output directory for QPX exchange format Parquet files. Writes per-input <basename>.psm.parquet and <basename>.pg.parquet, plus merged quantms.psm.parquet and quantms.pg.parquet. The pg files are written only when protein groups were actually inferred, which needs 'FDR:protein' > 0 together with decoys, or '-out_merged' with more than one input; an identification-only run produces no pg file rather than an empty one.", false, true);

      registerOutputDir_("out_parquet", "<dir>", "", "Output directory for OpenMS internal format Parquet files. Writes per-input <basename>.psm/proteins/pg/search_params.parquet, plus merged openms.* files.", false, true);

      registerOutputFile_("out_merged", "<file>", "", "Optional merged output file containing all PSMs pooled across input files with cross-file protein inference (BasicProteinInferenceAlgorithm) and optional picked-protein FDR (FDR:protein). Per-file outputs (-out / -out_qpx / -out_parquet) retain run-level information. Only useful with multiple -in files.", false);
      setValidFormats_("out_merged", ListUtils::create<std::string>("idXML,idparquet"));

      registerOutputFileList_("out_pin", "<files>", StringList(), "Output Percolator input (.pin/.tsv) file(s) for external rescoring. Must have the same number of entries as -in. Written independently of -rescore (i.e. you can produce .pin files without rescoring). With -rescore, the .pin is written from the rescored PSMs.", false, true);
      setValidFormats_("out_pin", ListUtils::create<std::string>("tsv"));

      registerOutputDir_("out_mod_analysis_dir", "<dir>", "", "Optional directory to write modification-analysis tables (delta-mass, PTM stats) when running in open-search mode. When set, per-file tables are written using each input file's basename and an additional aggregate table is written across all input files. Has no effect in closed-search mode.", false, true);

      registerOutputFile_("summary_out", "<file>", "", "Optional YAML file capturing the end-of-search report (per-file identification statistics, shared configuration/database/index facts, timing and output manifest) for pipeline ingestion. The same information is printed to the command line unless -no_summary is set.", false, true);
      setValidFormats_("summary_out", ListUtils::create<std::string>("yaml"));

      registerFlag_("no_summary", "Suppress the end-of-search report printed to the command line. The -summary_out YAML file (if requested) is still written.", true);

      registerFlag_("rescore", "Rescore the PSMs of each input file with Percolator before output, in-process (OpenMS' built-in Percolator implementation; no external executable needed). "
        "Rescored PSMs get Percolator q-values as main score. The PSM and protein FDR thresholds ('Search:FDR:*') are applied afterwards, on the rescored scores. "
        "A file with fewer than 100 PSMs or without decoys is not rescored and keeps its search engine scores (a warning is printed).");

      // put search algorithm parameters at Search: subtree of parameters
      Param search_algo_params_with_subsection;
      search_algo_params_with_subsection.insert("Search:", ProSEAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    /// Name of @p enzyme as PercolatorInfile computes the enzN/enzC/enzInt .pin features ("no_enzyme" if unknown).
    static std::string pinEnzymeName_(const std::string& enzyme)
    {
      std::string name;
      for (const char c : enzyme)
      {
        if (c != ' ') name += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      if (name == "trypsin/p") return "trypsinp";
      if (name == "pepsina") return "pepsin";
      static const std::set<std::string> pin_enzymes = {"elastase", "pepsin", "proteinasek", "thermolysin", "chymotrypsin",
                                                        "lys-n", "lys-c", "arg-c", "asp-n", "glu-c", "trypsin", "trypsinp"};
      return pin_enzymes.contains(name) ? name : "no_enzyme";
    }

    /**
      @brief Rescores the PSMs of one input file in-process with Percolator (see Percolator::rescorePSMs).

      Uses the settings ProSE used to run PercolatorAdapter with: PIN feature set with the
      search engine's 'extra_features', -train_best_positive, -post_processing_tdc, q-values as
      main score, PercolatorAdapter's defaults for everything else. Only the best hit per spectrum
      is kept afterwards.

      @return false (with a warning) if the file is not rescored: fewer than 100 PSMs, no decoys,
              or Percolator failed. The PSMs are left unchanged in that case.
    */
    bool rescoreWithPercolator_(const std::string& in_file,
                                std::vector<ProteinIdentification>& protein_ids,
                                PeptideIdentificationList& peptide_ids) const
    {
      if (protein_ids.empty() || peptide_ids.size() < 100)
      {
        OPENMS_LOG_WARN << "Skipping Percolator rescoring for " << in_file
                        << " (only " << peptide_ids.size() << " PSMs, need >= 100). Keeping the search engine scores." << endl;
        return false;
      }

      bool has_decoys = false;
      int min_charge = 10;
      int max_charge = 0;
      for (const PeptideIdentification& pid : peptide_ids)
      {
        for (const PeptideHit& hit : pid.getHits())
        {
          has_decoys = has_decoys || hit.isDecoy();
          min_charge = std::min(min_charge, hit.getCharge());
          max_charge = std::max(max_charge, hit.getCharge());
        }
      }
      if (!has_decoys)
      {
        OPENMS_LOG_WARN << "Skipping Percolator rescoring for " << in_file << ": no decoy PSMs for target/decoy competition. "
                        << "Use '-Search:decoys auto' / 'generate' or provide a FASTA with decoy proteins. Keeping the search engine scores." << endl;
        return false;
      }

      const ProteinIdentification::SearchParameters& sp = protein_ids.front().getSearchParameters();
      StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
      if (sp.metaValueExists("extra_features"))
      {
        StringList extra = ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString());
        feature_set.insert(feature_set.end(), extra.begin(), extra.end());
      }
      else
      {
        feature_set.push_back("score");
      }
      feature_set.push_back("Peptide");
      feature_set.push_back("Proteins");

      Percolator perc;
      Param pp = perc.getDefaults();
      pp.setValue("train_best_positive", "true");
      pp.setValue("post_processing_tdc", "true");
      pp.setValue("subset_max_train", 0);  // train on all PSMs
      pp.setValue("pep_method", "nonparametric");  // as the percolator executable
      // one model per cross-validation fold: more than 3 threads do not help
      const int threads = getIntOption_("threads");
      pp.setValue("num_threads", threads <= 0 ? 3 : std::min(threads, 3));
      perc.setParameters(pp);

      OPENMS_LOG_INFO << "[ProSE] Rescoring " << in_file << " with Percolator..." << endl;
      try
      {
        perc.rescorePSMs(peptide_ids, feature_set, pinEnzymeName_(sp.digestion_enzyme.getName()),
                         min_charge, max_charge, "q-value");
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Percolator rescoring failed for " << in_file << " (" << e.what()
                        << "). Keeping the search engine scores." << endl;
        return false;
      }
      IDFilter::keepNBestHits(peptide_ids, 1);

      // Mark the runs as post-processed by Percolator, as PercolatorAdapter does.
      for (ProteinIdentification& run : protein_ids)
      {
        run.setSearchEngine("Percolator");
        run.setSearchEngineVersion("3.08-vendored");
        run.setMetaValue("percolator", "ProSE");
        ProteinIdentification::SearchParameters run_sp = run.getSearchParameters();
        run_sp.setMetaValue("Percolator:peptide_level_fdrs", "false");
        run_sp.setMetaValue("Percolator:protein_level_fdrs", "false");
        run_sp.setMetaValue("Percolator:testFDR", static_cast<double>(pp.getValue("test_fdr")));
        run_sp.setMetaValue("Percolator:trainFDR", static_cast<double>(pp.getValue("train_fdr")));
        run_sp.setMetaValue("Percolator:maxiter", static_cast<int>(pp.getValue("num_iterations")));
        run_sp.setMetaValue("Percolator:subset_max_train", static_cast<int>(pp.getValue("subset_max_train")));
        run_sp.setMetaValue("Percolator:seed", static_cast<int>(pp.getValue("seed")));
        run_sp.setMetaValue("Percolator:post_processing_tdc", "true");
        run_sp.setMetaValue("Percolator:train_best_positive", "true");
        run.setSearchParameters(run_sp);
      }
      return true;
    }

    ExitCodes main_(int, const char**) override
    {
      const StringList in_list = getStringList_("in");
      const std::string database = getStringOption_("database");
      // -out_idxml is the deprecated, idXML-only predecessor of -out.
      const StringList out_new_list = getStringList_("out");
      const StringList out_idxml_old_list = getStringList_("out_idxml");
      if (!out_new_list.empty() && !out_idxml_old_list.empty())
      {
        OPENMS_LOG_ERROR << "-out_idxml is a deprecated alias of -out; give only one of them." << endl;
        return ILLEGAL_PARAMETERS;
      }
      if (!out_idxml_old_list.empty())
      {
        OPENMS_LOG_WARN << "-out_idxml is deprecated; use -out (accepts idXML and idparquet)." << endl;
      }
      const StringList& out_id_list = out_new_list.empty() ? out_idxml_old_list : out_new_list;
      const StringList out_pin_list = getStringList_("out_pin");
      const std::string out_merged = getStringOption_("out_merged");
      const std::string out_qpx_dir = getOutputDirOption("out_qpx");
      const std::string out_parquet_dir = getOutputDirOption("out_parquet");
      // getOutputDirOption auto-creates the directory if it doesn't exist (and returns "" if unset).
      const std::string out_mod_analysis_dir = getOutputDirOption("out_mod_analysis_dir");

      if (in_list.empty())
      {
        OPENMS_LOG_ERROR << "No input files provided (-in)." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // At least one output must be specified
      if (out_id_list.empty() && out_pin_list.empty() && out_qpx_dir.empty() && out_parquet_dir.empty() && out_merged.empty())
      {
        OPENMS_LOG_ERROR << "No output specified. Provide at least one of -out, -out_pin, -out_qpx, -out_parquet, or -out_merged." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // -out count must match -in count
      if (!out_id_list.empty() && in_list.size() != out_id_list.size())
      {
        OPENMS_LOG_ERROR << "Number of output files (-out, " << out_id_list.size()
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
      const bool rescore = getFlag_("rescore");
      const double user_protein_fdr = static_cast<double>(search_params.getValue("FDR:protein"));
      const double user_psm_fdr = static_cast<double>(search_params.getValue("FDR:PSM"));

      // When Percolator rescoring is enabled, defer BOTH PSM and protein FDR to
      // after rescoring. Applying FDR:PSM on raw HyperScores before Percolator
      // would (a) filter on the wrong score (raw HyperScore q-values, not
      // Percolator q-values) and (b) strip decoys, leaving Percolator with
      // nothing for target/decoy competition ("No decoys found").
      if (rescore && (user_protein_fdr > 0.0 || user_psm_fdr > 0.0))
      {
        OPENMS_LOG_INFO << "[ProSE] Percolator rescoring enabled: deferring PSM/protein FDR to post-rescoring." << endl;
        search_params.setValue("FDR:PSM", 0.0);
        search_params.setValue("FDR:protein", 0.0);
      }

      ProSEAlgorithm sse;
      sse.setParameters(search_params);

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
      // Timer spans search + Percolator + FDR + output writing (stopped at report time).
      StopWatch sw_total; sw_total.start();
      // Keep the algorithm's pooled aggregate only where the -out_merged block below consumes it verbatim; with Percolator that block re-merges the rescored per-file results, so the pre-rescoring pool would be waste.
      const bool build_pooled_aggregate = !out_merged.empty() && in_list.size() > 1 && !rescore;
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

      // Optional per-file Percolator rescoring (replaces HyperScore with
      // Percolator q-values for downstream FDR / protein inference).
      std::vector<bool> percolator_succeeded(in_list.size(), false);
      // Tracks whether PSM-level FDR filtering was actually applied per file in
      // this (Percolator) branch — false for files skipped for lack of decoys.
      std::vector<bool> psm_fdr_applied(in_list.size(), false);
      if (rescore)
      {
        Size n_rescorable = 0;
        for (Size i = 0; i < in_list.size(); ++i)
        {
          auto& result = mfres.per_file[i];
          if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
          ++n_rescorable;
          percolator_succeeded[i] = rescoreWithPercolator_(in_list[i], result.protein_ids, result.peptide_ids);
        }
        if (n_rescorable > 0 && std::none_of(percolator_succeeded.begin(), percolator_succeeded.end(), [](bool b) { return b; }))
        {
          OPENMS_LOG_WARN << "Percolator rescoring was requested (-rescore) but no input file could be rescored; "
                          << "all outputs carry the search engine scores." << endl;
        }

        // Apply deferred PSM-level FDR filtering. For files rescored by Percolator,
        // scores are already q-values (see rescoreWithPercolator_). For files that fell back to HyperScores (Percolator
        // skipped/failed), compute q-values via FalseDiscoveryRate first.
        // PSM-level FDR filters target+decoy PSMs alike by q-value; no decoy-specific stripping
        // here (categorical decoy removal is a protein-FDR finalization step, see below).
        if (user_psm_fdr > 0.0)
        {
          for (Size i = 0; i < in_list.size(); ++i)
          {
            auto& result = mfres.per_file[i];
            if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
            if (result.protein_ids.empty() || result.peptide_ids.empty()) continue;

            if (!percolator_succeeded[i])
            {
              bool file_has_decoys = false;
              for (const auto& ph : result.protein_ids[0].getHits())
              {
                if (ph.metaValueExists("target_decoy") && ph.getMetaValue("target_decoy").toString() == "decoy")
                {
                  file_has_decoys = true;
                  break;
                }
              }
              if (!file_has_decoys)
              {
                OPENMS_LOG_WARN << "FDR:PSM requested but no decoys available for " << in_list[i]
                                << " — skipping PSM FDR filtering for this file." << endl;
                continue;
              }
              sse.annotatePsmQValues(result.peptide_ids); // FDR:PSM_groups as in the algorithm's own FDR:PSM
            }

            IDFilter::filterHitsByScore(result.peptide_ids, user_psm_fdr);
            psm_fdr_applied[i] = true;
          }
        }
      }

      // Protein-level FDR. FDR does not compose across runs, so picked-protein FDR is
      // statistically valid only on a COMPLETE protein set: a single input file (handled
      // here) or the pooled -out_merged set (handled below). Filtering each file of a
      // multi-file run separately would not control the FDR of the combined protein list
      // (false positives accumulate across files), so that is deliberately never done.
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
              ProSEAlgorithm::applyCompleteSetProteinFDR(result.protein_ids, result.peptide_ids,
                                                         mfres.decoy_string, mfres.decoy_is_prefix, user_protein_fdr);
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
                          << "outputs (-out/-out_qpx/-out_parquet) retain target+decoy "
                          << "evidence as intermediates. Use -out_merged to obtain protein FDR "
                          << "over the aggregated set, or pool the per-file outputs and apply "
                          << "FDR downstream (e.g. IDMerger -> ProteinInference -> FalseDiscoveryRate)." << endl;
        }
      }

      // Percolator rescoring, deferred PSM FDR and protein-FDR finalization mutate the
      // per-file identifications after the search produced RunStatistics. Recompute the
      // result-level stats so the report reflects the final identifications the user
      // receives; the HyperScore distribution and per-phase timing captured during the
      // search are preserved by updateFinalStats().
      for (Size i = 0; i < in_list.size(); ++i)
      {
        auto& result = mfres.per_file[i];
        if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) { continue; }
        ProSEAlgorithm::updateFinalStats(result.stats, result.peptide_ids, mfres.shared.enzyme,
                                         /*fdr_applied=*/ result.stats.fdr_applied || psm_fdr_applied[i]);
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
      std::vector<std::string> written_idxml, written_idparquet, written_pin;
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
        bool input_failed = false;

        // -- identification output (idXML or idparquet) --
        if (!out_id_list.empty())
        {
          try
          {
            const FileTypes::Type id_type = identificationOutputType_(out_id_list[i]);
            FileHandler().storeIdentifications(out_id_list[i], result.protein_ids, result.peptide_ids, {id_type});
            if (id_type == FileTypes::IDPARQUET)
            {
              written_idparquet.push_back(out_id_list[i]);
            }
            else
            {
              written_idxml.push_back(out_id_list[i]);
            }
          }
          catch (const Exception::BaseException& e)
          {
            OPENMS_LOG_ERROR << "Failed to write identification output for " << in_file
                             << " -> " << out_id_list[i] << ": " << e.what() << endl;
            input_failed = true;
          }
        }

        // -- Percolator .pin output --
        if (!out_pin_list.empty())
        {
          try
          {
            if (result.protein_ids.empty())
            {
              OPENMS_LOG_ERROR << "Cannot write .pin for " << in_file
                               << ": no ProteinIdentification/search parameters available." << endl;
              input_failed = true;
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
              PercolatorInfile::store(out_pin_list[i], result.peptide_ids,
                                     feature_set, enz_str, min_charge, max_charge);
              written_pin.push_back(out_pin_list[i]);
            }
          }
          catch (const Exception::BaseException& e)
          {
            OPENMS_LOG_ERROR << "Failed to write .pin output for " << in_file
                             << " -> " << out_pin_list[i] << ": " << e.what() << endl;
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
            if (!QPXFile::exportToParquet(qpx_psm_table, qpx_psm_file, ParquetWriteConfig{}, psm_scan_format))
            {
              OPENMS_LOG_ERROR << "Failed to write QPX PSM parquet for " << in_file << " -> " << qpx_psm_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "QPX PSM table build returned null for " << in_file << " — skipping " << qpx_psm_file << endl;
          }

          // Protein groups — independent of PSM result.
          const std::string qpx_pg_file = out_qpx_dir + "/" + basename + ".pg.parquet";
          auto qpx_pg_table = ProteinGroupArrowExport::exportToArrow(result.protein_ids, result.peptide_ids);
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

        if (input_failed) { ++failed_count; }
      }

      // -- Write merged Parquet files --
      // Always write merged files (even for single input — provides a canonical name).
      bool merged_ok = true;
      if (!out_qpx_dir.empty())
      {
        // Concatenation inherits the first input's schema metadata, i.e. that per-run file's
        // uuid/creation_date. Stamp a fresh identity so the merged file is its own QPX file.
        // scan_format only describes the merged file if every contributing run agreed on a
        // convention; otherwise it is omitted rather than taken from an arbitrary run.
        std::map<std::string, std::string> psm_extra;
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
        if (pg_meta && !qpx_pg_tables.empty())
        {
          merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(
            qpx_pg_tables, out_qpx_dir + "/quantms.pg.parquet", ParquetWriteConfig{}, pg_meta) && merged_ok;
        }
      }

      if (!out_parquet_dir.empty())
      {
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_psm_tables,  out_parquet_dir + "/openms.psm.parquet")           && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_prot_tables, out_parquet_dir + "/openms.proteins.parquet")      && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_pg_tables,   out_parquet_dir + "/openms.pg.parquet")            && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_sp_tables,   out_parquet_dir + "/openms.search_params.parquet") && merged_ok;
      }

      // Write merged idXML (pool all per-file PSMs, possibly Percolator-rescored,
      // and run cross-file protein inference + optional protein FDR). Wrapped in
      // try/catch so any failure here does NOT discard per-file outputs that
      // already completed above. Per-file outputs (-out / -out_qpx /
      // -out_parquet) are not modified by this block — they retain run-level
      // information.
      bool merged_id_written = false;
      bool merged_id_failed = false;
      if (!out_merged.empty() && in_list.size() > 1)
      {
        try
        {
          // Use the effective decoy marker/position the search resolved (which
          // may be an auto-detected external marker, prefix or suffix), not the
          // raw decoy_prefix parameter.
          const std::string decoy_string = mfres.decoy_string;
          const bool decoy_is_prefix = mfres.decoy_is_prefix;

          vector<ProteinIdentification> merged_protein_ids;
          PeptideIdentificationList merged_peptides;

          // Without Percolator nothing has mutated the identifications in mfres.per_file since
          // the search, and the algorithm already produced exactly this merge: the multi-file
          // searchWithModificationAnalysis pools with the identical sequence of
          // insertRuns(pf.protein_ids, pf.peptide_ids) calls over the identical per-file results
          // and then setPrimaryMSRunPath(in_spectra_files) — so re-merging here would only build
          // a third copy of every PSM. Consume the aggregate by move instead. (An empty aggregate
          // means every file failed, or build_pooled_aggregate was false; fall through to the
          // merger, which reproduces the old result in both cases.)
          // The only run-level value IDMergerAlgorithm derives from its inputs is 'id_merge_index'
          // = the first-occurrence rank of a run's primary MS run path (the merged run's own path
          // is overwritten unconditionally right after this block). That rank is unchanged by the
          // 'file://<basename>' rewrite -test applies to the per-file runs above, because it maps
          // the -in paths one-to-one -- unless two -in files in different directories share a
          // basename, which only -test could ever collapse and which no ProSE test does.
          if (!rescore && !mfres.aggregate.protein_ids.empty())
          {
            merged_protein_ids.emplace_back(std::move(mfres.aggregate.protein_ids[0]));
            merged_peptides = std::move(mfres.aggregate.peptide_ids);
            mfres.aggregate.protein_ids.clear();
            mfres.aggregate.peptide_ids.clear();
          }
          else
          {
            // Percolator rescored the per-file PSMs after the search, so the algorithm's
            // pre-rescoring aggregate would be WRONG here (and was not built at all, see
            // build_pooled_aggregate above): merge the rescored results via IDMergerAlgorithm
            // (accession dedup + identifier remap). Moving out of per_file is safe: past this
            // block only pf.stats / pf.exit_code / pf.modification_analysis are read (the
            // search report and renderRunSummaryYaml), never the identifications.
            IDMergerAlgorithm merger;
            for (auto& pf : mfres.per_file)
            {
              if (pf.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
              merger.insertRuns(std::move(pf.protein_ids), std::move(pf.peptide_ids));
            }
            ProteinIdentification merged_proteins;
            merger.returnResultsAndClear(merged_proteins, merged_peptides);
            merged_protein_ids.emplace_back(std::move(merged_proteins));
          }

          merged_protein_ids[0].setPrimaryMSRunPath(in_list);
          // IDMergerAlgorithm carries over search engine/params but not a run date; set one so the
          // merged idXML is schema-valid (matches the single-file path in ProSEAlgorithm).
          merged_protein_ids[0].setDateTime(DateTime::now());

          // Protein inference: aggregate best PSM score per peptide per protein
          BasicProteinInferenceAlgorithm bpia;
          bpia.run(merged_peptides, merged_protein_ids);

          // Optional picked-protein FDR. Needs identified decoy proteins: gate on DB-level decoys
          // AND a result-level check (the merged inference could be target-only even with a decoy
          // database if no decoy survived in any file). Track whether FDR was actually applied so
          // the decoy-removal finalization below only runs when it was.
          bool merged_fdr_applied = false;
          if (user_protein_fdr > 0.0 && mfres.have_decoys && !decoy_string.empty())
          {
            const bool merged_has_decoy_proteins = std::any_of(
                merged_protein_ids[0].getHits().begin(), merged_protein_ids[0].getHits().end(),
                [&](const ProteinHit& ph) {
                  return decoy_is_prefix ? StringUtils::hasPrefix(ph.getAccession(), decoy_string)
                                         : StringUtils::hasSuffix(ph.getAccession(), decoy_string);
                });
            if (merged_has_decoy_proteins)
            {
              FalseDiscoveryRate fdr;
              fdr.applyPickedProteinFDR(merged_protein_ids[0], decoy_string, decoy_is_prefix);
              IDFilter::filterHitsByScore(merged_protein_ids, user_protein_fdr);
              merged_fdr_applied = true;

              OPENMS_LOG_INFO << "[ProSE] Merged protein inference + FDR: "
                              << merged_protein_ids[0].getHits().size() << " proteins at "
                              << user_protein_fdr * 100 << "% FDR" << endl;
            }
            else
            {
              OPENMS_LOG_WARN << "FDR:protein is set but the merged result has no decoy proteins "
                              << "(marker '" << decoy_string << "'); skipping merged protein FDR to "
                              << "avoid target-only q-values." << endl;
            }
          }

          // Decoy removal is a protein-FDR finalization step: strip decoys only when protein FDR was
          // actually applied. Otherwise the merged output is a pooled intermediate (cross-file
          // inference) that retains target+decoy evidence for downstream FDR.
          if (merged_fdr_applied)
          {
            IDFilter::removeDecoyHits(merged_peptides);
            IDFilter::removeEmptyIdentifications(merged_peptides);
            IDFilter::removeUnreferencedProteins(merged_protein_ids, merged_peptides);
            // Keep indistinguishable-protein and protein groups consistent with the filtered hit set.
            IDFilter::updateProteinGroups(merged_protein_ids[0].getIndistinguishableProteins(), merged_protein_ids[0].getHits());
            IDFilter::updateProteinGroups(merged_protein_ids[0].getProteinGroups(), merged_protein_ids[0].getHits());
            // Strip PeptideEvidence entries pointing to proteins removed by decoy/FDR
            // cleanup (target+decoy PSMs kept by removeDecoyHits would otherwise leave
            // dangling refs that IdXMLFile::store rejects).
            IDFilter::removeDanglingProteinReferences(merged_peptides, merged_protein_ids);
          }

          if (getFlag_("test") && !merged_protein_ids.empty())
          {
            StringList basenames;
            for (const auto& f : in_list) basenames.push_back("file://" + File::basename(f));
            merged_protein_ids[0].setPrimaryMSRunPath(basenames);
          }

          FileHandler().storeIdentifications(out_merged, merged_protein_ids, merged_peptides,
                                             {identificationOutputType_(out_merged)});
          merged_id_written = true;
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "Failed to write merged identification output -> " << out_merged
                           << ": " << e.what()
                           << ". Per-file outputs were written (check above for any per-file errors)." << endl;
          // Do NOT propagate; per-file outputs are already on disk. ProSE still fails at the end.
          merged_id_failed = true;
        }
      }
      else if (!out_merged.empty() && in_list.size() == 1)
      {
        OPENMS_LOG_WARN << "-out_merged is only useful with multiple -in files. Skipping merged output." << endl;
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
        add_manifest("idparquet", written_idparquet);
        add_manifest("pin", written_pin);
        if (merged_id_written)
        {
          add_manifest(identificationOutputType_(out_merged) == FileTypes::IDPARQUET ? "merged idparquet" : "merged idXML", {out_merged});
        }
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
            for (bool b : percolator_succeeded) { if (b) { ++n_ok; } }
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
      if (merged_id_failed)
      {
        OPENMS_LOG_ERROR << "ProSE failed to write the merged identification output " << out_merged << "." << endl;
        return INTERNAL_ERROR;
      }
      if (!merged_ok)
      {
        OPENMS_LOG_ERROR << "ProSE failed to write one or more merged Parquet files." << endl;
        return INTERNAL_ERROR;
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
