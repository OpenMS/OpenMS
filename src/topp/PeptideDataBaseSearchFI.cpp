// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// $Authors:  $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PeptideDataBaseSearchFI PeptideDataBaseSearchFI

@brief Identifies peptides in MS/MS spectra.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; PeptideDataBaseSearchFI &rarr;</td>
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
    @note Decoy handling: either enable '-Search:decoys' to generate decoys internally, or provide a FASTA database that already contains decoy proteins (e.g., from DecoyDatabase). In both cases, the decoy accession prefix must match '-Search:decoy_prefix' (default: "DECOY_").

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_PeptideDataBaseSearchFI.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_PeptideDataBaseSearchFI.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class PeptideDataBaseSearchFI :
    public TOPPBase
{
  public:
    PeptideDataBaseSearchFI() :
      TOPPBase("PeptideDataBaseSearchFI",
        "Annotates bottom-up MS/MS spectra using PeptideDataBaseSearchFI.",
        false)
    {
    }

  protected:
    void registerOptionsAndFlags_() override
    {
      registerInputFileList_("in", "<files>", StringList(), "Input spectrum file(s). Multiple files are searched against the same database; the fragment index is built once and reused.");
      setValidFormats_("in", ListUtils::create<String>("mzML,d"));

      registerInputFile_("database", "<file>", "", "Input protein sequence database in FASTA format.");
      setValidFormats_("database", ListUtils::create<String>("fasta"));

      registerOutputFileList_("out", "<files>", StringList(), "Output identification file(s). Must have the same number of entries as -in. Output format is auto-detected per file from the extension (.idXML or .parquet); a mix of formats within one run is allowed.");
      setValidFormats_("out", ListUtils::create<String>("idXML,parquet"));

      registerOutputDir_("out_mod_analysis_dir", "<dir>", "", "Optional directory to write modification-analysis tables (delta-mass, PTM stats) when running in open-search mode. When set, per-file tables are written using each input file's basename and an additional aggregate table is written across all input files. Has no effect in closed-search mode.", false, true);

      registerInputFile_("percolator_executable", "<executable>",
#ifdef OPENMS_WINDOWSPLATFORM
        "",
#else
        "",
#endif
        "Path to the Percolator executable. If set, per-file PSMs are rescored with Percolator (via PercolatorAdapter) before output. Leave empty to skip rescoring.", false, true, ListUtils::create<String>("skipexists"));

      // put search algorithm parameters at Search: subtree of parameters
      Param search_algo_params_with_subsection;
      search_algo_params_with_subsection.insert("Search:", PeptideSearchEngineFIAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    ExitCodes main_(int, const char**) override
    {
      const StringList in_list = getStringList_("in");
      const String database = getStringOption_("database");
      const StringList out_list = getStringList_("out");
      // getOutputDirOption auto-creates the directory if it doesn't exist (and returns "" if unset).
      const String out_mod_analysis_dir = getOutputDirOption("out_mod_analysis_dir");

      if (in_list.empty())
      {
        OPENMS_LOG_ERROR << "No input files provided (-in)." << endl;
        return ILLEGAL_PARAMETERS;
      }

      if (in_list.size() != out_list.size())
      {
        OPENMS_LOG_ERROR << "Number of output files (-out, " << out_list.size()
                         << ") must match number of input files (-in, " << in_list.size() << ")." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // Pre-flight validate every -out file extension before running the search,
      // so the user learns about a typo (e.g. ".idxml" -> "idxml") before paying the
      // index-build + search cost. registerOutputFileList_ does NOT validate per-entry
      // formats at the CLI layer, so we have to do it here.
      for (const String& out_file : out_list)
      {
        const FileTypes::Type out_type = FileHandler::getTypeByFileName(out_file);
        if (out_type != FileTypes::IDXML && out_type != FileTypes::PARQUET)
        {
          OPENMS_LOG_ERROR << "Unsupported output format for '" << out_file
                           << "' (expected .idXML or .parquet)." << endl;
          return ILLEGAL_PARAMETERS;
        }
      }

      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      const Param search_params = getParam_().copy("Search:", true);
      PeptideSearchEngineFIAlgorithm sse;
      sse.setParameters(search_params);

      // Determine open-search mode the same way the algorithm does, so we can give
      // the user actionable feedback if -out_mod_analysis_dir is set in closed-search.
      const double precursor_tol = static_cast<double>(search_params.getValue("precursor:mass_tolerance"));
      const String precursor_tol_unit = search_params.getValue("precursor:mass_tolerance_unit").toString();
      const bool open_search_mode = (precursor_tol_unit == "ppm")
        ? (precursor_tol > 1000.0)
        : (precursor_tol > 1.0);

      if (!out_mod_analysis_dir.empty() && !open_search_mode)
      {
        OPENMS_LOG_WARN << "-out_mod_analysis_dir was set but the search is in CLOSED-search mode "
                        << "(precursor tolerance " << precursor_tol << " " << precursor_tol_unit
                        << "). Modification-analysis tables will NOT be written." << endl;
      }

      // Build per-file modification-analysis base names if -out_mod_analysis_dir is set.
      // Each per-file base name is "<dir>/<input_basename_without_ext>" so the algorithm
      // appends suffixes like _ModificationAnalysis_DeltaMassStats.tsv.
      // Use a long, unlikely-to-collide aggregate stem so an input named "aggregate.mzML" (or "aggregate.d")
      // doesn't silently overwrite the aggregate output (or vice versa).
      static const String AGGREGATE_STEM = "_aggregate_across_files";
      std::vector<String> mod_analysis_base_names;
      String aggregate_base_name;
      if (!out_mod_analysis_dir.empty() && open_search_mode)
      {
        // Detect collisions: two inputs with the same File::stemName, OR an input whose
        // stem matches the aggregate stem. We refuse to run rather than silently overwrite.
        std::map<String, Size> stem_to_first_index;
        for (Size i = 0; i < in_list.size(); ++i)
        {
          const String stem = File::stemName(in_list[i]);
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
        for (const String& in_file : in_list)
        {
          mod_analysis_base_names.push_back(out_mod_analysis_dir + "/" + File::stemName(in_file));
        }
        aggregate_base_name = out_mod_analysis_dir + "/" + AGGREGATE_STEM;
      }

      // Single-shot multi-file search: builds the fragment index once and iterates over -in.
      PeptideSearchEngineFIAlgorithm::MultiFileSearchResult mfres =
        sse.searchWithModificationAnalysis(in_list, database, mod_analysis_base_names, aggregate_base_name);

      if (mfres.per_file.size() != in_list.size())
      {
        OPENMS_LOG_ERROR << "Internal error: per-file result count (" << mfres.per_file.size()
                         << ") does not match input file count (" << in_list.size() << ")." << endl;
        return INTERNAL_ERROR;
      }

      // Optional per-file Percolator rescoring (replaces HyperScore with
      // Percolator q-values for downstream FDR / protein inference).
      const String percolator_executable = getStringOption_("percolator_executable");
      if (!percolator_executable.empty())
      {
        // Percolator requires decoys for target/decoy competition.
        // Check per-file results for target_decoy annotations.
        bool has_decoys_for_percolator = false;
        for (const auto& pf : mfres.per_file)
        {
          if (pf.exit_code != PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK) continue;
          for (const auto& ph : pf.protein_ids[0].getHits())
          {
            if (ph.metaValueExists("target_decoy") && ph.getMetaValue("target_decoy").toString() == "decoy")
            {
              has_decoys_for_percolator = true;
              break;
            }
          }
          if (has_decoys_for_percolator) break;
        }

        if (!has_decoys_for_percolator)
        {
          OPENMS_LOG_WARN << "Percolator rescoring requires decoys but none found in results. "
                          << "Enable '-Search:decoys' or provide a FASTA with decoy proteins. "
                          << "Skipping rescoring." << endl;
        }
        else
        for (Size i = 0; i < in_list.size(); ++i)
        {
          auto& result = mfres.per_file[i];
          if (result.exit_code != PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK) continue;
          if (result.peptide_ids.size() < 100)
          {
            OPENMS_LOG_WARN << "Skipping Percolator rescoring for " << in_list[i]
                            << " (only " << result.peptide_ids.size() << " PSMs, need >= 100)." << endl;
            continue;
          }

          // Write intermediate idXML for PercolatorAdapter input
          String tmp_in = File::getTempDirectory() + "/" + File::stemName(in_list[i]) + "_perc_in.idXML";
          String tmp_out = File::getTempDirectory() + "/" + File::stemName(in_list[i]) + "_perc_out.idXML";
          String tmp_weights = File::getTempDirectory() + "/" + File::stemName(in_list[i]) + "_perc.weights";
          FileHandler().storeIdentifications(tmp_in, result.protein_ids, result.peptide_ids, {FileTypes::IDXML});

          std::vector<String> perc_params = {
            "-in", tmp_in,
            "-out", tmp_out,
            "-percolator_executable", percolator_executable,
            "-train_best_positive",
            "-score_type", "q-value",
            "-post_processing_tdc",
            "-weights", tmp_weights
          };

          OPENMS_LOG_INFO << "[PDBS-FI] Rescoring " << in_list[i] << " with Percolator..." << endl;
          TOPPBase::ExitCodes perc_ec = runExternalProcess_(String("PercolatorAdapter"), perc_params);

          if (perc_ec != EXECUTION_OK)
          {
            OPENMS_LOG_WARN << "Percolator rescoring failed for " << in_list[i]
                            << ". Using original HyperScore results." << endl;
          }
          else
          {
            result.protein_ids.clear();
            result.peptide_ids.clear();
            FileHandler().loadIdentifications(tmp_out, result.protein_ids, result.peptide_ids, {FileTypes::IDXML});
            IDFilter::keepNBestHits(result.peptide_ids, 1);
          }

          // Clean up temp files
          File::remove(tmp_in);
          File::remove(tmp_out);
          File::remove(tmp_weights);
        }
      }

      // Write per-file outputs and track failures.
      Size failed_count = 0;
      for (Size i = 0; i < in_list.size(); ++i)
      {
        const String& in_file = in_list[i];
        const String& out_file = out_list[i];
        PeptideSearchEngineFIAlgorithm::SearchResult& result = mfres.per_file[i];

        if (result.exit_code != PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_ERROR << "Search failed for " << in_file
                           << " (algorithm exit code " << static_cast<int>(result.exit_code) << "). Skipping output." << endl;
          ++failed_count;
          continue;
        }

        // In test mode, replace the absolute MS run path with file://<basename> for reproducible diffs.
        if (getFlag_("test") && !result.protein_ids.empty())
        {
          result.protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_file)});
        }

        // Dispatch on output extension (already validated above to be IDXML or PARQUET).
        const FileTypes::Type out_type = FileHandler::getTypeByFileName(out_file);
        if (out_type == FileTypes::IDXML)
        {
          FileHandler().storeIdentifications(out_file, result.protein_ids, result.peptide_ids, {FileTypes::IDXML});
        }
        else // FileTypes::PARQUET
        {
          if (!QPXFile::exportToParquet(result.protein_ids, result.peptide_ids, out_file, /*export_all_psms=*/false))
          {
            OPENMS_LOG_ERROR << "Failed to write parquet output for " << in_file << " -> " << out_file << endl;
            ++failed_count;
            continue;
          }
        }
      }

      if (failed_count > 0)
      {
        OPENMS_LOG_ERROR << "PeptideDataBaseSearchFI finished with " << failed_count
                         << " file(s) failing out of " << in_list.size() << "." << endl;
        return INTERNAL_ERROR;
      }
      return EXECUTION_OK;
    }
};

int main(int argc, const char** argv)
{
  PeptideDataBaseSearchFI tool;
  return tool.main(argc, argv);
}

///@endcond
