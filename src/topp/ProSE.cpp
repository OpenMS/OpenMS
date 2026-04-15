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
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/FORMAT/ProteinGroupArrowExport.h>
#include <OpenMS/FORMAT/ProteinIdentificationArrowIO.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <arrow/table.h>  // only for std::shared_ptr<arrow::Table> declarations

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
    @note Decoy handling: either enable '-Search:decoys' to generate decoys internally, or provide a FASTA database that already contains decoy proteins (e.g., from DecoyDatabase). In both cases, the decoy accession prefix must match '-Search:decoy_prefix' (default: "DECOY_").

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
        "Annotates bottom-up MS/MS spectra using ProSE.",
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

      registerOutputFileList_("out_idxml", "<files>", StringList(), "Output idXML identification file(s). Must have the same number of entries as -in.", false);
      setValidFormats_("out_idxml", ListUtils::create<String>("idXML"));

      registerOutputDir_("out_qpx", "<dir>", "", "Output directory for QPX exchange format Parquet files. Writes per-input <basename>.psm.parquet and <basename>.pg.parquet, plus merged quantms.psm.parquet and quantms.pg.parquet.", false, true);

      registerOutputDir_("out_parquet", "<dir>", "", "Output directory for OpenMS internal format Parquet files. Writes per-input <basename>.psm/proteins/pg/search_params.parquet, plus merged openms.* files.", false, true);

      registerOutputFile_("out_merged", "<file>", "", "Optional merged output file containing all PSMs pooled across input files with cross-file protein inference (BasicProteinInferenceAlgorithm) and optional picked-protein FDR (FDR:protein). Per-file outputs (-out_idxml / -out_qpx / -out_parquet) retain run-level information. Only useful with multiple -in files.", false);
      setValidFormats_("out_merged", ListUtils::create<String>("idXML"));

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
      search_algo_params_with_subsection.insert("Search:", ProSEAlgorithm().getDefaults());
      registerFullParam_(search_algo_params_with_subsection);
    }

    ExitCodes main_(int, const char**) override
    {
      const StringList in_list = getStringList_("in");
      const String database = getStringOption_("database");
      const StringList out_idxml_list = getStringList_("out_idxml");
      const String out_merged = getStringOption_("out_merged");
      const String out_qpx_dir = getOutputDirOption("out_qpx");
      const String out_parquet_dir = getOutputDirOption("out_parquet");
      // getOutputDirOption auto-creates the directory if it doesn't exist (and returns "" if unset).
      const String out_mod_analysis_dir = getOutputDirOption("out_mod_analysis_dir");

      if (in_list.empty())
      {
        OPENMS_LOG_ERROR << "No input files provided (-in)." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // At least one output must be specified
      if (out_idxml_list.empty() && out_qpx_dir.empty() && out_parquet_dir.empty())
      {
        OPENMS_LOG_ERROR << "No output specified. Provide at least one of -out_idxml, -out_qpx, or -out_parquet." << endl;
        return ILLEGAL_PARAMETERS;
      }

      // -out_idxml count must match -in count
      if (!out_idxml_list.empty() && in_list.size() != out_idxml_list.size())
      {
        OPENMS_LOG_ERROR << "Number of output files (-out_idxml, " << out_idxml_list.size()
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
      const String percolator_executable = getStringOption_("percolator_executable");
      const double user_protein_fdr = static_cast<double>(search_params.getValue("FDR:protein"));

      // When Percolator rescoring is enabled, defer protein FDR to after rescoring
      // so protein inference uses Percolator q-values, not raw HyperScores.
      if (!percolator_executable.empty() && user_protein_fdr > 0.0)
      {
        OPENMS_LOG_INFO << "[ProSE] Percolator rescoring enabled: deferring protein FDR to post-rescoring." << endl;
        search_params.setValue("FDR:protein", 0.0);
      }

      ProSEAlgorithm sse;
      sse.setParameters(search_params);

      // Determine open-search mode the same way the algorithm does, so we can give
      // the user actionable feedback if -out_mod_analysis_dir is set in closed-search.
      const double precursor_tol_lower = static_cast<double>(search_params.getValue("precursor:mass_tolerance_lower"));
      const double precursor_tol_upper = static_cast<double>(search_params.getValue("precursor:mass_tolerance_upper"));
      const String precursor_tol_unit = search_params.getValue("precursor:mass_tolerance_unit").toString();
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
      ProSEAlgorithm::MultiFileSearchResult mfres =
        sse.searchWithModificationAnalysis(in_list, database, mod_analysis_base_names, aggregate_base_name);

      if (mfres.per_file.size() != in_list.size())
      {
        OPENMS_LOG_ERROR << "Internal error: per-file result count (" << mfres.per_file.size()
                         << ") does not match input file count (" << in_list.size() << ")." << endl;
        return INTERNAL_ERROR;
      }

      // Optional per-file Percolator rescoring (replaces HyperScore with
      // Percolator q-values for downstream FDR / protein inference).
      if (!percolator_executable.empty())
      {
        // Percolator requires decoys for target/decoy competition.
        // Check per-file results for target_decoy annotations.
        bool has_decoys_for_percolator = false;
        for (const auto& pf : mfres.per_file)
        {
          if (pf.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
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
          if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
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

          OPENMS_LOG_INFO << "[ProSE] Rescoring " << in_list[i] << " with Percolator..." << endl;
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

      // Write merged output: pool all per-file PSMs (possibly Percolator-rescored),
      // run cross-file protein inference + optional protein FDR, write to -out_merged.
      // Per-file outputs (-out_idxml/-out_qpx/-out_parquet) are NOT modified — they retain run-level information.
      if (!out_merged.empty() && in_list.size() > 1)
      {
        const String decoy_prefix = search_params.getValue("decoy_prefix").toString();

        // Merge per-file results via IDMergerAlgorithm (accession dedup + identifier remap)
        IDMergerAlgorithm merger;
        for (const auto& pf : mfres.per_file)
        {
          if (pf.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK) continue;
          merger.insertRuns(pf.protein_ids, pf.peptide_ids);
        }
        ProteinIdentification merged_proteins;
        PeptideIdentificationList merged_peptides;
        merger.returnResultsAndClear(merged_proteins, merged_peptides);

        vector<ProteinIdentification> merged_protein_ids = {std::move(merged_proteins)};
        merged_protein_ids[0].setPrimaryMSRunPath(in_list);

        // Protein inference: aggregate best PSM score per peptide per protein
        BasicProteinInferenceAlgorithm bpia;
        bpia.run(merged_peptides, merged_protein_ids);

        // Optional picked-protein FDR
        if (user_protein_fdr > 0.0)
        {
          FalseDiscoveryRate fdr;
          fdr.applyPickedProteinFDR(merged_protein_ids[0], decoy_prefix, true);
          IDFilter::filterHitsByScore(merged_protein_ids, user_protein_fdr);

          OPENMS_LOG_INFO << "[ProSE] Merged protein inference + FDR: "
                          << merged_protein_ids[0].getHits().size() << " proteins at "
                          << user_protein_fdr * 100 << "% FDR" << endl;
        }

        // Clean up decoys in merged output
        IDFilter::removeDecoyHits(merged_peptides);
        IDFilter::removeEmptyIdentifications(merged_peptides);
        IDFilter::removeUnreferencedProteins(merged_protein_ids, merged_peptides);

        if (getFlag_("test") && !merged_protein_ids.empty())
        {
          StringList basenames;
          for (const auto& f : in_list) basenames.push_back("file://" + File::basename(f));
          merged_protein_ids[0].setPrimaryMSRunPath(basenames);
        }

        FileHandler().storeIdentifications(out_merged, merged_protein_ids, merged_peptides, {FileTypes::IDXML});
      }
      else if (!out_merged.empty() && in_list.size() == 1)
      {
        OPENMS_LOG_WARN << "-out_merged is only useful with multiple -in files. Skipping merged output." << endl;
      }

      // Write per-file outputs and track failures.
      // Accumulate Arrow tables for merged parquet output.
      std::vector<std::shared_ptr<arrow::Table>> qpx_psm_tables, qpx_pg_tables;
      std::vector<std::shared_ptr<arrow::Table>> oms_psm_tables, oms_prot_tables, oms_pg_tables, oms_sp_tables;

      Size failed_count = 0;
      for (Size i = 0; i < in_list.size(); ++i)
      {
        const String& in_file = in_list[i];
        ProSEAlgorithm::SearchResult& result = mfres.per_file[i];

        if (result.exit_code != ProSEAlgorithm::ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_ERROR << "Search failed for " << in_file
                           << " (algorithm exit code " << static_cast<int>(result.exit_code) << "). Skipping output." << endl;
          ++failed_count;
          continue;
        }

        const String basename = File::stemName(in_file);

        // In test mode, replace the absolute MS run path with file://<basename> for reproducible diffs.
        if (getFlag_("test") && !result.protein_ids.empty())
        {
          result.protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_file)});
        }

        // Output modes are independent: a failure in one mode must not skip the
        // sibling modes for the same input. Each mode contributes at most one
        // failure to input_failed; the input counts as failed if any mode failed.
        bool input_failed = false;

        // -- idXML output --
        if (!out_idxml_list.empty())
        {
          try
          {
            FileHandler().storeIdentifications(out_idxml_list[i], result.protein_ids, result.peptide_ids, {FileTypes::IDXML});
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
          // PSM — build table once, write to file, accumulate same table for merge.
          const String qpx_psm_file = out_qpx_dir + "/" + basename + ".psm.parquet";
          auto qpx_psm_table = QPXFile::exportPSMsToQPXArrow(result.protein_ids, result.peptide_ids, /*export_all_psms=*/false);
          if (qpx_psm_table)
          {
            qpx_psm_tables.push_back(qpx_psm_table);
            if (!QPXFile::exportToParquet(qpx_psm_table, qpx_psm_file))
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
          const String qpx_pg_file = out_qpx_dir + "/" + basename + ".pg.parquet";
          auto qpx_pg_table = ProteinGroupArrowExport::exportToArrow(result.protein_ids, result.peptide_ids);
          if (qpx_pg_table)
          {
            qpx_pg_tables.push_back(qpx_pg_table);
            if (!ProteinGroupArrowExport::exportToParquet(qpx_pg_table, qpx_pg_file))
            {
              OPENMS_LOG_ERROR << "Failed to write QPX PG parquet for " << in_file << " -> " << qpx_pg_file << endl;
              input_failed = true;
            }
          }
          else
          {
            OPENMS_LOG_WARN << "QPX PG table build returned null for " << in_file << " — skipping " << qpx_pg_file << endl;
          }
        }

        // -- OpenMS internal directory output --
        if (!out_parquet_dir.empty())
        {
          // PSM (internal PSMSchema) — build once, write, accumulate for merge.
          const String oms_psm_file = out_parquet_dir + "/" + basename + ".psm.parquet";
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
          const String oms_prot_file = out_parquet_dir + "/" + basename + ".proteins.parquet";
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
          const String oms_pg_file = out_parquet_dir + "/" + basename + ".pg.parquet";
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
          const String oms_sp_file = out_parquet_dir + "/" + basename + ".search_params.parquet";
          auto sp_table = ProteinIdentificationArrowIO::exportSearchParamsToArrow(result.protein_ids);
          if (sp_table)
          {
            oms_sp_tables.push_back(sp_table);
            if (!ProteinIdentificationArrowIO::exportSearchParamsToParquet(result.protein_ids, oms_sp_file))
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
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(qpx_psm_tables, out_qpx_dir + "/quantms.psm.parquet") && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(qpx_pg_tables,  out_qpx_dir + "/quantms.pg.parquet")  && merged_ok;
      }

      if (!out_parquet_dir.empty())
      {
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_psm_tables,  out_parquet_dir + "/openms.psm.parquet")           && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_prot_tables, out_parquet_dir + "/openms.proteins.parquet")      && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_pg_tables,   out_parquet_dir + "/openms.pg.parquet")            && merged_ok;
        merged_ok = ArrowIOHelpers::concatenateAndWriteToParquet(oms_sp_tables,   out_parquet_dir + "/openms.search_params.parquet") && merged_ok;
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
      return EXECUTION_OK;
    }
};

int main(int argc, const char** argv)
{
  ProSE tool;
  return tool.main(argc, argv);
}

///@endcond
