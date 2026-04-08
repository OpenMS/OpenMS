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
#include <OpenMS/FORMAT/QPXFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

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

      registerOutputPrefix_("out_mod_analysis_dir", "<dir>", "", "Optional directory to write modification-analysis tables (delta-mass, PTM stats) when running in open-search mode. When set, per-file tables are written using each input mzML's basename and an additional aggregate table is written across all input files. Has no effect in closed-search mode.", false, true);

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
      const String out_mod_analysis_dir = getStringOption_("out_mod_analysis_dir");

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

      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);

      PeptideSearchEngineFIAlgorithm sse;
      sse.setParameters(getParam_().copy("Search:", true));

      // Build per-file modification-analysis base names if -out_mod_analysis_dir is set.
      // Each per-file base name is "<dir>/<input_basename_without_ext>" so the algorithm
      // appends suffixes like _ModificationAnalysis_DeltaMassStats.tsv.
      std::vector<String> mod_analysis_base_names;
      String aggregate_base_name;
      if (!out_mod_analysis_dir.empty())
      {
        mod_analysis_base_names.reserve(in_list.size());
        for (const String& in_file : in_list)
        {
          String stem = File::stemName(in_file);
          mod_analysis_base_names.push_back(out_mod_analysis_dir + "/" + stem);
        }
        aggregate_base_name = out_mod_analysis_dir + "/aggregate";
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

      // Write per-file outputs and track failures.
      bool any_failure = false;
      for (Size i = 0; i < in_list.size(); ++i)
      {
        const String& in_file = in_list[i];
        const String& out_file = out_list[i];
        PeptideSearchEngineFIAlgorithm::SearchResult& result = mfres.per_file[i];

        if (result.exit_code != PeptideSearchEngineFIAlgorithm::ExitCodes::EXECUTION_OK)
        {
          OPENMS_LOG_ERROR << "Search failed for " << in_file
                           << " (algorithm exit code " << static_cast<int>(result.exit_code) << "). Skipping output." << endl;
          any_failure = true;
          continue;
        }

        // In test mode, replace the absolute MS run path with file://<basename> for reproducible diffs.
        if (getFlag_("test") && !result.protein_ids.empty())
        {
          result.protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_file)});
        }

        // Dispatch on output extension: .idXML -> idXML; .parquet -> QPX PSM parquet.
        const FileTypes::Type out_type = FileHandler::getTypeByFileName(out_file);
        if (out_type == FileTypes::IDXML)
        {
          FileHandler().storeIdentifications(out_file, result.protein_ids, result.peptide_ids, {FileTypes::IDXML});
        }
        else if (out_type == FileTypes::PARQUET)
        {
          if (!QPXFile::exportToParquet(result.protein_ids, result.peptide_ids, out_file, /*export_all_psms=*/false))
          {
            OPENMS_LOG_ERROR << "Failed to write parquet output for " << in_file << " -> " << out_file << endl;
            any_failure = true;
            continue;
          }
        }
        else
        {
          OPENMS_LOG_ERROR << "Unsupported output format for " << out_file
                           << " (expected .idXML or .parquet)." << endl;
          any_failure = true;
          continue;
        }
      }

      return any_failure ? INTERNAL_ERROR : EXECUTION_OK;
    }
};

int main(int argc, const char** argv)
{
  PeptideDataBaseSearchFI tool;
  return tool.main(argc, argv);
}

///@endcond
