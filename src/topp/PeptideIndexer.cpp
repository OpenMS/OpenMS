// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Knut Reinert $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;

/**
@page TOPP_PeptideIndexer PeptideIndexer

@brief Refreshes the protein references for all peptide hits from an idXML file and adds target/decoy information.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; PeptideIndexer &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
        </tr>
    </table>
</CENTER>

PeptideIndexer refreshes target/decoy information and mapping of peptides to proteins.
The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool. (For FDR calculations, peptides hitting both target and decoy proteins are counted as target hits.)

PeptideIndexer allows for ambiguous amino acids (B|J|Z|X) in the protein database and peptide sequence. 

Enzyme cutting rules and partial specificity are derived from input idXML automatically by default or can be specified explicitly by the user.

All peptide and protein hits are annotated with target/decoy information, using the meta value 'target_decoy'. 
For proteins the possible values are "target" and "decoy", depending on whether the protein accession contains the decoy pattern (parameter @p decoy_string) 
as a suffix or prefix, respectively (see parameter @p prefix). 
Resulting protein hits appear in the order of the FASTA file, except for orphaned proteins, which will appear first with an empty 'target_decoy' metavalue.
Duplicate protein accessions & sequences will not raise a warning, but create multiple hits (PeptideIndexer reads the FASTA file piecewise for efficiency
reasons, and thus might not see all accessions & sequences at once).

Peptide hits are annotated with metavalue 'protein_references', and if matched to at least one protein also with metavalue 'target_decoy'.
The possible values for 'target_decoy' in peptides are "target", "decoy" and "target+decoy", 
depending on whether the peptide sequence is found only in target proteins, only in decoy proteins, or in both. If the peptide is unmatched the metavalue is missing.

Runtime: PeptideIndexer is usually very fast (loading and storing the data takes the most time) and search speed can be further improved (linearly) by using more threads. 
Avoid allowing too many (>=4) ambiguous amino acids if your database contains long stretches of 'X' (exponential search space).

PeptideIndexer supports relative database filenames, which (when not found in the current working directory) are looked up in the directories specified
by @p OpenMS.ini:id_db_dir. The database is by default derived from the input idXML's metainformation ('auto' setting), but can be specified explicitly.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_PeptideIndexer.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_PeptideIndexer.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeptideIndexer :
  public TOPPBase
{
public:
  TOPPPeptideIndexer() :
    TOPPBase("PeptideIndexer",
             "Refreshes the protein references for all peptide hits.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input idXML file containing the identifications.");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerInputFile_("fasta", "<file>", "", "Input sequence database in FASTA format. "
                                              "Leave empty for using the same DB as used for the input idXML (this might fail). "
                                              "Non-existing relative filenames are looked up via 'OpenMS.ini:id_db_dir'", false, false, { "skipexists" });
    setValidFormats_("fasta", { "fasta" }, false);
    registerOutputFile_("out", "<file>", "", "Output idXML file.");
    setValidFormats_("out", {"idXML"});

    registerFullParam_(PeptideIndexing().getParameters());
   }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String db_name = getStringOption_("fasta"); // optional. Might be empty.

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // we stream the Fasta file
    std::vector<ProteinIdentification> prot_ids;
    std::vector<PeptideIdentification> pep_ids;

    FileHandler().loadIdentifications(in, prot_ids, pep_ids, {FileTypes::IDXML});

    if (db_name.empty())
    { // determine from metadata in idXML
      OPENMS_LOG_INFO << "Automatically deriving DB from meta data ...";
      for (const auto& pi : prot_ids)
      {
        if (!db_name.empty() && db_name != pi.getSearchParameters().db)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 
              "Different conflicting database names in idXML files are from multiple runs. Cannot automatically determine DB from these candidates:\n"
              "1) " + db_name + "\n"
              "2) " + pi.getSearchParameters().db);
        }
        db_name = pi.getSearchParameters().db;
      }
      OPENMS_LOG_INFO << "DB: " << db_name << std::endl;
    }
    
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    
    PeptideIndexing indexer;
    Param param = getParam_();
    Param param_pi = indexer.getParameters();
    param_pi.update(param, false, false, false, false, OpenMS_Log_debug); // suppress param. update message
    indexer.setParameters(param_pi);
    indexer.setLogType(this->log_type_);
    FASTAContainer<TFI_File> proteins(db_name);
    PeptideIndexing::ExitCodes indexer_exit = indexer.run(proteins, prot_ids, pep_ids);

    //-------------------------------------------------------------
    // calculate protein coverage
    //-------------------------------------------------------------

    if (param.getValue("write_protein_sequence").toBool())
    {
      for (Size i = 0; i < prot_ids.size(); ++i)
      {
        prot_ids[i].computeCoverage(pep_ids);
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    FileHandler().storeIdentifications(out, prot_ids, pep_ids, {FileTypes::IDXML});

    if (indexer_exit == PeptideIndexing::DATABASE_EMPTY)
    {
      return INPUT_FILE_EMPTY;
    }
    else if (indexer_exit == PeptideIndexing::UNEXPECTED_RESULT)
    {
      return UNEXPECTED_RESULT;
    }
    else if ((indexer_exit != PeptideIndexing::EXECUTION_OK) &&
             (indexer_exit != PeptideIndexing::PEPTIDE_IDS_EMPTY))
    {
      return UNKNOWN_ERROR;
    }
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPeptideIndexer tool;
  return tool.main(argc, argv);
}

/// @endcond
