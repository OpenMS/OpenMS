// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/XLMS/OpenPepXLLFAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/File.h>

using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenPepXLLF OpenPepXLLF

@brief Search for cross-linked peptide pairs in tandem MS spectra

This tool performs a search for cross-links in the given mass spectra.

It executes the following steps in order:
<ul>
  <li>Reading of MS2 spectra from the given mzML file, MS1 spectra are ignored for now</li>
  <li>Processing of spectra: deisotoping and filtering</li>
  <li>Digesting and preprocessing the protein database, building a peptide pair index dependent on the precursor masses of the MS2 spectra</li>
  <li>Generating theoretical spectra of cross-linked peptides and aligning the experimental spectra against those</li>
  <li>Scoring of cross-link spectrum matches</li>
  <li>Using PeptideIndexer to map the peptides to all possible source proteins</li>
  <li>Writing out the results in idXML,  mzid according to mzIdentML 1.2 specifications and/or in the xQuest output format</li>
</ul>

See below or have a look at the INI file (via "OpenPepXLLF -write_ini myini.ini") for available parameters and more functionality.

<h3>Input: MS2 spectra and fasta database of proteins expected to be cross-linked in the sample</h3>
The spectra should be provided as one mzML file. If you have multiple files, e.g. for multiple fractions, you should run this tool on each
file separately.
The database can either be provided as one merged file containing targets and decoys or as two separate files.

<h3>Parameters</h3>
The parameters for fixed and variable modifications refer to additional modifications beside the cross-linker.
The linker used in the experiment has to be described using the cross-linker specific parameters.
Only one mass is allowed for a cross-linker that links two peptides, while multiple masses are possible for mono-links of the same cross-linking reagent.
Mono-links are cross-linkers, that are linked to one peptide by one of their two reactive groups.
To search for isotopically labeled pairs of cross-linkers see the tool OpenPepXL.
The parameters -cross_linker:residue1 and -cross_linker:residue2 are used to enumerate the amino acids,
that each end of the linker can react with. This way any heterobifunctional cross-linker can be defined.
To define a homobifunctional cross-linker, these two parameters should have the same value.
The parameter -cross_linker:name is used to solve ambiguities caused by different cross-linkers with the same mass
after the linking reaction (see section on output for clarification).

<h3>Output: XL-MS Identifications with scores and linked positions in the proteins</h3>
There are three file formats for output of data possible. idXML is the internal format of OpenMS, and is recommended for post-processing using other TOPP tools like XFDR or TOPPView.
The second format is the output format of xQuest,which is a popular XL-MS ID tool.
This format is compatible with a number of post-processing and visulization tools,
like xProphet for FDR estimation (Leitner, A. et al., 2014, Nature protocols)
and through the xQuest Results Viewer also the XlinkAnalyzer for visualization and analysis using protein structures (Kosinski, J. et al., 2015, Journal of structural biology).
The third format is mzIdentML according to the specifications for XL-MS ID data in version 1.2 (Vizcaíno, J. A. et al., 2017, Mol Cell Proteomics).
This is a standardized format and will be compatible with complete submissions to the PRIDE database, which is part of the ProteomeXchange consortium.
The specification includes the XLMOD database of cross-linking reagents, and if the provided cross-link mass matches one from the
database, its accession and name are used. If the name is provided with the -cross_linker:name parameter, it is used
to solve ambiguities arising from different cross-linkers having the same mass after the linking reaction (e.g. DSS and BS3).
It is also used as the name of the linker, if no matching masses are found in the database.

<CENTER>
  <table>
      <tr>
          <th ALIGN = "center"> pot. predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=2> &rarr; OpenPepXLLF &rarr;</td>
          <th ALIGN = "center"> pot. successor tools </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
      </tr>
  </table>
</CENTER>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenPepXLLF.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenPepXLLF.html
*/

/// @cond TOPPCLASSES


class TOPPOpenPepXLLF :
  public TOPPBase
{
public:
  TOPPOpenPepXLLF() :
    TOPPBase("OpenPepXLLF", "Protein-protein cross linking with label-free linkers.", true)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // input files
    registerInputFile_("in", "<file>", "", "Input file containing the spectra.", true, false);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("database", "<file>", "", "Input file containing the protein database.", true, false);
    setValidFormats_("database", ListUtils::create<String>("fasta"));

    registerInputFile_("decoy_database", "<file>", "", "Input file containing the decoy protein database. Decoys can also be included in the normal database file instead (or additionally).", false, true);
    setValidFormats_("decoy_database", ListUtils::create<String>("fasta"));

    registerFullParam_(OpenPepXLLFAlgorithm().getDefaults());

    // output file
    registerOutputFile_("out_idXML", "<idXML_file>", "", "Results in idXML format (at least one of these output parameters should be set, otherwise you will not have any results).", false, false);
    setValidFormats_("out_idXML", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_mzIdentML", "<mzIdentML_file>","", "Results in mzIdentML (.mzid) format (at least one of these output parameters should be set, otherwise you will not have any results)", false, false);
    setValidFormats_("out_mzIdentML", ListUtils::create<String>("mzid"));

    registerOutputFile_("out_xquestxml", "<xQuestXML_file>", "", "Results in the xquest.xml format (at least one of these output parameters should be set, otherwise you will not have any results).", false, false);
    setValidFormats_("out_xquestxml", ListUtils::create<String>("xquest.xml"));

    registerOutputFile_("out_xquest_specxml", "<xQuestSpecXML_file>", "", "Matched spectra in the xQuest .spec.xml format for spectra visualization in the xQuest results manager.", false, false);
    setValidFormats_("out_xquest_specxml", ListUtils::create<String>("spec.xml"));
  }

  ExitCodes main_(int, const char**) override
  {
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);

    const string in_mzml(getStringOption_("in"));
    const string in_fasta(getStringOption_("database"));
    const string in_decoy_fasta(getStringOption_("decoy_database"));
    const string out_idXML(getStringOption_("out_idXML"));
    const string out_xquest = getStringOption_("out_xquestxml");
    const string out_xquest_specxml = getStringOption_("out_xquest_specxml");
    const string out_mzIdentML = getStringOption_("out_mzIdentML");

    OPENMS_LOG_INFO << "Analyzing file: " << endl;
    OPENMS_LOG_INFO << in_mzml << endl;

    // load MS2 map
    PeakMap unprocessed_spectra;
    FileHandler f;

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    options.addMSLevel(1);
    f.getOptions() = options;
    f.loadExperiment(in_mzml, unprocessed_spectra, {FileTypes::MZML}, log_type_);

    // load linked features
    // @FIXME Orphaned code
    //ConsensusMap cfeatures;
    //ConsensusXMLFile cf;

    // load fasta database
    progresslogger.startProgress(0, 1, "Load database from FASTA file...");
    FASTAFile fastaFile;
    vector<FASTAFile::FASTAEntry> fasta_db;
    fastaFile.load(in_fasta, fasta_db);

    if (!in_decoy_fasta.empty())
    {
      vector<FASTAFile::FASTAEntry> fasta_decoys;
      fastaFile.load(in_decoy_fasta, fasta_decoys);
      fasta_db.reserve(fasta_db.size() + fasta_decoys.size());
      fasta_db.insert(fasta_db.end(), fasta_decoys.begin(), fasta_decoys.end());
    }
    progresslogger.endProgress();


    // initialize solution vectors
    vector<ProteinIdentification> protein_ids(1);
    vector<PeptideIdentification> peptide_ids;

    // these are mainly necessary for writing out xQuest type spectrum files
    vector< vector< OPXLDataStructs::CrossLinkSpectrumMatch > > all_top_csms;
    PeakMap spectra;

    OpenPepXLLFAlgorithm search_algorithm;
    Param this_param = getParam_();
    Param algo_param = search_algorithm.getParameters();
    algo_param.update(this_param, false, false, false, false, OpenMS_Log_debug); // suppress param. update message
    search_algorithm.setParameters(algo_param);
    search_algorithm.setLogType(this->log_type_);

    ProteinIdentification::SearchParameters search_params;
    search_params.db = in_fasta;
    search_params.setMetaValue("input_mzML", in_mzml);
    search_params.setMetaValue("input_decoys", in_decoy_fasta);
    search_params.setMetaValue("out_xquest_specxml", out_xquest_specxml);
    protein_ids[0].setSearchParameters(search_params);

    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenPepXL");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    protein_ids[0].setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // crosslinking search = MS:1002494

    // run algorithm
    OpenPepXLLFAlgorithm::ExitCodes exit_code = search_algorithm.run(unprocessed_spectra, fasta_db, protein_ids, peptide_ids, all_top_csms, spectra);

    if (exit_code != OpenPepXLLFAlgorithm::EXECUTION_OK)
    {
      if (exit_code == OpenPepXLLFAlgorithm::ILLEGAL_PARAMETERS)
      {
        return ILLEGAL_PARAMETERS;
      }
    }

    // MS path already set in algorithm. Overwrite here so we get something testable
    if (getFlag_("test"))
    {
      // if test mode set, add file without path so we can compare it
      protein_ids[0].setPrimaryMSRunPath({"file://" + File::basename(in_mzml)});
    }

    // write output
    progresslogger.startProgress(0, 1, "Writing output...");
    if (!out_idXML.empty())
    {
      FileHandler().storeIdentifications(out_idXML, protein_ids, peptide_ids, {FileTypes::IDXML});
    }
    if (!out_mzIdentML.empty())
    {
      FileHandler().storeIdentifications(out_mzIdentML, protein_ids, peptide_ids, {FileTypes::MZIDENTML});
    }

    if (!out_xquest.empty() || !out_xquest_specxml.empty())
    {
      vector<String> input_split_dir;
      vector<String> input_split;
      getStringOption_("in").split(String("/"), input_split_dir);
      input_split_dir[input_split_dir.size()-1].split(String("."), input_split);
      String base_name = input_split[0];

      if (!out_xquest.empty())
      {
        FileHandler().storeIdentifications(out_xquest, protein_ids, peptide_ids, {FileTypes::XQUESTXML});
      }
      if (!out_xquest_specxml.empty())
      {
        XQuestResultXMLFile::writeXQuestXMLSpec(out_xquest_specxml, base_name, all_top_csms, spectra, test_mode_);
      }
    }
    progresslogger.endProgress();

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{

  TOPPOpenPepXLLF tool;

  return tool.main(argc, argv);
}

/// @endcond
