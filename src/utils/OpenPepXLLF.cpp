
// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/ANALYSIS/XLMS/XQuestScores.h>
#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/ANALYSIS/RNPXL/ModifiedPeptideGenerator.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>

#include <iostream>
#include <cmath>
#include <numeric>

using namespace std;
using namespace OpenMS;

#ifdef _OPENMP
#include <omp.h>
#define NUMBER_OF_THREADS (omp_get_num_threads())
#else
#define NUMBER_OF_THREADS (1)
#endif

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenPepXLLF OpenPepXLLF

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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ OpenPepXLLF \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> - </td>
        </tr>
    </table>
  </CENTER>

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenPepXLLF.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenPepXLLF.html
*/

class TOPPOpenPepXLLF :
  public TOPPBase
{
public:
  TOPPOpenPepXLLF() :
    TOPPBase("OpenPepXLLF", "Tool for protein-protein cross linking with label-free linkers.", false)
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

    registerStringOption_("decoy_string", "<string>", "decoy", "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.", false, false);
    registerFlag_("decoy_prefix", "Set flag, if the decoy_string is a prefix of accessions in the protein database. Otherwise it is a suffix.", false);

    registerTOPPSubsection_("precursor", "Precursor (Parent Ion) Options");
    registerDoubleOption_("precursor:mass_tolerance", "<tolerance>", 10.0, "Width of precursor mass tolerance window", false, false);
    StringList precursor_mass_tolerance_unit_valid_strings;
    precursor_mass_tolerance_unit_valid_strings.push_back("ppm");
    precursor_mass_tolerance_unit_valid_strings.push_back("Da");
    registerStringOption_("precursor:mass_tolerance_unit", "<unit>", "ppm", "Unit of precursor mass tolerance.",  false, false);
    setValidStrings_("precursor:mass_tolerance_unit", precursor_mass_tolerance_unit_valid_strings);

    registerIntOption_("precursor:min_charge", "<num>", 3, "Minimum precursor charge to be considered.", false, true);
    registerIntOption_("precursor:max_charge", "<num>", 7, "Maximum precursor charge to be considered.", false, true);
    registerIntList_("precursor:corrections", "<num>", ListUtils::create<int>("2, 1, 0"), "Monoisotopic peak correction. Matches candidates for possible monoisotopic precursor peaks for experimental mass m and given numbers n at masses (m - n * (C13-C12)). These should be ordered from more extreme to less extreme corrections. Numbers later in the list will be preferred in case of ambiguities.", false, false);


    registerTOPPSubsection_("fragment", "Fragments (Product Ion) Options");
    registerDoubleOption_("fragment:mass_tolerance", "<tolerance>", 20.0, "Fragment mass tolerance", false, false);
    registerDoubleOption_("fragment:mass_tolerance_xlinks", "<tolerance>", 20.0, "Fragment mass tolerance for cross-link ions", false, false);

    StringList fragment_mass_tolerance_unit_valid_strings;
    fragment_mass_tolerance_unit_valid_strings.push_back("ppm");
    fragment_mass_tolerance_unit_valid_strings.push_back("Da");

    registerStringOption_("fragment:mass_tolerance_unit", "<unit>", "ppm", "Unit of fragment m", false, false);
    setValidStrings_("fragment:mass_tolerance_unit", fragment_mass_tolerance_unit_valid_strings);

    registerTOPPSubsection_("modifications", "Modifications Options");
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("modifications:fixed", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false, false);
    setValidStrings_("modifications:fixed", all_mods);
    registerStringList_("modifications:variable", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Oxidation (M)'", false, false);
    setValidStrings_("modifications:variable", all_mods);
    registerIntOption_("modifications:variable_max_per_peptide", "<num>", 2, "Maximum number of residues carrying a variable modification per candidate peptide", false, false);

    registerTOPPSubsection_("peptide", "Peptide Options");
    registerIntOption_("peptide:min_size", "<num>", 5, "Minimum size a peptide must have after digestion to be considered in the search.", false, false);
    registerIntOption_("peptide:missed_cleavages", "<num>", 2, "Number of missed cleavages.", false, false);
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllNames(all_enzymes);
    registerStringOption_("peptide:enzyme", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("peptide:enzyme", all_enzymes);


    registerTOPPSubsection_("cross_linker", "Cross Linker Options");
    registerStringList_("cross_linker:residue1", "<one letter code>", ListUtils::create<String>("K,N-term"), "Comma separated residues, that the first side of a bifunctional cross-linker can attach to", false, false);
    registerStringList_("cross_linker:residue2", "<one letter code>", ListUtils::create<String>("K,N-term"), "Comma separated residues, that the second side of a bifunctional cross-linker can attach to", false, false);
    registerDoubleOption_("cross_linker:mass", "<mass>", 138.0680796, "Mass of the light cross-linker, linking two residues on one or two peptides", false, false);
    registerDoubleList_("cross_linker:mass_mono_link", "<mass>", ListUtils::create<double>("156.07864431, 155.094628715"), "Possible masses of the linker, when attached to only one peptide", false, false);
    registerStringOption_("cross_linker:name", "<string>", "DSS" ,  "Name of the searched cross-link, used to resolve ambiguity of equal masses (e.g. DSS or BS3)", false, false);

    registerTOPPSubsection_("algorithm", "Algorithm Options");
    registerFlag_("algorithm:pre_scoring", "Set flag, to use the pre-scoring heuristic. Otherwise a full enumeration will be performed.", false);
    registerIntOption_("algorithm:number_of_scored_candidates", "<num>", 10000, "Number of candidates that are scored using the full scoring function after ranking by a faster pre-scoring function.", false, false);
    registerIntOption_("algorithm:number_top_hits", "<num>", 5, "Number of top hits reported for each spectrum pair", false, false);
    StringList deisotope_strings = ListUtils::create<String>("true,false,auto");
    registerStringOption_("algorithm:deisotope", "<true/false/auto>", "auto", "Set to true, if the input spectra should be deisotoped before any other processing steps. If set to auto the spectra will be deisotoped, if the fragment mass tolerance is < 0.1 Da or < 100 ppm (0.1 Da at a mass of 1000)", false, true);
    setValidStrings_("algorithm:deisotope", deisotope_strings);

    StringList bool_strings = ListUtils::create<String>("true,false");
    registerTOPPSubsection_("ions", "Ion types to search for");
    registerStringOption_("ions:b_ions", "<true/false>", "true", "Search for peaks of b-ions.", false, true);
    setValidStrings_("ions:b_ions", bool_strings);
    registerStringOption_("ions:y_ions", "<true/false>", "true", "Search for peaks of y-ions.", false, true);
    setValidStrings_("ions:y_ions", bool_strings);
    registerStringOption_("ions:a_ions", "<true/false>", "false", "Search for peaks of a-ions.", false, true);
    setValidStrings_("ions:a_ions", bool_strings);
    registerStringOption_("ions:x_ions", "<true/false>", "false", "Search for peaks of x-ions.", false, true);
    setValidStrings_("ions:x_ions", bool_strings);
    registerStringOption_("ions:c_ions", "<true/false>", "false", "Search for peaks of c-ions.", false, true);
    setValidStrings_("ions:c_ions", bool_strings);
    registerStringOption_("ions:z_ions", "<true/false>", "false", "Search for peaks of z-ions.", false, true);
    setValidStrings_("ions:z_ions", bool_strings);
    registerStringOption_("ions:neutral_losses", "<true/false>", "true", "Search for neutral losses of H2O and H3N.", false, true);
    setValidStrings_("ions:neutral_losses", bool_strings);

    // output file
    registerOutputFile_("out_xquestxml", "<file>", "", "Results in the xquest.xml format (at least one of these output parameters should be set, otherwise you will not have any results).", false);
    setValidFormats_("out_xquestxml", ListUtils::create<String>("xml,xquest.xml"));

    registerOutputFile_("out_xquest_specxml", "<file>", "", "Matched spectra in the xQuest .spec.xml format for spectra visualization in the xQuest results manager.", false, false);
    setValidFormats_("out_xquest_specxml", ListUtils::create<String>("xml,spec.xml"));

    registerOutputFile_("out_idXML", "<file>", "", "Results in idXML format (at least one of these output parameters should be set, otherwise you will not have any results).", false);
    setValidFormats_("out_idXML", ListUtils::create<String>("idXML"));

    registerOutputFile_("out_mzIdentML", "<file>","", "Results in mzIdentML (.mzid) format (at least one of these output parameters should be set, otherwise you will not have any results)", false);
    setValidFormats_("out_mzIdentML", ListUtils::create<String>("mzid"));
  }

  vector<ResidueModification> getModifications_(StringList modNames)
  {
    vector<ResidueModification> modifications;

    // iterate over modification names and add to vector
    for (StringList::iterator mod_it = modNames.begin(); mod_it != modNames.end(); ++mod_it)
    {
      String modification(*mod_it);
      modifications.push_back(ModificationsDB::getInstance()->getModification(modification));
    }

    return modifications;
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

    const bool decoy_prefix(getFlag_("decoy_prefix"));
    const string decoy_string(getStringOption_("decoy_string"));

    Int min_precursor_charge = getIntOption_("precursor:min_charge");
    Int max_precursor_charge = getIntOption_("precursor:max_charge");
    double precursor_mass_tolerance = getDoubleOption_("precursor:mass_tolerance");
    bool precursor_mass_tolerance_unit_ppm = (getStringOption_("precursor:mass_tolerance_unit") == "ppm");
    IntList precursor_correction_steps = getIntList_("precursor:corrections");

    double fragment_mass_tolerance = getDoubleOption_("fragment:mass_tolerance");
    double fragment_mass_tolerance_xlinks = getDoubleOption_("fragment:mass_tolerance_xlinks");
    if (fragment_mass_tolerance_xlinks < fragment_mass_tolerance)
    {
      fragment_mass_tolerance_xlinks = fragment_mass_tolerance;
    }
    bool fragment_mass_tolerance_unit_ppm = (getStringOption_("fragment:mass_tolerance_unit") == "ppm");

    StringList cross_link_residue1 = getStringList_("cross_linker:residue1");
    StringList cross_link_residue2 = getStringList_("cross_linker:residue2");
    double cross_link_mass = getDoubleOption_("cross_linker:mass");
    DoubleList cross_link_mass_mono_link = getDoubleList_("cross_linker:mass_mono_link");
    std::sort(cross_link_mass_mono_link.begin(), cross_link_mass_mono_link.end(), std::greater< double >());
    String cross_link_name = getStringOption_("cross_linker:name");

    StringList fixedModNames = getStringList_("modifications:fixed");
    set<String> fixed_unique(fixedModNames.begin(), fixedModNames.end());

    Size peptide_min_size = getIntOption_("peptide:min_size");

    Int number_top_hits = getIntOption_("algorithm:number_top_hits");
    bool pre_scoring(getFlag_("algorithm:pre_scoring"));
    Int number_of_scored_candidates = getIntOption_("algorithm:number_of_scored_candidates");
    String deisotope_mode = getStringOption_("algorithm:deisotope");

    // deisotope if "true" or if "auto" and the tolerance is below the threshold (0.1 Da or 100 ppm)
    bool deisotope = (deisotope_mode == "true") ||
                      (deisotope_mode == "auto" &&
                      ((!fragment_mass_tolerance_unit_ppm && fragment_mass_tolerance < 0.1) ||
                      (fragment_mass_tolerance_unit_ppm && fragment_mass_tolerance < 100)));

    if (fixed_unique.size() != fixedModNames.size())
    {
      LOG_DEBUG << "duplicate fixed modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }

    StringList varModNames = getStringList_("modifications:variable");
    set<String> var_unique(varModNames.begin(), varModNames.end());
    if (var_unique.size() != varModNames.size())
    {
      LOG_DEBUG << "duplicate variable modification provided." << endl;
      return ILLEGAL_PARAMETERS;
    }
    vector<ResidueModification> fixed_modifications = getModifications_(fixedModNames);
    vector<ResidueModification> variable_modifications = getModifications_(varModNames);
    Size max_variable_mods_per_peptide = getIntOption_("modifications:variable_max_per_peptide");

    // load MS2 map
    PeakMap unprocessed_spectra;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in_mzml, unprocessed_spectra);

    // preprocess spectra (filter out 0 values, sort by position)
    progresslogger.startProgress(0, 1, "Filtering spectra...");
    PeakMap spectra = OPXLSpectrumProcessingAlgorithms::preprocessSpectra(unprocessed_spectra, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, peptide_min_size, min_precursor_charge, max_precursor_charge, deisotope, false);
    progresslogger.endProgress();

    // load linked features
    ConsensusMap cfeatures;
    ConsensusXMLFile cf;

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

    const Size missed_cleavages = getIntOption_("peptide:missed_cleavages");
    ProteaseDigestion digestor;
    String enzyme_name = getStringOption_("peptide:enzyme");
    digestor.setEnzyme(enzyme_name);
    digestor.setMissedCleavages(missed_cleavages);

    // set minimum size of peptide after digestion
    Size min_peptide_length = getIntOption_("peptide:min_size");

    // one identification run
    vector<ProteinIdentification> protein_ids(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("OpenXQuest");
    protein_ids[0].setSearchEngineVersion(VersionInfo::getVersion());
    StringList ms_runs;
    unprocessed_spectra.getPrimaryMSRunPath(ms_runs);
    protein_ids[0].setPrimaryMSRunPath(ms_runs);
    protein_ids[0].setMetaValue("SpectrumIdentificationProtocol", DataValue("MS:1002494")); // cross-linking search = MS:1002494

    ProteinIdentification::SearchParameters search_params;
    String searched_charges((String(min_precursor_charge)));
    for (int ch = min_precursor_charge+1; ch <= max_precursor_charge; ++ch)
    {
      searched_charges += "," + String(ch);
    }
    search_params.charges = searched_charges;
    search_params.db = in_fasta;
    search_params.digestion_enzyme = (*ProteaseDB::getInstance()->getEnzyme(enzyme_name));
    search_params.fixed_modifications = fixedModNames;
    search_params.variable_modifications = varModNames;
    search_params.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_params.missed_cleavages = missed_cleavages;
    search_params.fragment_mass_tolerance = fragment_mass_tolerance;
    search_params.fragment_mass_tolerance_ppm =  fragment_mass_tolerance_unit_ppm;
    search_params.precursor_mass_tolerance = precursor_mass_tolerance;
    search_params.precursor_mass_tolerance_ppm = precursor_mass_tolerance_unit_ppm;

    // As MetaValues
    search_params.setMetaValue("input_mzML", in_mzml);
    search_params.setMetaValue("input_decoys", in_decoy_fasta);
    search_params.setMetaValue("decoy_prefix", decoy_prefix);
    search_params.setMetaValue("decoy_string", decoy_string);
    search_params.setMetaValue("out_xquest_specxml", out_xquest_specxml);

    search_params.setMetaValue("precursor:min_charge", min_precursor_charge);
    search_params.setMetaValue("precursor:max_charge", max_precursor_charge);

    search_params.setMetaValue("fragment:mass_tolerance_xlinks", fragment_mass_tolerance_xlinks);
    search_params.setMetaValue("peptide:min_size", peptide_min_size);

    search_params.setMetaValue("cross_link:residue1", cross_link_residue1);
    search_params.setMetaValue("cross_link:residue2", cross_link_residue2);
    search_params.setMetaValue("cross_link:mass", cross_link_mass);
    search_params.setMetaValue("cross_link:mass_monolink", cross_link_mass_mono_link);
    search_params.setMetaValue("cross_link:name", cross_link_name);

    search_params.setMetaValue("modifications:variable_max_per_peptide", max_variable_mods_per_peptide);
    protein_ids[0].setSearchParameters(search_params);

    vector<PeptideIdentification> peptide_ids;

    // lookup for processed peptides. must be defined outside of omp section and synchronized
    vector<OPXLDataStructs::AASeqWithMass> peptide_masses;

    progresslogger.startProgress(0, 1, "Digesting peptides...");
    peptide_masses = OPXLHelper::digestDatabase(fasta_db, digestor, min_peptide_length, cross_link_residue1, cross_link_residue2, fixed_modifications,  variable_modifications, max_variable_mods_per_peptide);
    progresslogger.endProgress();

    // declare and set up spectrum generators
    TheoreticalSpectrumGeneratorXLMS specGen_fast;
    TheoreticalSpectrumGeneratorXLMS specGen_full;

    // Setting parameters for cross-link fragmentation
    // settings for pre-scoring, only a-, b- and y- ions without annotation etc.
    Param specGenParams_fast = specGen_fast.getParameters();
    specGenParams_fast.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams_fast.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    specGenParams_fast.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    specGenParams_fast.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    specGenParams_fast.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    specGenParams_fast.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    specGenParams_fast.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");

    specGenParams_fast.setValue("add_metainfo", "false");
    specGenParams_fast.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams_fast.setValue("max_isotope", 1, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");
    specGenParams_fast.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    specGenParams_fast.setValue("add_precursor_peaks", "false");
    specGenParams_fast.setValue("add_k_linked_ions", "false");
    specGen_fast.setParameters(specGenParams_fast);


    // settings fpr full-scoring, annotations, 2nd isotopic peaks, losses and precursors
    Param specGenParams_full = specGen_full.getParameters();
    specGenParams_full.setValue("add_b_ions", getStringOption_("ions:b_ions"), "Add peaks of y-ions to the spectrum");
    specGenParams_full.setValue("add_y_ions", getStringOption_("ions:y_ions"), "Add peaks of b-ions to the spectrum");
    specGenParams_full.setValue("add_a_ions", getStringOption_("ions:a_ions"), "Add peaks of a-ions to the spectrum");
    specGenParams_full.setValue("add_x_ions", getStringOption_("ions:x_ions"), "Add peaks of c-ions to the spectrum");
    specGenParams_full.setValue("add_c_ions", getStringOption_("ions:c_ions"), "Add peaks of x-ions to the spectrum");
    specGenParams_full.setValue("add_z_ions", getStringOption_("ions:z_ions"), "Add peaks of z-ions to the spectrum");
    specGenParams_full.setValue("add_losses", getStringOption_("ions:neutral_losses"), "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");

    specGenParams_full.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    specGenParams_full.setValue("add_metainfo", "true");
    specGenParams_full.setValue("add_isotopes", "true", "If set to 1 isotope peaks of the product ion peaks are added");
    specGenParams_full.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    specGenParams_full.setValue("add_precursor_peaks", "true");
    specGenParams_full.setValue("add_k_linked_ions", "true");
    specGen_full.setParameters(specGenParams_full);

    LOG_DEBUG << "Peptide candidates: " << peptide_masses.size() << endl;
    search_params = protein_ids[0].getSearchParameters();
    search_params.setMetaValue("MS:1001029", peptide_masses.size()); // number of sequences searched = MS:1001029
    protein_ids[0].setSearchParameters(search_params);

    // Collect precursor MZs for filtering enumerated peptide pairs
    vector< double > spectrum_precursors;
    for (Size i = 0; i < spectra.size(); i++)
    {
      double current_precursor_mz = spectra[i].getPrecursors()[0].getMZ();
      double current_precursor_charge = spectra[i].getPrecursors()[0].getCharge();
      double current_precursor_mass = (current_precursor_mz * current_precursor_charge) - (current_precursor_charge * Constants::PROTON_MASS_U);
      spectrum_precursors.push_back(current_precursor_mass);
    }
    sort(spectrum_precursors.begin(), spectrum_precursors.end());
    LOG_DEBUG << "Number of precursor masses in the spectra: " << spectrum_precursors.size() << endl;

    sort(peptide_masses.begin(), peptide_masses.end(), OPXLDataStructs::AASeqWithMassComparator());
    // The largest peptides given a fixed maximal precursor mass are possible with loop links
    // Filter peptides using maximal loop link mass first
    double max_precursor_mass = spectrum_precursors[spectrum_precursors.size()-1];

    // compute absolute tolerance from relative, if necessary
    double max_peptide_allowed_error = 0;
    if (precursor_mass_tolerance_unit_ppm) // ppm
    {
      max_peptide_allowed_error = max_precursor_mass * precursor_mass_tolerance * 1e-6;
    }
    else // Dalton
    {
      max_peptide_allowed_error= precursor_mass_tolerance;
    }

    double max_peptide_mass = max_precursor_mass - cross_link_mass + max_peptide_allowed_error;

    LOG_DEBUG << "Filtering peptides with precursors" << endl;

    // search for the first mass greater than the maximim, cut off everything larger
    vector<OPXLDataStructs::AASeqWithMass>::iterator last = upper_bound(peptide_masses.begin(), peptide_masses.end(), max_peptide_mass, OPXLDataStructs::AASeqWithMassComparator());
    vector<OPXLDataStructs::AASeqWithMass> filtered_peptide_masses;
    filtered_peptide_masses.assign(peptide_masses.begin(), last);

    // iterate over all spectra
    progresslogger.startProgress(0, 1, "Matching to theoretical spectra and scoring...");
    vector< vector< OPXLDataStructs::CrossLinkSpectrumMatch > > all_top_csms;

    Size spectrum_counter = 0;

    LOG_DEBUG << "Spectra left after preprocessing and filtering: " << spectra.size() << " of " << unprocessed_spectra.size() << endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(guided)
#endif
    for (SignedSize scan_index = 0; scan_index < static_cast<SignedSize>(spectra.size()); ++scan_index)
    {
      const PeakSpectrum& spectrum = spectra[scan_index];

#ifdef _OPENMP
#pragma omp critical
#endif
      {
        spectrum_counter++;
        cout << "Processing spectrum " << spectrum_counter << " / " << spectra.size() << " |\tSpectrum index: " << scan_index << "\t| at: " << DateTime::now().getTime() << endl;
      }

      const double precursor_charge = spectrum.getPrecursors()[0].getCharge();
      const double precursor_mz = spectrum.getPrecursors()[0].getMZ();
      const double precursor_mass = (precursor_mz * static_cast<double>(precursor_charge)) - (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U);

      vector< OPXLDataStructs::CrossLinkSpectrumMatch > top_csms_spectrum;

      // determine candidates
      vector< OPXLDataStructs::XLPrecursor > candidates;

      vector< double > spectrum_precursor_vector;
      vector< double > allowed_error_vector;

      for (int correction_mass : precursor_correction_steps)
      {
        double allowed_error = 0;

        double corrected_precursor_mass = precursor_mass - (static_cast<double>(correction_mass) * Constants::C13C12_MASSDIFF_U);

        if (precursor_mass_tolerance_unit_ppm) // ppm
        {
          allowed_error = corrected_precursor_mass * precursor_mass_tolerance * 1e-6;
        }
        else // Dalton
        {
          allowed_error = precursor_mass_tolerance;
        }

        spectrum_precursor_vector.push_back(corrected_precursor_mass);
        allowed_error_vector.push_back(allowed_error);

      } // end correction mass loop

      vector< int > precursor_correction_positions;
      candidates = OPXLHelper::enumerateCrossLinksAndMasses(filtered_peptide_masses, cross_link_mass, cross_link_mass_mono_link, cross_link_residue1, cross_link_residue2, spectrum_precursor_vector, precursor_correction_positions, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

      vector< int > precursor_corrections;
      for (Size pc = 0; pc < precursor_correction_positions.size(); ++pc)
      {
        precursor_corrections.push_back(precursor_correction_steps[precursor_correction_positions[pc]]);
      }

      vector <OPXLDataStructs::ProteinProteinCrossLink> cross_link_candidates = OPXLHelper::buildCandidates(candidates, precursor_corrections, precursor_correction_positions, filtered_peptide_masses, cross_link_residue1, cross_link_residue2, cross_link_mass, cross_link_mass_mono_link, spectrum_precursor_vector, allowed_error_vector, cross_link_name);

      LOG_DEBUG << "Size of enumerated candidates: " << double(cross_link_candidates.size()) * sizeof(OPXLDataStructs::ProteinProteinCrossLink) / 1024.0 / 1024.0 << " mb" << endl;

#ifdef _OPENMP
#pragma omp critical
#endif
      cout << "Spectrum number: " << spectrum_counter << " |\tNumber of peaks: " << spectrum.size() << " |\tNumber of candidates: " << cross_link_candidates.size() << endl;


      // lists for one spectrum, to determine best match to the spectrum
      vector< OPXLDataStructs::CrossLinkSpectrumMatch > prescore_csms_spectrum;
      vector< OPXLDataStructs::CrossLinkSpectrumMatch > all_csms_spectrum;

      if (pre_scoring)
      {
        progresslogger.startProgress(0, 1, "Start pre-scoring...");

        for (Size i = 0; i < cross_link_candidates.size(); ++i)
        {
          OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate = cross_link_candidates[i];

          PeakSpectrum theoretical_spec_linear_alpha;
          PeakSpectrum theoretical_spec_linear_beta;

          bool type_is_cross_link = cross_link_candidate.getType() == OPXLDataStructs::CROSS;
          bool type_is_loop = cross_link_candidate.getType() == OPXLDataStructs::LOOP;
          Size link_pos_B = 0;
          if (type_is_loop)
          {
            link_pos_B = cross_link_candidate.cross_link_position.second;
          }
          specGen_fast.getLinearIonSpectrum(theoretical_spec_linear_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, true, 1, link_pos_B);
          if (type_is_cross_link)
          {
            specGen_fast.getLinearIonSpectrum(theoretical_spec_linear_beta, cross_link_candidate.beta, cross_link_candidate.cross_link_position.second, false, 1);
          }

          // Something like this can happen, e.g. with a loop link connecting residues close to the ends of the peptide
          if ( (theoretical_spec_linear_alpha.size() < 3) )
          {
            continue;
          }

          PeakSpectrum theoretical_spec_linear;
          theoretical_spec_linear.insert(theoretical_spec_linear.end(), theoretical_spec_linear_alpha.begin(), theoretical_spec_linear_alpha.end());
          theoretical_spec_linear.insert(theoretical_spec_linear.end(), theoretical_spec_linear_beta.begin(), theoretical_spec_linear_beta.end());
          theoretical_spec_linear.sortByPosition();

          double xcorr_prescore = XQuestScores::xCorrelationPrescore(spectrum, theoretical_spec_linear, 0.1);

          if (xcorr_prescore == 0)
          {
            continue;
          }

          // store pre_score, take top 100 for further computations
          OPXLDataStructs::CrossLinkSpectrumMatch csm;
          csm.cross_link = cross_link_candidate;
          csm.score = xcorr_prescore;
          csm.pre_score = xcorr_prescore;
          csm.xcorrc_max = xcorr_prescore;

          prescore_csms_spectrum.push_back(csm);
        }
        progresslogger.endProgress();
        std::sort(prescore_csms_spectrum.rbegin(), prescore_csms_spectrum.rend(), OPXLDataStructs::CLSMScoreComparator());
      }
      else // if pre-scoring is turned off, then the pre-scored candidates will just be the list of all candidates
      {
        for (Size i = 0; i < cross_link_candidates.size(); ++i)
        {
          OPXLDataStructs::CrossLinkSpectrumMatch csm;
          csm.cross_link = cross_link_candidates[i];
          csm.pre_score = 0;
          prescore_csms_spectrum.push_back(csm);
        }
      }

      LOG_DEBUG << "Size of prescored candidates: " << double(prescore_csms_spectrum.size()) * sizeof(OPXLDataStructs::CrossLinkSpectrumMatch) / 1024.0 / 1024.0 << " mb" << endl;

      Size last_candidate_index = prescore_csms_spectrum.size();
      if (pre_scoring)
      {
        last_candidate_index = std::min(last_candidate_index, Size(number_of_scored_candidates));
      }

      for (Size i = 0; i < last_candidate_index ; ++i)
      {
        OPXLDataStructs::ProteinProteinCrossLink cross_link_candidate = prescore_csms_spectrum[i].cross_link;
        double candidate_mz = (cross_link_candidate.alpha.getMonoWeight() + cross_link_candidate.beta.getMonoWeight() +  cross_link_candidate.cross_linker_mass+ (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)) / precursor_charge;

        LOG_DEBUG << "Pair: " << cross_link_candidate.alpha.toString() << "-" << cross_link_candidate.beta.toString() << " matched to light spectrum " << scan_index << "\t and heavy spectrum " << scan_index
            << " with m/z: " << precursor_mz << "\t" << "and candidate m/z: " << candidate_mz << "\tK Positions: " << cross_link_candidate.cross_link_position.first << "\t" << cross_link_candidate.cross_link_position.second << endl;

        OPXLDataStructs::CrossLinkSpectrumMatch csm;
        csm.cross_link = cross_link_candidate;

        PeakSpectrum theoretical_spec_linear_alpha;
        PeakSpectrum theoretical_spec_linear_beta;
        PeakSpectrum theoretical_spec_xlinks_alpha;
        PeakSpectrum theoretical_spec_xlinks_beta;

        bool type_is_cross_link = cross_link_candidate.getType() == OPXLDataStructs::CROSS;
        bool type_is_loop = cross_link_candidate.getType() == OPXLDataStructs::LOOP;
        Size link_pos_B = 0;
        if (type_is_loop)
        {
          link_pos_B = cross_link_candidate.cross_link_position.second;
        }
        specGen_full.getLinearIonSpectrum(theoretical_spec_linear_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, true, 2, link_pos_B);
        if (type_is_cross_link)
        {
          specGen_full.getLinearIonSpectrum(theoretical_spec_linear_beta, cross_link_candidate.beta, cross_link_candidate.cross_link_position.second, false, 2);
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate, true, 1, precursor_charge);
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_beta, cross_link_candidate, false, 1, precursor_charge);
        }
        else
        {
          // Function for mono-links or loop-links
          specGen_full.getXLinkIonSpectrum(theoretical_spec_xlinks_alpha, cross_link_candidate.alpha, cross_link_candidate.cross_link_position.first, precursor_mass, true, 2, precursor_charge, link_pos_B);
        }

        // Something like this can happen, e.g. with a loop link connecting the first and last residue of a peptide
        if ( (theoretical_spec_linear_alpha.size() < 1) || (theoretical_spec_xlinks_alpha.size() < 1) )
        {
          continue;
        }

        vector< pair< Size, Size > > matched_spec_linear_alpha;
        vector< pair< Size, Size > > matched_spec_linear_beta;
        vector< pair< Size, Size > > matched_spec_xlinks_alpha;
        vector< pair< Size, Size > > matched_spec_xlinks_beta;

        DataArrays::FloatDataArray ppm_error_array_linear_alpha;
        DataArrays::FloatDataArray ppm_error_array_xlinks_alpha;
        DataArrays::FloatDataArray ppm_error_array_linear_beta;
        DataArrays::FloatDataArray ppm_error_array_xlinks_beta;

        PeakSpectrum::IntegerDataArray& theo_charges_la = theoretical_spec_linear_alpha.getIntegerDataArrays()[0];
        PeakSpectrum::IntegerDataArray theo_charges_xa;
        if (theoretical_spec_xlinks_alpha.getIntegerDataArrays().size() > 0)
        {
          theo_charges_xa = theoretical_spec_xlinks_alpha.getIntegerDataArrays()[0];
        }
        PeakSpectrum::IntegerDataArray theo_charges_lb;
        PeakSpectrum::IntegerDataArray theo_charges_xb;
        if (theoretical_spec_linear_beta.getIntegerDataArrays().size() > 0)
        {
          theo_charges_lb = theoretical_spec_linear_beta.getIntegerDataArrays()[0];
        }
        if (theoretical_spec_xlinks_beta.getIntegerDataArrays().size() > 0)
        {
          theo_charges_xb = theoretical_spec_xlinks_beta.getIntegerDataArrays()[0];
        }
        PeakSpectrum::IntegerDataArray exp_charges;
        if (spectrum.getIntegerDataArrays().size() > 0)
        {
          exp_charges = spectrum.getIntegerDataArrays()[0];
        }
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_alpha, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, theoretical_spec_linear_alpha, spectrum, theo_charges_la, exp_charges, ppm_error_array_linear_alpha);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_linear_beta, fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm, theoretical_spec_linear_beta, spectrum, theo_charges_lb, exp_charges, ppm_error_array_linear_beta);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_alpha, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, theoretical_spec_xlinks_alpha, spectrum, theo_charges_xa, exp_charges, ppm_error_array_xlinks_alpha);
        OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(matched_spec_xlinks_beta, fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, theoretical_spec_xlinks_beta, spectrum, theo_charges_xb, exp_charges, ppm_error_array_xlinks_beta);

        LOG_DEBUG << "Spectrum sizes: " << spectrum.size() << " || " << theoretical_spec_linear_alpha.size() <<  " | " << theoretical_spec_linear_beta.size()
                              <<  " | " << theoretical_spec_xlinks_alpha.size() <<  " | " << theoretical_spec_xlinks_beta.size() << endl;
        LOG_DEBUG << "Matched peaks: " << matched_spec_linear_alpha.size() << " | " << matched_spec_linear_beta.size()
                              <<  " | " << matched_spec_xlinks_alpha.size() <<  " | " << matched_spec_xlinks_beta.size() << endl;

        // TODO define good exclusion criteria for total crap
        Size matched_peaks = matched_spec_linear_alpha.size() + matched_spec_linear_beta.size() + matched_spec_xlinks_alpha.size() + matched_spec_xlinks_beta.size();
        if (matched_peaks < 1)
        {
          continue;
        }

        LOG_DEBUG << "Computing Intsum..." << endl;
        // compute intsum score
        double intsum = XQuestScores::totalMatchedCurrent(matched_spec_linear_alpha, matched_spec_linear_beta, matched_spec_xlinks_alpha, matched_spec_xlinks_beta, spectrum, spectrum);

        // Total ion intensity of light spectrum
        // sum over linear and xlink ion spectra instead of unfiltered
        double total_current = 0;
        for (SignedSize j = 0; j < static_cast<SignedSize>(spectrum.size()); ++j)
        {
          total_current += spectrum[j].getIntensity();
        }
        double TIC = intsum / total_current;

        // TIC_alpha and _beta (total ion current)
        double intsum_alpha = XQuestScores::matchedCurrentChain(matched_spec_linear_alpha, matched_spec_xlinks_alpha, spectrum, spectrum);
        double intsum_beta = 0;
        if (type_is_cross_link)
        {
          intsum_beta = XQuestScores::matchedCurrentChain(matched_spec_linear_beta, matched_spec_xlinks_beta, spectrum, spectrum);
        }

        // normalize TIC_alpha and  _beta
        if ((intsum_alpha + intsum_beta) > 0.0)
        {
          intsum_alpha = intsum_alpha * intsum / (intsum_alpha + intsum_beta);
          intsum_beta = intsum_beta *  intsum / (intsum_alpha + intsum_beta);
        }

        LOG_DEBUG << "Computing TIC..." << endl;

        // compute weighted TIC
        double wTIC = XQuestScores::weightedTICScore(cross_link_candidate.alpha.size(), cross_link_candidate.beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);
        double wTICold = XQuestScores::weightedTICScoreXQuest(cross_link_candidate.alpha.size(), cross_link_candidate.beta.size(), intsum_alpha, intsum_beta, total_current, type_is_cross_link);

        // maximal xlink ion charge = (Precursor charge - 1), minimal xlink ion charge: 2
        Size n_xlink_charges = (precursor_charge - 1) - 2;
        if (n_xlink_charges < 1) n_xlink_charges = 1;

        LOG_DEBUG << "Computing Match-Odds..." << endl;
        // compute match odds (unweighted), the 3 is the number of charge states in the theoretical spectra
        double match_odds_c_alpha = XQuestScores::matchOddsScore(theoretical_spec_linear_alpha, matched_spec_linear_alpha.size(), fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
        double match_odds_x_alpha = XQuestScores::matchOddsScore(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha.size(), fragment_mass_tolerance_xlinks , fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
        double log_occu_c_alpha = XQuestScores::logOccupancyProb(theoretical_spec_linear_alpha, matched_spec_linear_alpha.size(), fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
        double log_occu_x_alpha = XQuestScores::logOccupancyProb(theoretical_spec_xlinks_alpha, matched_spec_xlinks_alpha.size(), fragment_mass_tolerance_xlinks , fragment_mass_tolerance_unit_ppm);
        double match_odds = 0;
        double log_occu = 0;

        double match_odds_alpha = 0;
        double match_odds_beta = 0;
        double log_occu_alpha = 0;
        double log_occu_beta = 0;

        if (type_is_cross_link)
        {
          double match_odds_c_beta = XQuestScores::matchOddsScore(theoretical_spec_linear_beta, matched_spec_linear_beta.size(), fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
          double match_odds_x_beta = XQuestScores::matchOddsScore(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta.size(), fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm, true, n_xlink_charges);
          double log_occu_c_beta = XQuestScores::logOccupancyProb(theoretical_spec_linear_beta, matched_spec_linear_beta.size(), fragment_mass_tolerance, fragment_mass_tolerance_unit_ppm);
          double log_occu_x_beta = XQuestScores::logOccupancyProb(theoretical_spec_xlinks_beta, matched_spec_xlinks_beta.size(), fragment_mass_tolerance_xlinks, fragment_mass_tolerance_unit_ppm);
          match_odds = (match_odds_c_alpha + match_odds_x_alpha + match_odds_c_beta + match_odds_x_beta) / 4;
          match_odds_alpha = (match_odds_c_alpha + match_odds_x_alpha) / 2;
          match_odds_beta = (match_odds_c_beta + match_odds_x_beta) / 2;
          log_occu = (log_occu_c_alpha + log_occu_x_alpha + log_occu_c_beta + log_occu_x_beta) / 4;
          log_occu_alpha = (log_occu_c_alpha + log_occu_x_alpha) / 2;
          log_occu_beta = (log_occu_c_beta + log_occu_x_beta) / 2;
        }
        else
        {
          match_odds = (match_odds_c_alpha + match_odds_x_alpha) / 2;
          match_odds_alpha = match_odds;
          log_occu = (log_occu_c_alpha + log_occu_x_alpha) / 2;
          log_occu_alpha = log_occu;
        }

        //Cross-correlation
        PeakSpectrum theoretical_spec_linear;
        PeakSpectrum theoretical_spec_xlinks;
        if (type_is_cross_link)
        {
          theoretical_spec_linear = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear_alpha, theoretical_spec_linear_beta);
          theoretical_spec_xlinks = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_xlinks_alpha, theoretical_spec_xlinks_beta);
        }
        else
        {
          theoretical_spec_linear = theoretical_spec_linear_alpha;
          theoretical_spec_xlinks = theoretical_spec_xlinks_alpha;
        }
        //
        PeakSpectrum theoretical_spec = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theoretical_spec_linear, theoretical_spec_xlinks);
        double xcorrx_max = XQuestScores::xCorrelationPrescore(spectrum, theoretical_spec_xlinks, 0.1);
        double xcorrc_max = XQuestScores::xCorrelationPrescore(spectrum, theoretical_spec_linear, 0.1);

        // Compute score from the 4 scores and 4 weights
        // The weights are adapted from the xQuest algorithm (O. Rinner et al., 2008, "Identification of cross-linked peptides from large sequence databases"),
        // they were determined by an Linear Discriminant Analysis on CID fragmentation data.
        // The match-odds score does not work very well on HCD data and label-free cross-linkers (has the maximal possible value very often), so its weight was drastically reduced here.
        double xcorrx_weight = 2.488;
        double xcorrc_weight = 21.279;
        double match_odds_weight = 1.973;
        double wTIC_weight = 12.829;
        double intsum_weight = 1.8;

        double xquest_score = xcorrx_weight * xcorrx_max + xcorrc_weight * xcorrc_max + match_odds_weight * match_odds + wTIC_weight * wTICold + intsum_weight * intsum;
        csm.xquest_score = xquest_score;

        // Error calculation
        double weight = cross_link_candidate.alpha.getMonoWeight();
        if (cross_link_candidate.getType() == OPXLDataStructs::CROSS)
        {
          weight += cross_link_candidate.beta.getMonoWeight() + cross_link_candidate.cross_linker_mass;
        }
        else
        {
          weight += cross_link_candidate.cross_linker_mass;
        }
        double precursor_mass = (precursor_mz * static_cast<double>(precursor_charge)) - (static_cast<double>(precursor_charge) * Constants::PROTON_MASS_U)
                                  - (static_cast<double>(cross_link_candidate.precursor_correction) * Constants::C13C12_MASSDIFF_U);
        double error = precursor_mass - weight;
        double rel_error = (error / precursor_mass) / 1e-6;

        double new_match_odds_weight = 0.2;
        double new_rel_error_weight = -0.03;
        double new_score = new_match_odds_weight * std::log(1e-7 + match_odds) + new_rel_error_weight * abs(rel_error);

        csm.score = new_score;
        csm.pre_score = prescore_csms_spectrum[i].pre_score;
        csm.percTIC = TIC;
        csm.wTIC = wTIC;
        csm.wTICold = wTICold;
        csm.int_sum = intsum;
        csm.intsum_alpha = intsum_alpha;
        csm.intsum_beta = intsum_beta;
        csm.total_current = total_current;
        csm.precursor_error_ppm = rel_error;

        csm.match_odds = match_odds;
        csm.match_odds_alpha = match_odds_alpha;
        csm.match_odds_beta = match_odds_beta;
        csm.log_occupancy = log_occu;
        csm.log_occupancy_alpha = log_occu_alpha;
        csm.log_occupancy_beta = log_occu_beta;

        csm.xcorrx_max = xcorrx_max;
        csm.xcorrc_max = xcorrc_max;

        csm.matched_linear_alpha = matched_spec_linear_alpha.size();
        csm.matched_linear_beta = matched_spec_linear_beta.size();
        csm.matched_xlink_alpha = matched_spec_xlinks_alpha.size();
        csm.matched_xlink_beta = matched_spec_xlinks_beta.size();
        csm.scan_index_light = scan_index;
        csm.scan_index_heavy = -1;

        // default values to avoid too many "else {var = 0}" in the code below
        csm.num_iso_peaks_mean = 0;
        csm.num_iso_peaks_mean_linear_alpha = 0;
        csm.num_iso_peaks_mean_linear_beta = 0;
        csm.num_iso_peaks_mean_xlinks_alpha = 0;
        csm.num_iso_peaks_mean_xlinks_beta = 0;

        csm.precursor_correction = cross_link_candidate.precursor_correction;

        // num_iso_peaks array from deisotoping
        if (deisotope)
        {
          LOG_DEBUG << "Computing Iso Peak summeries..." << endl;

          DataArrays::IntegerDataArray num_iso_peaks_array;
          auto num_iso_peaks_array_it = getDataArrayByName(spectrum.getIntegerDataArrays(), "NumIsoPeaks");
          num_iso_peaks_array = *num_iso_peaks_array_it;

          csm.num_iso_peaks_mean = Math::mean(num_iso_peaks_array.begin(), num_iso_peaks_array.end());

          vector< double > iso_peaks_linear_alpha;
          vector< double > iso_peaks_linear_beta;
          vector< double > iso_peaks_xlinks_alpha;
          vector< double > iso_peaks_xlinks_beta;

          if (!matched_spec_linear_alpha.empty())
          {
            for (auto match : matched_spec_linear_alpha)
            {
              iso_peaks_linear_alpha.push_back(num_iso_peaks_array[match.second]);
            }
            csm.num_iso_peaks_mean_linear_alpha = Math::mean(iso_peaks_linear_alpha.begin(), iso_peaks_linear_alpha.end());
          }

          if (!matched_spec_linear_beta.empty())
          {
            for (auto match : matched_spec_linear_beta)
            {
              iso_peaks_linear_beta.push_back(num_iso_peaks_array[match.second]);
            }
            csm.num_iso_peaks_mean_linear_beta = Math::mean(iso_peaks_linear_beta.begin(), iso_peaks_linear_beta.end());
          }

          if (!matched_spec_xlinks_alpha.empty())
          {
            for (auto match : matched_spec_xlinks_alpha)
            {
              iso_peaks_xlinks_alpha.push_back(num_iso_peaks_array[match.second]);
            }
            csm.num_iso_peaks_mean_xlinks_alpha = Math::mean(iso_peaks_xlinks_alpha.begin(), iso_peaks_xlinks_alpha.end());
          }

          if (!matched_spec_xlinks_beta.empty())
          {
            for (auto match : matched_spec_xlinks_beta)
            {
              iso_peaks_xlinks_beta.push_back(num_iso_peaks_array[match.second]);
            }
            csm.num_iso_peaks_mean_xlinks_beta = Math::mean(iso_peaks_xlinks_beta.begin(), iso_peaks_xlinks_beta.end());
          }
        }
        csm.ppm_error_abs_sum_linear_alpha = 0;
        csm.ppm_error_abs_sum_linear_beta = 0;
        csm.ppm_error_abs_sum_xlinks_alpha = 0;
        csm.ppm_error_abs_sum_xlinks_beta = 0;
        csm.ppm_error_abs_sum_linear = 0;
        csm.ppm_error_abs_sum_xlinks = 0;
        csm.ppm_error_abs_sum_alpha = 0;
        csm.ppm_error_abs_sum_beta = 0;
        csm.ppm_error_abs_sum = 0;

        LOG_DEBUG << "Computing ppm error summeries..." << endl;

        // TODO find a better way to compute the absolute sum
        if (ppm_error_array_linear_alpha.size() > 0)
        {
          for (Size k = 0; k < ppm_error_array_linear_alpha.size(); ++k)
          {
            csm.ppm_error_abs_sum_linear_alpha += abs(ppm_error_array_linear_alpha[k]);
          }
          csm.ppm_error_abs_sum_linear_alpha = csm.ppm_error_abs_sum_linear_alpha / ppm_error_array_linear_alpha.size();
        }

        if (ppm_error_array_linear_beta.size() > 0)
        {
          for (Size k = 0; k < ppm_error_array_linear_beta.size(); ++k)
          {
            csm.ppm_error_abs_sum_linear_beta += abs(ppm_error_array_linear_beta[k]);
          }
          csm.ppm_error_abs_sum_linear_beta = csm.ppm_error_abs_sum_linear_beta / ppm_error_array_linear_beta.size();
        }

        if (ppm_error_array_xlinks_alpha.size() > 0)
        {
          for (Size k = 0; k < ppm_error_array_xlinks_alpha.size(); ++k)
          {
            csm.ppm_error_abs_sum_xlinks_alpha += abs(ppm_error_array_xlinks_alpha[k]);
          }
          csm.ppm_error_abs_sum_xlinks_alpha = csm.ppm_error_abs_sum_xlinks_alpha / ppm_error_array_xlinks_alpha.size();
        }

        if (ppm_error_array_xlinks_beta.size() > 0)
        {
          for (Size k = 0; k < ppm_error_array_xlinks_beta.size(); ++k)
          {
            csm.ppm_error_abs_sum_xlinks_beta += abs(ppm_error_array_xlinks_beta[k]);
          }
          csm.ppm_error_abs_sum_xlinks_beta = csm.ppm_error_abs_sum_xlinks_beta / ppm_error_array_xlinks_beta.size();
        }

        DataArrays::FloatDataArray ppm_error_array_linear;
        DataArrays::FloatDataArray ppm_error_array_xlinks;
        DataArrays::FloatDataArray ppm_error_array_alpha;
        DataArrays::FloatDataArray ppm_error_array_beta;
        DataArrays::FloatDataArray ppm_error_array;
        ppm_error_array_linear.insert(ppm_error_array_linear.end(), ppm_error_array_linear_alpha.begin(), ppm_error_array_linear_alpha.end());
        ppm_error_array_linear.insert(ppm_error_array_linear.end(), ppm_error_array_linear_beta.begin(), ppm_error_array_linear_beta.end());
        ppm_error_array_xlinks.insert(ppm_error_array_xlinks.end(), ppm_error_array_xlinks_alpha.begin(), ppm_error_array_xlinks_alpha.end());
        ppm_error_array_xlinks.insert(ppm_error_array_xlinks.end(), ppm_error_array_xlinks_beta.begin(), ppm_error_array_xlinks_beta.end());
        ppm_error_array_alpha.insert(ppm_error_array_alpha.end(), ppm_error_array_linear_alpha.begin(), ppm_error_array_linear_alpha.end());
        ppm_error_array_alpha.insert(ppm_error_array_alpha.end(), ppm_error_array_xlinks_alpha.begin(), ppm_error_array_xlinks_alpha.end());
        ppm_error_array_beta.insert(ppm_error_array_beta.end(), ppm_error_array_linear_beta.begin(), ppm_error_array_linear_beta.end());
        ppm_error_array_beta.insert(ppm_error_array_beta.end(), ppm_error_array_xlinks_beta.begin(), ppm_error_array_xlinks_beta.end());
        ppm_error_array.insert(ppm_error_array.end(), ppm_error_array_linear.begin(), ppm_error_array_linear.end());
        ppm_error_array.insert(ppm_error_array.end(), ppm_error_array_xlinks.begin(), ppm_error_array_xlinks.end());

        if (ppm_error_array_linear.size() > 0)
        {
          for (double ppm_error : ppm_error_array_linear)
          {
            csm.ppm_error_abs_sum_linear += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_linear = csm.ppm_error_abs_sum_linear / ppm_error_array_linear.size();
        }

        if (ppm_error_array_xlinks.size() > 0)
        {
          for (double ppm_error : ppm_error_array_xlinks)
          {
            csm.ppm_error_abs_sum_xlinks += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_xlinks = csm.ppm_error_abs_sum_xlinks / ppm_error_array_xlinks.size();
        }

        if (ppm_error_array_alpha.size() > 0)
        {
          for (double ppm_error : ppm_error_array_alpha)
          {
            csm.ppm_error_abs_sum_alpha += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_alpha = csm.ppm_error_abs_sum_alpha / ppm_error_array_alpha.size();
        }

        if (ppm_error_array_beta.size() > 0)
        {
          for (double ppm_error : ppm_error_array_beta)
          {
            csm.ppm_error_abs_sum_beta += abs(ppm_error);
          }
          csm.ppm_error_abs_sum_beta = csm.ppm_error_abs_sum_beta / ppm_error_array_beta.size();
        }

        if (ppm_error_array.size() > 0)
        {
          for (double ppm_error : ppm_error_array)
          {
            csm.ppm_error_abs_sum += abs(ppm_error);
          }
          csm.ppm_error_abs_sum = csm.ppm_error_abs_sum / ppm_error_array.size();
        }

        // write fragment annotations
        LOG_DEBUG << "Start writing annotations" << endl;
        vector<PeptideHit::PeakAnnotation> frag_annotations;

        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_alpha, theoretical_spec_linear_alpha, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_linear_beta, theoretical_spec_linear_beta, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_alpha, theoretical_spec_xlinks_alpha, spectrum);
        OPXLHelper::buildFragmentAnnotations(frag_annotations, matched_spec_xlinks_beta, theoretical_spec_xlinks_beta, spectrum);
        LOG_DEBUG << "End writing annotations, size: " << frag_annotations.size() << endl;

        // make annotations unique
        sort(frag_annotations.begin(), frag_annotations.end());
        vector<PeptideHit::PeakAnnotation>::iterator last_unique_anno = unique(frag_annotations.begin(), frag_annotations.end());
        if (last_unique_anno != frag_annotations.end())
        {
          frag_annotations.erase(last_unique_anno, frag_annotations.end());
        }

        csm.frag_annotations = frag_annotations;

        all_csms_spectrum.push_back(csm);
      } // candidates for peak finished, determine best matching candidate

      // collect top n matches to spectrum
      sort(all_csms_spectrum.rbegin(), all_csms_spectrum.rend(), OPXLDataStructs::CLSMScoreComparator());
      Size max_hit = min(all_csms_spectrum.size(), static_cast<Size>(number_top_hits));

      for (Size top = 0; top < max_hit; top++)
      {
        all_csms_spectrum[top].rank = top+1;
        top_csms_spectrum.push_back(all_csms_spectrum[top]);
      }

      Size all_top_csms_current_index = 0;
#ifdef _OPENMP
#pragma omp critical (all_top_csms_access)
#endif
      {
        if (!top_csms_spectrum.empty())
        {
          all_top_csms.push_back(top_csms_spectrum);
          all_top_csms_current_index = all_top_csms.size()-1;
        }
      }

      // Write PeptideIdentifications and PeptideHits for n top hits
      if (!top_csms_spectrum.empty())
      {
        OPXLHelper::buildPeptideIDs(peptide_ids, top_csms_spectrum, all_top_csms, all_top_csms_current_index, spectra, scan_index, scan_index);
      }

      LOG_DEBUG << "Next Spectrum ##################################" << endl;
    }

    // end of matching / scoring
    progresslogger.endProgress();

    // Add protein identifications
    PeptideIndexing pep_indexing;
    Param indexing_param = pep_indexing.getParameters();

    String d_prefix = decoy_prefix ? "prefix" : "suffix";
    indexing_param.setValue("decoy_string_position", d_prefix, "If set, protein accessions in the database contain 'decoy_string' as prefix.");
    indexing_param.setValue("decoy_string", decoy_string, "String that was appended (or prefixed - see 'prefix' flag below) to the accessions in the protein database to indicate decoy proteins.");
    indexing_param.setValue("missing_decoy_action", "warn");
    indexing_param.setValue("enzyme:name", enzyme_name);
    pep_indexing.setParameters(indexing_param);

    pep_indexing.run(fasta_db, protein_ids, peptide_ids);

    OPXLHelper::addProteinPositionMetaValues(peptide_ids);

    // write output
    progresslogger.startProgress(0, 1, "Writing output...");
    if (out_idXML.size() > 0)
    {
      IdXMLFile().store(out_idXML, protein_ids, peptide_ids);
    }
    if (out_mzIdentML.size() > 0)
    {
      MzIdentMLFile().store(out_mzIdentML, protein_ids, peptide_ids);
    }

    if (out_xquest.size() > 0 || out_xquest_specxml.size() > 0)
    {
      vector<String> input_split_dir;
      vector<String> input_split;
      getStringOption_("in").split(String("/"), input_split_dir);
      input_split_dir[input_split_dir.size()-1].split(String("."), input_split);
      String base_name = input_split[0];

      if (out_xquest.size() > 0)
      {
        // String precursor_mass_tolerance_unit_string = precursor_mass_tolerance_unit_ppm ? "ppm" : "Da";
        // String fragment_mass_tolerance_unit_string = fragment_mass_tolerance_unit_ppm ? "ppm" : "Da";
        XQuestResultXMLFile().store(out_xquest, protein_ids, peptide_ids);
      }
      if (out_xquest_specxml.size() > 0)
      {
        XQuestResultXMLFile::writeXQuestXMLSpec(out_xquest_specxml, base_name, all_top_csms, spectra);
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
