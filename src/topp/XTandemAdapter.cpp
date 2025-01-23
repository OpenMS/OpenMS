// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/SearchEngineBase.h>

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/ANALYSIS/ID/PeptideIndexing.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>

#include <fstream>
#include <QStringList>


using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_XTandemAdapter XTandemAdapter

@brief Identifies peptides in MS/MS spectra via the search engine X! Tandem.

@note This adapter is deprecated (and likely removed in OpenMS 3.3). Contact us, if you feel it's irreplaceable.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; XTandemAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em X!Tandem must be installed before this wrapper can be used.
    This wrapper has been successfully tested with several versions of @em X!Tandem.
    The earliest version known to work is "PILEDRIVER" (2015-04-01). The latest is "ALANINE" (2017-02-01).

    @note @em X!Tandem only support <b>uncompressed mzML files</b> (e.g. no zlib compression or other fancy things like numpress) may be used internally!
    This converter only forwards the mzML filename and you will get an error like 'Fatal error: unsupported CODEC used for mzML peak data (CODEC type=zlib compression)'.
    If this happens, preprocess the mzML files using OpenMS' @ref TOPP_FileConverter to write a plain mzML which @em X!Tandem understands.

    @em X!Tandem has a build-in adventitious cleavage rule for Asp|Pro (Aspartate/D | Proline/P), which it allows as cutting site for all enzymes.
    Furthermore, it treats any occurence of 'X' as stop codon (and thus as cleavage site). The resulting peptide will be non- or semi-tryptic.

    To speed up computations, FASTA databases can be compressed using the fasta_pro.exe tool of @em X!Tandem.
    It is contained in the "bin" folder of the @em X!Tandem installation.
    Refer to the documentation of @em X!Tandem for further information about settings.

    This adapter supports relative database filenames.
    If a database is not found in the current working directory, it is looked up in the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    @em X!Tandem settings not exposed by this adapter (especially refinement settings) can be directly adjusted using an XML configuration file.
    By default, all (!) parameters available explicitly via this wrapper take precedence over the XML configuration file.
    The parameter @p default_config_file can be used to specify such a custom configuration.
    An example of a configuration file (named "default_input.xml") is contained in the "bin" folder of the @em X!Tandem installation and in the %OpenMS installation
    under <code>OpenMS/share/CHEMISTRY/XTandem_default_config.xml</code>.
    If you want to use the XML configuration file and @em ignore most of the parameters set via this adapter, use the @p ignore_adapter_param flag.
    Then, the config given via @p default_config_file is used exclusively and only the values for the parameters @p in, @p out, @p database and @p xtandem_executable are taken from this adapter.

    @note This adapter supports <b>15N labeling</b> by using the <tt>XTandem_residue_mass.bioml.xml</tt> file (which defines modified AA masses) as provided in
          <code>OpenMS/share/OpenMS/CHEMISTRY/</code>. To use it, specify the full path (which will depend on your system!) to this bioml.xml file 
          within the <tt>XTandem_default_config.xml</tt> config file (see above).
    Within this config file, modify the path in the following line to match your system's configuration.
@code
<note type="input" label="protein, modified residue mass file">/path/to/XTandem_residue_mass.bioml.xml</note>
@endcode
    and pass the config file's filename via the <tt>default_config_file</tt> parameter to XTandemAdapter.
    For more details, see https://www.thegpm.org/TANDEM/api/pmrmf.html.
    <br>Warning: If the path to XTandem_residue_mass.bioml.xml is invalid, @em X!Tandem will simply ignore the setting without feedback!
    <br>Warning: The resulting peptide sequences in the idXML file will not contain any N15 labeling information because 
                 X!Tandem simply received modified AA masses without further information what they mean. 
                 Add the 15N modification information by calling the @ref TOPP_StaticModification tool on the idXML file created by this adapter.

    <br>

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.


    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_XTandemAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_XTandemAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPXTandemAdapter :
  public SearchEngineBase
{
public:
  TOPPXTandemAdapter() :
    SearchEngineBase("XTandemAdapter", "Annotates MS/MS spectra using X! Tandem.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {

    registerInputFile_("in", "<file>", "", "Input file containing MS2 spectra");
    setValidFormats_("in", {"mzML"});
    registerOutputFile_("out", "<file>", "", "Output file containing search results", false);
    setValidFormats_("out", {"idXML"});
    registerOutputFile_("xml_out", "<file>", "", "Raw output file directly from X! Tandem. Either 'out' or 'xml_out' are required. They can be used together.", false);
    setValidFormats_("xml_out", {"xml"});
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, {"skipexists"});
    setValidFormats_("database", {"FASTA"});
    registerInputFile_("xtandem_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      // X! Tandem compiles as tandem on OSX and tandem.exe on any other platform
#if  defined(__APPLE__)
      "tandem",
#else
      "tandem.exe",
#endif
      "X! Tandem executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", true, false, {"is_executable"});
    registerInputFile_("default_config_file", "<file>", "", "Default X! Tandem configuration file. All parameters of this adapter take precedence over the file - use it for parameters not available here. A template file can be found at 'OpenMS/share/CHEMISTRY/XTandem_default_config.xml'.", false, false, {"skipexists"});
    setValidFormats_("default_config_file", {"xml"});
    registerFlag_("ignore_adapter_param", "Set this to use the configuration given in 'default_config_file' exclusively, ignoring other parameters (apart from 'in', 'out', 'database', 'xtandem_executable') set via this adapter.");

    addEmptyLine_();
    //
    // Optional parameters (if '-ignore_adapter_param' is set)
    //
    registerDoubleOption_("precursor_mass_tolerance", "<value>", 10.0, "Precursor mass tolerance", false);
    registerDoubleOption_("fragment_mass_tolerance", "<value>", 0.3, "Fragment mass error", false);

    registerStringOption_("precursor_error_units", "<unit>", "ppm", "Parent monoisotopic mass error units", false);
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    const vector<String> valid_strings = {"ppm", "Da"};
    setValidStrings_("precursor_error_units", valid_strings);
    setValidStrings_("fragment_error_units", valid_strings);

    registerIntOption_("max_precursor_charge", "<number>", 4, "Maximum precursor charge ('0' to use X! Tandem default)", false);
    setMinInt_("max_precursor_charge", 0);

    registerFlag_("no_isotope_error", "By default, misassignment to the first and second isotopic 13C peak are also considered. Set this flag to disable.", false);

    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("fixed_modifications", "<mods>", {"Carbamidomethyl (C)"}, "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", {"Oxidation (M)"}, "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    registerDoubleOption_("minimum_fragment_mz", "<number>", 150.0, "Minimum fragment m/z", false);

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllXTandemNames(all_enzymes);
    registerStringOption_("enzyme", "<choice>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("enzyme", all_enzymes);
    registerIntOption_("missed_cleavages", "<number>", 1, "Number of possible cleavage sites missed by the enzyme", false);
    registerFlag_("semi_cleavage", "Require only peptide end to have a valid cleavage site, not both.");

    registerStringOption_("output_results", "<choice>", "all", "Which hits should be reported. All, valid ones (passing the E-Value threshold), or stochastic (failing the threshold)", false);
    setValidStrings_("output_results", { "all", "valid", "stochastic" });

    registerDoubleOption_("max_valid_expect", "<value>", 0.1, "Maximal E-Value of a hit to be reported (only evaluated if 'output_result' is 'valid' or 'stochastic')", false);

    // register peptide indexing parameter (with defaults for this search engine) TODO: check if search engine defaults are needed
    registerPeptideIndexingParameter_(PeptideIndexing().getParameters()); 
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in = getRawfileName();
    String out = getStringOption_("out");
    String xml_out = getStringOption_("xml_out");
    if (xml_out.empty() && out.empty())
    {
      writeLogError_("Fatal error: no output file given (parameter 'out' or 'xml_out')");
      return ILLEGAL_PARAMETERS;
    }

    // write input xml file
    File::TempDir dir(debug_level_ >= 2);
    String input_filename = dir.getPath() + "tandem_input.xml";
    String tandem_input_filename = in;
    String tandem_output_filename = dir.getPath() + "tandem_output.xml";
    String tandem_taxonomy_filename = dir.getPath() + "tandem_taxonomy.xml";

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    String db_name = getDBFilename();

    /// check if X!Tandem is available (warn early, since loading/storing of mzML below will delay the error -- which is not user friendly)
    String xtandem_executable = getStringOption_("xtandem_executable");

    PeakMap exp;
    FileHandler mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load MS level 2
    mzml_file.getOptions().setFillData(false); // do not fill the actual spectra. We only need RT and mz info for mapping
    mzml_file.loadExperiment(in, exp, {FileTypes::MZML});

    ofstream tax_out(tandem_taxonomy_filename.c_str());
    tax_out << "<?xml version=\"1.0\"?>" << "\n";
    tax_out << "\t<bioml label=\"x! taxon-to-file matching list\">" << "\n";
    tax_out << "\t\t<taxon label=\"OpenMS_dummy_taxonomy\">" << "\n";
    tax_out << "\t\t\t<file format=\"peptide\" URL=\"" << db_name << "\" />" << "\n";
    tax_out << "\t</taxon>" << "\n";
    tax_out << "</bioml>" << "\n";
    tax_out.close();

    //
    //  Prepare the XML configuration file
    //
    XTandemInfile infile;
    infile.setInputFilename(tandem_input_filename);
    infile.setOutputFilename(tandem_output_filename);
    infile.setTaxonomyFilename(tandem_taxonomy_filename); // contains the FASTA name

    if (getStringOption_("precursor_error_units") == "Da")
    {
      infile.setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
    }
    else
    {
      infile.setPrecursorMassErrorUnit(XTandemInfile::PPM);
    }

    if (getStringOption_("fragment_error_units") == "Da")
    {
      infile.setFragmentMassErrorUnit(XTandemInfile::DALTONS);
    }
    else
    {
      infile.setFragmentMassErrorUnit(XTandemInfile::PPM);
    }

    double precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
    infile.setPrecursorMassTolerancePlus(precursor_mass_tolerance);
    infile.setPrecursorMassToleranceMinus(precursor_mass_tolerance);
    infile.setFragmentMassTolerance(getDoubleOption_("fragment_mass_tolerance"));
    infile.setMaxPrecursorCharge(getIntOption_("max_precursor_charge"));
    infile.setNumberOfThreads(getIntOption_("threads"));
    infile.setModifications(ModificationDefinitionsSet(getStringList_("fixed_modifications"), getStringList_("variable_modifications")));
    infile.setTaxon("OpenMS_dummy_taxonomy");
    String output_results = getStringOption_("output_results");
    infile.setOutputResults(output_results);
    double max_evalue = getDoubleOption_("max_valid_expect");
    infile.setMaxValidEValue(max_evalue);
    String enzyme_name = getStringOption_("enzyme");
    infile.setCleavageSite(ProteaseDB::getInstance()->getEnzyme(enzyme_name)->getXTandemID());
    infile.setNumberOfMissedCleavages(getIntOption_("missed_cleavages"));
    infile.setSemiCleavage(getFlag_("semi_cleavage"));
    infile.setAllowIsotopeError(!getFlag_("no_isotope_error"));

    String default_XML_config = getStringOption_("default_config_file");
    if (!default_XML_config.empty())
    {
      // augment with absolute path. If absolute filename is already given, this is a no-op.
      default_XML_config = File::find(default_XML_config);
      infile.setDefaultParametersFilename(default_XML_config);
    }

    infile.write(input_filename, getFlag_("ignore_adapter_param"), getFlag_("force"));

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    TOPPBase::ExitCodes exit_code = runExternalProcess_(xtandem_executable.toQString(), QStringList(input_filename.toQString())); // does automatic escaping etc...
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    vector<ProteinIdentification> protein_ids;
    ProteinIdentification protein_id;
    protein_id.setPrimaryMSRunPath({in}, exp);
    vector<PeptideIdentification> peptide_ids;

    // read the output of X! Tandem and write it to idXML
    XTandemXMLFile tandem_output;
    ModificationDefinitionsSet mod_def_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
    tandem_output.load(tandem_output_filename, protein_id, peptide_ids, mod_def_set);

    // add RT and precursor m/z to the peptide IDs (look them up in the spectra):
    SpectrumLookup lookup;
    lookup.readSpectra(exp);

    for (PeptideIdentification& pep : peptide_ids)
    {
      String ref = pep.getSpectrumReference();
      Size index = lookup.findByNativeID(ref);
      if (index < exp.size())
      {
        pep.setRT(exp[index].getRT());
        if (!exp[index].getPrecursors().empty())
        {
          pep.setMZ(exp[index].getPrecursors()[0].getMZ());
        }
      }
      else
      {
        OPENMS_LOG_ERROR << "Error: spectrum with ID '" << ref << "' not found in input data! RT and precursor m/z values could not be looked up." << endl;
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if (!xml_out.empty())
    { // move the temporary file to the actual destination:
      if (!File::rename(tandem_output_filename, xml_out))
      {
        return CANNOT_WRITE_OUTPUT_FILE;
      }
    }
    
    if (!out.empty())
    {
      // handle the search parameters
      ProteinIdentification::SearchParameters search_parameters;
      search_parameters.db = getStringOption_("database");

      ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
      search_parameters.mass_type = mass_type;
      set<String> mods = mod_def_set.getFixedModificationNames();
      search_parameters.fixed_modifications.reserve(mods.size());
      search_parameters.fixed_modifications.insert(
        search_parameters.fixed_modifications.end(), mods.begin(), mods.end());
      mods = mod_def_set.getVariableModificationNames();
      search_parameters.variable_modifications.reserve(mods.size());
      search_parameters.variable_modifications.insert(
        search_parameters.variable_modifications.end(), mods.begin(),
        mods.end());
      search_parameters.missed_cleavages = getIntOption_("missed_cleavages");
      search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
      search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
      search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_error_units") == "ppm" ? true : false;
      search_parameters.fragment_mass_tolerance_ppm = getStringOption_("fragment_error_units") == "ppm" ? true : false;
      search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme_name));
      protein_id.setSearchParameters(search_parameters);
      protein_id.setSearchEngineVersion("");
      protein_id.setSearchEngine("XTandem");

      // write all (!) parameters as metavalues to the search parameters
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_id.getSearchParameters(), this->getToolPrefix());

      protein_ids.push_back(protein_id);

    // if "reindex" parameter is set to true will perform reindexing
      if (auto ret = reindex_(protein_ids, peptide_ids); ret != EXECUTION_OK) return ret;

      StringList feature_set;
      PercolatorFeatureSetHelper::addXTANDEMFeatures(peptide_ids, feature_set);
      protein_ids.front().getSearchParameters().setMetaValue("extra_features", ListUtils::concatenate(feature_set, ","));

      FileHandler().storeIdentifications(out, protein_ids, peptide_ids, {FileTypes::IDXML});
    }

    // some stats (note that only MS2 spectra were loaded into "exp"):
    Int percent = peptide_ids.size() * 100.0 / exp.size();
    OPENMS_LOG_INFO << "Statistics:\n"
             << "- identified MS2 spectra: " << peptide_ids.size() << " / "
             << exp.size() << " = " << percent << "%";
    if (output_results != "all")
    {
      OPENMS_LOG_INFO << " (with E-value " << (output_results == "valid" ? "< " : "> ")
               << String(max_evalue) << ")";
    }
    OPENMS_LOG_INFO << std::endl;

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPXTandemAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
