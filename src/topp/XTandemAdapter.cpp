// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_XTandemAdapter XTandemAdapter

    @brief Identifies peptides in MS/MS spectra via the search engine X! Tandem.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ XTandemAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em X! Tandem must be installed before this wrapper can be used.
    This wrapper has been successfully tested with several versions of X! Tandem.
    The earliest version known to work is "PILEDRIVER" (2015-04-01). The latest is "ALANINE" (2017-02-01).

    To speed up computations, FASTA databases can be compressed using the fasta_pro.exe tool of @em X! Tandem.
    It is contained in the "bin" folder of the @em X! Tandem installation.
    Refer to the documentation of @em X! Tandem for further information about settings.

    This adapter supports relative database filenames.
    If a database is not found in the current working directory, it is looked up in the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    @em X! Tandem settings not exposed by this adapter (especially refinement settings) can be directly adjusted using an XML configuration file.
    By default, all (!) parameters available explicitly via this wrapper take precedence over the XML configuration file.
    The parameter @p default_config_file can be used to specify such a custom configuration.
    An example of a configuration file (named "default_input.xml") is contained in the "bin" folder of the @em X! Tandem installation and in the OpenMS installation under OpenMS/share/CHEMISTRY/XTandem_default_input.xml.
    If you want to use the XML configuration file and @em ignore most of the parameters set via this adapter, use the @p ignore_adapter_param flag.
    Then, the config given via @p default_config_file is used exclusively and only the values for the paramters @p in, @p out, @p database and @p xtandem_executable are taken from this adapter.

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_XTandemAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_XTandemAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPXTandemAdapter :
  public TOPPBase
{
public:
  TOPPXTandemAdapter() :
    TOPPBase("XTandemAdapter", "Annotates MS/MS spectra using X! Tandem.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {

    registerInputFile_("in", "<file>", "", "Input file containing MS2 spectra");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file containing search results", false);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("xml_out", "<file>", "", "Raw output file directly from X! Tandem. Either 'out' or 'xml_out' are required. They can be used together.", false);
    setValidFormats_("xml_out", ListUtils::create<String>("xml"));
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    registerInputFile_("xtandem_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      // X! Tandem compiles as tandem on OSX and tandem.exe on any other platform
#if  defined(__APPLE__)
      "tandem",
#else
      "tandem.exe",
#endif
      "X! Tandem executable of the installation e.g. 'tandem.exe'", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("default_config_file", "<file>", "", "Default X! Tandem configuration file. All parameters of this adapter take precedence over the file - use it for parameters not available here. A template file can be found at 'OpenMS/share/CHEMISTRY/XTandem_default_input.xml'.", false, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("default_config_file", ListUtils::create<String>("xml"));
    registerFlag_("ignore_adapter_param", "Set this to use the configuration given in 'default_config_file' exclusively, ignoring other parameters (apart from 'in', 'out', 'database', 'xtandem_executable') set via this adapter.");

    addEmptyLine_();
    //
    // Optional parameters (if '-ignore_adapter_param' is set)
    //
    registerDoubleOption_("precursor_mass_tolerance", "<value>", 10.0, "Precursor mass tolerance", false);
    registerDoubleOption_("fragment_mass_tolerance", "<value>", 0.3, "Fragment mass error", false);

    registerStringOption_("precursor_error_units", "<unit>", "ppm", "Parent monoisotopic mass error units", false);
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    vector<String> valid_strings = ListUtils::create<String>("ppm,Da");
    setValidStrings_("precursor_error_units", valid_strings);
    setValidStrings_("fragment_error_units", valid_strings);

    registerIntOption_("max_precursor_charge", "<number>", 4, "Maximum precursor charge ('0' to use X! Tandem default)", false);
    setMinInt_("max_precursor_charge", 0);

    registerFlag_("no_isotope_error", "By default, misassignment to the first and second isotopic 13C peak are also considered. Set this flag to disable.", false);

    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using Unimod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);

    registerDoubleOption_("minimum_fragment_mz", "<number>", 150.0, "Minimum fragment mz", false);

    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllXTandemNames(all_enzymes);
    registerStringOption_("enzyme", "<choice>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("enzyme", all_enzymes);
    registerIntOption_("missed_cleavages", "<number>", 1, "Number of possible cleavage sites missed by the enzyme", false);
    registerFlag_("semi_cleavage", "Require only peptide end to have a valid cleavage site, not both.");

    registerStringOption_("output_results", "<choice>", "all", "Which hits should be reported. All, valid ones (passing the E-Value threshold), or stochastic (failing the threshold)", false);
    valid_strings = ListUtils::create<String>("all,valid,stochastic", ',');
    setValidStrings_("output_results", valid_strings);

    registerDoubleOption_("max_valid_expect", "<value>", 0.1, "Maximal E-Value of a hit to be reported (only evaluated if 'output_result' is 'valid' or 'stochastic')", false);
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String xml_out = getStringOption_("xml_out");
    if (xml_out.empty() && out.empty())
    {
      writeLog_("Fatal error: no output file given (parameter 'out' or 'xml_out')");
      return ILLEGAL_PARAMETERS;
    }

    // write input xml file
    String temp_directory = makeTempDirectory_();
    String input_filename = temp_directory + "tandem_input.xml";
    String tandem_input_filename = in;
    String tandem_output_filename = temp_directory + "tandem_output.xml";
    String tandem_taxonomy_filename = temp_directory + "tandem_taxonomy.xml";

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    String db_name(getStringOption_("database"));
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

    PeakMap exp;
    MzMLFile mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load MS level 2
    mzml_file.setLogType(log_type_);
    mzml_file.load(in, exp);

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS2 spectra in input file.");
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();

    if (spectrum_type == SpectrumSettings::RAWDATA)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the 'force' flag.");
      }
    }

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

    infile.write(input_filename, getFlag_("ignore_adapter_param"),
                 getFlag_("force"));

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    String xtandem_executable(getStringOption_("xtandem_executable"));
    int status = QProcess::execute(xtandem_executable.toQString(), QStringList(input_filename.toQString())); // does automatic escaping etc...
    if (status != 0)
    {
      writeLog_("X! Tandem problem. Aborting! Calling command was: '" + xtandem_executable + " \"" + input_filename + "\"'.\nDoes the X! Tandem executable exist?");
      // clean temporary files
      removeTempDirectory_(temp_directory);
      return EXTERNAL_PROGRAM_ERROR;
    }

    vector<ProteinIdentification> protein_ids;
    ProteinIdentification protein_id;
    StringList ms_runs;
    exp.getPrimaryMSRunPath(ms_runs);
    protein_id.setPrimaryMSRunPath(ms_runs);
    vector<PeptideIdentification> peptide_ids;

    // read the output of X! Tandem and write it to idXML
    XTandemXMLFile tandem_output;
    ModificationDefinitionsSet mod_def_set(getStringList_("fixed_modifications"), getStringList_("variable_modifications"));
    tandem_output.load(tandem_output_filename, protein_id, peptide_ids, mod_def_set);

    // add RT and precursor m/z to the peptide IDs (look them up in the spectra):
    SpectrumLookup lookup;
    lookup.readSpectra(exp);

    for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
    {
      String ref = it->getMetaValue("spectrum_reference");
      Size index = lookup.findByNativeID(ref);
      if (index < exp.size())
      {
        it->setRT(exp[index].getRT());
        if (!exp[index].getPrecursors().empty())
        {
          it->setMZ(exp[index].getPrecursors()[0].getMZ());
        }
      }
      else
      {
        LOG_ERROR << "Error: spectrum with ID '" << ref << "' not found in input data! RT and precursor m/z values could not be looked up." << endl;
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if (!xml_out.empty())
    {
      // existing file? Qt won't overwrite, so try to remove it:
      if (QFile::exists(xml_out.toQString()) &&
          !QFile::remove(xml_out.toQString()))
      {
        writeLog_("Fatal error: Could not overwrite existing file '" + xml_out + "'");
        return CANNOT_WRITE_OUTPUT_FILE;
      }
      // move the temporary file to the actual destination:
      if (!QFile::rename(tandem_output_filename.toQString(), xml_out.toQString()))
      {
        writeLog_("Fatal error: Could not move temporary X! Tandem XML file to '" + xml_out + "'");
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

      protein_ids.push_back(protein_id);

      IdXMLFile().store(out, protein_ids, peptide_ids);
    }

    /// Deletion of temporary files
    removeTempDirectory_(temp_directory);

    // some stats (note that only MS2 spectra were loaded into "exp"):
    Int percent = peptide_ids.size() * 100.0 / exp.size();
    LOG_INFO << "Statistics:\n"
             << "- identified MS2 spectra: " << peptide_ids.size() << " / "
             << exp.size() << " = " << percent << "%";
    if (output_results != "all")
    {
      LOG_INFO << " (with E-value " << (output_results == "valid" ? "< " : "> ")
               << String(max_evalue) << ")";
    }
    LOG_INFO << std::endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPXTandemAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
