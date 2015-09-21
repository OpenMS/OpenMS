// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/XTandemXMLFile.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/EnzymesDB.h>

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

    @brief Identifies peptides in MS/MS spectra via XTandem.

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

    @em X!Tandem must be installed before this wrapper can be used. This wrapper
    has been successfully tested with several versions of X!Tandem.
  The last known version to work is 2009-04-01. We encountered problems with
  later versions (namely 2010-01-01).

    To speed up computations, FASTA databases can be compressed using the fasta_pro.exe
    tool of @em X!Tandem. It is contained in the "bin" folder of the @em X!Tandem installation.
    Refer to the docu of @em X!Tandem for further information about settings.

  This adapter supports relative database filenames, which (when not found in the current working directory) is looked up in
  the directories specified by 'OpenMS.ini:id_db_dir' (see @subpage TOPP_advanced).

    The major part of the setting can be directly adjusted using the "default_input.xml" of
    @em X!Tandem. Parameters set by this wrapper overwrite the default settings given in the
    "default_input.xml", even those parameters not set explicitly, but defaulting to a value.
    An example of such a "default_input.xml" is contained in the "bin" folder of the
    @em X!Tandem installation. The parameter "default_input_file" must point to a valid
    file. "Masterfiles" for "default_input.xml" parameter importing other xml input files
    are not recommended, use at own risk.

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
    TOPPBase("XTandemAdapter", "Annotates MS/MS spectra using XTandem.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {

    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "Precursor mass tolerance", false);
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error", false);

    addEmptyLine_();
    registerStringOption_("precursor_error_units", "<unit>", "ppm", "Parent monoisotopic mass error units", false);
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    vector<String> valid_strings;
    valid_strings.push_back("ppm");
    valid_strings.push_back("Da");
    setValidStrings_("precursor_error_units", valid_strings);
    setValidStrings_("fragment_error_units", valid_strings);
    registerIntOption_("min_precursor_charge", "<charge>", 1, "Minimum precursor charge", false);
    registerIntOption_("max_precursor_charge", "<charge>", 4, "Maximum precursor charge", false);

    registerStringOption_("allow_isotope_error", "<error>", "yes", "If set, misassignment to the first and second isotopic 13C peak are also considered.", false);
    valid_strings.clear();
    valid_strings.push_back("yes");
    valid_strings.push_back("no");
    setValidStrings_("allow_isotope_error", valid_strings);

    registerStringList_("fixed_modifications", "<mods>", ListUtils::create<String>(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", ListUtils::create<String>(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    setValidStrings_("variable_modifications", all_mods);
    registerIntOption_("missed_cleavages", "<num>", 1, "Number of possible cleavage sites missed by the enzyme", false);

    addEmptyLine_();
    registerInputFile_("xtandem_executable", "<executable>",
// choose the default value according to the platform where it will be executed
// xtandem compiles as tandem on osx and tandem.exe on any other platform
#if  defined(__APPLE__)
                       "tandem",
#else
                       "tandem.exe",
#endif
                       "X!Tandem executable of the installation e.g. 'tandem.exe'", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("default_input_file", "<file>", "", "Default parameters input file, if not given default parameters are used", false);
    registerDoubleOption_("minimum_fragment_mz", "<num>", 150.0, "Minimum fragment mz", false);
    vector<String> all_enzymes;
    EnzymesDB::getInstance()->getAllXTandemNames(all_enzymes);
    registerStringOption_("cleavage_site", "<cleavage site>", "Trypsin", "The enzyme used for peptide digestion.", false);
    setValidStrings_("cleavage_site", all_enzymes);
    registerStringOption_("output_results", "<result reporting>", "all", "Which hits should be reported. All, valid ones (passing the E-Ealue threshold), or stochastic (failing the threshold)", false);
    valid_strings.clear();
    valid_strings.push_back("all");
    valid_strings.push_back("valid");
    valid_strings.push_back("stochastic");
    setValidStrings_("output_results", valid_strings);

    registerDoubleOption_("max_valid_expect", "<E-Value>", 0.1, "Maximal E-Value of a hit to be reported (only evaluated if 'output_result' is 'valid' or 'stochastic'", false);
    registerFlag_("refinement", "Enable the refinement. For most applications (especially when using FDR, PEP approaches) it is NOT recommended to set this flag.");
    registerFlag_("use_noise_suppression", "Enable the use of the noise suppression routines.");
    registerFlag_("semi_cleavage", "If set, both termini must NOT follow the cutting rule. For most applications it is NOT recommended to set this flag.");
  }

  ExitCodes main_(int, const char**)
  {
    // instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
    String ini_location;
    // path to the log file
    String logfile(getStringOption_("log"));
    String xtandem_executable(getStringOption_("xtandem_executable"));
    String inputfile_name;
    String outputfile_name;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name == "")
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    outputfile_name = getStringOption_("out");
    writeDebug_(String("Output file: ") + outputfile_name, 1);
    if (outputfile_name == "")
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    // write input xml file
    String temp_directory = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }

    String input_filename(temp_directory + "_tandem_input_file.xml");
    String tandem_input_filename(temp_directory + "_tandem_input_file.mzData");
    String tandem_output_filename(temp_directory + "_tandem_output_file.xml");
    String tandem_taxonomy_filename(temp_directory + "_tandem_taxonomy_file.xml");

    //-------------------------------------------------------------
    // Validate user parameters
    //-------------------------------------------------------------
    if (getIntOption_("min_precursor_charge") > getIntOption_("max_precursor_charge"))
    {
      LOG_ERROR << "Given charge range is invalid: max_precursor_charge needs to be >= min_precursor_charge." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

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
    mzml_file.getOptions().addMSLevel(2); // only load msLevel 2
    mzml_file.setLogType(log_type_);
    mzml_file.load(inputfile_name, exp);

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
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    // we need to replace the native id with a simple numbering schema, to be able to
    // map the IDs back to the spectra (RT, and MZ information)
    Size native_id(0);
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      it->setNativeID(++native_id);
    }

    // We store the file in mzData file format, because MGF files somehow produce in most
    // of the cases IDs with charge 2+. We do not use the input file directly
    // because XTandem sometimes stumbles over misleading substrings in the filename,
    // e.g. mzXML ...
    MzDataFile mzdata_outfile;
    mzdata_outfile.store(tandem_input_filename, exp);

    XTandemInfile infile;
    infile.setInputFilename(tandem_input_filename);
    infile.setOutputFilename(tandem_output_filename);

    ofstream tax_out(tandem_taxonomy_filename.c_str());
    tax_out << "<?xml version=\"1.0\"?>" << "\n";
    tax_out << "\t<bioml label=\"x! taxon-to-file matching list\">" << "\n";
    tax_out << "\t\t<taxon label=\"OpenMS_dummy_taxonomy\">" << "\n";
    tax_out << "\t\t\t<file format=\"peptide\" URL=\"" << db_name << "\" />" << "\n";
    tax_out << "\t</taxon>" << "\n";
    tax_out << "</bioml>" << "\n";
    tax_out.close();

    infile.setTaxonomyFilename(tandem_taxonomy_filename);

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

    if (getStringOption_("default_input_file") != "")
    {
      infile.load(getStringOption_("default_input_file"));
      infile.setDefaultParametersFilename(getStringOption_("default_input_file"));
    }
    else
    {
      String default_file = File::find("CHEMISTRY/XTandem_default_input.xml");
      infile.load(default_file);
      infile.setDefaultParametersFilename(default_file);
    }

    infile.setPrecursorMassTolerancePlus(getDoubleOption_("precursor_mass_tolerance"));
    infile.setPrecursorMassToleranceMinus(getDoubleOption_("precursor_mass_tolerance"));
    infile.setFragmentMassTolerance(getDoubleOption_("fragment_mass_tolerance"));
    infile.setMaxPrecursorCharge(getIntOption_("max_precursor_charge"));
    infile.setNumberOfThreads(getIntOption_("threads"));
    infile.setModifications(ModificationDefinitionsSet(getStringList_("fixed_modifications"), getStringList_("variable_modifications")));
    infile.setTaxon("OpenMS_dummy_taxonomy");
    infile.setOutputResults(getStringOption_("output_results"));
    infile.setMaxValidEValue(getDoubleOption_("max_valid_expect"));
    String enzyme_name = getStringOption_("cleavage_site");
    infile.setCleavageSite(EnzymesDB::getInstance()->getEnzyme(enzyme_name)->getXTANDEMid());
    infile.setNumberOfMissedCleavages(getIntOption_("missed_cleavages"));
    infile.setRefine(getFlag_("refinement"));
    infile.setNoiseSuppression(getFlag_("use_noise_suppression"));
    infile.setSemiCleavage(getFlag_("semi_cleavage"));
    bool allow_isotope_error = getStringOption_("allow_isotope_error") == "yes" ? true : false;
    infile.setAllowIsotopeError(allow_isotope_error);

    infile.write(input_filename);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    int status = QProcess::execute(xtandem_executable.toQString(), QStringList(input_filename.toQString())); // does automatic escaping etc...
    if (status != 0)
    {
      writeLog_("XTandem problem. Aborting! Calling command was: '" + xtandem_executable + " \"" + input_filename + "\"'.\nDoes the !XTandem executable exist?");
      // clean temporary files
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
      }
      return EXTERNAL_PROGRAM_ERROR;
    }

    vector<ProteinIdentification> protein_ids;
    ProteinIdentification protein_id;
    vector<PeptideIdentification> peptide_ids;

    // read the output of X!Tandem and write it to idXML
    XTandemXMLFile tandem_output;
    tandem_output.setModificationDefinitionsSet(ModificationDefinitionsSet(getStringList_("fixed_modifications"), getStringList_("variable_modifications")));
    // find the file, because XTandem extends the filename with a timestamp we do not know (exactly)
    StringList files;
    File::fileList(temp_directory, "_tandem_output_file*.xml", files);
    if (files.size() != 1)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, tandem_output_filename);
    }
    tandem_output.load(temp_directory + files[0], protein_id, peptide_ids);

    // now put the RTs into the peptide_ids from the spectrum ids
    for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
    {
      UInt id = (Int)it->getMetaValue("spectrum_id");
      --id; // native IDs were written 1-based
      if (id < exp.size())
      {
        it->setRT(exp[id].getRT());
        double pre_mz(0.0);
        if (!exp[id].getPrecursors().empty()) pre_mz = exp[id].getPrecursors()[0].getMZ();
        it->setMZ(pre_mz);
        //it->removeMetaValue("spectrum_id");
      }
      else
      {
        LOG_ERROR << "XTandemAdapter: Error: id '" << id << "' not found in peak map!" << endl;
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // handle the search parameters
    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.db = getStringOption_("database");
    search_parameters.charges = "+" + String(getIntOption_("min_precursor_charge")) + "-+" + String(getIntOption_("max_precursor_charge"));

    ProteinIdentification::PeakMassType mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.mass_type = mass_type;
    search_parameters.fixed_modifications = getStringList_("fixed_modifications");
    search_parameters.variable_modifications = getStringList_("variable_modifications");
    search_parameters.missed_cleavages = getIntOption_("missed_cleavages");
    search_parameters.peak_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
    search_parameters.precursor_tolerance = getDoubleOption_("precursor_mass_tolerance");
    search_parameters.digestion_enzyme = *EnzymesDB::getInstance()->getEnzyme(enzyme_name);
    protein_id.setSearchParameters(search_parameters);
    protein_id.setSearchEngineVersion("");
    protein_id.setSearchEngine("XTandem");

    protein_ids.push_back(protein_id);

    IdXMLFile().store(outputfile_name, protein_ids, peptide_ids);

    /// Deletion of temporary files
    if (this->debug_level_ < 2)
    {
      File::removeDirRecursively(temp_directory);
      LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
    }

    // some stats
    LOG_INFO << "Statistics:\n"
             << "  identified MS2 spectra: " << peptide_ids.size() << " / " << exp.size() << " = " << int(peptide_ids.size() * 100.0 / exp.size()) << "% (with e-value < " << String(getDoubleOption_("max_valid_expect")) << ")" << std::endl;

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPXTandemAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
