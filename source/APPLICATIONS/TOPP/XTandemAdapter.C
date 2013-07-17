// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

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
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", StringList::create("idXML"));
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 1.5, "Precursor mass tolerance", false);
    registerDoubleOption_("fragment_mass_tolerance", "<tolerance>", 0.3, "Fragment mass error", false);

    addEmptyLine_();
    registerStringOption_("precursor_error_units", "<unit>", "ppm", "Parent monoisotopic mass error units", false);
    registerStringOption_("fragment_error_units", "<unit>", "Da", "Fragment monoisotopic mass error units", false);
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", true, false, StringList::create("skipexists"));
    setValidFormats_("database", StringList::create("FASTA"));
    vector<String> valid_strings;
    valid_strings.push_back("ppm");
    valid_strings.push_back("Da");
    setValidStrings_("precursor_error_units", valid_strings);
    setValidStrings_("fragment_error_units", valid_strings);
    registerIntOption_("min_precursor_charge", "<charge>", 1, "Minimum precursor charge", false);
    registerIntOption_("max_precursor_charge", "<charge>", 4, "Maximum precursor charge", false);

    registerStringList_("fixed_modifications", "<mods>", StringList::create(""), "Fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    setValidStrings_("fixed_modifications", all_mods);
    registerStringList_("variable_modifications", "<mods>", StringList::create(""), "Variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
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
                       "X!Tandem executable of the installation e.g. 'tandem.exe'", true, false, StringList::create("skipexists"));
    registerInputFile_("default_input_file", "<file>", "", "Default parameters input file, if not given default parameters are used", false);
    registerDoubleOption_("minimum_fragment_mz", "<num>", 150.0, "Minimum fragment mz", false);
    registerStringOption_("cleavage_site", "<cleavage site>", "[RK]|{P}", "Cleavage site of the used enzyme as regular expression ([RK]|{P} (i.e. tryptic clevage) is default, [X]|[X] (i.e. every site) would be best for peptide input or unspecific digestion).", false);
    registerDoubleOption_("max_valid_expect", "<E-Value>", 0.1, "Maximal E-Value of a hit to be reported", false);
    registerFlag_("refinement", "Enable the refinement. For most applications (especially when using FDR, PEP approaches) it is NOT recommended to set this flag.");
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
    mzml_file.getOptions().addMSLevel(2);     // only load msLevel 2
    mzml_file.setLogType(log_type_);
    mzml_file.load(inputfile_name, exp);

    // we need to replace the native id with a simple numbering schema, to be able to
    // map the IDs back to the spectra (RT, and MZ information)
    Size native_id(0);
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      it->setNativeID(++native_id);
    }

    // We store the file in mzData file format, because mgf file somehow produce in most
    // of the cases ids with charge 2+. We do not use the input file of this TOPP-tools
    // because XTandem sometimes stumbles over misleading substrings in the filename,
    // e.g. mzXML ...
    MzDataFile mzdata_outfile;
    mzdata_outfile.store(tandem_input_filename, exp);

    XTandemInfile infile;
    infile.setInputFilename(tandem_input_filename);
    infile.setOutputFilename(tandem_output_filename);

    ofstream tax_out(tandem_taxonomy_filename.c_str());
    tax_out << "<?xml version=\"1.0\"?>" << endl;
    tax_out << "\t<bioml label=\"x! taxon-to-file matching list\">" << endl;
    tax_out << "\t\t<taxon label=\"OpenMS_dummy_taxonomy\">" << endl;
    tax_out << "\t\t\t<file format=\"peptide\" URL=\"" << db_name << "\" />" << endl;
    tax_out << "\t</taxon>" << endl;
    tax_out << "</bioml>" << endl;
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
    infile.setMaxValidEValue(getDoubleOption_("max_valid_expect"));
    infile.setNumberOfMissedCleavages(getIntOption_("missed_cleavages"));
    infile.setRefine(getFlag_("refinement"));
    infile.setSemiCleavage(getFlag_("semi_cleavage"));

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
      id -= 1; // native IDs were written 1-based
      if (id < exp.size())
      {
        it->setMetaValue("RT", exp[id].getRT());
        DoubleReal pre_mz = 0.0;
        if (!exp[id].getPrecursors().empty()) pre_mz = exp[id].getPrecursors()[0].getMZ();
        it->setMetaValue("MZ", pre_mz);
        it->removeMetaValue("spectrum_id");
      }
      else
      {
        cerr << "XTandemAdapter: Error: id '" << id << "' not found in peak map!" << endl;
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

    protein_id.setSearchParameters(search_parameters);
    protein_id.setSearchEngineVersion("");
    protein_id.setSearchEngine("XTandem");

    protein_ids.push_back(protein_id);

    IdXMLFile id_output;
    id_output.store(outputfile_name, protein_ids, peptide_ids);

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

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPXTandemAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
