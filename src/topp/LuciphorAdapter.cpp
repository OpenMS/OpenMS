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
// $Maintainer: Petra Gutenbrunner $
// $Authors: Petra Gutenbrunner $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

#include <cstddef>
#include <fstream>
#include <map>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_LuciphorAdapter LuciphorAdapter

   @brief Adapter for the LuciPHOr2: a site localisation tool of generic post-translational modifications from tandem mass spectrometry data.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ LuciphorAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    LuciPHOr2 must be installed before this wrapper can be used. Please make sure that Java and LuciPHOr2 are working.@n
    At the time of writing, it could be downloaded from http://luciphor2.sourceforge.net.

    Input spectra for LuciPHOr2 have to be in pepXML file format. The input mzML file must be the same as the one used to create the pepXML input file.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_LuciphorAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_LuciphorAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

using namespace OpenMS;
using namespace std;

class LuciphorAdapter :
  public TOPPBase
{
public:
  LuciphorAdapter() :
    TOPPBase("LuciphorAdapter", "Modification site localisation using LuciPHOr2.", true),
    // parameter choices (the order of the values must be the same as in the LuciPHOr2 parameters!):
    fragment_methods_(ListUtils::create<String>("CID,HCD")),
    fragment_error_units_(ListUtils::create<String>("Daltons,PPM")),
    input_types_(ListUtils::create<String>("pepXML"))
  {
  }

protected:
  // lists of allowed parameter values:
  vector<String> fragment_methods_, fragment_error_units_, input_types_;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("spectrum_in", "<file>", "", "Input spectrum file");	
    setValidFormats_("spectrum_in", ListUtils::create<String>("mzML,mgf,mzXML"));
    
    registerInputFile_("in", "<file>", "", "Input file");	
    setValidFormats_("in", input_types_);
	
    registerOutputFile_("out", "<file>", "", "Output file", false);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
	
    registerInputFile_("executable", "<file>", "luciphor2.jar", "LuciPHOr2 .jar file, e.g. 'c:\\program files\\luciphor2.jar'", true, false, ListUtils::create<String>("skipexists"));

    registerStringOption_("fragment_method", "<choice>", fragment_methods_[0], "Fragmentation method", false);
    setValidStrings_("fragment_method", fragment_methods_);
    
    registerDoubleOption_("fragment_mass_tol", "<value>", 0.5, "Tolerance of the peaks in the fragment spectrum", false);
    registerStringOption_("fragment_error_units", "<choice>", fragment_error_units_[1], "Unit of fragment mass tolerance", false);
    setValidStrings_("fragment_error_units", fragment_error_units_);
	
    registerDoubleOption_("min_mz", "<value>", 100.0, "Do not consider peaks below this value for matching fragment ions", false);
    
    vector<String> all_mods;
    ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
    registerStringList_("target_modifications", "<mods>", vector<String>(), "List the amino acids to be searched for and their mass modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)'", false);
    setValidStrings_("target_modifications", all_mods);
    
    registerStringList_("neutral_losses", "<value>", vector<String>(), "List the types of neutral losses that you want to consider. The residue field is case sensitive. For example: lower case 'sty' implies that the neutral loss can only occur if the specified modification is present. Syntax: NL = <RESDIUES> -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>", false);
    
    registerDoubleOption_("decoy_mass", "<value>", 79.966331, "How much to add to an amino acid to make it a decoy", false);
    setMinFloat_("decoy_mass", 1.0);
    registerStringList_("decoy_neutral_losses", "<value>", vector<String>(), "For handling the neutral loss from a decoy sequence. The syntax for this is identical to that of the normal neutral losses given above except that the residue is always 'X'. Syntax: DECOY_NL = X -<NEUTRAL_LOSS_MOLECULAR_FORMULA> <MASS_LOST>", false);
    
    registerIntOption_("max_charge_state", "<num>", 5, "Do not consider PSMs with a charge state above this value", false);
    setMinInt_("max_charge_state", 1);
	
    registerIntOption_("max_peptide_length", "<num>", 40, "Restrict scoring to peptides with a length shorter than this value", false);
    setMinInt_("max_peptide_length", 1);
	
    registerIntOption_("max_num_perm", "<num>", 16384, "Maximum number of permutations a sequence can have", false);
    setMinInt_("max_num_perm", 1);

    registerStringOption_("selection_method", "<choice>", "0", "Score selection method: 0 = Peptide Prophet probability (default), 1 = Mascot Ion Score, 2 = -log(E-value), 3 = X!Tandem Hyperscore, 4 = Sequest Xcorr", false);
    setValidStrings_("selection_method", ListUtils::create<String>("0,1,2,3,4"));
	
    registerDoubleOption_("modeling_score_threshold", "<value>", 0.95, "Minimum score a PSM needs to be considered for modeling", false);
    setMinFloat_("modeling_score_threshold", 0.0);
    
    registerDoubleOption_("scoring_threshold", "<value>", 0.0, "PSMs below this value will be discarded", false);
    setMinFloat_("scoring_threshold", 0.0);
    
    registerIntOption_("min_num_psms_model", "<num>", 50, "The minimum number of PSMs you need for any charge state in order to build a model for it", false);
    setMinInt_("min_num_psms_model", 1);

    registerIntOption_("num_threads", "<num>", 6, "For multi-threading, zero = use all CPU found by JAVA", false);
    setMinInt_("num_threads", 0);

    registerStringOption_("run_mode", "<choice>", "0", "Determines how Luciphor will run: 0 = calculate FLR then rerun scoring without decoys (two iterations), 1 = Report Decoys: calculate FLR but don't rescore PSMs, all decoy hits will be reported", false);
    setValidStrings_("run_mode", ListUtils::create<String>("0,1")); 
    
    //registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);
    //registerIntOption_("java_permgen", "<num>", 0, "Maximum Java permanent generation space (in MB); only for Java 7 and below", false, true);
  }
  
  String makeModString_(const String& mod_name)
  {
    ResidueModification mod = ModificationsDB::getInstance()->getModification(mod_name);
    String residue = mod.getOrigin();
    
    return String(residue +  " " + mod.getDiffMonoMass());    
  }
 
  ExitCodes parseParameters_(map<String, vector<String> >& config_map)
  {
    FileHandler fh;
    
    String spectrum_in = getStringOption_("spectrum_in");    
    config_map["SPECTRUM_PATH"].push_back(File::path(File::absolutePath(spectrum_in)));
    config_map["SPECTRUM_SUFFIX"].push_back(FileTypes::typeToName(fh.getTypeByFileName(spectrum_in)));
    
    String in = getStringOption_("in");
    config_map["INPUT_DATA"].push_back(getStringOption_("in"));
    
    String type = FileTypes::typeToName(fh.getTypeByFileName(in));
    config_map["INPUT_TYPE"].push_back(ListUtils::getIndex<String>(input_types_, type));
    
    config_map["ALGORITHM"].push_back(ListUtils::getIndex<String>(fragment_methods_, getStringOption_("fragment_method")));
    config_map["MS2_TOL"].push_back(getDoubleOption_("fragment_mass_tol"));
    config_map["MS2_TOL_UNITS"].push_back(ListUtils::getIndex<String>(fragment_error_units_, getStringOption_("fragment_error_units")));
    config_map["MIN_MZ"].push_back(getDoubleOption_("min_mz"));
    config_map["OUTPUT_FILE"].push_back(getStringOption_("out"));
    config_map["DECOY_MASS"].push_back(getDoubleOption_("decoy_mass"));
    config_map["MAX_CHARGE_STATE"].push_back(getIntOption_("max_charge_state"));
    config_map["MAX_PEP_LEN"].push_back(getIntOption_("max_peptide_length"));
    config_map["MAX_NUM_PERM"].push_back(getIntOption_("max_num_perm"));
    config_map["SELECTION_METHOD"].push_back(getStringOption_("selection_method"));
    config_map["MODELING_SCORE_THRESHOLD"].push_back(getDoubleOption_("modeling_score_threshold"));
    config_map["SCORING_THRESHOLD"].push_back(getDoubleOption_("scoring_threshold"));
    config_map["MIN_NUM_PSMS_MODEL"].push_back(getIntOption_("min_num_psms_model"));
    config_map["NUM_THREADS"].push_back(getIntOption_("num_threads"));
    config_map["RUN_MODE"].push_back(getStringOption_("run_mode"));
    
    // list values
    vector<String> target_mods = getStringList_("target_modifications");
    if (target_mods.empty())
    {
      writeLog_("Error: No target modification existing.");
      return ILLEGAL_PARAMETERS;
    }
    
    for (vector<String>::const_iterator it = target_mods.begin(); it != target_mods.end(); ++it)
    {
      config_map["TARGET_MOD"].push_back(makeModString_(*it));
    }
    
    vector<String> neutral_losses = getStringList_("neutral_losses");
    for (vector<String>::const_iterator it = neutral_losses.begin(); it != neutral_losses.end(); ++it)
    {
      config_map["NL"].push_back(*it);
    }
    
    vector<String> dcy_neutral_losses = getStringList_("decoy_neutral_losses");
    for (vector<String>::const_iterator it = dcy_neutral_losses.begin(); it != dcy_neutral_losses.end(); ++it)
    {
      config_map["DECOY_NL"].push_back(*it);
    }
    
    return EXECUTION_OK;
  }
  
  void writeConfigurationFile_(const String& out_path, map<String, vector<String> >& config_map)
  {
    ofstream output(out_path.c_str());
    output << "## Input file for Luciphor2 (aka: LucXor). (part of OpenMS)\n\n";
	
    for (std::map<String, vector<String> >::iterator it = config_map.begin(); it != config_map.end(); ++it)
    {
      String key = it->first;
      if (!key.empty())
      {
        for (vector<String>::iterator it_val = it->second.begin(); it_val != it->second.end(); ++it_val)
        {
          output << key << " = " << *it_val << "\n";
        }        
      }
    }
	
    //------------------------------------------------------------------
    // static parameter definition
    //------------------------------------------------------------------
    // Generate a tab-delimited file of all the matched peaks
    output << "WRITE_MATCHED_PEAKS_FILE = 0\n";
    
    output << "MOD_PEP_REP = 0 ## 0 = show single character modifications, 1 = show TPP-formatted modifications\n";
    
    output << "## This option can be used to help diagnose problems with Luciphor. Multi-threading is disabled in debug mode.\n";
    output << "DEBUG_MODE = 0 ## 0 = default: turn off debugging\n";
    output << "               ## 1 = write peaks selected for modeling to disk\n";
    output << "               ## 2 = write the scores of all permutations for each PSM to disk\n";
    output << "               ## 3 = write the matched peaks for the top scoring permutation to disk\n";
    output << "               ## 4 = write HCD non-parametric models to disk (HCD-mode only option)\n";
  }
  
  void removeTempDir_(const String& temp_dir)
  {
    if (temp_dir.empty()) return; // no temp. dir. created

    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + temp_dir + "'. Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (debug_level_ == 1) writeDebug_("Deleting temporary directory '" + temp_dir + "'. Set debug level to 2 or higher to keep it.", 1);
      File::removeDirRecursively(temp_dir);
    }
  }
  
  ExitCodes main_(int, const char**)
  {
    if (!getFlag_("force"))
    {
      if (!JavaInfo::canRun("java"))
      {
        writeLog_("Fatal error: Java not found, or the Java process timed out. Java is needed to run LuciPHOr2. Make sure that it can be executed by calling 'java', e.g. add the directory containing the Java binary to your PATH variable. If you are certain java is installed, please set the 'force' flag in order to avoid this error message.");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
    else
    {
      writeLog_("The installation of Java was not checked.");
    }
	
    // create temporary directory
    String temp_dir, conf_file;
    temp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString());
    writeDebug_("Creating temporary directory '" + temp_dir + "'", 1);
    QDir d;
    d.mkpath(temp_dir.toQString());
	
    // create a temporary config file for LuciPHOr2 parameters
    conf_file = temp_dir + "luciphor2_input_template.txt";
	
    // initialize map
    map<String, vector<String> > config_map;	
    ExitCodes ret = parseParameters_(config_map);
    if (ret != EXECUTION_OK)
    {
      return ret;
    }
    
    writeConfigurationFile_(conf_file, config_map);    
    QString executable = getStringOption_("executable").toQString();
    
    // Hack for KNIME. Looks for LUCIPHOR_PATH in the environment which is set in binaries.ini
    QProcessEnvironment env;
    String luciphorpath = "LUCIPHOR_PATH";
    QString qluciphorpath = env.systemEnvironment().value(luciphorpath.toQString());

    if (!qluciphorpath.isEmpty())
    {
      executable = qluciphorpath;
    }

    QStringList process_params; // the actual process is Java, not LuciPHOr2!
    process_params << "-jar" << executable << conf_file.toQString();

    // execute LuciPHOr2    
    int status = QProcess::execute("java", process_params);
    if (status != 0)
    {
      writeLog_("Fatal error: Running LuciPHOr2 returned an error code. Does the LuciPHOr2 executable (.jar file) exist?");
      return EXTERNAL_PROGRAM_ERROR;
    }

    removeTempDir_(temp_dir);

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  LuciphorAdapter tool;
  return tool.main(argc, argv);
}
