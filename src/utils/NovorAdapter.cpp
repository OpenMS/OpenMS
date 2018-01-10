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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <QDir>
#include <QProcess>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_NovorAdapter NovorAdapter

    @brief Novoradapter does ...

    This tool can be used for scientific stuff.

    And more scientific applications.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_NovorAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_NovorAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNovorAdapter :
  public TOPPBase
{
public:
  TOPPNovorAdapter() :
    TOPPBase("NovorAdapter", "Template for Tool creation", false)
  {

  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    // thirdparty executable 
    registerInputFile_("executable", "<executable>", "", "novor executable", false, false, ListUtils::create<String>("skipexists"));
    // input and output
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));
    registerOutputFile_("out", "<file>", "", "Novor idXML output");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    // enzyme
    registerStringOption_("enzyme", "<choice>", "Trypsin", "Digestion enzyme - currently only Trypsin is supported ", false);
    setValidStrings_("enzyme", ListUtils::create<String>("Trypsin"));
    // instrument
    registerStringOption_("fragmentation", "<choice>", "CID", "Fragmentation method", false);
    setValidStrings_("fragmentation", ListUtils::create<String>("CID,HCD"));
    registerStringOption_("massAnalyzer", "<choice>" , "Trap", "MassAnalyzer e.g. (Oritrap CID-Trap, CID-FT, HCD-FT; QTof CID-TOF)", false);
    setValidStrings_("massAnalyzer", ListUtils::create<String>("Trap,TOF,FT"));
    // mass error tolerance
    registerStringOption_("fragmentIon_error_tolerance", "<string>", "0.5", "Fragmentation error tolerance  (Da)", false);
    registerStringOption_("precursor_error_tolerance", "<string>" , "15", "Precursor error tolerance  (ppm or Da)", false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("Da,ppm"));
    // post-translational-modification
    registerStringList_("variable_modifications", "<mods>", vector<String>(), "Variable modifications", false);
    setValidStrings_("variable_modifications", ListUtils::create<String>("Acetyl (K),Acetyl (N-term),Amidated (C-term),Ammonia-loss (N-term C),Biotin (K),Biotin (N-term),Carbamidomethyl (C),Carbamyl (K),Carbamyl (N-term),Carboxymethyl (C),Deamidated (NQ),Dehydrated (N-term C),Dioxidation (M),Methyl (C-term),Methyl (DE),Oxidation (M),Oxidation (HW),Phospho (ST),Phospho (Y),Pyro-carbamidomethyl (N-term C),Pyro-Glu (E),Pyro-Glu (Q),Sodium (C-term),Sodium (DE),Sulfo (STY),Trimethyl (RK)"));
    registerStringList_("fixed_modifications", "<mods>", vector<String>(), "Fixed modifications", false);
    setValidStrings_("fixed_modifications", ListUtils::create<String>("Acetyl (K),Acetyl (N-term),Amidated (C-term),Ammonia-loss (N-term C),Biotin (K),Biotin (N-term),Carbamidomethyl (C),Carbamyl (K),Carbamyl (N-term),Carboxymethyl (C),Deamidated (NQ),Dehydrated (N-term C),Dioxidation (M),Methyl (C-term),Methyl (DE),Oxidation (M),Oxidation (HW),Phospho (ST),Phospho (Y),Pyro-carbamidomethyl (N-term C),Pyro-Glu (E),Pyro-Glu (Q),Sodium (C-term),Sodium (DE),Sulfo (STY),Trimethyl (RK)"));
   // forbidden residues
   registerStringList_("forbiddenResidues", "<mods>", vector<String>(), "Forbidden Resiudes", false);
   setValidStrings_("forbiddenResidues", ListUtils::create<String>("I,U"));
 
   // parameter novorFile will not be wrapped here
  }

  // remove temporary folder 
  void removeTempDir_(const String& tmp_dir)
  {
    if (tmp_dir.empty()) {return;} // no temporary directory created

    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + tmp_dir + "'. Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (debug_level_ == 1) 
      {
        writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 1);
      }
      File::removeDirRecursively(tmp_dir);
    }
  }

  void createParamFile_(ostream& os)
  {
  vector<String> variable_mods = getStringList_("variable__modifications");
  vector<String> fixed_mods = getStringList_("fixed_modifications");
  vector<String> forbidden_residues = getStringList_("forbiddenResidues");
  bool no_mods = fixed_mods.empty() && variable_mods.empty();
  if (!no_mods)
    {
      writeLog_("Warning: Modifications are defined ('fixed_modifications'/'variable_modifications'), but the number of modifications is zero. Is that intended?");
    }
  
  String variable_mod = ListUtils::concatenate(variable_mods, ',');
  String fixed_mod = ListUtils::concatenate(fixed_mods, ',');
  String forbidden_res = ListUtils::concatenate(forbidden_residues, ',');

  os << "enzyme = " << getStringOption_("enzyme") << "\n"
     << "fragmentation = " << getStringOption_("fragmentation") << "\n"
     << "massAnalyzer = " << getStringOption_("massAnalyzer") << "\n"
     << "fragmentIonError = " << getStringOption_("fragmentIon_error_tolerance") << "Da" << "\n"
     << "precursorErrorTol = " << getStringOption_("precursor_error_tolerance") << getStringOption_("precursor_error_units") << "\n"
     << "variableModifications = " << variable_mod << "\n"
     << "fixedModifications = "    << fixed_mod << "\n"
     << "forbiddenResidues = " << forbidden_res << "\n";

  std::cout << os << std::endl;
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");     
    
    if (out.empty())
    {
      writeLog_("Fatal error: no output file given");
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // determining the executable
    //-------------------------------------------------------------
    
    QString executable = getStringOptions_("exectuable").toQString();   

    if (executable.isEmpty())
    {
      const QProcessEnvironment env;
      const QString & qnovorpathenv = env.systemEnvironment().value("NOVOR_PATH");
      if (novorpathenv.isEmpty())
      {
        writeLog_( "FATAL: Executable of Novor could not be found. Please either use NOVOR_PATH env variable or provide with -executable");
        return MISSING_PARAMETERS;
      }
      executable = qnovorpathenv;
    }

    //Normalize file path
    QFileInfo file_info(executable);
    executable = file_info.canonicalFilePath();

    writeLog_("Executable is: " + executable);
    const QString & path_to_executable = File::path(executable).toQString();
                
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    
    //tmp_dir
    const String tmp_dir = makeTempDirectory_();
    writeDebug_("Creating temporary directory '" + tmp_dir + "'", 1);

    // parameter file
    String tmp_param = tmp_dir + "param.txt";    
    ofstream os(tmp_param.c_str());
    createParamFile_(os);
    os.close();

    // convert mzML to mgf format
    MzMLFile f;
    MSExperiment exp;
    f.setLogType(log_type_);
    f.load(in, exp);
 
    String tmp_mgf = tmp_dir + "tmp_mgf.mgf";
    MascotGenericFile mgf;
    mgf.setLogType(log_type_);
    mgf.store(tmp_mgf,exp);

    //-------------------------------------------------------------
    // process
    //-------------------------------------------------------------

    String tmp_out = tmp_dir + "novor_output";
    QStringList process_params;
    process_params << "-f" 
                   << "-o" << tmp_out.toQString()
                   << "-p" << tmp_param.toQString(;)
    
    //TODO: How does ist work with a .bat (Batchfile - Win & sh file linux/mac) 
    QProcess qp;
    qp.setWorkingDirectory(path_to_executable);
    qp.start("/bin/sh",executable, process_params);
    const bool success = qp.waitForFinished(-1));

    if (success == false || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_( "FATAL: External invocation of Novor failed. Standard output and error were:");
      const QString nr_stdout(qp.readAllStandardOutput());
      const QString nr_stderr(qp.readAllErrorOutput());
      writeLog_(nr_stdout);
      writeLog_(nr_stderr);
      writeLog_(String(qp.exitCode()));

      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // TODO: parse novor output into internal data structure
    // and write idXML to output 
    // delete tmp dir 

  }
  return EXECUTION_OK;
};


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPNovorAdapter tool;
  return tool.main(argc, argv);
}
/// @endcond
