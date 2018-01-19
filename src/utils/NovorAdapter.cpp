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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

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

    @brief Novoradapter for de novo sequencing from tandem mass spectrometry data.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> MS2-Filtering with @ref TOPP_FileFilter</td>
            <td VALIGN = "middle" ROWSPAN=2> \f$ \longrightarrow \f$ NovorAdapter \f$</td>
        </tr>
    </table>
</CENTER>


    This tool can be used for de novo sequencing of peptides from MS/MS data.
    Please use MS2-Spectra only. If filtering is needed please use the @ref TOPP_FileFilter.

    Novor must be installed before this wrapper can be used. This wrapper was successfully tested with version v1.06.0634 (stable).
    
    Novor settings can be either used via command line or directly using a param file (param.txt).

    Parameter names have been changed to match names found in other search engine adapters. For further information  check the Novor wiki (http://wiki.rapidnovor.com/wiki/Main_Page) and the official tool website (https://www.rapidnovor.com/). 

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
    TOPPBase("NovorAdapter", "Template for Tool creation", false, 
    {
      { "Ma Bin", "Novor: Real-Time Peptide de Novo Sequencing Software", "Journal of The American Society for Mass Spectrometry; 30 June 2015", "0.1007/s13361-015-1204-0"        
      }
    })
    {}

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_()
  {
    // thirdparty executable 
    registerInputFile_("executable", "<jar>", "novor.jar", "novor.jar", false, false, ListUtils::create<String>("skipexists"));
    // input, output and parameter file 
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
    registerDoubleOption_("fragment_mass_tolerance", "<double>", 0.5, "Fragmentation error tolerance  (Da)", false);
    registerDoubleOption_("precursor_mass_tolerance", "<double>" , 15.0, "Precursor error tolerance  (ppm or Da)", false);
    registerStringOption_("precursor_error_units", "<choice>", "Da", "Unit of precursor mass tolerance", false);
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
   registerInputFile_("novorFile", "<file>", "", "File to introduce customized algorithm parameters for advanced users (otional .novor file)", false);
   setValidFormats_("novorFile", ListUtils::create<String>("novor"));
   
   registerInputFile_("java_executable", "<file>", "java", "The Java executable. Usually Java is on the system PATH. If Java is not found, use this parameter to specify the full path to Java", false, false, ListUtils::create<String>("skipexists"));
   registerIntOption_("java_memory", "<num>", 3500, "Maximum Java heap size (in MB)", false);

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
    vector<String> variable_mods = getStringList_("variable_modifications");
    vector<String> fixed_mods = getStringList_("fixed_modifications");
    vector<String> forbidden_residues = getStringList_("forbiddenResidues");
  
    String variable_mod = ListUtils::concatenate(variable_mods, ',');
    String fixed_mod = ListUtils::concatenate(fixed_mods, ',');
    String forbidden_res = ListUtils::concatenate(forbidden_residues, ',');

    os << "enzyme = " << getStringOption_("enzyme") << "\n"
       << "fragmentation = " << getStringOption_("fragmentation") << "\n"
       << "massAnalyzer = " << getStringOption_("massAnalyzer") << "\n"
       << "fragmentIonErrorTol = " << getDoubleOption_("fragment_mass_tolerance") << "Da" << "\n"
       << "precursorErrorTol = " << getDoubleOption_("precursor_mass_tolerance") << getStringOption_("precursor_error_units") << "\n"
       << "variableModifications = " << variable_mod << "\n"
       << "fixedModifications = "    << fixed_mod << "\n"
       << "forbiddenResidues = " << forbidden_res << "\n";
  
    // novorFile for custom alogrithm parameters of nova
    String cparamfile = getStringOption_("novorFile");
    ifstream cpfile(cparamfile);
    if (!cpfile)
    {
      os << "novorFile" << cparamfile << "\n";
    }
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
    String java_executable = getStringOption_("java_executable");
    QString java_memory = "-Xmx" + QString::number(getIntOption_("java_memory")) + "m";

    QString executable = getStringOption_("executable").toQString();   

    if (executable.isEmpty())
    {
      const QProcessEnvironment env;
      const QString & qnovorpathenv = env.systemEnvironment().value("NOVOR_PATH");
      if (qnovorpathenv.isEmpty())
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

    String tmp_out = tmp_dir + "tmp_out_novor.csv";
  
    QStringList process_params;
    process_params << java_memory
                   << "-jar" << executable
                   << "-f" 
                   << "-o" << tmp_out.toQString()               
                   << "-p" << tmp_param.toQString()
                   << tmp_mgf.toQString();
    
    QProcess qp;
    qp.setWorkingDirectory(path_to_executable);
    qp.start(java_executable.toQString(), process_params);
 
    // novor command line
    std::stringstream ss;
    ss << "COMMAND: " << executable.toStdString();
    for (QStringList::const_iterator it = process_params.begin(); it != process_params.end(); ++it)
    {
        ss << " " << it->toStdString();
    }
    LOG_DEBUG << ss.str() << endl;

    // see if process was successfull
    const bool success = qp.waitForFinished(-1);

    if (success == false || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_( "FATAL: External invocation of Novor failed. Standard output and error were:");
      const QString nr_stdout(qp.readAllStandardOutput());
      const QString nr_stderr(qp.readAllStandardError());
      writeLog_(nr_stdout);
      writeLog_(nr_stderr);
      writeLog_(String(qp.exitCode()));

      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

   ifstream file(tmp_out);
   if (file) 
   {
     CsvFile csv(tmp_out, ',');
          
     vector<PeptideIdentification> peptide_ids;
     for (Size i = 0; i != csv.rowCount(); ++i)
     {
       StringList sl;
       csv.getRow(i, sl);
       
       if (sl.empty() || sl[0][0] == '#') { continue; }
       
       PeptideIdentification pi;
       pi.setMetaValue("scan_index", sl[1].toDouble());
       pi.setScoreType("novorscore");
       pi.setHigherScoreBetter(true);
       pi.setRT(sl[2].toDouble());
       pi.setMZ(sl[3].toDouble());

       PeptideHit ph;
       ph.setCharge(sl[4].toInt());
       ph.setScore(sl[8].toDouble());

       // replace PTM name (see http://wiki.rapidnovor.com/wiki/Built-in_PTMs)
       String sequence = sl[9];
       sequence.substitute("(Cam)", "(Carbamidomethyl)");
       sequence.substitute("(O)","(Oxidation)");
       sequence.substitute("(PyroCam)", "(Pyro-carbamidomethyl)");
       
       ph.setSequence(AASequence::fromString(sequence));      
       ph.setMetaValue("pepMass(denovo)", sl[5].toDouble());
       ph.setMetaValue("err(data-denovo)", sl[6].toDouble());
       ph.setMetaValue("ppm(1e6*err/(mz*z))", sl[7].toDouble());
       ph.setMetaValue("aaScore", sl[10].toQString());

       pi.getHits().push_back(ph);   
       peptide_ids.emplace_back(pi);
       
     } 

     // extract version from comment 
     // #              v1.06.0634 (stable)
     vector<ProteinIdentification> protein_ids;
     StringList versionrow;
     csv.getRow(2, versionrow);
     versionrow[0].suffix('#').trim();
        
     protein_ids = vector<ProteinIdentification>(1);
     protein_ids[0].setDateTime(DateTime::now());
     protein_ids[0].setSearchEngine("Novor");
     protein_ids[0].setSearchEngineVersion(versionrow[0]);

     ProteinIdentification::SearchParameters search_parameters;
     search_parameters.db = "denovo";
     search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
    
     // if a parameter file is used the modifications need to be parsed from the novor output csv
     search_parameters.fixed_modifications = getStringList_("fixed_modifications");
     search_parameters.variable_modifications = getStringList_("variable_modifications");
     search_parameters.fragment_mass_tolerance = getDoubleOption_("fragment_mass_tolerance");
     search_parameters.precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
     search_parameters.precursor_mass_tolerance_ppm = getStringOption_("precursor_error_units") == "ppm" ? true : false;
     search_parameters.fragment_mass_tolerance_ppm = false;
     search_parameters.digestion_enzyme = *ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme"));
     
     //StringList inputFile;
     //inputFile[0] = in;
     //protein_ids[0].setPrimaryMSRunPath(inputFile); 
     protein_ids[0].setSearchParameters(search_parameters);
     
     IdXMLFile().store(out, protein_ids, peptide_ids);

   }
   else
   {
     writeLog_("Novor output is empty! No IdXML output was generated.");
   } 

   // remove tempdir
  removeTempDir_(tmp_dir);

  return EXECUTION_OK;
  }
};


// the actual main function needed to create an executable
int main(int argc, const char ** argv)
{
  TOPPNovorAdapter tool;
  return tool.main(argc, argv);
}
/// @endcond
