// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/METADATA/SpectrumLookup.h>
// TOOD Remove once we've moved transform to handler
#include <OpenMS/FORMAT/MzMLFile.h>
// TODO remove once we have store spectrum in handler
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>

#include <OpenMS/SYSTEM/JavaInfo.h>

#include <QFileInfo>

#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_NovorAdapter NovorAdapter

@brief Novoradapter for de novo sequencing from tandem mass spectrometry data.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; NovorAdapter &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> signal-/preprocessing e.g. @ref TOPP_FileFilter  @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>


This tool can be used for de novo sequencing of peptides from MS/MS data.

Novor must be installed before this wrapper can be used. This wrapper was successfully tested with version v1.06.0634 (stable).

Novor settings can be either used via command line or directly using a param file (param.txt).

Parameter names have been changed to match names found in other search engine adapters. For further information  check the Novor wiki (http://wiki.rapidnovor.com/wiki/Main_Page) and the official tool website (https://www.rapidnovor.com/). 

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_NovorAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_NovorAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNovorAdapter :
  public TOPPBase
{
public:
  TOPPNovorAdapter() :
    TOPPBase("NovorAdapter", "Performs de novo sequencing of peptides from MS/MS data with Novor.", true, 
    {
      Citation{"Ma Bin",
               "Novor: Real-Time Peptide de Novo Sequencing Software",
               "Journal of The American Society for Mass Spectrometry; 30 June 2015",
               "10.1007/s13361-015-1204-0"}
    })
    {}

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    // thirdparty executable 
    registerInputFile_("executable", "<jar>", "novor.jar", "novor.jar", false, false, ListUtils::create<String>("skipexists"));
    // input, output and parameter file 
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
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
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "Unit of precursor mass tolerance", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("ppm,Da"));
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

  void createParamFile_(ofstream& os)
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
    if (cpfile)
    {
      os << "novorFile" << cparamfile << "\n";
    }
    os.close();
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // determine the executable
    //-------------------------------------------------------------
    const String java_executable = getStringOption_("java_executable");
    QString java_memory = "-Xmx" + QString::number(getIntOption_("java_memory")) + "m";

    QString executable = getStringOption_("executable").toQString();   

    if (executable.isEmpty())
    {
      const char* novor_path_env = getenv("NOVOR_PATH");
      if (novor_path_env == nullptr || strlen(novor_path_env) == 0)
      {
        writeLogError_("FATAL: Executable of Novor could not be found. Please either use NOVOR_PATH env variable or provide via '-executable'!");
        return MISSING_PARAMETERS;
      }
      executable = novor_path_env;
    }

    // Normalize file path
    QFileInfo file_info(executable);
    executable = file_info.canonicalFilePath();

    writeLogInfo_("Executable is: " + executable);
    const QString & path_to_executable = File::path(executable).toQString();
    
    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    
    // tmp_dir
    File::TempDir tmp_dir(debug_level_ >= 2);

    // parameter file
    String tmp_param = tmp_dir.getPath() + "param.txt";    
    ofstream os(tmp_param.c_str());
    createParamFile_(os);

    // convert mzML to mgf format
    Size count_written{ 0 };
    String tmp_mgf = tmp_dir.getPath() + "tmp_mgf.mgf";
    {
      std::ofstream ofs(tmp_mgf, std::ofstream::out);
      MascotGenericFile mgf;

      auto f = [&ofs, &in, &mgf, &count_written](const MSSpectrum& spec)
      {
        if (spec.getType(true) == MSSpectrum::SpectrumType::CENTROID)
        {
          mgf.writeSpectrum(ofs, spec, in, "UNKNOWN");
          ++count_written;
        }
      };
      MSDataTransformingConsumer c;
      c.setSpectraProcessingFunc(f);
      MzMLFile mzml_file;
      mzml_file.getOptions().setMSLevels({2});
      mzml_file.transform(in, &c, true);
      OPENMS_LOG_INFO << "Info: " << count_written << " MS2 spectra will be passed to NOVOR\n";
      if (count_written == 0)
      {
        if (getFlag_("force"))
        {
          OPENMS_LOG_WARN << "Warning: no MS2 spectra were passed to NOVOR. The result will be empty! To make this an error, remove the '-force' flag!" << std::endl;
        }
        else
        {
          OPENMS_LOG_ERROR << "Error: no MS2 spectra were passed to NOVOR. This is probably not intentional. Please check the input data '" << in << "'. To make this a warning, set the '-force' flag." << std::endl;
          return ExitCodes::INPUT_FILE_EMPTY;
        }
      }
    }

    //-------------------------------------------------------------
    // process
    //-------------------------------------------------------------

    String tmp_out = tmp_dir.getPath() + "tmp_out_novor.csv";

    QStringList process_params;
    process_params << java_memory
                   << "-jar" << executable
                   << "-f" 
                   << "-o" << tmp_out.toQString()               
                   << "-p" << tmp_param.toQString()
                   << tmp_mgf.toQString();


    // print novor command line
    TOPPBase::ExitCodes exit_code = runExternalProcess_(java_executable.toQString(), process_params, path_to_executable);
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // Build mapping to get correct spectrum reference later
    MSExperiment exp;
    MzMLFile m;
    PeakFileOptions op;
    op.setFillData(false); // no actual peak data
    op.setMSLevels({ 2 }); //only MS2
    m.setOptions(op);
    m.load(in, exp);
    SpectrumLookup mapping;
    mapping.readSpectra(exp);

    if (!File::exists(tmp_out))
    {
      OPENMS_LOG_ERROR << "Error: NOVOR did not write the output file '" << tmp_out << "' as requested!\n";
      return ExitCodes::EXTERNAL_PROGRAM_ERROR;
    }
    CsvFile csv(tmp_out, ',');
        
    vector<PeptideIdentification> peptide_ids;
    for (Size i = 0; i != csv.rowCount(); ++i)
    {
      StringList sl;
      csv.getRow(i, sl);
        
      if (sl.empty() || sl[0][0] == '#') { continue; }
        
      PeptideIdentification pi;
      pi.setSpectrumReference( exp[mapping.findByScanNumber(sl[1].toInt())].getNativeID());
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

      pi.getHits().push_back(std::move(ph));   
      peptide_ids.push_back(std::move(pi));
    } 

    // extract version from comment 
    // #              v1.06.0634 (stable)
    // v1.06.0634 (stable)
    vector<ProteinIdentification> protein_ids;
    StringList versionrow;
    csv.getRow(2, versionrow);
    String version = versionrow[0].substr(versionrow[0].find("v."));
      
    protein_ids = vector<ProteinIdentification>(1);
    protein_ids[0].setDateTime(DateTime::now());
    protein_ids[0].setSearchEngine("Novor");
    protein_ids[0].setSearchEngineVersion(version);

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
     
    OPENMS_LOG_INFO << "NOVOR created " << peptide_ids.size() << " PSMs from " << count_written << " MS2 spectra (" << (peptide_ids.size() * 100 / count_written) << "% annotated)\n";

    FileHandler().storeIdentifications(out, protein_ids, peptide_ids, {FileTypes::IDXML});

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
