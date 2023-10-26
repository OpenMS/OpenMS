// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>
#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedAssay.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusFragmentAnnotation.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDebug>
#include <QDir>
#include <QDirIterator>
#include <QtCore/QProcess>
#include <fstream>

using namespace OpenMS;
using namespace std;

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page TOPP_SiriusAdapter SiriusAdapter

  @brief De novo metabolite identification.

  CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

  To use this feature, the Sirius command line tool as well as a java installation is needed.

  Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/ 

  Please use Sirius Version 4.0.1

  If you want to use the software with the Gurobi solver or CPLEX instead of GLPK, please follow the instructions in the sirius manual.

  Internal procedure in SiriusAdapter \n
  1. Input mzML (and optional featureXML) \n
  2. Preprocessing (see below)\n
  3. Parsed by SiriusMSConverter into (sirius internal) .ms format \n
  4. Submission of .ms and additional parameters to wrapped SIRIUS CLI \n
  5. Sirius output saved in internal temporary folder structure \n
  6. Sirius output is parsed (SiriusMzTabWriter/CsiFingerIDMzTabWriter) \n
  7. Merge corresponding output in one mzTab (out_sirius/out_fingerid) \n

  Preprocessing (featureXML): 
  By providing a featureXML, the feature information can be used for feature mapping.
  Sirius will then process the internally merged MS2 spectra allocated to one feature (instead of all available MS2).
  To reduce the feature space even further a masstrace filter can be set. 
  Additional adduct information can be provided using a featureXML from the MetaboliteAdductDecharger or AccurateMassSearch.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SiriusAdapter.html
 */

/// @cond TOPPCLASSES

class TOPPSiriusAdapter :
 public TOPPBase
{
 public:
  TOPPSiriusAdapter() :
    TOPPBase("SiriusAdapter", "Metabolite identification using single and tandem mass spectrometry", false,
      {
        {"Kai Duehrkop and Sebastian Boecker",
         "Fragmentation trees reloaded",
         "J Cheminform; 2016",
         "10.1186/s13321-016-0116-8"},
        {"Kai Duehrkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Boecker",
         "Searching molecular structure databases with tandem mass spectra using CSI:FingerID",
         "Proceedings of the National Academy of Sciences; 2015",
         "10.1073/pnas.1509788112"}
      })
    {}

private:
  SiriusAdapterAlgorithm algorithm;

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("sirius_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
#ifdef OPENMS_WINDOWSPLATFORM
      "sirius.bat",
#else
      "sirius",
#endif
      "The Sirius executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", false, false, {"is_executable"});

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("in_featureinfo", "<file>", "", "FeatureXML input with feature and adduct information", false);
    setValidFormats_("in_featureinfo", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab output file for SIRIUS results", false);
    setValidFormats_("out_sirius", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_fingerid","<file>", "", "MzTab output file for CSI:FingerID, if this parameter is given, SIRIUS will search for a molecular structure using CSI:FingerID after determining the sum formula", false);
    setValidFormats_("out_fingerid", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_ms","<file>", "", "Internal SIRIUS .ms format after OpenMS preprocessing", false);
    setValidFormats_("out_ms", ListUtils::create<String>("ms"));

    registerOutputFile_("out_annotated_spectra","<file>", "", "Export spectra with fragment annotations from SIRIUS", false);
    setValidFormats_("out_annotated_spectra", ListUtils::create<String>("mzML"));

    registerStringOption_("out_project_space", "<directory>", "", "Output directory for SIRIUS project space", false);

    registerStringOption_("sirius_user_email", "<string>", "", "E-mail for your SIRIUS account.", false);
    
    registerStringOption_("sirius_user_password", "<string>", "", "Password for your SIRIUS account.", false);

    registerFlag_("converter_mode", "Use this flag in combination with the out_ms file to convert the input mzML and featureXML to a .ms file. Without further SIRIUS processing.", true);

    addEmptyLine_();

    auto defaults = algorithm.getDefaults();
    defaults.remove("project:processors");

    registerFullParam_(defaults);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------
    String sirius_executable = getStringOption_("sirius_executable");
    String in = getStringOption_("in");
    String out_sirius = getStringOption_("out_sirius");
    String out_csifingerid = getStringOption_("out_fingerid");
    String featureinfo = getStringOption_("in_featureinfo");
    String out_ms = getStringOption_("out_ms");
    String out_ann_spectra = getStringOption_("out_annotated_spectra");
    String sirius_workspace_directory = getStringOption_("out_project_space");
    String sirius_user_email = getStringOption_("sirius_user_email");
    String sirius_user_password = getStringOption_("sirius_user_password");
    bool converter_mode = getFlag_("converter_mode");

    auto params = getParam_();
    if (debug_level_ > 3)
    {
      params.setValue("read_sirius_stdout", "true");
    }
    params.setValue("project:processors", params.getValue("threads"));
    algorithm.updateExistingParameter(params);

    writeDebug_("Parameters passed to SiriusAdapterAlgorithm", algorithm.getParameters(), 3);

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------
    MSExperiment spectra;
    FileHandler().loadExperiment(in, spectra, {FileTypes::MZML}, log_type_);

    // make temporary files
    SiriusAdapterAlgorithm::SiriusTemporaryFileSystemObjects sirius_tmp(debug_level_);

    // run masstrace filter and feature mapping
    FeatureMapping::FeatureMappingInfo fm_info;
    FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in Feature Maps stored in fm_info using a KDTree
    algorithm.preprocessingSirius(featureinfo,
                                  spectra,
                                  fm_info,
                                  feature_mapping);

    // returns Log of feature and/or spectra number
    algorithm.logFeatureSpectraNumber(featureinfo, feature_mapping, spectra);

    // write msfile and store the compound information in CompoundInfo Object
    vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
    SiriusMSFile::store(spectra,
                        sirius_tmp.getTmpMsFile(),
                        feature_mapping,
                        algorithm.isFeatureOnly(),
                        algorithm.getIsotopePatternIterations(),
                        algorithm.isNoMasstraceInfoIsotopePattern(),
                        v_cmpinfo);

    // converter_mode enabled (only needed for SiriusAdapter)
    if (!out_ms.empty() && converter_mode)
    {
      QFile::copy(sirius_tmp.getTmpMsFile().toQString(), out_ms.toQString());
      
      OPENMS_LOG_WARN << "SiriusAdapter was used in converter mode and is terminated after OpenMS preprocessing. \n"
                         "If you would like to run SIRIUS internally please disable the converter mode." << std::endl;
      
      return EXECUTION_OK;
    }

    // log in to Sirius account if email and password are specified
    if (!sirius_user_email.empty() && !sirius_user_password.empty())
    {
      algorithm.logInSiriusAccount(sirius_executable, sirius_user_email, sirius_user_password);
    }
    else OPENMS_LOG_WARN << "No Sirius user account login information specified!" << std::endl;

    // calls Sirius and returns vector of paths to sirius folder structure
    std::vector<String> subdirs = algorithm.callSiriusQProcess(sirius_tmp.getTmpMsFile(),
                                           sirius_tmp.getTmpOutDir(),
                                           sirius_executable,
                                           out_csifingerid,
                                           false);
    
    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if (subdirs.empty())
    {
        throw OpenMS::Exception::Postcondition(__FILE__,__LINE__, OPENMS_PRETTY_FUNCTION, "Sirius was executed, but an empty output was generated");
    }

    // sort vector path list
    SiriusAdapterAlgorithm::sortSiriusWorkspacePathsByScanIndex(subdirs);
    // TODO make parameter(s)?
    double score_threshold = 0.0;
    bool use_exact_mass = false;
    if (!out_ann_spectra.empty())
    {
      MSExperiment annotations{};
      // extract Sirius FragmentAnnotations from subdirs
      // do not resolve for concat native IDs (i.e., per feature) since we want to write the annotated spectrum for every candidate
      // use 0.0 to not have a score_threshold
      annotations.setSpectra(SiriusFragmentAnnotation::extractSiriusAnnotationsTgtOnly(subdirs, score_threshold, use_exact_mass, false));
      // TODO check if we have duplicate native IDs without resolution and if this is a problem
      FileHandler().storeExperiment(out_ann_spectra, annotations, {FileTypes::MZML}, log_type_);

      // TODO remove the following or use it to add more info to the spectra
      // combine compound information (SiriusMSFile) with annotated spectra (SiriusFragmentAnnotation)
      //vector<MetaboTargetedAssay::CompoundSpectrumPair> v_cmp_spec = MetaboTargetedAssay::pairCompoundWithAnnotatedSpectra(v_cmpinfo, annotated_spectra);
    }

    // convert sirius_output to mztab and store file
    const int candidates = algorithm.getNumberOfSiriusCandidates();
    MzTab sirius_result;
    MzTabFile siriusfile;
    SiriusMzTabWriter::read(subdirs, in, candidates, sirius_result);
    siriusfile.store(out_sirius, sirius_result);

    // convert sirius_output to mztab and store file
    if (!out_csifingerid.empty())
    {
      MzTab csi_result;
      MzTabFile csifile;
      CsiFingerIdMzTabWriter::read(subdirs, in, candidates, csi_result);
      csifile.store(out_csifingerid, csi_result);
    }

    // should the sirius workspace be retained
    if (!sirius_workspace_directory.empty())
    {
      // convert path to absolute path
      QDir sw_dir(sirius_workspace_directory.toQString());
      sirius_workspace_directory = String(sw_dir.absolutePath());
      
      // move tmp folder to new location
      bool copy_status = File::copyDirRecursively(sirius_tmp.getTmpDir().toQString(), sirius_workspace_directory.toQString());
      if (copy_status)
      { 
        OPENMS_LOG_INFO << "Sirius workspace was successfully copied to " << sirius_workspace_directory << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "Sirius workspace could not be copied to " << sirius_workspace_directory << ". Please run SiriusAdapter with debug >= 2." << std::endl;
      }
    }
   
    // should the ms file be retained (non-converter mode)
    if (!out_ms.empty())
    {  
      QFile::copy(sirius_tmp.getTmpMsFile().toQString(), out_ms.toQString());
      OPENMS_LOG_INFO << "Preprocessed .ms files were moved to " << out_ms << std::endl;
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPSiriusAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond
