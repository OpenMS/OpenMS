// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/SiriusExportAlgorithm.h>
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
  @page TOPP_SiriusExport SiriusExport

  @brief De novo metabolite identification.

  CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

  To use this feature, the Sirius command line tool as well as a java installation is needed.

  Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/ 

  Please use Sirius Version 4.0.1

  If you want to use the software with the Gurobi solver or CPLEX instead of GLPK, please follow the instructions in the sirius manual.

  Internal procedure in SiriusExport \n
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
  @verbinclude TOPP_SiriusExport.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_SiriusExport.html
 */

/// @cond TOPPCLASSES

class TOPPSiriusExport :
 public TOPPBase
{
 public:
  TOPPSiriusExport() :
    TOPPBase("SiriusExport", "Metabolite identification using single and tandem mass spectrometry", false,
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
  SiriusExportAlgorithm algorithm;

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("in_featureinfo", "<file>", "", "FeatureXML input with feature and adduct information", false);
    setValidFormats_("in_featureinfo", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_ms","<file>", "", "Internal SIRIUS .ms format after OpenMS preprocessing", false);
    setValidFormats_("out_ms", ListUtils::create<String>("ms"));

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
    String in = getStringOption_("in");
    String featureinfo = getStringOption_("in_featureinfo");
    String out_ms = getStringOption_("out_ms");

    auto params = getParam_();
    if (debug_level_ > 3)
    {
      params.setValue("read_sirius_stdout", "true");
    }
    params.setValue("project:processors", params.getValue("threads"));
    algorithm.updateExistingParameter(params);

    writeDebug_("Parameters passed to SiriusExportAlgorithm", algorithm.getParameters(), 3);

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------
    MSExperiment spectra;
    FileHandler().loadExperiment(in, spectra, {FileTypes::MZML}, log_type_);

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
                        out_ms.toQString(),
                        feature_mapping,
                        algorithm.isFeatureOnly(),
                        algorithm.getIsotopePatternIterations(),
                        algorithm.isNoMasstraceInfoIsotopePattern(),
                        v_cmpinfo);

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPSiriusExport tool;
  return tool.main(argc, argv);
}

/// @endcond