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
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <OpenMS/ANALYSIS/ID/SiriusAdapterAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h> // store
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h> // read
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h> // read

#include <QtCore/QProcess>
#include <QDir>
#include <QDebug>
#include <QDirIterator>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <fstream>

using namespace OpenMS;
using namespace std;

//----------------------------------------------------------
// Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_SiriusAdapter SiriusAdapter

  @brief De novo metabolite identification.

  CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

  To use this feature, the Sirius command line tool as well as a java installation is needed.

  Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/ 

  Please use Sirius Version 4.0.1

  If you want to use the software with the Gurobi solver or CPLEX instead of GLPK, please follow the instructions in the sirius manual.

  Internal procedure in SiriusAdpater \n
  1. Input mzML (and optional featureXML) \n
  2. Preprocessing (see below)\n
  3. Parsed by SiriusMSConverter into (sirius internal) .ms format \n
  4. Submission of .ms and additional parameters to wrapped SIRIUS.jar \n
  5. Sirius output saved in interal temporary folder structure \n
  6. Sirius output is parsed (SiriusMzTabWriter/CsiFingerIDMzTabWriter) \n
  7. Merge corresponding output in one mzTab (out_sirius/out_fingerid) \n

  Preprocessing (featureXML): 
  By providing a featureXML, the feature information can be used for feature mapping.
  Sirius will then process the internally merged MS2 spectra allocated to one feature (instead of all available MS2).
  To reduce the feature space even further a masstrace filter can be set. 
  Additional adduct information can be provided using a featureXML from the MetaboliteAdductDecharger or AccurateMassSearch.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES

class TOPPSiriusAdapter :
 public TOPPBase
{
 public:
  TOPPSiriusAdapter() :
    TOPPBase("SiriusAdapter", "Tool for metabolite identification using single and tandem mass spectrometry", false,
      {
        {"Kai Dührkop and Sebastian Böcker",
         "Fragmentation trees reloaded",
         "J Cheminform; 2016",
         "10.1186/s13321-016-0116-8"},
        {"Kai Dührkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Böcker",
         "Searching molecular structure databases with tandem mass spectra using CSI:FingerID",
         "Proceedings of the National Academy of Sciences; 2015",
         "10.1073/pnas.1509788112"}
      })
    {}


protected:

  static bool extractAndCompareScanIndexLess_(const String& i, const String& j)
  {
    return (SiriusMzTabWriter::extract_scan_index(i) < SiriusMzTabWriter::extract_scan_index(j));
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("executable", "<executable>", "", "sirius executable e.g. sirius", false, false, ListUtils::create<String>("skipexists"));

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("in_featureinfo", "<file>", "", "FeatureXML input with feature and adduct information", false);
    setValidFormats_("in_featureinfo", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab Output file for SiriusAdapter results", false);
    setValidFormats_("out_sirius", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_fingerid","<file>", "", "MzTab output file for CSI:FingerID, if this parameter is given, SIRIUS will search for a molecular structure using CSI:FingerID after determining the sum formula", false);
    setValidFormats_("out_fingerid", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_ms","<file>", "", "Internal SIRIUS .ms format after OpenMS preprocessing", false);
    setValidFormats_("out_ms", ListUtils::create<String>("ms"));

    registerStringOption_("out_workspace_directory", "<directory>", "", "Output directory for SIRIUS workspace", false);

    registerFlag_("converter_mode", "Use this flag in combination with the out_ms file to only convert the input mzML and featureXML to an .ms file. Without further SIRIUS processing.", true);

    addEmptyLine_();
    registerFullParam_(SiriusAdapterAlgorithm().getDefaults());
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    // param SiriusAdapter
    String executable = getStringOption_("executable");
    String in = getStringOption_("in");
    String out_sirius = getStringOption_("out_sirius");
    String out_csifingerid = getStringOption_("out_fingerid");
    String featureinfo = getStringOption_("in_featureinfo");
    String out_ms = getStringOption_("out_ms");
    String sirius_workspace_directory = getStringOption_("out_workspace_directory");
    bool converter_mode = getFlag_("converter_mode");

    // param SiriusAdapterAlgorithm
    Param combined; 
    SiriusAdapterAlgorithm sirius_algo;
    Param preprocessing = getParam_().copy("preprocessing", false);
    Param sirius = getParam_().copy("sirius", false);
    combined.insert("", preprocessing);
    combined.insert("", sirius);
    sirius_algo.setParameters(combined);

    writeDebug_("Parameters passed to SiriusAdapterAlgorithm", combined, 3);

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------
    MSExperiment spectra;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, spectra);

    // make temporary files
    SiriusAdapterAlgorithm::SiriusTmpStruct sirius_tmp = SiriusAdapterAlgorithm::constructSiriusTmpStruct();
    String tmp_dir = sirius_tmp.tmp_dir;
    String tmp_ms_file = sirius_tmp.tmp_ms_file;
    String tmp_out_dir = sirius_tmp.tmp_out_dir;

    // run masstrace filter and feature mapping
    vector<FeatureMap> v_fp; // copy FeatureMap via push_back
    KDTreeFeatureMaps fp_map_kd; // reference to *basefeature in vector<FeatureMap>
    FeatureMapping::FeatureToMs2Indices feature_mapping; // reference to *basefeature in vector<FeatureMap>
    SiriusAdapterAlgorithm::preprocessingSirius(featureinfo,
                                                spectra,
                                                v_fp,
                                                fp_map_kd,
                                                sirius_algo,
                                                feature_mapping);

    // returns Log of feature and/or spectra number
    SiriusAdapterAlgorithm::checkFeatureSpectraNumber(featureinfo,
                                                      feature_mapping,
                                                      spectra,
                                                      sirius_algo);

    // write msfile and store the compound information in CompoundInfo Object
    vector<SiriusMSFile::CompoundInfo> v_cmpinfo;
    bool feature_only = (sirius_algo.getFeatureOnly() == "true") ? true : false;
    bool no_mt_info = (sirius_algo.getNoMasstraceInfoIsotopePattern() == "true") ? true : false;
    int isotope_pattern_iterations = sirius_algo.getIsotopePatternIterations();
    SiriusMSFile::store(spectra,
                        tmp_ms_file,
                        feature_mapping,
                        feature_only,
                        isotope_pattern_iterations,
                        no_mt_info,
                        v_cmpinfo);

    // converter_mode enabled (only needed for SiriusAdapter)
    if (!out_ms.empty() && converter_mode)
    {
      QFile::copy(tmp_ms_file.toQString(), out_ms.toQString());
      OPENMS_LOG_WARN << "SiriusAdapter was used in converter mode and is terminated after openms preprocessing. \n"
                  "If you would like to run SIRIUS internally please disable the converter mode." << std::endl; 
      return EXECUTION_OK;
    }

    // calls Sirius and returns vector of paths to sirius folder structure
    std::vector<String> subdirs;
    subdirs = SiriusAdapterAlgorithm::callSiriusQProcess(tmp_ms_file,
                                                         tmp_out_dir,
                                                         executable,
                                                         out_csifingerid,
                                                         sirius_algo);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // sort vector path list
    std::sort(subdirs.begin(), subdirs.end(), extractAndCompareScanIndexLess_);

    // convert sirius_output to mztab and store file
    int candidates = sirius_algo.getCandidates();
    MzTab sirius_result;
    MzTabFile siriusfile;
    SiriusMzTabWriter::read(subdirs, in, candidates, sirius_result);
    siriusfile.store(out_sirius, sirius_result);

    // convert sirius_output to mztab and store file
    if (!out_csifingerid.empty())
    {
      int top_n_hits = sirius_algo.getTopNHits();
      MzTab csi_result;
      MzTabFile csifile;
      CsiFingerIdMzTabWriter::read(subdirs, in, top_n_hits, csi_result);
      csifile.store(out_csifingerid, csi_result);
    }

    // should the sirius workspace be retained
    if (!sirius_workspace_directory.empty())
    {
      // convert path to absolute path
      QDir sw_dir(sirius_workspace_directory.toQString());
      sirius_workspace_directory = String(sw_dir.absolutePath());
      
      // move tmp folder to new location
      bool copy_status = File::copyDirRecursively(tmp_dir.toQString(), sirius_workspace_directory.toQString());
      if (copy_status)
      { 
        OPENMS_LOG_INFO << "Sirius Workspace was successfully copied to " << sirius_workspace_directory << std::endl;
      }
      else
      {
        OPENMS_LOG_INFO << "Sirius Workspace could not be copied to " << sirius_workspace_directory << ". Please run SiriusAdapter with debug >= 2." << std::endl;
      }
    }
   
    // should the ms file be retained (non-converter mode)
    if (!out_ms.empty())
    {  
      QFile::copy(tmp_ms_file.toQString(), out_ms.toQString());
      OPENMS_LOG_INFO << "Preprocessed .ms files was moved to " << out_ms << std::endl; 
    }

    // clean tmp directory if debug level < 2 
    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + tmp_dir + " and msfile at this location "+ tmp_ms_file + ". Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (tmp_dir.empty() == false)
      {
        writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 0);
        File::removeDir(tmp_dir.toQString());
      }
      if (tmp_ms_file.empty() == false)
      {
        writeDebug_("Deleting temporary msfile '" + tmp_ms_file + "'. Set debug level to 2 or higher to keep it.", 0);
        File::remove(tmp_ms_file); // remove msfile
      }
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
