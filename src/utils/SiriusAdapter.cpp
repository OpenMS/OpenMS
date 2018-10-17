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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/SYSTEM/File.h>

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

  Please use Sirius Version 4.0.

  If you want to use the software with the Gurobi solver or CPLEX instead of GLPK, please follow the instructions in the sirius manual.

  Internal procedure in SiriusAdpater
  1. Input mzML (and optional featureXML)
  2. Parsed by SiriusMSConverter into (sirius internal) .ms format
  3. Submission of .ms and additional parameters to wrapped SIRIUS.jar
  4. Sirius output saved in interal temporary folder structure
  5. Sirius output is parsed (SiriusMzTabWriter/CsiFingerIDMzTabWriter)
  6. Merge corresponding output in one mzTab (out_sirius/out_fingerid)

  By providing a featureXML, the feature information can be used for feature mapping.
  Sirius will then process the mappend MS2 spectra (instead of all available MS2).
  If the featureXML provides additional adduct information (e.g. from the MetaboliteAdductDecharger)
  this can be used to speed the Sirius calculation.

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
    registerInputFile_("executable", "<executable>", "",
                       "sirius executable e.g. sirius", false, false, ListUtils::create<String>("skipexists"));

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerInputFile_("in_featureinfo", "<file>", "", "FeatureXML input with feature and adduct information", false);
    setValidFormats_("in_featureinfo", ListUtils::create<String>("featurexml"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab Output file for SiriusAdapter results", false);
    setValidFormats_("out_sirius", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_fingerid","<file>", "", "MzTab output file for CSI:FingerID, if this parameter is given, SIRIUS will search for a molecular structure using CSI:FingerID after determining the sum formula", false);
    setValidFormats_("out_fingerid", ListUtils::create<String>("mzTab"));

    registerOutputFile_("out_ms","<file>", "", "Internal SIRIUS .ms format after OpenMS preprocessing", false);
    registerStringOption_("sirius_workspace_directory","<directory>", "", "Output directory with SIRIUS workspace", false);


    // adapter parameters
    registerIntOption_("filter_by_num_masstraces", "<num>", 1, "Features have to have at least x MassTraces. To use this parameter feature_only is neccessary", false);
    setMinInt_("filter_by_num_masstraces", 1);
    registerFlag_("feature_only", "Uses the feature information from in_featureinfo to reduce the search space to only MS2 associated with a feature", false);
    registerDoubleOption_("precursor_mz_tolerance", "<num>", 0.005, "Tolerance window for precursor selection (Feature selection in regard to the precursor)", false);
    registerStringOption_("precursor_mz_tolerance_unit", "<choice>", "Da", "Unit of the precursor_mz_tolerance", false);
    setValidStrings_("precursor_mz_tolerance_unit", ListUtils::create<String>("Da,ppm"));
    registerDoubleOption_("precursor_rt_tolerance", "<num>", 5, "Tolerance window (left and right) for precursor selection [seconds]", false);
    registerIntOption_("isotope_pattern_iterations", "<num>", 3, "Number of iterations that should be performed to extract the C13 isotope pattern. If no peak is found (C13 distance) the function will abort. Be careful with noisy data - since this can lead to wrong isotope patterns.", false, true);
    registerFlag_("no_masstrace_info_isotope_pattern", "Use this flag if the masstrace information from a feature should be discarded and the isotope_pattern_iterations should be used instead.", true);
    registerFlag_("converter_mode", "Use this flag in combination with the out_ms file to only convert the input mzML and featureXML to an .ms file. Without further SIRIUS processing.", true);

    // internal sirius parameters
    registerStringOption_("profile", "<choice>", "qtof", "Specify the used analysis profile", false);
    setValidStrings_("profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
    registerIntOption_("candidates", "<num>", 5, "The number of candidates in the SIRIUS output.", false);
    registerStringOption_("database", "<choice>", "all", "search formulas in given database", false);
    setValidStrings_("database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));
    registerIntOption_("noise", "<num>", 0, "median intensity of noise peaks", false);
    registerIntOption_("ppm_max", "<num>", 10, "allowed ppm for decomposing masses", false);
    registerStringOption_("isotope", "<choice>", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring. Use 'omit' to ignore isotope pattern.", false);
    setValidStrings_("isotope", ListUtils::create<String>("score,filter,both,omit"));
    registerStringOption_("elements", "<choice>", "CHNOP[5]S[8]Cl[1]", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1].", false);
    registerIntOption_("compound_timeout", "<num>", 10, "Time out in seconds per compound. To disable the timeout set the value to 0", false);
    registerIntOption_("tree_timeout", "<num>", 0, "Time out in seconds per fragmentation tree computation.", false);
    registerIntOption_("top_n_hits", "<num>", 10, "The number of top hits for each compound written to the CSI:FingerID output", false);

    registerFlag_("auto_charge", "Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.", false);
    registerFlag_("ion_tree", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", false);
    registerFlag_("no_recalibration", "If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", false);
    registerFlag_("most_intense_ms2", "SIRIUS uses the fragmentation spectrum with the most intense precursor peak (for each spectrum)", false);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out_sirius = getStringOption_("out_sirius");
    String out_csifingerid = getStringOption_("out_fingerid");
    String featureinfo = getStringOption_("in_featureinfo");

    String out_ms = getStringOption_("out_ms");
    String sirius_workspace_directory = getStringOption_("sirius_workspace_directory");

    bool converter_mode = getFlag_("converter_mode");

    // parameter for SiriusAdapter
    bool feature_only = getFlag_("feature_only");
    unsigned int num_masstrace_filter = getIntOption_("filter_by_num_masstraces");
    if (num_masstrace_filter != 1 && !feature_only)
    {
      num_masstrace_filter = 1;
      LOG_WARN << "Parameter: filter_by_num_masstraces, was set to 1 to retain the adduct information for all MS2 spectra, if available. Please use the masstrace filter in combination with feature_only." << endl;
    }

    double precursor_mz_tol = getDoubleOption_("precursor_mz_tolerance");
    String unit_prec = getStringOption_("precursor_mz_tolerance_unit");
    bool ppm_prec = unit_prec == "ppm" ? true : false;
    double precursor_rt_tol = getDoubleOption_("precursor_rt_tolerance");
    int isotope_pattern_iterations = getIntOption_("isotope_pattern_iterations");
    bool no_mt_info = getFlag_("no_masstrace_info_isotope_pattern");

    // needed for counting
    int top_n_hits = getIntOption_("top_n_hits");

    // parameter for Sirius
    QString executable = getStringOption_("executable").toQString();
    const QString profile = getStringOption_("profile").toQString();
    const QString elements = getStringOption_("elements").toQString();
    const QString database = getStringOption_("database").toQString();
    const QString isotope = getStringOption_("isotope").toQString();
    const QString noise = QString::number(getIntOption_("noise"));
    const QString ppm_max = QString::number(getIntOption_("ppm_max"));
    const Size candidates = getIntOption_("candidates");
    const QString compound_timeout = QString::number(getIntOption_("compound_timeout"));
    const QString tree_timeout = QString::number(getIntOption_("tree_timeout"));

    bool auto_charge = getFlag_("auto_charge");
    bool no_recalibration = getFlag_("no_recalibration");
    bool ion_tree = getFlag_("ion_tree");
    bool most_intense_ms2 = getFlag_("most_intense_ms2");

    //-------------------------------------------------------------
    // Determination of the Executable
    //-------------------------------------------------------------

    // parameter executable not provided
    if (executable.isEmpty())
    {
      const QProcessEnvironment env;
      const QString & qsiriuspathenv = env.systemEnvironment().value("SIRIUS_PATH");
      if (qsiriuspathenv.isEmpty())
      {
        writeLog_( "FATAL: Executable of Sirius could not be found. Please either use SIRIUS_PATH env variable or provide with -executable");
        return MISSING_PARAMETERS;
      }
      executable = qsiriuspathenv;
    }
    // normalize file path
    QFileInfo file_info(executable);
    executable = file_info.canonicalFilePath();

    writeLog_("Executable is: " + executable);
    const QString & path_to_executable = File::path(executable).toQString();

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------
    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, spectra);
    std::vector<String> subdirs;

    QString tmp_base_dir = File::getTempDirectory().toQString();
    QString tmp_dir = QDir(tmp_base_dir).filePath(File::getUniqueName().toQString());

    String tmp_ms_file = QDir(tmp_base_dir).filePath((File::getUniqueName() + ".ms").toQString());
    String out_dir = QDir(tmp_dir).filePath("sirius_out");

    FeatureMapping::FeatureToMs2Indices feature_mapping;
    FeatureMap feature_map;
    KDTreeFeatureMaps fp_map_kd;
    vector<FeatureMap> v_fp;

    // if fileparameter is given and should be not empty
    if (!featureinfo.empty())
    {
      if (File::exists(featureinfo) && !File::empty(featureinfo))
      {
        // read featureXML
        FeatureXMLFile fxml;
        fxml.load(featureinfo, feature_map);

        // filter feature by number of masstraces
        auto map_it = remove_if(feature_map.begin(), feature_map.end(),
                                [&num_masstrace_filter](const Feature &feat) -> bool
                                {
                                  unsigned int n_masstraces = feat.getMetaValue("num_of_masstraces");
                                  return n_masstraces < num_masstrace_filter;
                                });
        feature_map.erase(map_it, feature_map.end());

        v_fp.push_back(feature_map);
        fp_map_kd.addMaps(v_fp);

        // mapping of MS2 spectra to features
        feature_mapping = FeatureMapping::assignMS2IndexToFeature(spectra,
                                                                  fp_map_kd,
                                                                  precursor_mz_tol,
                                                                  precursor_rt_tol,
                                                                  ppm_prec);
      }
      else
      {
        throw OpenMS::Exception::FileEmpty(__FILE__,
                                           __LINE__,
                                           __FUNCTION__,
                                           "Error: FeatureXML was empty, please provide a valid file.");
      }
    }

    // write msfile
    SiriusMSFile::store(spectra, tmp_ms_file, feature_mapping, feature_only, isotope_pattern_iterations, no_mt_info);

    // converter_mode enabled 
    if (!out_ms.empty() && converter_mode)
    {
      QFile::copy(tmp_ms_file.toQString(), out_ms.toQString());
      LOG_WARN << "SiriusAdapter was used in converter mode and is terminated after openms preprocessing. \n"
                  "If you would like to run SIRIUS internally please disable the converter mode." << std::endl; 
      return EXECUTION_OK;
    }

    // assemble SIRIUS parameters
    QStringList process_params;
    process_params << "-p" << profile
                   << "-e" << elements
                   << "-d" << database
                   << "-s" << isotope
                   << "--noise" << noise
                   << "--candidates" << QString::number(candidates)
                   << "--ppm-max" << ppm_max
                   << "--compound-timeout" << compound_timeout
                   << "--tree-timeout" << tree_timeout 
                   << "--quiet"
                   << "--output" << out_dir.toQString(); //internal output folder for temporary SIRIUS output file storage

    // add flags
    if (no_recalibration)
    {
      process_params << "--no-recalibration";
    }
    if (!out_csifingerid.empty())
    {
      process_params << "--fingerid";
    }
    if (ion_tree)
    {
      process_params << "--iontree";
    }
    if (auto_charge)
    {
      process_params << "--auto-charge";
    }
    if (most_intense_ms2)
    {
      process_params << "--mostintense-ms2";
    }

    process_params << tmp_ms_file.toQString();

    // the actual process
    QProcess qp;
    qp.setWorkingDirectory(path_to_executable); //since library paths are relative to sirius executable path
    qp.start(executable, process_params); // does automatic escaping etc... start
    std::stringstream ss;
    ss << "COMMAND: " << executable.toStdString();
    for (QStringList::const_iterator it = process_params.begin(); it != process_params.end(); ++it)
    {
        ss << " " << it->toStdString();
    }
    LOG_DEBUG << ss.str() << endl;
    writeLog_("Executing: " + String(executable));
    writeLog_("Working Dir is: " + path_to_executable);
    const bool success = qp.waitForFinished(-1); // wait till job is finished
    qp.close();

    if (!success || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_( "FATAL: External invocation of Sirius failed. Standard output and error were:");
      const QString sirius_stdout(qp.readAllStandardOutput());
      const QString sirius_stderr(qp.readAllStandardError());
      writeLog_(sirius_stdout);
      writeLog_(sirius_stderr);
      writeLog_(String(qp.exitCode()));

      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    // extract path to subfolders (sirius internal folder structure)
    QDirIterator it(out_dir.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
    while (it.hasNext())
    {
      subdirs.push_back(it.next());
    }

    // sort vector path list
    std::sort(subdirs.begin(), subdirs.end(), extractAndCompareScanIndexLess_);

    // convert sirius_output to mztab and store file
    MzTab sirius_result;
    MzTabFile siriusfile;
    SiriusMzTabWriter::read(subdirs, in, candidates, sirius_result);
    siriusfile.store(out_sirius, sirius_result);

    // convert sirius_output to mztab and store file
    if (!out_csifingerid.empty())
    {
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
      
      // try to create directory if not present
      if (!sw_dir.exists())
      {
        sw_dir.mkpath(sirius_workspace_directory.toQString());
      }
      
      // move tmp folder to new location
      std::rename(tmp_dir.toStdString().c_str(), sirius_workspace_directory.c_str());
      LOG_WARN << "Sirius Workspace was moved to " << sirius_workspace_directory << std::endl;
    }
   
    // should the ms file be retained (non-converter mode)
    if (!out_ms.empty())
    {  
      QFile::copy(tmp_ms_file.toQString(), out_ms.toQString());
      LOG_WARN << "Preprocessed .ms files was moved to " << out_ms << std::endl; 
    }


    // clean tmp directory if debug level < 2 
    // if out_ms and sirius_workspace_directoy is set - the files/folders have already be moved to 
    // the designated location
    if (debug_level_ >= 2 && out_ms.empty() && sirius_workspace_directory.empty())
    {
      writeDebug_("Keeping temporary files in directory '" + String(tmp_dir) + " and msfile at this location "+ tmp_ms_file + ". Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if (tmp_dir.isEmpty() == false)
      {
        writeDebug_("Deleting temporary directory '" + String(tmp_dir) + "'. Set debug level to 2 or higher to keep it.", 0);
        File::removeDir(tmp_dir);
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
