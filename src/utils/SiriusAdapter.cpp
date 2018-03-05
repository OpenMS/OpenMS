// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <QtCore/QProcess>
#include <QDir>
#include <QDebug>
#include <QDirIterator>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <fstream>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_SiriusAdapter SiriusAdapter

  @brief De novo metabolite identification.

  CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

  To use this feature, the Sirius command line tool as well as a java installation is needed.

  Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/ 

  Please use Sirius Version 4.0.

  If you want to use the software with the Gurobi solver (free academic license) instead of GLPK, please follow the instructions in the sirius manual.

  Internal procedure in SiriusAdpater
  1. Input mzML
  2. Parsed by SiriusMSConverter into (sirius internal) .ms format
  3. Submission of .ms and additional parameters to wrapped SIRIUS.jar
  4. Sirius output saved in interal temporary folder structure
  5. Sirius output is parsed (SiriusMzTabWriter/CsiFingerIDMzTabWriter)
  6. Merge corresponding output in one mzTab (out_sirius/out_fingerid)

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
    }){}

protected:

  static bool sortByScanIndex(const String & i, const String & j)
  {
    return (atoi(SiriusMzTabWriter::extract_scan_index(i).c_str()) < atoi(SiriusMzTabWriter::extract_scan_index(j).c_str()));
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("executable", "<executable>", "",
                       "sirius executable e.g. sirius", false, false, ListUtils::create<String>("skipexists"));

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerInputFile_("in_adductinfo", "<file>", "", "FeatureXML Input with adduct information", false);
    setValidFormats_("in_adductinfo", ListUtils::create<String>("featurexml"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab Output file for SiriusAdapter results");
    setValidFormats_("out_sirius", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_fingerid","<file>", "", "MzTab output file for CSI:FingerID, if this parameter is given, SIRIUS will search for a molecular structure using CSI:FingerID after determining the sum formula", false);
    setValidFormats_("out_fingerid", ListUtils::create<String>("tsv"));

    registerStringOption_("profile", "<choice>", "qtof", "Specify the used analysis profile", false);
    setValidStrings_("profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
    registerIntOption_("candidates", "<num>", 5, "The number of candidates in the output.", false);
    registerStringOption_("database", "<choice>", "all", "search formulas in given database", false);
    setValidStrings_("database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));    
    registerIntOption_("noise", "<num>", 0, "median intensity of noise peaks", false);
    registerIntOption_("ppm_max", "<num>", 10, "allowed ppm for decomposing masses", false);
    registerStringOption_("isotope", "<choice>", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring. Use 'omit' to ignore isotope pattern.", false);
    setValidStrings_("isotope", ListUtils::create<String>("score,filter,both,omit"));
    registerStringOption_("elements", "<choice>", "CHNOP[5]S", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1].", false);
    registerIntOption_("tree_timeout", "<num>", 10, "Time out in seconds per fragmenation tree computation. To disable the tree timout set the value to 0", false);
    registerIntOption_("top_n_hits", "<num>", 10, "The top_n_hit for each compound written to the output", false);

    registerFlag_("auto_charge", "Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.", false);
    registerFlag_("ion_tree", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", false);
    registerFlag_("no_recalibration", "If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", false);
    registerFlag_("most_intense_ms2", "Sirius uses the fragmentation sepctrum with the most intense precursor peak (for each spectrum)", false);
  }

  // extract adduct information from featureXML (MetaboliteAdductDecharger)
  void extractAdductInformation(const PeakMap & spectra, const FeatureMap & feature_map, map<size_t, StringList> & map_precursor_to_adducts)
  {
    KDTreeFeatureMaps adduct_map_kd;
    vector<FeatureMap> adduct_map;
    adduct_map.push_back(feature_map);
    adduct_map_kd.addMaps(adduct_map);
    
    // map precursors to closest feature and retrieve annotated adducts
    for (size_t index = 0; index != spectra.size(); ++index)
    {
      if (spectra[index].getMSLevel() != 2) { continue; }
      
      // get precursor meta data (m/z, rt)
      const vector<Precursor> & pcs = spectra[index].getPrecursors();

      if (!pcs.empty())
      {
        const double mz = pcs[0].getMZ();
        const double rt = spectra[index].getRT();

        // query features in tolerance window
        vector<Size> matches;
        adduct_map_kd.queryRegion(rt - 5.0, rt + 5.0, mz - 0.2, mz + 0.2, matches, true);

        // no adduct information found - will use defaults in SIRIUS 
        if (matches.empty()) { continue; }
        
        // in the case of multiple features in tolerance window, select the one closest in m/z to the precursor
        Size min_distance_feature_index(0);
        double min_distance(1e11);
        for (auto const & k_idx : matches)
        {
          const double f_mz = adduct_map_kd.mz(k_idx);
          const double distance = fabs(f_mz - mz);
          if (distance < min_distance)
          {
            min_distance = distance;
            min_distance_feature_index = k_idx;
          }
        }
        const BaseFeature * min_distance_feature = adduct_map_kd.feature(min_distance_feature_index);
        
        // extract adducts from featureXML and associate with precursor
        if (min_distance_feature->metaValueExists("adducts"))
        {
          StringList adducts = min_distance_feature->getMetaValue("adducts");
          map_precursor_to_adducts[index] = adducts;
        }
      }
    }
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out_sirius = getStringOption_("out_sirius");
    String out_csifingerid = getStringOption_("out_fingerid");
    String adductinfo = getStringOption_("in_adductinfo");

    // needed for counting
    int top_n_hits = getIntOption_("top_n_hits"); 

    // Parameter for Sirius3
    QString executable = getStringOption_("executable").toQString();
    const QString profile = getStringOption_("profile").toQString();
    const QString elements = getStringOption_("elements").toQString();
    const QString database = getStringOption_("database").toQString();
    const QString isotope = getStringOption_("isotope").toQString();
    const QString noise = QString::number(getIntOption_("noise"));
    const QString ppm_max = QString::number(getIntOption_("ppm_max"));
    const QString candidates = QString::number(getIntOption_("candidates"));
    const QString tree_timeout = QString::number(getIntOption_("tree_timeout"));

    bool auto_charge = getFlag_("auto_charge");
    bool no_recalibration = getFlag_("no_recalibration");
    bool ion_tree = getFlag_("ion_tree");
    bool most_intense_ms2 = getFlag_("most_intense_ms2");

    //-------------------------------------------------------------
    // Determination of the Executable
    //-------------------------------------------------------------

    // Parameter executable not provided
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
    // Normalize file path
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

    // Read FeatureXML in KDTree for range query
    map<size_t, StringList> map_precursor_to_adducts;
    std::ifstream afile(adductinfo);    
    if (afile)
    {
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(adductinfo, feature_map);
      extractAdductInformation(spectra, feature_map, map_precursor_to_adducts);
    }   

    // Write msfile
    SiriusMSFile::store(spectra, tmp_ms_file, map_precursor_to_adducts);

    // Assemble SIRIUS parameters
    QStringList process_params;
    process_params << "-p" << profile
                   << "-e" << elements
                   << "-d" << database
                   << "-s" << isotope
                   << "--noise" << noise
                   << "--candidates" << candidates
                   << "--ppm-max" << ppm_max
                   << "--tree-timeout" << tree_timeout 
                   << "--quiet"
                   << "--output" << out_dir.toQString(); //internal output folder for temporary SIRIUS output file storage

    // Add flags
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

    // The actual process
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

    if (success == false || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_( "FATAL: External invocation of Sirius failed. Standard output and error were:");
      const QString sirius_stdout(qp.readAllStandardOutput());
      const QString sirius_stderr(qp.readAllStandardOutput());
      writeLog_(sirius_stdout);
      writeLog_(sirius_stderr);
      writeLog_(String(qp.exitCode()));

      return EXTERNAL_PROGRAM_ERROR;
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //Extract path to subfolders (sirius internal folder structure)
    QDirIterator it(out_dir.toQString(), QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::NoIteratorFlags);
    while (it.hasNext())
    {
      subdirs.push_back(it.next());
    }

    //sort vector path list
    std::sort(subdirs.begin(), subdirs.end(), sortByScanIndex);

    //Convert sirius_output to mztab and store file
    MzTab sirius_result;
    MzTabFile siriusfile;
    SiriusMzTabWriter::read(subdirs, in, top_n_hits, sirius_result);
    siriusfile.store(out_sirius, sirius_result);

    //Convert sirius_output to mztab and store file
    if (out_csifingerid.empty() == false)
    {
      MzTab csi_result;
      MzTabFile csifile;
      CsiFingerIdMzTabWriter::read(subdirs, in, top_n_hits, csi_result);
      csifile.store(out_csifingerid, csi_result);
    }

    //clean tmp directory if debug level < 2
    if (debug_level_ >= 2)
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
