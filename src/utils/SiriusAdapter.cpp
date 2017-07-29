// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

//not sure if more #include directives are needed

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/JavaInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>
#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>
#include <OpenMS/FORMAT/DATAACCESS/CsiFingerIdMzTabWriter.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>
#include <QDebug>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>

#include <QDirIterator>
#include <regex>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_SiriusAdapter

  @brief De novo metabolite identification.

  CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

  To use this feature, the Sirius command line tool as well as a java installation is needed.

  Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/

  If you want to use the software with the Gurobi solver (free academic license) instead of GLPK, please follow the instructions in the sirius manual.

  Please see the following publications:

  Kai Dührkop and Sebastian Böcker. Fragmentation trees reloaded.  J Cheminform, 8:5, 2016. (Cite this for fragmentation pattern analysis and fragmentation tree computation)

  Kai Dührkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Böcker. Searching molecular structure databases with tandem mass spectra using CSI:FingerID. Proc Natl Acad Sci U S A, 112(41):12580-12585, 2015. (Cite this when using CSI:FingerID)

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
    TOPPBase("SiriusAdapter", "Tool for metabolite identification using single and tandem mass spectrometry", false)
  {
  }

protected:

  static bool sortByScanIndex(const String & i, const String & j)
  {
    return (atoi(SiriusMzTabWriter::extract_scan_index(i).c_str()) < atoi(SiriusMzTabWriter::extract_scan_index(j).c_str()));
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("executable", "<executable>",
#if  defined(__APPLE__)
                       "sirius",
#else
                       "sirius-console-64.exe",
#endif
                       "sirius executable e.g. sirius", true, false, ListUtils::create<String>("skipexists"));

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab Output file for SiriusAdapter results");
    setValidFormats_("out_sirius", ListUtils::create<String>("csv"));

    registerOutputFile_("out_CSIFingerID","<file>", "", "MzTab ouput file for CSI:FingerID", false);
    setValidFormats_("out_CSIFingerID", ListUtils::create<String>("csv"));

    registerStringOption_("profile", "<choice>", "qtof", "Specify the used analysis profile", false);
    setValidStrings_("profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
    registerIntOption_("candidates", "<num>", 5, "The number of candidates in the output. Default 5 best candidates", false);
    registerStringOption_("database", "<choice>", "all", "search formulas in given database", false);
    setValidStrings_("database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));
    registerIntOption_("noise", "<num>", 0, "median intensity of noise peaks", false);
    registerIntOption_("ppm_max", "<num>", 10, "allowed ppm for decomposing masses", false);
    registerStringOption_("isotope", "<choice>", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring (default). Use 'omit' to ignore isotope pattern.", false);
    setValidStrings_("isotope", ListUtils::create<String>("score,filter,both,omit"));
    registerStringOption_("elements", "<choice>", "CHNOP[5]S", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurence of these elements: CHNOP[5]S[8]Cl[1]. By default CHNOP[5]S is used.", false);

    registerIntOption_("number", "<num>", 10, "The number of compounds used in the output", false);

    registerFlag_("auto_charge", "Use this option if the charge of your compounds is unknown and you do not want to assume [M+H]+ as default. With the auto charge option SIRIUS will not care about charges and allow arbitrary adducts for the precursor peak.", false);
    registerFlag_("iontree", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", false);
    registerFlag_("no_recalibration", "If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", false);
    registerFlag_("fingerid", "If this option is set, SIRIUS will search for a molecular structure using CSI:FingerID after determining the sum formula", false);
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out_sirius = getStringOption_("out_sirius");
    String out_csifingerid = getStringOption_("out_CSIFingerID");

    // needed for counting
    int number_compounds = getIntOption_("number") + 1;  // +1 needed to write the correct number of compounds

    // Parameter for Sirius3
    QString executable = getStringOption_("executable").toQString();
    QString profile = getStringOption_("profile").toQString();
    QString elements = getStringOption_("elements").toQString();
    QString database = getStringOption_("database").toQString();
    QString isotope = getStringOption_("isotope").toQString();
    QString noise = QString::number(getIntOption_("noise"));
    QString ppm_max = QString::number(getIntOption_("ppm_max"));
    QString candidates = QString::number(getIntOption_("candidates"));

    QString path_to_executable = File::path(getStringOption_("executable")).toQString();

    bool auto_charge = getFlag_("auto_charge");
    bool no_recalibration = getFlag_("no_recalibration");
    bool fingerid = getFlag_("fingerid");
    bool iontree = getFlag_("iontree");

    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------

    PeakMap spectra;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, spectra);
    std::vector<String> subdirs;

    //Write msfile
    String tmp_dir = QDir::toNativeSeparators(String(File::getTempDirectory()).toQString()) + "/" + String(File::getUniqueName()).toQString();
    String tmp_ms_file = tmp_dir + "/" + "msfile";
    String out_dir = tmp_dir + "/" + "sirius_out";
    SiriusMSFile::store(spectra, tmp_ms_file);

    //Knime hack
    QProcessEnvironment env;
    String siriuspath = "SIRIUS_PATH";
    QString qsiriuspath = env.systemEnvironment().value(siriuspath.toQString());

    if (!qsiriuspath.isEmpty())
    {
      executable = qsiriuspath;
    }

    //Start Sirius
    QStringList process_params; // the actual process
    process_params << "-p" << profile
                   << "-e" << elements
                   << "-d" << database
                   << "-s" << isotope
                   << "--noise" << noise
                   << "--candidates" << candidates
                   << "--ppm-max" << ppm_max
                   << "--quiet"
                   << "--output" << out_dir.toQString(); //internal output folder for temporary

    if (no_recalibration)
    {
      process_params << "--no-recalibration";
    }
    if (fingerid)
    {
      process_params << "--fingerid";
    }
    if (iontree)
    {
      process_params << "--iontree";
    }
    if (auto_charge)
    {
      process_params << "--auto-charge";
    }

    process_params << tmp_ms_file.toQString();

    QProcess qp;
    qp.workingDirectory();
    qp.setWorkingDirectory(path_to_executable); //since library paths are relative to sirius executable path
    qp.start(executable, process_params); // does automatic escaping etc... start
    bool success = qp.waitForFinished(-1); // wait till job is finished
    qp.close();

    if ( ! success || qp.exitStatus() != 0 || qp.exitCode() != 0)
    {
      writeLog_( "Fatal error: Running SiriusAdapter returned an error code or could no compute the input" );
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
    SiriusMzTabWriter::read(subdirs, number_compounds, sirius_result);
    siriusfile.store(out_sirius, sirius_result);

    //Convert sirius_output to mztab and store file
    if (!out_csifingerid.empty() && fingerid)
    {
      MzTab csi_result;
      MzTabFile csifile;
      CsiFingerIdMzTabWriter::read(subdirs, number_compounds, csi_result);
      csifile.store(out_csifingerid, csi_result);
    }

    //clean tmp directory if debug level < 2
    if (debug_level_ >= 2)
    {
      writeDebug_("Keeping temporary files in directory '" + tmp_dir + " and msfile at this location "+ tmp_ms_file + ". Set debug level to 1 or lower to remove them.", 2);
    }
    else
    {
      if ( ! tmp_dir.empty())
      {
        writeDebug_("Deleting temporary directory '" + tmp_dir + "'. Set debug level to 2 or higher to keep it.", 0);
        File::removeDir(tmp_dir.toQString());
      }
      if ( ! tmp_ms_file.empty() )
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
