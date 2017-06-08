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
//#include <OpenMS/?/SiriusMzTabWriter.h>
//#include <OpenMS/?/CsiFingerIDMzTabWriter.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>
#include <OpenMS/ANALYSIS/ID/SiriusMSConverter.h>


using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
@page UTILS_SiriusAdapter

@brief Metabolite identification using single and tandem mass spectrometry.

CSI:FingerID (Compound Structure Identification: FingerID) is a method for searching a tandem mass spectrum of a small molecule (metabolite) in a database of molecular structures.

To use this feature, the Sirius command line tool as well as a java installation is needed.

Sirius can be found on https://bio.informatik.uni-jena.de/software/sirius/

If you want to use the software with the Gurobi solver (free academic license) instead of GLPK, please follow the instructions in the sirius manual.

Please see the following publications:

Kai Dührkop and Sebastian Böcker. Fragmentation trees reloaded.  J Cheminform, 8:5, 2016. (Cite this for fragmentation pattern analysis and fragmentation tree computation)

Kai Dührkop, Huibin Shen, Marvin Meusel, Juho Rousu, and Sebastian Böcker. Searching molecular structure databases with tandem mass spectra using CSI:FingerID. Proc Natl Acad Sci U S A, 112(41):12580-12585, 2015. (cite this when using CSI:FingerID)

Sebastian Böcker, Matthias C. Letzel, Zsuzsanna Lipták and Anton Pervukhin. SIRIUS: decomposing isotope patterns for metabolite identification. Bioinformatics (2009) 25 (2): 218-224. (Cite this for isotope pattern analysis)

Florian Rasche, Aleš Svatoš, Ravi Kumar Maddula, Christoph Böttcher, and Sebastian Böcker. Computing Fragmentation Trees from Tandem Mass Spectrometry Data. Analytical Chemistry (2011) 83 (4): 1243–1251. (Cite this for introduction of fragmentation trees as used by SIRIUS)

Sebastian Böcker and Florian Rasche. Towards de novo identification of metabolites by analyzing tandem mass spectra. Bioinformatics (2008) 24 (16): i49-i55. (The very first paper to mention fragmentation trees as used by SIRIUS)

<B>The command line parameters of this tool are:</B>
@verbinclude UTILS_SiriusAdapter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude UTILS_SiriusAdapter.html
*/

// We do not want this class to show up in the docu:
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

  struct SiriusAdapterHit
  {
    String inchikey2D;
    String inchi;
    unsigned int rank;
    double score;
    String name;
    String smiles;
    vector<String> pubchemids;
    String links;
  };

  struct SiriusAdapterIdentification
  {
    String id;
    unsigned int scan_index;
    int charge;
    int mz;
    int rt;
    vector<SiriusAdapterHit> hits;
  };

  struct SiriusAdapterRun
  {
    vector<SiriusAdapterIdentification> identifications;
  };

  void registerOptionsAndFlags_()
  {
    registerInputFile_("executable", "<file>", "", "sirius file e.g. /bin/sirius3", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerOutputFile_("out_sirius", "<file>", "", "MzTab Output file for SiriusAdapter results");
    setValidFormats_("out_sirius", ListUtils::create<String>("csv"));

    registerOutputFile_("out_CsiFingerID","<file>", "", "MzTab ouput file for CsiFingerID", false);
    setValidFormats_("out_CsiFingerID", ListUtils::create<String>("csv"));

    registerStringOption_("analysis_profile", "<choice>", "qtof", "Specify the used analysis profile", false);
    setValidStrings_("analysis_profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
    registerIntOption_("candidates", "<num>", 5, "The number of candidates in the output. Default 5 best candidates", false);
    registerStringOption_("database", "<choice>", "all", "search formulas in given database", false);
    setValidStrings_("database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));
    registerIntOption_("noise", "<num>", 0, "median intensity of noise peaks", false);
    registerIntOption_("ppm_max", "<num>", 10, "allowed ppm for decomposing masses", false);
    registerStringOption_("isotope", "<choice>", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring (default). Use 'omit' to ignore isotope pattern.", false);
    setValidStrings_("isotope", ListUtils::create<String>("score,filter,both,omit"));
    registerStringOption_("elements", "<choice>", "CHNOP[5]S", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurence of these elements: CHNOP[5]S[8]Cl[1]. By default CHNOP[5]S is used.", false);
    registerIntOption_("mass_deviation", "<num>", 5, "Specify the allowed mass deviation of the fragment peak in ppm.", false);

    registerIntOption_("batch_size", "<num>", 0, "Number of files in one .ms file (Only needed if Rest-query is used)", false);

    registerIntOption_("number", "<num>", 10, "The number of compounds used in the output", false);

    registerFlag_("iontree", "Print molecular formulas and node labels with the ion formula instead of the neutral formula", false);
    registerFlag_("no_recalibration", "If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", false);
    registerFlag_("fingerid", "If this option is set, SIRIUS will search for molecular structure using CSI:FingerId after determining the molecIf this option is set, SIRIUS will search for molecular structure using CSI:FingerId after determining the molecular formulular formula", false);
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String out1 = getStringOption_("out_sirius");
    String out2 = getStringOption_("out_CsiFingerID");

    // needed for counting
    int number = getIntOption_("number");
    number = number + 1;

    // Parameter for Sirius3
    QString executable = getStringOption_("executable").toQString();
    QString analysis_profile = getStringOption_("analysis_profile").toQString();
    QString elements = getStringOption_("elements").toQString();
    QString database = getStringOption_("database").toQString();
    QString isotope = getStringOption_("isotope").toQString();
    QString noise = QString::number(getIntOption_("noise"));
    QString ppm_max = QString::number(getIntOption_("ppm_max"));
    QString candidates = QString::number(getIntOption_("candidates"));

    Int batch_size = getIntOption_("batch_size");
    if(batch_size == 0)
    {
      batch_size = static_cast<int>(in.size());
    }

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

    std::vector<String> ms_files = SiriusMSFile::store(spectra, batch_size);

    for(int i = 0; i < ms_files.size(); ++i)
    {

      QStringList process_params; // the actual process
      process_params << "-p" << analysis_profile
                     << "-e" << elements
                     << "-d" << database
                     << "-s" << isotope
                     << "--noise" << noise
                     << "--candidates" << candidates
                     << "--ppm-max" << ppm_max
                     << "--output" << File::removeExtension(ms_files[i]).toQString(); //internal output folder for temporary files

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

      process_params << ms_files[i].toQString();

      QProcess qp;
      qp.start(executable, process_params); // does automatic escaping etc... start
      bool success = qp.waitForFinished(30000); // exits job after 30 seconds

      if (!success || qp.exitStatus() != 0 || qp.exitCode() != 0)
      {
        qp.close();
        writeLog_( "Fatal error: Running SiriusAdapter returned an error code or could no compute the input within 30 seconds" );
      }

      //close the process
      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      //TODO: output from SirusMzTabWriter
      //TODO: output from CsiFingerIDMzTabWriter

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
