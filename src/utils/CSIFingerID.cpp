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
   #include <OpenMS/FORMAT/MzMLFile.h>
   #include <OpenMS/FORMAT/FileHandler.h>
   #include <OpenMS/DATASTRUCTURES/ListUtils.h>
   #include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
   #include <OpenMS/FORMAT/MzTab.h>
   #include <OpenMS/FORMAT/CsvFile.h>

   #include <QtCore/QFile>
   #include <QtCore/QProcess>
   #include <QDir>


   using namespace OpenMS;
   using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_CSIFingerID

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
    @verbinclude UTILS_CSIFingerID.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_CSIFingerID.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPCSIFingerID :
  public TOPPBase
{
public:
  TOPPCSIFingerID() :
    TOPPBase("CSIFingerID", "Tool for metabolite identification using single and tandem mass spectrometry", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
   registerInputFile_("executable", "<file>", "", "sirius file e.g. /bin/sirius3", true, false, ListUtils::create<String>("skipexists"));
   registerInputFile_("in", "<file>", "", "mzML File");
   setValidFormats_("in", ListUtils::create<String>("mzml"));

   //registerOutputFile_("out", "<file>", "", "mzTab File");
   //setValidFormats_("out", ListUtils::create<String>("csv"));

   registerStringOption_("analysis_profile", "<choice>", "qtof", "Specify the used analysis profile", false);
   setValidStrings_("analysis_profile", ListUtils::create<String>("qtof,orbitrap,fticr"));
   registerIntOption_("candidates", "<num>", 5, "The number of candidates in the output. Default 5 best candidates", false);
   registerStringOption_("database", "<choice>", "all", "search formulas in given database", false);
   setValidStrings_("database", ListUtils::create<String>("all,chebi,custom,kegg,bio,natural products,pubmed,hmdb,biocyc,hsdb,knapsack,biological,zinc bio,gnps,pubchem,mesh,maconda"));
   registerIntOption_("noise", "<num>", 0, "median intensity of noise peaks", false);
   registerIntOption_("ppm_max", "<num>", 10, "allowed ppm for decomposing masses", false);
   registerStringOption_("formula", "<choice>", "","specify the neutral molecular formula of the measured compound to compute its tree or a list of candidate formulas the method should discriminate. Omit this option if you want to consider all possible molecular formulas", false);
   registerStringOption_("isotope", "<choice>", "both", "how to handle isotope pattern data. Use 'score' to use them for ranking or 'filter' if you just want to remove candidates with bad isotope pattern. With 'both' you can use isotopes for filtering and scoring (default). Use 'omit' to ignore isotope pattern.", false);
   setValidStrings_("isotope", ListUtils::create<String>("score,filter,both,omit"));
   registerStringOption_("elements", "<choice>", "CHNOP[5]S", "The allowed elements. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the maximal allowed occurence of these elements: CHNOP[5]S[8]Cl[1]. By default CHNOP[5]S is used.", false);
   registerIntOption_("mass_deviation", "<num>", 5, "Specify the allowed mass deviation of the fragment peak in ppm.", false);
   registerFlag_("iontree", "'--iontree' Print molecular formulas and node labels with the ion formula instead of the neutral formula", false);
   registerFlag_("no_recalibration", "'--no-recalibration' If this option is set, SIRIUS will not recalibrate the spectrum during the analysis.", false);
   registerFlag_("fingerid", "'--fingerid' If this option is set, SIRIUS will search for molecular structure using CSI:FingerId after determining the molecIf this option is set, SIRIUS will search for molecular structure using CSI:FingerId after determining the molecular formulular formula", false);
  }

  ExitCodes main_(int, const char **)
  {
   //-------------------------------------------------------------
   // Parsing parameters
   //-------------------------------------------------------------

  String in = getStringOption_("in");
  //String out = getStringOption_("out");

  // Parameter for Sirius3
  QString executable = getStringOption_("executable").toQString();

  QString analysis_profile = getStringOption_("analysis_profile").toQString();
  QString elements = getStringOption_("elements").toQString();
  QString database = getStringOption_("database").toQString();
  //QString formula = getStringOption_("formula").toQString();
  QString isotope = getStringOption_("isotope").toQString();
  QString noise = QString::number(getIntOption_("noise"));
  QString ppm_max = QString::number(getIntOption_("ppm_max"));
  QString candidates = QString::number(getIntOption_("candidates"));

  bool no_recalibration = getFlag_("no_recalibration");
  bool fingerid = getFlag_("fingerid");
  bool iontree = getFlag_("iontree");

   // laod spectra

   PeakMap spectra;
   MzMLFile f;
   f.setLogType(log_type_);
   f.load(in, spectra);

   // loop over all spectra
   for (PeakMap::ConstIterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
   {
     //process only MS2 spectra
     if (s_it->getMSLevel() != 2) continue;

     const MSSpectrum<Peak1D>& spectrum = *s_it;

     int scan_index = s_it - spectra.begin();

     const vector<Precursor>& precursor = spectrum.getPrecursors();

     IonSource::Polarity p = spectrum.getInstrumentSettings().getPolarity(); //charge

     // needed later for writing in ms file
     int int_charge(1);
     if (p == IonSource::Polarity::NEGATIVE)
     {
        int_charge = -1;
     }
     else if (p == IonSource::Polarity::POSITIVE)
     {
        int_charge = +1;
     }

     //there should be only one precursor and MS2 should contain peaks to be considered
     if (precursor.size() == 1 && !spectrum.empty())
     {
       //read charge annotated to MS2
       int precursor_charge = precursor[0].getCharge();

       //sirius only supports +1 / -1 charge so far
       if (precursor_charge > 1 || precursor_charge < -1)
       {
         LOG_WARN << "Sirius only support mono charges analytes." << endl;
       }

       double precursor_mz = precursor[0].getMZ();

       //store temporary data
       String query_id = String("unknown") + String(scan_index);

       String unique_name =  String(File::getUniqueName()).toQString(); //if not done this way - always new "unique name"
       String tmp_dir = QDir::toNativeSeparators(String(File::getTempDirectory()).toQString()) + "/" + unique_name.toQString() + "_out";
       String tmp_filename = QDir::toNativeSeparators(String(File::getTempDirectory()).toQString()) + "/" + unique_name.toQString() + "_" + query_id.toQString() + ".ms";

       //to get the path and filename
       cout << "\n" << "Temp_output_folder: " << tmp_dir << "\n" << endl;
       cout << "Temp_filename: " << tmp_filename << "\n" << endl;

       // create temporary input file (.ms)
       ofstream os(tmp_filename.c_str());
       if (!os)
       {
         throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tmp_filename);
       }

       // create temporary output folder
       if (!QDir().mkdir(tmp_dir.toQString()))
       {
         throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tmp_dir);
       }

       //TODO: MS1 data m/z and intensity of precursor and precursor isotope pattern
       //TODO: IF no MS1 information is present -  run without MS1 information- might lead to an incorrect identification
       //TODO: Collision energy optional for MS2

       //write internal unique .ms data as sirius input
       os.precision(12);
       os << fixed;
       os << ">compound " << query_id << "\n"
          << ">parentmass " << precursor_mz << fixed << "\n"
          << ">charge " << int_charge << "\n\n"
          << ">ms2" << "\n";

       //single spectrum peaks
       for (Size i = 0; i != spectrum.size(); ++i)
       {
         const Peak1D& peak = spectrum[i];
         double mz = peak.getMZ();
         float intensity = peak.getIntensity();

          os << mz << " " << intensity << "\n";
        }
       os.close();

       QStringList process_params; // the actual process
       process_params << "-p" << analysis_profile
                      << "-e" << elements
                      << "-d" << database
                      //<< "-f" << formula
                      << "-s" << isotope
                      << "--noise" << noise
                      << "-c" << candidates
                      << "--ppm-max" << ppm_max
                      << "--output" << tmp_dir.toQString(); //internal output folder for temporary files

                      if (no_recalibration) process_params << "--no-recalibration";
                      if (fingerid) process_params << "--fingerid";
                      if (iontree) process_params << "--iontree";

       process_params << tmp_filename.toQString();

       //int status = QProcess::startDetached((executable, process_params)); //instead of QProcess::execute
       int status = QProcess::execute(executable, process_params);
       if (status != 0)
       {
         writeLog_("Fatal error: Running sirius returned an error code");
         return EXTERNAL_PROGRAM_ERROR;
       }

       // read results from sirius output files
       CsvFile compounds(tmp_dir + "/" + unique_name + "_" + query_id +".csv", '\t');



       for (Size j = 0; j != compounds.rowCount(); ++j)
       {
         StringList sl;
         compounds.getRow(j, sl);
         for (Size k = 0; k != sl.size(); ++k)
         {
           cout << sl[k];
           if (k != sl.size() - 1)
           {
             cout << "\t";
           } else
           {
            cout << "\n";
           }
         }
       }

       // clean up temporary input files and output folder
       if (getIntOption_("debug") < 2)
       {
         writeDebug_("Removing temporary files", 1);

         // remove temporary input file
         if (File::exists(tmp_filename) && !File::remove(tmp_filename))
         {
           LOG_WARN << "Unable to remove temporary file: " << tmp_filename << endl;
         }

         // remove temporary output folder
         if (File::removeDirRecursively(tmp_dir))
         {
           LOG_WARN << "Unable to remove temporary folder: " << tmp_dir << endl;
         }
       }
       else
       {
         writeDebug_(String("Input to sirius kept for inspection at ") + tmp_filename + "\n", 2);
         writeDebug_(String("Output folder kept for inspection at ") + tmp_dir + "\n", 2);
       }
     }
   }
   return EXECUTION_OK;
 }
};

int main(int argc, const char ** argv)
{
  TOPPCSIFingerID tool;
  return tool.main(argc, argv);
}

/// @endcond
