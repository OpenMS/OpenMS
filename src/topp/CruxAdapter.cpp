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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>
#include <QDebug>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_CruxAdapter CruxAdapter

    @brief Identifies peptides in MS/MS spectra via Crux.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ CruxAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    @em Crux must be installed before this wrapper can be used. 

    The default parameters are set for a high resolution instrument.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_CruxAdapter.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_CruxAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPCruxAdapter :
  public TOPPBase
{
public:
  TOPPCruxAdapter() :
    TOPPBase("CruxAdapter", "Annotates MS/MS spectra using Crux.")
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("database", "<file>", "", "FASTA file", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    registerInputFile_("crux_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      "crux.exe",
      "Crux executable of the installation e.g. 'Crux.exe'", true, false, ListUtils::create<String>("skipexists"));
    //
    // Optional parameters //
    //

    //Files
    registerOutputFile_("pin_out", "<file>", "", "Output file - for Percolator input", false);
    setValidFormats_("pin_out", ListUtils::create<String>("csv"));

    //Masses
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance (Crux parameter: peptide_mass_tolerance)", false, false);
    registerStringOption_("precursor_error_units", "<choice>", "ppm", "peptide_mass_units 0=amu, 1=mmu, 2=ppm", false, false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("amu,ppm,Da"));
    //registerIntOption_("mass_type_parent", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("mass_type_fragment", "<num>", 1, "0=average masses, 1=monoisotopic masses", false, true);
    //registerIntOption_("precursor_tolerance_type", "<num>", 0, "0=average masses, 1=monoisotopic masses", false, false);
    registerStringOption_("isotope_error", "<choice>", "off", "0=off, 1=on -1/0/1/2/3 (standard C13 error), 2= -8/-4/0/4/8 (for +4/+8 labeling)", false, false);
    setValidStrings_("isotope_error", ListUtils::create<String>("off,-1/0/1/2/3,-8/-4/0/4/8"));
  }

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

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name.empty())
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    String out = getStringOption_("out");
    writeDebug_(String("Output file___real one: ") + out, 1);
    if (out.empty())
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    String db_name(getStringOption_("database"));
    if (!File::readable(db_name))
    {
      String full_db_name;
      try
      {
        full_db_name = File::findDatabase(db_name);
      }
      catch (...)
      {
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }
      db_name = full_db_name;
    }

    //tmp_dir
    //const String tmp_dir = QDir::toNativeSeparators((File::getTempDirectory() + "/").toQString());
    const String tmp_dir = makeTempDirectory_(); //OpenMS::File::getTempDirectory() + "/";

    String output_dir = tmp_dir + "crux-output";
    String out_dir_q = QDir::toNativeSeparators((output_dir + "/").toQString());
    String concat = " --concat T"; // concat target and decoy
    String parser = " --spectrum-parser mstoolkit "; // only this parser correctly parses our .mzML files

    writeDebug_("Creating temporary directory '" + tmp_dir + "'", 1);
    String tmp_xml = tmp_dir + "input.mzML";
    String tmp_pepxml = tmp_dir + "result.pep.xml";
    String tmp_pin = tmp_dir + "result.pin";
    String default_params = getStringOption_("default_params_file");
    String tmp_file;

    PeakMap exp;
    MzMLFile mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load msLevel 2 //TO DO: setMSLevels or clearMSLevels
    mzml_file.setLogType(log_type_);
    mzml_file.load(inputfile_name, exp);

    // TODO: check if we need to convert for TPP compatibility
    {
      mzml_file.getOptions().setForceTPPCompatability(true);
      mzml_file.store(tmp_xml, exp);
    }

    //
    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS2 spectra in input file.");
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();

    if (spectrum_type == SpectrumSettings::RAWDATA)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__,
            "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String crux_executable = getStringOption_("crux_executable");
    String idx_name = tmp_dir + "tmp_idx";

    // create index
    {
      // $ ~/projects/crux-toolkit/src/crux tide-index --overwrite T --peptide-list T spyo_combined_withiRT.fasta spyo_idx
      String tool = "tide-index";
      String params = "--overwrite T --peptide-list T";
      params.simplify();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (auto s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << db_name.toQString() << idx_name.toQString();

      qDebug() << process_params;

      int status = QProcess::execute(crux_executable.toQString(), process_params); // does automatic escaping etc...
      if (status != 0)
      {
        writeLog_("Crux problem. Aborting! Calling command was: '" + crux_executable + " \"" + inputfile_name + "\"'.\nDoes the Crux executable exist?");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }

    // run crux tide-search
    {
      // $ ~/projects/crux-toolkit/src/crux tide-search --overwrite T --file-column F --spectrum-charge  all /tmp/out.mzML spyo_idx --precursor-window 3  --precursor-window-type mass --spectrum-parser mstoolkit
      // String params = "--overwrite T --peptide-list T" + tmp_file;
      String tool = "tide-search";
      String output = " --output-dir " + output_dir;
      String debug_args = "--verbosity 30";
      if (debug_level_ > 0) debug_args = " --verbosity 60";

      String extra_args = " --spectrum-charge all --precursor-window 3    --precursor-window-type mass  " + debug_args;
      extra_args += " --mzid-output T" + concat;

      String params = "--overwrite T --file-column F  " + parser + " " + extra_args + " " + output;
      params.simplify();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (auto s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << tmp_xml.toQString() << idx_name.toQString();
      qDebug() << process_params;

      int status = QProcess::execute(crux_executable.toQString(), process_params); // does automatic escaping etc...
      if (status != 0)
      {
        writeLog_("Crux problem. Aborting! Calling command was: '" + crux_executable + " \"" + inputfile_name + "\"'.\nDoes the Crux executable exist?");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }

    // run crux percolator
    bool run_percolator = true;
    if (run_percolator)
    {
      // $ ~/bin/crux-3.1.Linux.x86_64/bin/crux  percolator --overwrite T --train-fdr 0.05 --test-fdr 0.05 crux-output/tide-search.target.txt 
      // String params = "--overwrite T --peptide-list T" + tmp_file;
      String tool = "percolator";
      String output = " --output-dir " + output_dir;
      // output = "";
      String debug_args = "--verbosity 30";
      String input = out_dir_q + "tide-search.txt";
      if (debug_level_ > 0) debug_args = " --verbosity 60";
      String extra_args = "" + debug_args + concat;
      // extra_args += " --mzid-output T";

      String params = "--overwrite T --train-fdr 0.05 --test-fdr 0.05 " + extra_args + output;
      params.simplify();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (auto s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << input.toQString();
      qDebug() << process_params;

      int status = QProcess::execute(crux_executable.toQString(), process_params); // does automatic escaping etc...
      if (status != 0)
      {
        writeLog_("Crux problem. Aborting! Calling command was: '" + crux_executable + " \"" + inputfile_name + "\"'.\nDoes the Crux executable exist?");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the mzIdentML output of Crux and write it to idXML
    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    String mzid = out_dir_q + "tide-search.mzid";
    writeDebug_("load mzIdentml", 1);
    Identification target_id;
    Identification decoy_id;
    MzIdentMLFile().load(mzid, protein_identifications, peptide_identifications);
    writeDebug_("write idXMLFile", 1);
    writeDebug_(out, 1);
    IdXMLFile().store(out, protein_identifications, peptide_identifications);

    // remove tempdir
    if (this->debug_level_ == 0)
    {
        removeTempDir_(tmp_dir);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << tmp_dir << "'" << std::endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << tmp_pepxml << "'. Set debug level to 0 to remove them." << std::endl;
    }

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPCruxAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
