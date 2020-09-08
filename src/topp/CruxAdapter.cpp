// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
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

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>

#include <QtCore/QDir>

#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_CruxAdapter CruxAdapter

    @brief Identifies peptides in MS/MS spectra via Crux and tide-search.

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

    @em Crux must be installed before this wrapper can be used. This should be the case for the installers.

    The default parameters are set for a high resolution instrument.
    See the following publication: <br>
    Christopher Y. Park, Aaron A. Klammer, Lukas Käll, Michael J. MacCoss and William Stafford Noble.
    "Rapid and accurate peptide identification from tandem mass spectra." Journal of Proteome Research. 7(7):3022-3027, 2008. doi: 10.1021/pr800127y

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
    TOPPBase("CruxAdapter", "Identifies MS/MS spectra using Crux.", true,
      {
        { "Park CI, Klammer AA, Käll L, MacCoss MJ, Noble WS", "Rapid and accurate peptide identification from tandem mass spectra", "J Proteome Res 7(7):3022-3027, 2008.", "10.1021/pr800127y" }
      }
    )
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerInputFile_("database", "<file>", "", "FASTA file", true, false, ListUtils::create<String>("skipexists"));
    setValidFormats_("database", ListUtils::create<String>("FASTA"));
    registerInputFile_("crux_executable", "<executable>",
      // choose the default value according to the platform where it will be executed
      #ifdef OPENMS_WINDOWSPLATFORM
                     "crux.exe",
      #else
                     "crux",
      #endif
      "Crux executable of the installation e.g. 'crux.exe'", true, false, ListUtils::create<String>("skipexists"));

    //
    // Optional parameters //
    //
    registerStringOption_("extra_index_args", "<choice>", "", "Extra arguments to be passed to tide-index", false, false);
    registerStringOption_("extra_search_args", "<choice>", "", "Extra arguments to be passed to tide-search", false, false);
    registerStringOption_("extra_percolator_args", "<choice>", "", "Extra arguments to be passed to percolator", false, false);

    //Masses
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance (Crux parameter: peptide_mass_tolerance)", false, false);
    registerStringOption_("precursor_mass_units", "<choice>", "ppm", "Unit of precursor mass tolerance (amu, m/z or ppm)", false, false);
    setValidStrings_("precursor_mass_units", ListUtils::create<String>("mass,mz,ppm"));
    registerDoubleOption_("fragment_bin_offset", "<offset>", 0.0, "In the discretization of the m/z axes of the observed and theoretical spectra, this parameter specifies the location of the left edge of the first bin, relative to mass = 0 (i.e., mz-bin-offset = 0.xx means the left edge of the first bin will be located at +0.xx Da).", false, false);
    registerDoubleOption_("fragment_bin_width", "<width>", 0.02, "Before calculation of the XCorr score, the m/z axes of the observed and theoretical spectra are discretized. This parameter specifies the size of each bin. The exact formula for computing the discretized m/z value is floor((x/mz-bin-width) + 1.0 - mz-bin-offset), where x is the observed m/z value. For low resolution ion trap ms/ms data 1.0005079 and for high resolution ms/ms 0.02 is recommended.", false, false);
    registerStringOption_(Constants::UserParam::ISOTOPE_ERROR, "<choice>", "", "List of positive, non-zero integers.", false, false);

    registerStringOption_("run_percolator", "<true/false>", "true", "Whether to run percolator after tide-search", false, false);
    setValidStrings_("run_percolator", ListUtils::create<String>("true,false"));

    //Search Enzyme
    vector<String> all_enzymes;
    ProteaseDB::getInstance()->getAllCruxNames(all_enzymes);
    registerStringOption_("enzyme", "<cleavage site>", "trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", all_enzymes);
    registerStringOption_("digestion", "<choice>", "full-digest", "Full, partial or non specific digestion", false, false);
    setValidStrings_("digestion", ListUtils::create<String>("full-digest,partial-digest,non-specific-digest"));
    registerIntOption_("allowed_missed_cleavages", "<num>", 0, "Number of possible cleavage sites missed by the enzyme, maximum value is 5; for enzyme search", false, false);
    registerStringOption_("custom_enzyme", "<enzyme description>", "", "Specify rules for in silico digestion of protein sequences. Overrides the enzyme option. Two lists of residues are given enclosed in square brackets or curly braces and separated by a |. The first list contains residues required/prohibited before the cleavage site and the second list is residues after the cleavage site.  ", false, true);
    registerStringOption_("decoy_prefix", "<decoy_prefix>", "decoy_", "Specifies the prefix of the protein names that indicate a decoy", false, true);

    registerStringOption_("decoy_format", "<choice>", "shuffle", "Decoy generation method either by reversing the sequence or shuffling it.", false, false);
    setValidStrings_("decoy_format", ListUtils::create<String>("none,shuffle,peptide-reverse,protein-reverse"));
    registerStringOption_("keep_terminal_aminos", "<choice>", "NC", "Whether to keep N and C terminal in place or also shuffled / reversed.", false, false);
    setValidStrings_("keep_terminal_aminos", ListUtils::create<String>("N,C,NC,none"));

    //Modifications
    registerStringOption_("cterm_modifications", "<mods>", "", "Specifies C-terminal static and variable mass modifications on peptides.  Specify a comma-separated list of C-terminal modification sequences of the form: X+21.9819 Default = <empty>.", false, false);
    registerStringOption_("nterm_modifications", "<mods>", "", "Specifies N-terminal static and variable mass modifications on peptides.  Specify a comma-separated list of N-terminal modification sequences of the form: 1E-18.0106,C-17.0265 Default = <empty>.", false, false);
    registerStringOption_("modifications", "<mods>", "", "Expression for static and variable mass modifications to include. Specify a comma-separated list of modification sequences of the form: C+57.02146,2M+15.9949,1STY+79.966331,... Default = C+57.02146.", false, false);

    // Percolator
    registerDoubleOption_("test_fdr", "<fdr>", 0.01, "False discovery rate threshold used in selecting hyperparameters during internal cross-validation and for reporting the final results.", false, false);
    registerDoubleOption_("train_fdr", "<fdr>", 0.01, "False discovery rate threshold to define positive examples in training.", false, false);

    registerFlag_("deisotope", "Deisotope spectra before searching", true);
    registerFlag_("report_decoys", "Include decoys in the final reported dataset", true);
  }

  String argumentPassthrough(const String& arg_)
  {
    // get arguments that are passed to the tools directly (first escape the argument)
    String arg = arg_;
    if (arg.hasPrefix('\\'))
    {
      arg = arg.substr(1, arg.size());
    }
    return arg;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    bool deisotope = getFlag_("deisotope");
    bool report_decoys = getFlag_("report_decoys");
    bool run_percolator = getStringOption_("run_percolator") == "true";

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
    File::TempDir tmp_dir(debug_level_ >= 2);

    String output_dir = tmp_dir.getPath() + "crux-output";
    String out_dir_q = QDir::toNativeSeparators((output_dir + "/").toQString());
    String concat = " --concat T"; // concat target and decoy
    String parser = " --spectrum-parser mstoolkit "; // only this parser correctly parses our .mzML files

    String tmp_mzml = tmp_dir.getPath() + "input.mzML";

    // Low memory conversion
    {
      MzMLFile mzml_file;

      PlainMSDataWritingConsumer consumer(tmp_mzml);
      consumer.getOptions().setForceTPPCompatability(true);
      consumer.getOptions().addMSLevel(2); // only load msLevel 2
      bool skip_full_count = true;
      mzml_file.transform(inputfile_name, &consumer, skip_full_count);
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    String crux_executable = getStringOption_("crux_executable");
    String idx_name = tmp_dir.getPath() + "tmp_idx";

    // create index
    {
      String tool = "tide-index";
      String params = "--overwrite T --peptide-list T --num-threads " + String(getIntOption_("threads"));
      params += " --missed-cleavages " + String(getIntOption_("allowed_missed_cleavages"));
      params += " --digestion " + getStringOption_("digestion");
      params += " --decoy-format " + getStringOption_("decoy_format");
      params += " --keep-terminal-aminos " + getStringOption_("keep_terminal_aminos");
      if (!getStringOption_("enzyme").empty()) params += " --enzyme " + getStringOption_("enzyme");
      if (!getStringOption_("custom_enzyme").empty()) params += " --custom-enzyme " + getStringOption_("custom_enzyme");
      if (!getStringOption_("modifications").empty()) params += " --mods-spec " + getStringOption_("modifications");
      if (!getStringOption_("cterm_modifications").empty()) params += " --cterm-peptide-mods-spec " + getStringOption_("cterm_modifications");
      if (!getStringOption_("nterm_modifications").empty()) params += " --nterm-peptide-mods-spec " + getStringOption_("nterm_modifications");

      // add extra arguments passed on the command-line (pass through args)
      if (!getStringOption_("extra_index_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_index_args"));

      params.trim();
      params.simplify();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (const auto& s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << db_name.toQString() << idx_name.toQString();

      //-------------------------------------------------------------
      // run tide-index
      //-------------------------------------------------------------
      writeLog_("Executing Crux (tide-index)...");
      TOPPBase::ExitCodes exit_code = runExternalProcess_(crux_executable.toQString(), process_params);
      if (exit_code != EXECUTION_OK)
      {
        return exit_code;
      }
    }

    // run crux tide-search
    {
      String tool = "tide-search";
      String params = "--overwrite T --file-column F --num-threads " + String(getIntOption_("threads"));
      params += " --output-dir " + output_dir;
      String debug_args = " --verbosity 30 ";
      if (debug_level_ > 5) debug_args = " --verbosity 60 ";
      params += debug_args;

      String extra_args;
      if (!run_percolator) extra_args += " --mzid-output T"; // not recommended, too slow

      params += concat;
      params += extra_args;
      params += parser;

      params += " --precursor-window " + String(getDoubleOption_("precursor_mass_tolerance"));
      params += " --precursor-window-type " + getStringOption_("precursor_mass_units");
      params += " --mz-bin-offset " + String(getDoubleOption_("fragment_bin_offset"));
      params += " --mz-bin-width " + String(getDoubleOption_("fragment_bin_width"));
      if (deisotope) params += " --deisotope ";
      if (!getStringOption_(Constants::UserParam::ISOTOPE_ERROR).empty()) params += " --isotope-error " + getStringOption_(Constants::UserParam::ISOTOPE_ERROR);

      // add extra arguments passed on the command-line (pass through args)
      if (!getStringOption_("extra_search_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_search_args"));

      params.simplify();
      params.trim();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (const auto& s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << tmp_mzml.toQString() << idx_name.toQString();

      //-------------------------------------------------------------
      // run tide-index
      //-------------------------------------------------------------
      writeLog_("Executed Crux (tide-search)...");
      TOPPBase::ExitCodes exit_code = runExternalProcess_(crux_executable.toQString(), process_params);
      if (exit_code != EXECUTION_OK)
      {
        return exit_code;
      }
    }

    // run crux percolator (currently we dont have much choice in the matter)
    if (run_percolator)
    {
      String tool = "percolator";
      String params = " --output-dir " + output_dir;
      String input = out_dir_q + "tide-search.txt";
      String debug_args = " --verbosity 30 ";
      if (debug_level_ > 5) debug_args = " --verbosity 60 ";
      params += debug_args;

      String extra_args = concat;
      params += extra_args;

      params += " --mzid-output T --decoy-xml-output T ";
      params += " --test-fdr " + String(getDoubleOption_("test_fdr"));
      params += " --train-fdr " + String(getDoubleOption_("train_fdr"));
      params += " --decoy-prefix " + getStringOption_("decoy_prefix");
      params += " --overwrite T ";

      // add extra arguments passed on the command-line (pass through args)
      if (!getStringOption_("extra_percolator_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_percolator_args"));

      params.simplify();
      params.trim();

      vector<String> substrings;
      QStringList process_params;
      process_params << tool.toQString();
      params.split(' ', substrings);
      for (const auto& s : substrings)
      {
        process_params << s.toQString();
      }
      process_params << input.toQString();

      //-------------------------------------------------------------
      // run percolator
      //-------------------------------------------------------------
      writeLog_("Executing Crux (percolator)...");
      TOPPBase::ExitCodes exit_code = runExternalProcess_(crux_executable.toQString(), process_params);
      if (exit_code != EXECUTION_OK)
      {
        return exit_code;
      }
    }

    //-------------------------------------------------------------
    // writing IdXML output
    //-------------------------------------------------------------

    // read the mzIdentML output of Crux and write it to idXML
    vector<PeptideIdentification> peptide_identifications;
    vector<ProteinIdentification> protein_identifications;

    // fill search parameters
    ProteinIdentification::SearchParameters sp;
    sp.db = getStringOption_("database");
    //sp.charges = getIntList_("charge"); //dont know. Seems like tide doesnt support ranges and usually searches all?
    //TODO input options do not follow our standard so we cant just copy here
    sp.fixed_modifications = {}; //getStringList_("fixed_modifications");
    sp.variable_modifications = {}; //getStringList_("variable_modifications");
    sp.missed_cleavages = getIntOption_("allowed_missed_cleavages");
    sp.fragment_mass_tolerance = getDoubleOption_("fragment_bin_width");
    sp.fragment_mass_tolerance_ppm = "Da";
    sp.precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
    sp.precursor_mass_tolerance_ppm = getStringOption_("precursor_mass_units") == "ppm";
    sp.digestion_enzyme = *static_cast<const DigestionEnzymeProtein*>(ProteaseDB::getInstance()->getEnzyme(getStringOption_("enzyme")));

    std::cout << " will load file now " << std::endl;
    if (run_percolator)
    {
      String mzid = out_dir_q + "percolator.target.mzid";
      String mzid_decoy = out_dir_q + "percolator.decoy.mzid";
      MzIdentMLFile().load(mzid, protein_identifications, peptide_identifications);
      for (auto& protID : protein_identifications)
      {
        protID.setSearchEngine("Percolator");
        String SE = "tide_search";
        protID.setMetaValue("SE:"+SE,"");
        protID.setMetaValue(SE+":db",sp.db);
        protID.setMetaValue(SE+":db_version",sp.db_version);
        protID.setMetaValue(SE+":taxonomy",sp.taxonomy);
        protID.setMetaValue(SE+":charges",sp.charges);
        protID.setMetaValue(SE+":fixed_modifications",ListUtils::concatenate(sp.fixed_modifications, ","));
        protID.setMetaValue(SE+":variable_modifications",ListUtils::concatenate(sp.variable_modifications, ","));
        protID.setMetaValue(SE+":missed_cleavages",sp.missed_cleavages);
        protID.setMetaValue(SE+":fragment_mass_tolerance",sp.fragment_mass_tolerance);
        protID.setMetaValue(SE+":fragment_mass_tolerance_ppm",sp.fragment_mass_tolerance_ppm);
        protID.setMetaValue(SE+":precursor_mass_tolerance",sp.precursor_mass_tolerance);
        protID.setMetaValue(SE+":precursor_mass_tolerance_ppm",sp.precursor_mass_tolerance_ppm);
        protID.setMetaValue(SE+":digestion_enzyme",sp.digestion_enzyme.getName());
      }

      // also load the decoys
      if (report_decoys)
      {
        MzIdentMLFile().load(mzid_decoy, protein_identifications, peptide_identifications);
      }
    }
    else
    {
      String mzid = out_dir_q + "tide-search.mzid";
      MzIdentMLFile().load(mzid, protein_identifications, peptide_identifications);
      for (auto& protID : protein_identifications)
      {
        protID.setSearchEngine("tide-search");
        protID.setSearchParameters(sp);
      }
    }

    // write all (!) parameters as metavalues to the search parameters
    if (!protein_identifications.empty())
    {
      DefaultParamHandler::writeParametersToMetaValues(this->getParam_(), protein_identifications[0].getSearchParameters(), this->getToolPrefix());
    }

    IdXMLFile().store(out, protein_identifications, peptide_identifications);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPCruxAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
