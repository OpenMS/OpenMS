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
    registerStringOption_("extra_index_args", "<choice>", "", "Extra arguments to be passed to tide-index", false, false);
    registerStringOption_("extra_search_args", "<choice>", "", "Extra arguments to be passed to tide-search", false, false);
    registerStringOption_("extra_percolator_args", "<choice>", "", "Extra arguments to be passed to percolator", false, false);

    //Files
    registerOutputFile_("pin_out", "<file>", "", "Output file - for Percolator input", false);
    setValidFormats_("pin_out", ListUtils::create<String>("csv"));

    //Masses
    registerDoubleOption_("precursor_mass_tolerance", "<tolerance>", 10.0, "Precursor monoisotopic mass tolerance (Crux parameter: peptide_mass_tolerance)", false, false);
    registerStringOption_("precursor_mass_units", "<choice>", "ppm", "Unit of precursor mass tolerance (amu, m/z or ppm)", false, false);
    setValidStrings_("precursor_mass_units", ListUtils::create<String>("mass,mz,ppm"));
    registerDoubleOption_("fragment_bin_offset", "<offset>", 0.0, "In the discretization of the m/z axes of the observed and theoretical spectra, this parameter specifies the location of the left edge of the first bin, relative to mass = 0 (i.e., mz-bin-offset = 0.xx means the left edge of the first bin will be located at +0.xx Da).", false, false);
    registerDoubleOption_("fragment_bin_width", "<width>", 0.02, "Before calculation of the XCorr score, the m/z axes of the observed and theoretical spectra are discretized. This parameter specifies the size of each bin. The exact formula for computing the discretized m/z value is floor((x/mz-bin-width) + 1.0 - mz-bin-offset), where x is the observed m/z value. For low resolution ion trap ms/ms data 1.0005079 and for high resolution ms/ms 0.02 is recommended.", false, false);
    registerStringOption_("isotope_error", "<choice>", "", "List of positive, non-zero integers.", false, false);

    //Search Enzyme
    registerStringOption_("enzyme", "<cleavage site>", "trypsin", "The enzyme used for peptide digestion.", false, false);
    setValidStrings_("enzyme", ListUtils::create<String>("no-enzyme,trypsin,trypsin/p,chymotrypsin,elastase,clostripain,cyanogen-bromide,iodosobenzoate,proline-endopeptidase,staph-protease,asp-n,lys-c,lys-n,arg-c,glu-c,pepsin-a,elastase-trypsin-chymotrypsin,custom-enzyme"));
    registerStringOption_("digestion", "<choice>", "full-digest", "Full, partial or non specific digestion", false, false);
    setValidStrings_("digestion", ListUtils::create<String>("full-digest,partial-digest,non-specific-digest"));
    registerIntOption_("allowed_missed_cleavages", "<num>", 0, "Number of possible cleavage sites missed by the enzyme, maximum value is 5; for enzyme search", false, false);
    registerStringOption_("custom_enzyme", "<enzyme description>", "", "Specify rules for in silico digestion of protein sequences. Overrides the enzyme option. Two lists of residues are given enclosed in square brackets or curly braces and separated by a |. The first list contains residues required/prohibited before the cleavage site and the second list is residues after the cleavage site.  ", false, false); 
    registerStringOption_("decoy_prefix", "<decoy_prefix>", "decoy_", "Specifies the prefix of the protein names that indicate a decoy", false, false); 

    //Modifications
    registerStringOption_("cterm_modifications", "<mods>", "", "Specifies C-terminal static and variable mass modifications on peptides.  Specify a comma-separated list of C-terminal modification sequences of the form: X+21.9819 Default = <empty>.", false, false);
    registerStringOption_("nterm_modifications", "<mods>", "", "Specifies N-terminal static and variable mass modifications on peptides.  Specify a comma-separated list of N-terminal modification sequences of the form: 1E-18.0106,C-17.0265 Default = <empty>.", false, false);
    registerStringOption_("modifications", "<mods>", "", "Expression for static and variable mass modifications to include. Specify a comma-separated list of modification sequences of the form: C+57.02146,2M+15.9949,1STY+79.966331,... Default = C+57.02146.", false, false);
     
    // Percolator
    registerDoubleOption_("test_fdr", "<fdr>", 0.01, "False discovery rate threshold used in selecting hyperparameters during internal cross-validation and for reporting the final results.", false, false);
    registerDoubleOption_("train_fdr", "<fdr>", 0.01, "False discovery rate threshold to define positive examples in training.", false, false);

    registerFlag_("deisotope", "Deisotope spectra before searching", true);
  }

/*
  
FATAL: Expected 1 arguments, but found 0

USAGE:

  crux percolator [options] <peptide-spectrum matches>

REQUIRED ARGUMENTS:

  <peptide-spectrum matches> A collection of target and decoy peptide-spectrum
  matches (PSMs). Input may be in one of five formats: PIN, SQT, pepXML, Crux
  tab-delimited text, or a list of files (when list-of-files=T). Note that if
  the input is provided as SQT, pepXML, or Crux tab-delimited text, then a PIN
  file will be generated in the output directory prior to execution.Crux
  determines the format of the input file by examining its filename extension.  

OPTIONAL ARGUMENTS:

  [--fileroot <string>]
     The fileroot string will be added as a prefix to all output file names.
     Default = <empty>.
  [--output-dir <string>]
     The name of the directory where output files will be created. Default =
     crux-output.
  [--overwrite T|F]
     Replace existing files if true or fail when trying to overwrite a file if
     false. Default = false.
  [--txt-output T|F]
     Output a tab-delimited results file to the output directory. Default =
     true.
  [--pout-output T|F]
     Output a Percolator pout.xml format results file to the output directory.
     Default = false.
  [--mzid-output T|F]
     Output an mzIdentML results file to the output directory. Default = false.
  [--pepxml-output T|F]
     Output a pepXML results file to the output directory. Default = false.
  [--feature-file-out T|F]
     Output the computed features in tab-delimited text format. Default = false.
  [--list-of-files T|F]
     Specify that the search results are provided as lists of files, rather than
     as individual files. Default = false.
  [--parameter-file <string>]
     A file containing parameters.  Default = <empty>.
  [--feature-file-in T|F]
     When set to T, interpret the input file as a PIN file. Default = false.
  [--picked-protein <string>]
     Use the picked protein-level FDR to infer protein probabilities, provide
     the fasta file as the argument to this flag. Default = <empty>.
  [--protein-enzyme no_enzyme|elastase|pepsin|proteinasek|thermolysin|trypsinp|chymotrypsin|lys-n|lys-c|arg-c|asp-n|glu-c|trypsin]
     Type of enzyme Default = trypsin.
  [--protein-report-fragments T|F]
     By default, if the peptides associated with protein A are a proper subset
     of the peptides associated with protein B, then protein A is eliminated and
     all the peptides are considered as evidence for protein B. Note that this
     filtering is done based on the complete set of peptides in the database,
     not based on the identified peptides in the search results. Alternatively,
     if this option is set and if all of the identified peptides associated with
     protein B are also associated with protein A, then Percolator will report a
     comma-separated list of protein IDs, where the full-length protein B is
     first in the list and the fragment protein A is listed second. Not
     available for Fido. Default = false.
  [--protein-report-duplicates T|F]
     If multiple database proteins contain exactly the same set of peptides,
     then Percolator will randomly discard all but one of the proteins. If this
     option is set, then the IDs of these duplicated proteins will be reported
     as a comma-separated list. Not available for Fido. Default = false.
  [--protein T|F]
     Use the Fido algorithm to infer protein probabilities. Must be true to use
     any of the Fido options. Default = false.
  [--decoy-xml-output T|F]
     Include decoys (PSMs, peptides, and/or proteins) in the XML output. Default
     = false.
  [--decoy-prefix <string>]
     Specifies the prefix of the protein names that indicate a decoy. Default =
     decoy_.
  [--subset-max-train <integer>]
     Only train Percolator on a subset of PSMs, and use the resulting score
     vector to evaluate the other PSMs. Recommended when analyzing huge numbers
     (>1 million) of PSMs. When set to 0, all PSMs are used for training as
     normal. Default = 0.
  [--c-pos <float>]
     Penalty for mistakes made on positive examples. If this value is set to 0,
     then it is set via cross validation over the values {0.1, 1, 10}, selecting
     the value that yields the largest number of PSMs identified at the q-value
     threshold set via the --test-fdr parameter. Default = 0.
  [--c-neg <float>]
     Penalty for mistake made on negative examples. If not specified, then this
     value is set by cross validation over {0.1, 1, 10}. Default = 0.
  [--train-fdr <float>]
     False discovery rate threshold to define positive examples in training.
     Default = 0.01.
  [--test-fdr <float>]
     False discovery rate threshold used in selecting hyperparameters during
     internal cross-validation and for reporting the final results. Default =
     0.01.
  [--maxiter <integer>]
     Maximum number of iterations for training. Default = 10.
  [--quick-validation T|F]
     Quicker execution by reduced internal cross-validation. Default = false.
  [--output-weights T|F]
     Output final weights to a file named "percolator.weights.txt". Default =
     false.
  [--init-weights <string>]
     Read initial weights from the given file (one per line). Default = <empty>.
  [--default-direction <string>]
     In its initial round of training, Percolator uses one feature to induce a
     ranking of PSMs. By default, Percolator will select the feature that
     produces the largest set of target PSMs at a specified FDR threshold (cf.
     --train-fdr). This option allows the user to specify which feature is used
     for the initial ranking, using the name as a string. The name can be
     preceded by a hyphen (e.g. "-XCorr") to indicate that a lower value is
     better. Default = <empty>.
  [--unitnorm T|F]
     Use unit normalization (i.e., linearly rescale each PSM's feature vector to
     have a Euclidean length of 1), instead of standard deviation normalization.
     Default = false.
  [--fido-alpha <float>]
     Specify the probability with which a present protein emits an associated
     peptide. Set by grid search (see --fido-gridsearch-depth parameter) if not
     specified. Default = 0.
  [--fido-beta <float>]
     Specify the probability of the creation of a peptide from noise. Set by
     grid search (see --fido-gridsearch-depth parameter) if not specified.
     Default = 0.
  [--fido-gamma <float>]
     Specify the prior probability that a protein is present in the sample. Set
     by grid search (see --fido-gridsearch-depth parameter) if not specified.
     Default = 0.
  [--test-each-iteration T|F]
     Measure performance on test set each iteration. Default = false.
  [--override T|F]
     By default, Percolator will examine the learned weights for each feature,
     and if the weight appears to be problematic, then percolator will discard
     the learned weights and instead employ a previously trained, static score
     vector. This switch allows this error checking to be overriden. Default =
     false.
  [--percolator-seed <string>]
     When given a unsigned integer value seeds the random number generator with
     that value. When given the string "time" seeds the random number generator
     with the system time. Default = 1.
  [--klammer T|F]
     Use retention time features calculated as in "Improving tandem mass
     spectrum identification using peptide retention time prediction across
     diverse chromatography conditions" by Klammer AA, Yi X, MacCoss MJ and
     Noble WS. (Analytical Chemistry. 2007 Aug 15;79(16):6111-8.). Default =
     false.
  [--only-psms T|F]
     Do not remove redundant peptides; keep all PSMs and exclude peptide level
     probability. Default = false.
  [--fido-empirical-protein-q T|F]
     Estimate empirical p-values and q-values for proteins using target-decoy
     analysis. Default = false.
  [--fido-gridsearch-depth <integer>]
     Set depth of the grid search for alpha, beta and gamma estimation. Default
     = 0.
  [--fido-gridsearch-mse-threshold <float>]
     Q-value threshold that will be used in the computation of the MSE and ROC
     AUC score in the grid search. Default = 0.05.
  [--fido-fast-gridsearch <float>]
     Apply the specified threshold to PSM, peptide and protein probabilities to
     obtain a faster estimate of the alpha, beta and gamma parameters. Default =
     0.
  [--fido-protein-truncation-threshold <float>]
     To speed up inference, proteins for which none of the associated peptides
     has a probability exceeding the specified threshold will be assigned
     probability = 0. Default = 0.01.
  [--fido-no-split-large-components T|F]
     Do not approximate the posterior distribution by allowing large graph
     components to be split into subgraphs. The splitting is done by duplicating
     peptides with low probabilities. Splitting continues until the number of
     possible configurations of each subgraph is below 2^18 Default = false.
  [--tdc T|F]
     Use target-decoy competition to assign q-values and PEPs. When set to F,
     the mix-max method, which estimates the proportion pi0 of incorrect target
     PSMs, is used instead. Default = true.
  [--verbosity <integer>]
     Specify the verbosity of the current processes. Each level prints the
     following messages, including all those at lower verbosity levels: 0-fatal
     errors, 10-non-fatal errors, 20-warnings, 30-information on the progress of
     execution, 40-more progress information, 50-debug info, 60-detailed debug
     info. Default = 30.
  [--top-match <integer>]
     Specify the number of matches to report for each spectrum. Default = 5.
  [--search-input auto|separate|concatenated]
     Specify the type of target-decoy search. Using 'auto', percolator attempts
     to detect the search type automatically.  Using 'separate' specifies two
     searches: one against target and one against decoy protein db. Using
     'concatenated' specifies a single search on concatenated target-decoy
     protein db. Default = auto.

*/

  String argumentPassthrough(const String& arg_)
  {
    String arg = arg_;
    if (arg.hasPrefix('\\'))
    {
      arg = arg.substr(1, arg.size());
    }
    return arg;
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

    bool deisotope = getFlag_("deisotope");
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
    String tmp_file;

    PeakMap exp;
    MzMLFile mzml_file;
    mzml_file.getOptions().addMSLevel(2); // only load msLevel 2 //TO DO: setMSLevels or clearMSLevels
    mzml_file.setLogType(log_type_);
    mzml_file.load(inputfile_name, exp);

#if 0
    // Low memory conversion
    {
      PlainMSDataWritingConsumer consumer(tmp_xml);
      // consumer.getOptions().setWriteIndex(write_scan_index);
      consumer.getOptions().setForceTPPCompatability(true);
      bool skip_full_count = true;
      mzml_file.setLogType(log_type_);
      mzml_file.transform(inputfile_name, &consumer, skip_full_count);
    }
#endif 
    // TODO: check if we really need to convert for TPP compatibility
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
      String tool = "tide-index";
      String params = "--overwrite T --peptide-list T --num-threads " + String(getIntOption_("threads"));
      params += " --missed-cleavages " + String(getIntOption_("allowed_missed_cleavages"));
      params += " --digestion " + getStringOption_("digestion");
			if (!getStringOption_("custom_enzyme").empty()) params += " --custom-enzyme " + getStringOption_("custom_enzyme");
			if (!getStringOption_("modifications").empty()) params += " --mods-spec " + getStringOption_("modifications");
			if (!getStringOption_("cterm_modifications").empty()) params += " --cterm-peptide-mods-spec " + getStringOption_("cterm_modifications");
			if (!getStringOption_("nterm_modifications").empty()) params += " --nterm-peptide-mods-spec " + getStringOption_("nterm_modifications");

			if (!getStringOption_("extra_index_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_index_args"));

      params.trim();
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

    std::cout << " Done running tide-index ... " << std::endl;

    // run crux tide-search
    {
      String tool = "tide-search";
      String params = "--overwrite T --file-column F   --num-threads " + String(getIntOption_("threads"));
      params += " --output-dir " + output_dir;
      String debug_args = " --verbosity 30 ";
      if (debug_level_ > 5) debug_args = " --verbosity 60 ";
      params += debug_args;

      String extra_args;
      // extra_args += " --mzid-output T";
      params += concat;
      params += extra_args;
      params += parser;

			params += " --precursor-window " + String(getDoubleOption_("precursor_mass_tolerance"));
			params += " --precursor-window-type " + getStringOption_("precursor_mass_units");
			params += " --mz-bin-offset " + String(getDoubleOption_("fragment_bin_offset"));
			params += " --mz-bin-width " + String(getDoubleOption_("fragment_bin_width"));
      if (deisotope) params += " --deisotope ";
			if (!getStringOption_("isotope_error").empty()) params += " --isotope-error " + getStringOption_("isotope_error");

			if (!getStringOption_("extra_search_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_search_args"));

      params.simplify();
      params.trim();

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
        writeLog_("Crux problem. Aborting! Calling command was: '" + 
            crux_executable + " \"" + params + " " + tmp_xml + " " + idx_name + "\"'.\nDoes the Crux executable exist?");
        return EXTERNAL_PROGRAM_ERROR;
      }
    }
    std::cout << " Done running tide-search ... " << std::endl;

    // run crux percolator
    bool run_percolator = true;
    //run_percolator = false;
    if (run_percolator)
    {
      // $ ~/bin/crux-3.1.Linux.x86_64/bin/crux  percolator --overwrite T --train-fdr 0.05 --test-fdr 0.05 crux-output/tide-search.target.txt 
      // String params = "--overwrite T --peptide-list T" + tmp_file;
      String tool = "percolator";
      String params = " --output-dir " + output_dir;
      String input = out_dir_q + "tide-search.txt";
      String debug_args = " --verbosity 30 ";
      if (debug_level_ > 5) debug_args = " --verbosity 60 ";
      params += debug_args;

      String extra_args = concat;
      params += extra_args;

      params += " --mzid-output T";
			params += " --test-fdr " + String(getDoubleOption_("test_fdr"));
			params += " --train-fdr " + String(getDoubleOption_("train_fdr"));
      params += " --overwrite T ";

			if (!getStringOption_("extra_percolator_args").empty()) params += " " + argumentPassthrough(getStringOption_("extra_percolator_args"));

      params.simplify();
      params.trim();

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

    std::cout << " Done running percolator ... " << std::endl;

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
