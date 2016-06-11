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
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer $
// --------------------------------------------------------------------------
#include <OpenMS/ANALYSIS/ID/TopPerc.h>
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/ANALYSIS/ID/TopPerc.h>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QProcess>

#include <iostream>
#include <cmath>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_TopPerc TopPerc

  @brief TopPerc facilitates the input to, the call of and output integration of Percolator.
  Percolator (http://per-colator.com/) is a tool to apply semi-supervised learning for peptide
  identification from shotgun proteomics datasets.

  @experimental This tool is work in progress and usage and input requirements might change.

  <center>
  <table>
  <tr>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
  <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ MSGF+\f$ \longrightarrow \f$</td>
  <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDFilter</td>
  <td VALIGN="middle" ALIGN ="center" ROWSPAN=1> @ref TOPP_IDMapper</td>
  </tr>
  </table>
  </center>

  <p>Percolator is search engine sensitive, i.e. it's input features vary, depending on the search engine.</p>

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TopPerc.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_TopPerc.html

  Percolator is written by Lukas Käll (http://per-colator.com/ Copyright Lukas Käll <lukas.kall@scilifelab.se>)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPercolator :
  public TOPPBase
{
public:
  TOPPPercolator() :
    TOPPBase("TopPerc", "Facilitate input to Percolator and reintegrate.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", true);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerInputFile_("in_decoy", "<file>", "", "Input decoy file", false);
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid"));
    registerOutputFile_("out", "<file>", "", "Output file", true);
    std::string enzs = "no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin";
    registerStringOption_("enzyme", "<enzyme>", "trypsin", "Type of enzyme: "+enzs , false);
    setValidStrings_("enzyme", ListUtils::create<String>(enzs));
    registerInputFile_("percolator_executable", "<executable>",
        // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "Percolator executable of the installation e.g. 'percolator.exe'", true, false, ListUtils::create<String>("skipexists")
    );

    //Advanced parameters
//  //registerOutputFile_("r", "<file>", "out", "Output tab delimited results to a file instead of stdout", false, true);
    registerOutputFile_("X", "<file>", "", "path to file in xml-output format (pout). Default is: pout.tab", false, true);
    registerFlag_("e", "read xml-input format (pin) from standard input", true);
    registerFlag_("Z", "Include decoys (PSMs, peptides and/or proteins) in the xml-output. Only available if -X is used.", true);
    registerDoubleOption_("p", "<value>", 0.0, "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.", false, true);
    registerDoubleOption_("n", "<value>", 0.0, "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.", false, true);
    registerDoubleOption_("F", "<value>", 0.01, "False discovery rate threshold to define positive examples in training. Set by cross validation if 0. Default is 0.01.", false, true);
    registerDoubleOption_("t", "<value>", 0.01, "False discovery rate threshold for evaluating best cross validation result and the reported end result. Default is 0.01.", false, true);
    registerIntOption_("i", "<number>", 0, "Maximal number of iterations", false, true);
    registerFlag_("x", "Quicker execution by reduced internal cross-validation.", true);
    registerDoubleOption_("f", "<value>", 0.6, "Fraction of the negative data set to be used as train set when only providing one negative set, remaining examples will be used as test set. Set to 0.6 by default.", false, true);
    registerOutputFile_("J", "<file>", "", "Output the computed features to the given file in tab-delimited format. A file with the features with the given file name will be created", false, true);
    registerInputFile_("k", "<file>", "", "Input file given in the deprecated pin-xml format generated by e.g. sqt2pin with the -k option", false, true);
    registerOutputFile_("w", "<file>", "", "Output final weights to the given file", false, true);
    registerInputFile_("W", "<file>", "", "Read initial weights to the given file", false, true);
    registerStringOption_("V", "<featurename>", "", "The most informative feature given as the feature name, can be negated to indicate that a lower value is better.", false, true);
    registerIntOption_("v", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all, default is 2", false, true);
    registerFlag_("u", "Use unit normalization [0-1] instead of standard deviation normalization", true);
    registerFlag_("R", "Measure performance on test set each iteration", true);
    registerFlag_("O", "Override error check and do not fall back on default score vector in case of suspect score vector", true);
    registerIntOption_("S", "<value>", 1, "Setting seed of the random number generator. Default value is 1", false, true);
    registerFlag_("K", "Retention time features calculated as in Klammer et al.", true);
    registerFlag_("D", "Include description of correct features", true);
    registerOutputFile_("B", "<file>", "", "Output tab delimited results for decoys into a file", false, true);
    registerFlag_("U", "Do not remove redundant peptides, keep all PSMS and exclude peptide level probabilities.", true);
    registerFlag_("s", "skip validation of input file against xml schema", true);
    registerFlag_("A", "output protein level probabilities", true);
    registerDoubleOption_("a", "<value>", 0.0, "Probability with which a present protein emits an associated peptide (to be used jointly with the -A option). Set by grid search if not specified.", false, true);
    registerDoubleOption_("b", "<value>", 0.0, "Probability of the creation of a peptide from noise (to be used jointly with the -A option). Set by grid search if not specified", false, true);
    registerDoubleOption_("G", "<value>", 0.0, "Prior probability of that a protein is present in the sample ( to be used with the -A option). Set by grid search if not specified", false, true);
    registerFlag_("g", "treat ties as if it were one protein (Only valid if option -A is active).", true);
    registerFlag_("I", "use pi_0 value when calculating empirical q-values (no effect if option Q is activated) (Only valid if option -A is active).", true);
    registerFlag_("q", "output empirical q-values and p-values (from target-decoy analysis) (Only valid if option -A is active).", true);
    registerFlag_("N", "disactivates the grouping of proteins with similar connectivity, for example if proteins P1 and P2 have the same peptides matching both of them, P1 and P2 will not be grouped as one protein (Only valid if option -A is active).", true);
    registerFlag_("E", "Proteins graph will not be separated in sub-graphs (Only valid if option -A is active).", true);
    registerFlag_("C", "it does not prune peptides with a very low score (~0.0) which means that if a peptide with a very low score is matching two proteins, when we prune the peptide,it will be duplicated to generate two new protein groups (Only valid if option -A is active).", true);
    registerIntOption_("d", "<value>", 0, "Setting depth 0 or 1 or 2 from low depth to high depth(less computational time) of the grid search for the estimation Alpha,Beta and Gamma parameters for fido(Only valid if option -A is active). Default value is 0", false, true);
    registerStringOption_("P", "<value>", "random", "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (by default : random) (Only valid if option -A  is active).", false, true);
    registerFlag_("T", "Reduce the tree of proteins (removing low scored proteins) in order to estimate alpha,beta and gamma faster.(Only valid if option -A is active).", true);
    registerFlag_("Y", "Use target decoy competition to compute peptide probabilities.(recommended when using -A).", true);
    registerFlag_("H", "Q-value threshold that will be used in the computation of the MSE and ROC AUC score in the grid search (recommended 0.05 for normal size datasets and 0.1 for big size datasets).(Only valid if option -A is active).", true);
    registerFlag_("fido-truncation", "Proteins with a very low score (< 0.001) will be truncated (assigned 0.0 probability).(Only valid if option -A is active)", true);
    registerFlag_("Q", "Uses protein group level inference, each cluster of proteins is either present or not, therefore when grouping proteins discard all possible combinations for each group.(Only valid if option -A is active and -N is inactive).", true);
    registerFlag_("MHC", "Add a feature for MHC ligand properties to the specific PSM.", true);
    registerFlag_("same_search_db", "Manual override to ckeck if same settings for multiple search engines were applied.", true);
    registerFlag_("concat", "Manual override to concatenate multiple search results instead of merging on scan level.", true);
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // general variables and data to perform TopPerc
    //-------------------------------------------------------------
    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList in_list = getStringList_("in");
    const String in_decoy = getStringOption_("in_decoy");
    LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << " & " << in_decoy << " (decoy)" << endl;

    const String percolator_executable(getStringOption_("percolator_executable"));
    writeDebug_(String("Path to the percolator: ") + percolator_executable, 2);
    if (percolator_executable.empty()) //TODO? - TOPPBase::findExecutable after registerInputFile_("percolator_executable"... ???
    {
      writeLog_("No percolator executable specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------
    vector<vector<PeptideIdentification> > peptide_ids_list;
    vector<vector<ProteinIdentification> > protein_ids_list;
    for (size_t i = 0; i < in_list.size(); ++i)
    {
      String in = in_list[i];
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);
      if (in_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in, protein_ids, peptide_ids);
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
        MzIdentMLFile().load(in, protein_ids, peptide_ids);
      }
      //else catched by TOPPBase:registerInput being mandatory mzid or idxml

      if (peptide_ids.empty())
      {
        writeLog_("No or empty input file specified. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      //being paranoid about the presence of target decoy denominations, which are crucial to the percolator process
      for (std::vector<PeptideIdentification>::iterator pit = peptide_ids.begin(); pit != peptide_ids.end(); ++pit)
      {
        for (vector<PeptideHit>::iterator pht = pit->getHits().begin(); pht != pit->getHits().end(); ++pht)
        {
          // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (!pht->metaValueExists("target_decoy"))
          {
            if (!in_decoy.empty())
            {
              pht->setMetaValue("target_decoy", "target");
            }
            else
            {
              writeLog_("No target decoy search results discrimination possible. Aborting!");
              printUsage_();
              return ILLEGAL_PARAMETERS;
            }
          }
        }
      }
      peptide_ids_list.push_back(peptide_ids);
      protein_ids_list.push_back(protein_ids);
    }

    //-------------------------------------------------------------
    // read more input if necessary
    //-------------------------------------------------------------
    if (!in_decoy.empty() && in_list.size() == 1)
    {
      vector<PeptideIdentification> decoy_peptide_ids;
      vector<ProteinIdentification> decoy_protein_ids;
      FileHandler fh;
      FileTypes::Type in_decoy_type = fh.getType(in_decoy);
      if (in_decoy_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in_decoy, decoy_protein_ids, decoy_peptide_ids);
      }
      else if (in_decoy_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
        MzIdentMLFile().load(in_decoy, decoy_protein_ids, decoy_peptide_ids);
      }

      //paranoia check if this comes from the same search engine! (only in the first proteinidentification of the first proteinidentifications vector vector)
      {
        if (decoy_protein_ids.front().getSearchEngine()                             != protein_ids_list.front().front().getSearchEngine()                             )
        {
          LOG_WARN << "Warning about differing SearchEngine between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getScoreType()                                != protein_ids_list.front().front().getScoreType()                                )
        {
          LOG_WARN << "Warning about differing SoreType between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getPrimaryMSRunPath()                         != protein_ids_list.front().front().getPrimaryMSRunPath()                         )
        {
          LOG_WARN << "Warning about differing SearchInput between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().digestion_enzyme        != protein_ids_list.front().front().getSearchParameters().digestion_enzyme        )
        {
          LOG_WARN << "Warning about differing DigestionEnzyme between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().variable_modifications  != protein_ids_list.front().front().getSearchParameters().variable_modifications  )
        {
          LOG_WARN << "Warning about differing VarMods between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().fixed_modifications     != protein_ids_list.front().front().getSearchParameters().fixed_modifications     )
        {
          LOG_WARN << "Warning about differing FixMods between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().charges                 != protein_ids_list.front().front().getSearchParameters().charges                 )
        {
          LOG_WARN << "Warning about differing SearchCharges between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().fragment_mass_tolerance != protein_ids_list.front().front().getSearchParameters().fragment_mass_tolerance )
        {
          LOG_WARN << "Warning about differing FragTol between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getSearchParameters().precursor_tolerance     != protein_ids_list.front().front().getSearchParameters().precursor_tolerance     )
        {
          LOG_WARN << "Warning about differing PrecTol between target and decoy run" << endl;
        }
      }

      //being paranoid about the presence of target decoy denominations, which are crucial to the percolator process
      for (std::vector<PeptideIdentification>::iterator pit = decoy_peptide_ids.begin(); pit != decoy_peptide_ids.end(); ++pit)
      {
        for (std::vector<PeptideHit>::iterator pht = pit->getHits().begin(); pht != pit->getHits().end(); ++pht)
        {
          pht->setMetaValue("target_decoy", "decoy");
          //TODO what about proteins - internal target decoy handling is shitty - rework pls
        }
      }
      //TODO check overlap of ids in terms of spectrum id/reference
      peptide_ids_list.front().insert( peptide_ids.end(), decoy_peptide_ids.begin(), decoy_peptide_ids.end() );
      protein_ids_list.front().insert( protein_ids.end(), decoy_protein_ids.begin(), decoy_protein_ids.end() );
      writeLog_("Using decoy hits from separate file.");
    }
    else
    {
      writeLog_("Using decoy hits from input id file. You did you use a target decoy search, did you?");
    }


    //-------------------------------------------------------------
    // extract search engine and prepare pin
    //-------------------------------------------------------------
    String se = protein_ids_list.front().front().getSearchEngine();
    for (vector<vector<ProteinIdentification> >::iterator pilit = protein_ids_list.begin(); pilit != protein_ids_list.end(); ++pilit)
    {
      if (se != pilit->front().getSearchEngine())
      {
        se = "multiple";
        break;
      }
    }
    LOG_DEBUG << "Registered search engine: " << se << endl;
    TextFile txt;

    //TODO introduce min/max charge to parameters for now take available range
    int max_charge = 0;
    int min_charge = 10;
    for (vector<vector<PeptideIdentification> >::iterator pilit = peptide_ids_list.begin(); pilit != peptide_ids_list.end(); ++pilit)
    {
      for (vector<PeptideIdentification>::iterator it = pilit->begin(); it != pilit->end(); ++it)
      {
        for (vector<PeptideHit>::const_iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          if (hit->getCharge() > max_charge)
          {
            max_charge = hit->getCharge();
          }
          if (hit->getCharge() < min_charge)
          {
            min_charge = hit->getCharge();
          }
        }
      }
    }
    LOG_DEBUG << "Using min/max charges of " << min_charge << "/" << max_charge << endl;

    string enz_str = getStringOption_("enzyme");

    //ignore all but first input if NOT multiple for now
    if (se == "multiple")
    {
      if (getFlag_("concat"))
      {
        LOG_DEBUG << "Concatenating " << protein_ids_list.size() << " and " << peptide_ids_list.size() << endl;
        TopPerc::prepareCONCATpin(peptide_ids_list, protein_ids_list, enz_str, txt, min_charge, max_charge);
      }
      else
      {
        TopPerc::mergeMULTIids(protein_ids_list,peptide_ids_list, getFlag_("same_search_db"));  // will collapse the list (reference)
        LOG_DEBUG << "Merged to sizes " << protein_ids_list.size() << " and " << peptide_ids_list.size() << endl;
        TopPerc::prepareMULTIpin(peptide_ids_list.front(), protein_ids_list.front().front(), enz_str, txt, min_charge, max_charge);
      }
    }
    //TODO introduce custom feature selection from TopPerc::prepareCUSTOMpin to parameters
    else if (se == "MS-GF+") TopPerc::prepareMSGFpin(peptide_ids_list.front(), enz_str, txt, min_charge, max_charge, getFlag_("MHC"));
    else if (se == "Mascot") TopPerc::prepareMASCOTpin(peptide_ids_list.front(), enz_str, txt, min_charge, max_charge);
    else if (se == "XTandem") TopPerc::prepareXTANDEMpin(peptide_ids_list.front(), enz_str, txt, min_charge, max_charge);
    else
    {
      writeLog_("No known input to create percolator features from. Aborting");
      return INCOMPATIBLE_INPUT_DATA;
    }

    // create temp directory to store percolator in file pin.tab temporarily
    String temp_directory_body = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory_body.toQString());
    }
    String txt_designator = File::getUniqueName();
    String pin_file(temp_directory_body + txt_designator + "_pin.tab");
    String pout_file(temp_directory_body + txt_designator + "_pout.tab");
    txt.store(pin_file);

    QStringList arguments;
    // Check all set parameters and get them into arguments StringList
    {
      arguments << "-r" << pout_file.toQString();
      if (getFlag_("e")) arguments << "-e";
      if (getFlag_("Z")) arguments << "-Z";
      if (getDoubleOption_("p") != 0.0) arguments << "-p" << String(getDoubleOption_("p")).toQString();
      if (getDoubleOption_("n") != 0.0) arguments << "-n" << String(getDoubleOption_("n")).toQString();
      if (getDoubleOption_("F") != 0.01) arguments << "-F" << String(getDoubleOption_("F")).toQString();
      if (getDoubleOption_("t") != 0.01) arguments << "-t" << String(getDoubleOption_("t")).toQString();
      if (getIntOption_("i") != 0) arguments << "-i" << String(getIntOption_("i")).toQString();
      if (getFlag_("x")) arguments << "-x";
      if (getDoubleOption_("f") != 0.6) arguments << "-f" << String(getDoubleOption_("f")).toQString();
      if (getStringOption_("J") != "") arguments << "-J" << getStringOption_("J").toQString();
      if (getStringOption_("k") != "") arguments << "-k" << getStringOption_("k").toQString();
      if (getStringOption_("w") != "") arguments << "-w" << getStringOption_("w").toQString();
      if (getStringOption_("W") != "") arguments << "-W" << getStringOption_("W").toQString();
      if (getStringOption_("V") != "") arguments << "-V" << getStringOption_("V").toQString();
      if (getIntOption_("v") != 2) arguments << "-v" << String(getIntOption_("v")).toQString();
      if (getFlag_("u")) arguments << "-u";
      if (getFlag_("R")) arguments << "-R";
      if (getFlag_("O")) arguments << "-O";
      if (getIntOption_("S") != 1) arguments << "-S" << String(getDoubleOption_("S")).toQString();
      if (getFlag_("K")) arguments << "-K";
      if (getFlag_("D")) arguments << "-D";
      if (getStringOption_("B") != "") arguments << "-B" << getStringOption_("B").toQString();
      if (getFlag_("U")) arguments << "-U";
      if (getFlag_("s")) arguments << "-s";
      if (getFlag_("A")) arguments << "-A";
      if (getDoubleOption_("a") != 0.0) arguments << "-a" << String(getDoubleOption_("a")).toQString();
      if (getDoubleOption_("b") != 0.0) arguments << "-b" << String(getDoubleOption_("b")).toQString();
      if (getDoubleOption_("G") != 0.0) arguments << "-G" << String(getDoubleOption_("G")).toQString();
      if (getFlag_("g")) arguments << "-g";
      if (getFlag_("I")) arguments << "-I";
      if (getFlag_("q")) arguments << "-q";
      if (getFlag_("N")) arguments << "-N";
      if (getFlag_("E")) arguments << "-E";
      if (getFlag_("C")) arguments << "-C";
      if (getIntOption_("d") != 0) arguments << "-d" << String(getIntOption_("d")).toQString();
      if (getStringOption_("P") != "random") arguments << "-P" << getStringOption_("P").toQString();
      if (getFlag_("T")) arguments << "-T";
      if (getFlag_("Y")) arguments << "-Y";
      if (getFlag_("H")) arguments << "-H";
      if (getFlag_("fido-truncation")) arguments << "--fido-truncation";
      if (getFlag_("Q")) arguments << "-Q";
      arguments << "-U";
      arguments << pin_file.toQString();
    }
    writeLog_("Prepared percolator input.");

    //-------------------------------------------------------------
    // run percolator
    //-------------------------------------------------------------
    // Percolator execution with the executable ant the arguments StringList
    int status = QProcess::execute(percolator_executable.toQString(), arguments); // does automatic escaping etc...
    if (status != 0)
    {
      writeLog_("Percolator problem. Aborting! Calling command was: '" + percolator_executable + " \"" + arguments.join("-").toStdString() + "\".");
      // clean temporary files
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory_body);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory_body << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory_body << "'. Set debug level to <2 to remove them." << std::endl;
      }
      return EXTERNAL_PROGRAM_ERROR;
    }
    writeLog_("Executed percolator!");


    //-------------------------------------------------------------
    // reintegrate pout results
    //-------------------------------------------------------------
    // when percolator finished calculation, it stores the results -r option (with or without -U) or -m (which seems to be not working)
    //  WARNING: The -r option cannot be used in conjunction with -U: no peptide level statistics are calculated, redirecting PSM level statistics to provided file instead.

    CsvFile csv_file(pout_file, '\t');

    map<String, vector<TopPerc::PercolatorResult> > pep_map;
    StringList row;

    for (size_t i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      TopPerc::PercolatorResult res(row);
      StringList spl;
      res.PSMId.split("_",spl);
      String spec_ref = spl.front();
      if (pep_map.find(spec_ref) == pep_map.end())
      {
        pep_map[spec_ref] = vector<TopPerc::PercolatorResult>();
      }
      pep_map[spec_ref].push_back(res);
    }

    // As the percolator output file is not needed anymore, the temporary directory is going to be deleted
    if (this->debug_level_ < 99)
    {
      File::removeDirRecursively(temp_directory_body);
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << temp_directory_body << "'. Set debug level to <2 to remove them." << std::endl;
    }

    // Add the percolator results to the peptide vector of the original input file
    size_t c_debug = 0;
    size_t cnt = 0;
    for (vector<PeptideIdentification>::iterator it = peptide_ids_list.front().begin(); it != peptide_ids_list.front().end(); ++it)
    {
      String sid = it->getMetaValue("spectrum_reference");
      sid = sid.removeWhitespaces();
      if (pep_map.find(sid) == pep_map.end())
      {
        String sid_ = sid;
        vector<String> sr;
        sid.split('=', sr);
        sid = sr.back();
        if (pep_map.find(sid) == pep_map.end())
        {
          ++c_debug;
          LOG_DEBUG << "No suitable PeptideIdentification entry found for .pout entry " << sid << " or " << sid_ << endl;
          continue;
        }
      }

      //check each PeptideHit for compliance with one of the PercolatorResults (by sequence)
      for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
      {
        String pis = hit->getSequence().toUnmodifiedString();
        for (vector<TopPerc::PercolatorResult>::iterator pr = pep_map.find(sid)->second.begin(); pr != pep_map.find(sid)->second.end(); ++pr)
        {
          if (pis == pr->peptide &&
                    pr->preAA == hit->getPeptideEvidences().front().getAABefore() &&
                    pr->postAA == hit->getPeptideEvidences().front().getAAAfter())
          {
            hit->setMetaValue("MS:1001492", pr->score);  // svm score
            hit->setMetaValue("MS:1001491", pr->qvalue);  // percolator q value
            hit->setMetaValue("MS:1001493", pr->posterior_error_prob);  // percolator pep
            ++cnt;
          }
        }
      }
    }
    LOG_INFO << "No suitable PeptideIdentification for " << c_debug << " out of " << peptide_ids_list.front().size() << endl;
    LOG_INFO << "No suitable PeptideHits for " << cnt << " found." << endl;

    for (vector<ProteinIdentification>::iterator it = protein_ids_list.front().begin(); it != protein_ids_list.front().end(); ++it)
    {
      //will not be set because ALL decoy hits got no new score
      //it->setSearchEngine("Percolator");
      //it->setScoreType("q-value");
      //it->setHigherScoreBetter(false);

      //TODO add software percolator and topperc
      it->setMetaValue("percolator", "TopPerc");
      ProteinIdentification::SearchParameters sp = it->getSearchParameters();
      //TODO write all percolator parameters as set here in sp
      it->setSearchParameters(sp);
    }

    // Storing the PeptideHits with calculated q-value, pep and svm score
    MzIdentMLFile().store(getStringOption_("out").toQString().toStdString(), protein_ids_list.front(), peptide_ids_list.front());

    writeLog_("TopPerc finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPPercolator tool;

  return tool.main(argc, argv);
}

/// @endcond
