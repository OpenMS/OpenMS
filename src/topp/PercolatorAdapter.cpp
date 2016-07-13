// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
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
  @page TOPP_PercolatorAdapter PercolatorAdapter

  @brief PercolatorAdapter facilitates the input to, the call of and output integration of Percolator.
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
  @verbinclude TOPP_PercolatorAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_PercolatorAdapter.html

  Percolator is written by Lukas Käll (http://per-colator.com/ Copyright Lukas Käll <lukas.kall@scilifelab.se>)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class PercolatorAdapter :
  public TOPPBase
{
public:
  PercolatorAdapter() :
    TOPPBase("PercolatorAdapter", "Facilitate input to Percolator and reintegrate.", false)
  {
  }

protected:
  struct PercolatorResult
    {
      String PSMId;
      double score;
      double qvalue;
      double posterior_error_prob;
      String peptide;
      char preAA;
      char postAA;
      StringList proteinIds;

      PercolatorResult(const String& pid, const double s, const double q, const String& p, const char pre, const char pos, const StringList& pl):
          PSMId (pid),
          score (s),
          qvalue (q),
          peptide (p),
          preAA (pre),
          postAA (pos),
          proteinIds (pl)
      {
      }

      PercolatorResult(StringList& row):
      proteinIds()
      {
        // peptide sequence
        StringList pep;
        row[4].split(".", pep);
        //TODO test pep size 3
        peptide = pep[1];
        preAA = pep[0]=="-"?'[':pep[0].c_str()[0];  // const char PeptideEvidence::N_TERMINAL_AA = '[';
        postAA = pep[2]=="-"?']':pep[2].c_str()[0]; // const char PeptideEvidence::C_TERMINAL_AA = ']';
        // SVM-score
        score = row[1].toDouble();
        // q-Value
        qvalue = row[2].toDouble();
        // PEP
        posterior_error_prob = row[3].toDouble();
        // scannr. as written in preparePIN
        PSMId = row[0];
        proteinIds = vector<String>(row.begin()+5,row.end());
      }

      bool operator!=(const PercolatorResult& rhs) const
      {
        if (PSMId != rhs.PSMId || score != rhs.score || qvalue != rhs.qvalue ||
            posterior_error_prob != rhs.posterior_error_prob || peptide != rhs.peptide ||
            proteinIds != rhs.proteinIds)
          return true;
        return false;
      }

      bool operator==(const PercolatorResult& rhs) const
      {
        return !(operator !=(rhs));
      }
  };
  
  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", true);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerInputFileList_("in_decoy", "<files>", StringList(), "Input decoy file(s) in case of separate searches", false);
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid,idXML"));
    registerOutputFile_("out", "<file>", "", "Output file in idXML format", false);
    registerOutputFile_("mzid_out", "<file>", "", "Output file in mzid format", false);
    String enzs = "no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin";
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
    registerDoubleOption_("cpos", "<value>", 0.0, "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.", false, true);
    registerDoubleOption_("cneg", "<value>", 0.0, "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.", false, true);
    registerDoubleOption_("testFDR", "<value>", 0.01, "False discovery rate threshold for evaluating best cross validation result and the reported end result.", false, true);
    registerDoubleOption_("trainFDR", "<value>", 0.01, "False discovery rate threshold to define positive examples in training. Set to testFDR if 0.", false, true);
    registerIntOption_("maxiter", "<number>", 10, "Maximal number of iterations", false, true);
    registerFlag_("quick-validation", "Quicker execution by reduced internal cross-validation.", true);
    registerOutputFile_("weights", "<file>", "", "Output final weights to the given file", false, true);
    registerInputFile_("init-weights", "<file>", "", "Read initial weights to the given file", false, true);
    registerStringOption_("default-direction", "<featurename>", "", "The most informative feature given as the feature name, can be negated to indicate that a lower value is better.", false, true);
    registerIntOption_("verbose", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all.", false, true);
    registerFlag_("unitnorm", "Use unit normalization [0-1] instead of standard deviation normalization", true);
    registerFlag_("test-each-iteration", "Measure performance on test set each iteration", true);
    registerFlag_("override", "Override error check and do not fall back on default score vector in case of suspect score vector", true);
    registerIntOption_("seed", "<value>", 1, "Setting seed of the random number generator.", false, true);
    registerIntOption_("doc", "<value>", 0, "Include description of correct features", false, true);
    registerFlag_("klammer", "Retention time features calculated as in Klammer et al. Only available if -doc is set", true);
    registerFlag_("picked-protein", "Use the picked protein-level FDR to infer protein probabilities.", true);
    registerInputFile_("fasta", "<file>", "", "Provide the fasta file as the argument to this flag, which will be used for protein grouping based on an in-silico digest (only valid if option -picked-protein is active).", false, true);
    setValidFormats_("fasta", ListUtils::create<String>("FASTA"));
    registerStringOption_("decoy-pattern", "<value>", "random", "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (Only valid if option -picked-protein is active).", false, true);
    registerFlag_("post-processing-tdc", "Use target-decoy competition to assign q-values and PEPs.", true);
  }
  
  String getScanIdentifier_(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
  {
    String scan_identifier = it->getMetaValue("spectrum_reference");
    if (scan_identifier.empty())
    {
      scan_identifier = String(it->getMetaValue("spectrum_id"));
      if (scan_identifier.empty())
      {
        scan_identifier = String(it - start + 1);
        LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
      }
    }
    return scan_identifier.removeWhitespaces();
  }
  
  Int getScanNumber_(String scan_identifier)
  {
    Size idx = 0;
    if ((idx = scan_identifier.find("index=")) != string::npos) 
    {
      scan_identifier = scan_identifier.substr(idx + 6);
    }
    else if ((idx = scan_identifier.find("scan=")) != string::npos) 
    {
      scan_identifier = scan_identifier.substr(idx + 5);
    }
    return scan_identifier.toInt();
  }
    
  void readPoutAsMap_(String pout_file, Map<String, vector<PercolatorResult> >& pep_map)
  {
    CsvFile csv_file(pout_file, '\t');
    StringList row;

    for (Size i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      PercolatorResult res(row);
      String spec_ref = res.PSMId;
      if (pep_map.find(spec_ref) == pep_map.end())
      {
        pep_map[spec_ref] = vector<PercolatorResult>();
      }
      pep_map[spec_ref].push_back(res);
    }
  }
  
  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // general variables and data to perform PercolatorAdapter
    //-------------------------------------------------------------
    vector<PeptideIdentification> peptide_ids;
    vector<ProteinIdentification> protein_ids;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList in_list = getStringList_("in");
    const StringList in_decoy = getStringList_("in_decoy");
    LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << " & " << ListUtils::concatenate(in_decoy, ",") << " (decoy)" << endl;

    const String percolator_executable(getStringOption_("percolator_executable"));
    writeDebug_(String("Path to the percolator: ") + percolator_executable, 2);
    if (percolator_executable.empty()) //TODO? - TOPPBase::findExecutable after registerInputFile_("percolator_executable"... ???
    {
      writeLog_("No percolator executable specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    const String mzid_out(getStringOption_("mzid_out"));
    const String out(getStringOption_("out"));
    if (mzid_out.empty() && out.empty())
    {
      writeLog_("Fatal error: no output file given (parameter 'out' or 'mzid_out')");
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
      for (vector<PeptideIdentification>::iterator pit = peptide_ids.begin(); pit != peptide_ids.end(); ++pit)
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
      FileTypes::Type in_decoy_type = fh.getType(in_decoy.front());
      if (in_decoy_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in_decoy.front(), decoy_protein_ids, decoy_peptide_ids);
      }
      else if (in_decoy_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: you might experience loss of information depending on the capabilities of the target format." << endl;
        MzIdentMLFile().load(in_decoy.front(), decoy_protein_ids, decoy_peptide_ids);
      }

      //paranoia check if this comes from the same search engine! (only in the first proteinidentification of the first proteinidentifications vector vector)
      {
        if (decoy_protein_ids.front().getSearchEngine()                             != protein_ids_list.front().front().getSearchEngine()                         )
        {
          LOG_WARN << "Warning about differing SearchEngine between target and decoy run" << endl;
        }
        if (decoy_protein_ids.front().getScoreType()                                != protein_ids_list.front().front().getScoreType()                         )
        {
          LOG_WARN << "Warning about differing ScoreType between target and decoy run" << endl;
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
      for (vector<PeptideIdentification>::iterator pit = decoy_peptide_ids.begin(); pit != decoy_peptide_ids.end(); ++pit)
      {
        for (vector<PeptideHit>::iterator pht = pit->getHits().begin(); pht != pit->getHits().end(); ++pht)
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

    
    StringList feature_set;
    feature_set.push_back("SpecId");
    feature_set.push_back("Label");
    feature_set.push_back("ScanNr");
    feature_set.push_back("ExpMass");
    feature_set.push_back("CalcMass");
    feature_set.push_back("mass");
    feature_set.push_back("score");
    feature_set.push_back("peplen");
    for (int i = min_charge; i <= max_charge; ++i)
    {
       feature_set.push_back("charge" + String(i));
    }
    feature_set.push_back("enzN");
    feature_set.push_back("enzC");
    feature_set.push_back("enzInt");
    feature_set.push_back("dm");
    feature_set.push_back("absdm");
    feature_set.push_back("Peptide");
    feature_set.push_back("Proteins");
    
    string enz_str = getStringOption_("enzyme");     
    txt.addLine(ListUtils::concatenate(feature_set, '\t'));
    TopPerc::preparePin(peptide_ids_list.front(), feature_set, enz_str, txt, min_charge, max_charge);

    // create temp directory to store percolator in file pin.tab temporarily
    String temp_directory_body = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory_body.toQString());
    }
    String txt_designator = File::getUniqueName();
    String pin_file(temp_directory_body + txt_designator + "_pin.tab");
    String pout_target_file(temp_directory_body + txt_designator + "_target_pout.tab");
    String pout_decoy_file(temp_directory_body + txt_designator + "_decoy_pout.tab");
    txt.store(pin_file);

    QStringList arguments;
    // Check all set parameters and get them into arguments StringList
    {
      arguments << "-U";
      arguments << "-m" << pout_target_file.toQString();
      arguments << "-M" << pout_decoy_file.toQString();
      
      double cpos = getDoubleOption_("cpos");
      double cneg = getDoubleOption_("cneg");
      if (cpos != 0.0) arguments << "-p" << String(cpos).toQString();
      if (cneg != 0.0) arguments << "-n" << String(cneg).toQString();
      
      double train_FDR = getDoubleOption_("trainFDR");
      double test_FDR = getDoubleOption_("testFDR");
      if (train_FDR != 0.01) arguments << "-F" << String(train_FDR).toQString();
      if (test_FDR != 0.01) arguments << "-t" << String(test_FDR).toQString();
      
      Int max_iter = getIntOption_("maxiter");
      if (max_iter != 10) arguments << "-i" << String(max_iter).toQString();
      if (getFlag_("quick-validation")) arguments << "-x";
      if (getFlag_("post-processing-tdc")) arguments << "-Y";
      
      String weights_file = getStringOption_("weights");
      String init_weights_file = getStringOption_("init-weights");
      String default_search_direction = getStringOption_("default-direction");
      if (!weights_file.empty()) arguments << "-w" << weights_file.toQString();
      if (!init_weights_file.empty()) arguments << "-W" << init_weights_file.toQString();
      if (!default_search_direction.empty()) arguments << "-V" << default_search_direction.toQString();
      
      Int verbose_level = getIntOption_("verbose");
      if (verbose_level != 2) arguments << "-v" << String(verbose_level).toQString();
      if (getFlag_("unitnorm")) arguments << "-u";
      if (getFlag_("test-each-iteration")) arguments << "-R";
      if (getFlag_("override")) arguments << "-O";
      
      Int seed = getIntOption_("seed");
      if (seed != 1) arguments << "-S" << String(seed).toQString();
      if (getFlag_("klammer")) arguments << "-K";
      
      Int description_of_correct = getIntOption_("doc");
      if (description_of_correct != 0) arguments << "-D" << String(description_of_correct).toQString();

      String decoy_pattern = getStringOption_("decoy-pattern");
      if (decoy_pattern != "random") arguments << "-P" << decoy_pattern.toQString();
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
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory_body << "'" << endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory_body << "'. Set debug level to <2 to remove them." << endl;
      }
      return EXTERNAL_PROGRAM_ERROR;
    }
    writeLog_("Executed percolator!");


    //-------------------------------------------------------------
    // reintegrate pout results
    //-------------------------------------------------------------
    // when percolator finished calculation, it stores the results -r option (with or without -U) or -m (which seems to be not working)
    //  WARNING: The -r option cannot be used in conjunction with -U: no peptide level statistics are calculated, redirecting PSM level statistics to provided file instead.
    Map<String, vector<PercolatorResult> > pep_map;
    readPoutAsMap_(pout_target_file, pep_map);
    readPoutAsMap_(pout_decoy_file, pep_map);

    // As the percolator output file is not needed anymore, the temporary directory is going to be deleted
    if (this->debug_level_ < 99)
    {
      File::removeDirRecursively(temp_directory_body);
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << temp_directory_body << "'. Set debug level to <99 to remove them." << endl;
    }

    // Add the percolator results to the peptide vector of the original input file
    size_t c_debug = 0;
    size_t cnt = 0;
    for (vector<PeptideIdentification>::iterator it = peptide_ids_list.front().begin(); it != peptide_ids_list.front().end(); ++it)
    {
      String scan_identifier = getScanIdentifier_(it, peptide_ids_list.front().begin());
      if (pep_map.find(scan_identifier) == pep_map.end())
      {
        ++c_debug;
        LOG_DEBUG << "No suitable PeptideIdentification entry found for .pout entry " << scan_identifier << endl;
        continue;
      }

      //check each PeptideHit for compliance with one of the PercolatorResults (by sequence)
      for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
      {
        String pis = hit->getSequence().toString();
        for (vector<PercolatorResult>::iterator pr = pep_map.find(scan_identifier)->second.begin(); pr != pep_map.find(scan_identifier)->second.end(); ++pr)
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
    LOG_INFO << "Suitable PeptideHits for " << cnt << " found." << endl;

    for (vector<ProteinIdentification>::iterator it = protein_ids_list.front().begin(); it != protein_ids_list.front().end(); ++it)
    {
      //will not be set because ALL decoy hits got no new score
      //it->setSearchEngine("Percolator");
      //it->setScoreType("q-value");
      //it->setHigherScoreBetter(false);

      //TODO add software percolator and PercolatorAdapter
      it->setMetaValue("percolator", "PercolatorAdapter");
      ProteinIdentification::SearchParameters sp = it->getSearchParameters();
      //TODO write all percolator parameters as set here in sp
      it->setSearchParameters(sp);
    }
    
    // Storing the PeptideHits with calculated q-value, pep and svm score
    if (!mzid_out.empty())
    {
      MzIdentMLFile().store(mzid_out.toQString().toStdString(), protein_ids_list.front(), peptide_ids_list.front());
    }
    if (!out.empty())
    {
      IdXMLFile().store(out.toQString().toStdString(), protein_ids_list.front(), peptide_ids_list.front());
    }

    writeLog_("PercolatorAdapter finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  PercolatorAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
