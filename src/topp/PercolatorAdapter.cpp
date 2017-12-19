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
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QProcess>

#include <iostream>
#include <cmath>
#include <string>
#include <set>

#include <boost/algorithm/clamp.hpp>
#include <typeinfo>

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
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ PercolatorAdapter \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_PSMFeatureExtractor </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
        </tr>
    </table>
  </center>
  <p>Percolator is search engine sensitive, i.e. it's input features vary,
depending on the search engine. Must be prepared beforehand. If you do not want
to use the specific features, use the generic-feature-set flag. Will incorporate
the score attribute of a PSM, so be sure, the score you want is set as main
score with @ref TOPP_IDScoreSwitcher . Be aware, that you might very well
experience a perfomance loss compared to the search engine specific features.</p>

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
    TOPPBase("PercolatorAdapter", "Facilitate input to Percolator and reintegrate.", true)
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
          posterior_error_prob (0.0),
          peptide (p),
          preAA (pre),
          postAA (pos),
          proteinIds (pl)
      {
      }
      
      explicit PercolatorResult(StringList& row):
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
  
  struct PercolatorProteinResult
  {
    String protein_accession;
    double qvalue;
    double posterior_error_prob;

    PercolatorProteinResult(const String& pid, const double q, const double pep):
        protein_accession (pid),
        qvalue (q),
        posterior_error_prob (pep)
    {
    }

    bool operator!=(const PercolatorProteinResult& rhs) const
    {
      if (protein_accession != rhs.protein_accession || qvalue != rhs.qvalue ||
          posterior_error_prob != rhs.posterior_error_prob)
      {
        return true;
      }
      return false;
    }

    bool operator==(const PercolatorProteinResult& rhs) const
    {
      return !(operator !=(rhs));
    }
  };
  
  void registerOptionsAndFlags_() override
  {
    static const bool is_required(true);
    static const bool is_advanced_option(true);
    
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", !is_required);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerInputFileList_("in_decoy", "<files>", StringList(), "Input decoy file(s) in case of separate searches", !is_required);
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid,idXML"));
    registerInputFile_("in_osw", "<file>", "", "Input file in OSW format", !is_required);
    setValidFormats_("in_osw", ListUtils::create<String>("OSW"));
    registerOutputFile_("out", "<file>", "", "Output file in idXML format", !is_required);
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerOutputFile_("mzid_out", "<file>", "", "Output file in mzid format", !is_required);
    setValidFormats_("mzid_out", ListUtils::create<String>("mzid"));
    registerOutputFile_("osw_out", "<file>", "", "Output file in OSW format", !is_required);
    setValidFormats_("osw_out", ListUtils::create<String>("OSW"));
    String enzs = "no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin";
    registerStringOption_("enzyme", "<enzyme>", "trypsin", "Type of enzyme: "+enzs , !is_required);
    setValidStrings_("enzyme", ListUtils::create<String>(enzs));
    registerInputFile_("percolator_executable", "<executable>",
        // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "Percolator executable of the installation e.g. 'percolator.exe'", is_required, !is_advanced_option, ListUtils::create<String>("skipexists")
    );
    registerFlag_("peptide-level-fdrs", "Calculate peptide-level FDRs instead of PSM-level FDRs.");
    registerFlag_("protein-level-fdrs", "Use the picked protein-level FDR to infer protein probabilities. Use the -fasta option and -decoy-pattern to set the Fasta file and decoy pattern.");
    registerStringOption_("osw_level", "<osw_level>", "ms2", "OSW: Either \"ms1\", \"ms2\" or \"transition\"; the data level selected for scoring.", !is_required);

    //Advanced parameters
    registerFlag_("generic-feature-set", "Use only generic (i.e. not search engine specific) features. Generating search engine specific features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.", is_advanced_option);
    registerIntOption_("subset-max-train", "<number>", 0, "Only train an SVM on a subset of <x> PSMs, and use the resulting score vector to evaluate the other PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.", !is_required, is_advanced_option);
    registerDoubleOption_("cpos", "<value>", 0.0, "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.", !is_required, is_advanced_option);
    registerDoubleOption_("cneg", "<value>", 0.0, "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.", !is_required, is_advanced_option);
    registerDoubleOption_("testFDR", "<value>", 0.01, "False discovery rate threshold for evaluating best cross validation result and the reported end result.", !is_required, is_advanced_option);
    registerDoubleOption_("trainFDR", "<value>", 0.01, "False discovery rate threshold to define positive examples in training. Set to testFDR if 0.", !is_required, is_advanced_option);
    registerIntOption_("maxiter", "<number>", 10, "Maximal number of iterations", !is_required, is_advanced_option);
    registerFlag_("quick-validation", "Quicker execution by reduced internal cross-validation.", is_advanced_option);
    registerOutputFile_("weights", "<file>", "", "Output final weights to the given file", !is_required, is_advanced_option);
    registerInputFile_("init-weights", "<file>", "", "Read initial weights to the given file", !is_required, is_advanced_option);
    registerStringOption_("default-direction", "<featurename>", "", "The most informative feature given as the feature name, can be negated to indicate that a lower value is better.", !is_required, is_advanced_option);
    registerIntOption_("verbose", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all.", !is_required, is_advanced_option);
    registerFlag_("unitnorm", "Use unit normalization [0-1] instead of standard deviation normalization", is_advanced_option);
    registerFlag_("test-each-iteration", "Measure performance on test set each iteration", is_advanced_option);
    registerFlag_("override", "Override error check and do not fall back on default score vector in case of suspect score vector", is_advanced_option);
    registerIntOption_("seed", "<value>", 1, "Setting seed of the random number generator.", !is_required, is_advanced_option);
    registerIntOption_("doc", "<value>", 0, "Include description of correct features", !is_required, is_advanced_option);
    registerFlag_("klammer", "Retention time features calculated as in Klammer et al. Only available if -doc is set", is_advanced_option);
    registerInputFile_("fasta", "<file>", "", "Provide the fasta file as the argument to this flag, which will be used for protein grouping based on an in-silico digest (only valid if option -protein-level-fdrs is active).", !is_required, is_advanced_option);
    setValidFormats_("fasta", ListUtils::create<String>("FASTA"));
    registerStringOption_("decoy-pattern", "<value>", "random", "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (Only valid if option -protein-level-fdrs is active).", !is_required, is_advanced_option);
    registerFlag_("post-processing-tdc", "Use target-decoy competition to assign q-values and PEPs.", is_advanced_option);

    //OSW/IPF parameters
    registerDoubleOption_("ipf_max_peakgroup_pep", "<value>", 0.7, "OSW/IPF: Assess transitions only for candidate peak groups until maximum posterior error probability.", !is_required, is_advanced_option);
    registerDoubleOption_("ipf_max_transition_isotope_overlap", "<value>", 0.5, "OSW/IPF: Maximum isotope overlap to consider transitions in IPF.", !is_required, is_advanced_option);
    registerDoubleOption_("ipf_min_transition_sn", "<value>", 0, "OSW/IPF: Minimum log signal-to-noise level to consider transitions in IPF. Set -1 to disable this filter.", !is_required, is_advanced_option);

  }
  
  // TODO replace with TopPerc::getScanMergeKey
  String getScanIdentifier_(vector<PeptideIdentification>::iterator it, vector<PeptideIdentification>::iterator start)
  {
    // MSGF+ uses this field, is empty if not specified
    String scan_identifier = it->getMetaValue("spectrum_reference");
    if (scan_identifier.empty())
    {
      // XTandem uses this (integer) field
      // these ids are 1-based in contrast to the index which is 0-based. This might be problematic to use for merging
      if (it->metaValueExists("spectrum_id") && !it->getMetaValue("spectrum_id").toString().empty())
      {
        scan_identifier = "scan=" + it->getMetaValue("spectrum_id").toString();
      }
      else
      {
        scan_identifier = "index=" + String(it - start + 1);
        LOG_WARN << "no known spectrum identifiers, using index [1,n] - use at own risk." << endl;
      }
    }
    return scan_identifier.removeWhitespaces();
  }
  
  // TODO replace with TopPerc::getScanMergeKey
  Int getScanNumber_(String scan_identifier)
  {
    Int scan_number = 0;
    StringList fields = ListUtils::create<String>(scan_identifier);
    for (StringList::const_iterator it = fields.begin(); it != fields.end(); ++it)
    {
      // if scan number is not available, use the scan index
      Size idx = 0;
      if ((idx = it->find("scan=")) != string::npos)
      {
        scan_number = it->substr(idx + 5).toInt();
        break;
      }
      else if ((idx = it->find("index=")) != string::npos)
      {
        scan_number = it->substr(idx + 6).toInt();
      }
    }
    return scan_number;
  }
  
  // Function adapted from Enzyme.h in Percolator converter
  // TODO: adapt to OpenMS enzymes. Use existing functionality in EnzymaticDigestion.
  bool isEnz_(const char& n, const char& c, string& enz)
  {
    if (enz == "trypsin")
    {
      return ((n == 'K' || n == 'R') && c != 'P') || n == '-' || c == '-';
    }
    else if (enz == "chymotrypsin")
    {
      return ((n == 'F' || n == 'W' || n == 'Y' || n == 'L') && c != 'P') || n == '-' || c == '-';
    }
    else if (enz == "thermolysin")
    {
      return ((c == 'A' || c == 'F' || c == 'I' || c == 'L' || c == 'M'
              || c == 'V' || (n == 'R' && c == 'G')) && n != 'D' && n != 'E') || n == '-' || c == '-';
    }
    else if (enz == "proteinasek")
    {
      return (n == 'A' || n == 'E' || n == 'F' || n == 'I' || n == 'L'
             || n == 'T' || n == 'V' || n == 'W' || n == 'Y') || n == '-' || c == '-';
    }
    else if (enz == "pepsin")
    {
      return ((c == 'F' || c == 'L' || c == 'W' || c == 'Y' || n == 'F'
              || n == 'L' || n == 'W' || n == 'Y') && n != 'R') || n == '-' || c == '-';
    }
    else if (enz == "elastase")
    {
      return ((n == 'L' || n == 'V' || n == 'A' || n == 'G') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "lys-n")
    {
      return (c == 'K')
             || n == '-' || c == '-';
    }
    else if (enz == "lys-c")
    {
      return ((n == 'K') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "arg-c")
    {
      return ((n == 'R') && c != 'P')
             || n == '-' || c == '-';
    }
    else if (enz == "asp-n")
    {
      return (c == 'D')
             || n == '-' || c == '-';
    }
    else if (enz == "glu-c")
    {
      return ((n == 'E') && (c != 'P'))
             || n == '-' || c == '-';
    }
    else
    {
      return true;
    }
  }

  // Function adapted from Enzyme.h in Percolator converter
  // TODO: Use existing OpenMS functionality.
  Size countEnzymatic_(String peptide, string& enz)
  {
    Size count = 0;
    for (Size ix = 1; ix < peptide.size(); ++ix)
    {
      if (isEnz_(peptide[ix - 1], peptide[ix], enz))
      {
        ++count;
      }
    }
    return count;
  }

  //id <tab> label <tab> scannr <tab> calcmass <tab> expmass <tab> feature1 <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM
  void preparePin_(vector<PeptideIdentification>& peptide_ids, StringList& feature_set, std::string& enz, TextFile& txt, int min_charge, int max_charge)
  {
    for (vector<PeptideIdentification>::iterator it = peptide_ids.begin(); it != peptide_ids.end(); ++it)
    {
      String scan_identifier = getScanIdentifier_(it, peptide_ids.begin());
      Int scan_number = getScanNumber_(scan_identifier);
      
      double exp_mass = it->getMZ();
      for (vector<PeptideHit>::const_iterator jt = it->getHits().begin(); jt != it->getHits().end(); ++jt)
      {
        PeptideHit hit(*jt); // make a copy of the hit to store temporary features
        hit.setMetaValue("SpecId", scan_identifier);
        hit.setMetaValue("ScanNr", scan_number);
        
        if (!hit.metaValueExists("target_decoy") || hit.getMetaValue("target_decoy").toString().empty()) continue;
        
        int label = 1;
        if (String(hit.getMetaValue("target_decoy")).hasSubstring("decoy"))
        {
          label = -1;
        }
        hit.setMetaValue("Label", label);
        
        int charge = hit.getCharge();
        String unmodified_sequence = hit.getSequence().toUnmodifiedString();
        
        double calc_mass = hit.getSequence().getMonoWeight(Residue::Full, charge)/charge;
        hit.setMetaValue("CalcMass", calc_mass);
        
        
        hit.setMetaValue("ExpMass", exp_mass);
        hit.setMetaValue("mass", exp_mass);
        
        double score = hit.getScore();
        hit.setMetaValue("score", score);
        
        int peptide_length = unmodified_sequence.size();
        hit.setMetaValue("peplen", peptide_length);
        
        for (int i = min_charge; i <= max_charge; ++i)
        {
           hit.setMetaValue("charge" + String(i), charge == i);
        }
        
        bool enzN = isEnz_(hit.getPeptideEvidences().front().getAABefore(), unmodified_sequence.prefix(1)[0], enz);
        hit.setMetaValue("enzN", enzN);
        bool enzC = isEnz_(unmodified_sequence.suffix(1)[0], hit.getPeptideEvidences().front().getAAAfter(), enz);
        hit.setMetaValue("enzC", enzC);
        int enzInt = countEnzymatic_(unmodified_sequence, enz);
        hit.setMetaValue("enzInt", enzInt);
        
        double delta_mass = exp_mass - calc_mass;
        hit.setMetaValue("dm", delta_mass);
        
        double abs_delta_mass = abs(delta_mass);
        hit.setMetaValue("absdm", abs_delta_mass);
        
        //peptide
        String sequence = "";
        // just first peptide evidence
        String aa_before(hit.getPeptideEvidences().front().getAABefore());
        String aa_after(hit.getPeptideEvidences().front().getAAAfter());
        aa_before = aa_before=="["?'-':aa_before;
        aa_after = aa_after=="]"?'-':aa_after;
        sequence += aa_before; 
        sequence += "." + hit.getSequence().toString() + ".";
        sequence += aa_after;
        hit.setMetaValue("Peptide", sequence);
        
        //proteinId1
        StringList proteins;
        for (vector<PeptideEvidence>::const_iterator kt = hit.getPeptideEvidences().begin(); kt != hit.getPeptideEvidences().end(); ++kt)
        {
          proteins.push_back(kt->getProteinAccession());
        }
        hit.setMetaValue("Proteins", ListUtils::concatenate(proteins, '\t'));
        
        StringList feats;
        for (vector<String>::const_iterator feat = feature_set.begin(); feat != feature_set.end(); ++feat)
        {
        // Some Hits have no NumMatchedMainIons, and MeanError, etc. values. Have to ignore them!
          if (hit.metaValueExists(*feat))
          {
            feats.push_back(hit.getMetaValue(*feat).toString());
          }
        }
        if (feats.size() == feature_set.size())
        { // only if all feats were present add
          txt.addLine(ListUtils::concatenate(feats, '\t'));
        }
      }
    }
  }
  
  void readPoutAsMap_(String pout_file, std::map<String, PercolatorResult>& pep_map)
  {
    CsvFile csv_file(pout_file, '\t');
    StringList row;

    for (Size i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      PercolatorResult res(row);
      String spec_ref = res.PSMId + res.peptide;
      // retain only the best result in the unlikely case that a PSMId+peptide combination occurs multiple times
      if (pep_map.find(spec_ref) == pep_map.end())
      {
        pep_map.insert( map<String, PercolatorResult>::value_type ( spec_ref, res ) );
      }
    }
  }
  
  void readProteinPoutAsMap_(String pout_protein_file, std::map<String, PercolatorProteinResult>& protein_map)
  {
    CsvFile csv_file(pout_protein_file, '\t');
    StringList row;

    for (Size i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      StringList protein_accessions;
      row[0].split(",", protein_accessions);
      double qvalue = row[2].toDouble();
      double posterior_error_prob = row[3].toDouble();
      for (StringList::iterator it = protein_accessions.begin(); it != protein_accessions.end(); ++it) 
      {
        protein_map.insert( map<String, PercolatorProteinResult>::value_type ( *it, PercolatorProteinResult(*it, qvalue, posterior_error_prob ) ) );
      }
    }
  }
  
  ExitCodes readInputFiles_(StringList in_list, vector<PeptideIdentification>& all_peptide_ids, vector<ProteinIdentification>& all_protein_ids, bool isDecoy, bool& found_decoys, int& min_charge, int& max_charge)
  {
    for (StringList::iterator fit = in_list.begin(); fit != in_list.end(); ++fit)
    {
      String file_idx(distance(in_list.begin(), fit));
      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      String in = *fit;
      FileHandler fh;
      FileTypes::Type in_type = fh.getType(in);
      LOG_INFO << "Loading input file: " << in << endl;
      if (in_type == FileTypes::IDXML)
      {
        IdXMLFile().load(in, protein_ids, peptide_ids);
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        LOG_WARN << "Converting from mzid: possible loss of information depending on target format." << endl;
        MzIdentMLFile().load(in, protein_ids, peptide_ids);
      }
      //else catched by TOPPBase:registerInput being mandatory mzid or idxml

      //being paranoid about the presence of target decoy denominations, which are crucial to the percolator process
      for (vector<PeptideIdentification>::iterator pit = peptide_ids.begin(); pit != peptide_ids.end(); ++pit)
      {
        if (in_list.size() > 1)
        {
          String scan_identifier = getScanIdentifier_(pit, peptide_ids.begin());
          scan_identifier = "file=" + file_idx + "," + scan_identifier;
          pit->setMetaValue("spectrum_reference", scan_identifier);
        }
        for (vector<PeptideHit>::iterator pht = pit->getHits().begin(); pht != pit->getHits().end(); ++pht)
        {
          if (!pht->metaValueExists("target_decoy"))
          {
            if (isDecoy)
            {
              pht->setMetaValue("target_decoy", "decoy");
              found_decoys = true;
            }
            else
            {
              pht->setMetaValue("target_decoy", "target");
            }
          }
          else if (pht->getMetaValue("target_decoy").toString().hasSubstring("decoy"))
          {
            found_decoys = true;
          }
          
          if (pht->getCharge() > max_charge)
          {
            max_charge = pht->getCharge();
          }
          if (pht->getCharge() < min_charge)
          {
            min_charge = pht->getCharge();
          }

          // TODO: set min/max scores?
        }
      }
      
      //paranoia check if this comes from the same search engine! (only in the first proteinidentification of the first proteinidentifications vector vector)
      if (!all_protein_ids.empty()) 
      {
        if (protein_ids.front().getSearchEngine() != all_protein_ids.front().getSearchEngine())
        {
          writeLog_("Input files are not all from the same search engine: " + protein_ids.front().getSearchEngine() + " and " + all_protein_ids.front().getSearchEngine() + ". Use TOPP_PSMFeatureExtractor to merge results from different search engines if desired. Aborting!");
          return INCOMPATIBLE_INPUT_DATA;
        }
        
        bool identical_extra_features = true;
        ProteinIdentification::SearchParameters all_search_parameters = all_protein_ids.front().getSearchParameters();
        ProteinIdentification::SearchParameters search_parameters = protein_ids.front().getSearchParameters();
        if (all_search_parameters.metaValueExists("extra_features"))
        {
          StringList all_search_feature_list = ListUtils::create<String>(all_search_parameters.getMetaValue("extra_features").toString());
          set<String> all_search_feature_set(all_search_feature_list.begin(),all_search_feature_list.end());
          if (search_parameters.metaValueExists("extra_features"))
          {
            StringList search_feature_list = ListUtils::create<String>(search_parameters.getMetaValue("extra_features").toString());
            set<String> search_feature_set(search_feature_list.begin(), search_feature_list.end());
            identical_extra_features = (search_feature_set == all_search_feature_set);
          }
          else
          {
            identical_extra_features = false;
          }
        }
        if (!identical_extra_features) 
        {
          writeLog_("Input files do not have the same set of extra features from TOPP_PSMFeatureExtractor. Aborting!");
          return INCOMPATIBLE_INPUT_DATA;
        }
        
        if (protein_ids.front().getScoreType() != all_protein_ids.front().getScoreType())
        {
          LOG_WARN << "Warning: differing ScoreType between input files" << endl;
        }
        if (search_parameters.digestion_enzyme != all_search_parameters.digestion_enzyme)
        {
          LOG_WARN << "Warning: differing DigestionEnzyme between input files" << endl;
        }
        if (search_parameters.variable_modifications != all_search_parameters.variable_modifications)
        {
          LOG_WARN << "Warning: differing VarMods between input files" << endl;
        }
        if (search_parameters.fixed_modifications != all_search_parameters.fixed_modifications)
        {
          LOG_WARN << "Warning: differing FixMods between input files" << endl;
        }
        if (search_parameters.charges != all_search_parameters.charges)
        {
          LOG_WARN << "Warning: differing SearchCharges between input files" << endl;
        }
        if (search_parameters.fragment_mass_tolerance != all_search_parameters.fragment_mass_tolerance)
        {
          LOG_WARN << "Warning: differing FragTol between input files" << endl;
        }
        if (search_parameters.precursor_mass_tolerance != all_search_parameters.precursor_mass_tolerance)
        {
          LOG_WARN << "Warning: differing PrecTol between input files" << endl;
        }
      }
      LOG_INFO << "Merging peptide ids." << endl;
      all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
      LOG_INFO << "Merging protein ids." << endl;
      PercolatorFeatureSetHelper::mergeMULTISEProteinIds(all_protein_ids, protein_ids);
    }
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // general variables and data to perform PercolatorAdapter
    //-------------------------------------------------------------
    vector<PeptideIdentification> all_peptide_ids;
    vector<ProteinIdentification> all_protein_ids;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    const StringList in_list = getStringList_("in");
    const StringList in_decoy = getStringList_("in_decoy");
    LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << " & " << ListUtils::concatenate(in_decoy, ",") << " (decoy)" << endl;
    const String in_osw = getStringOption_("in_osw");
    const String osw_level = getStringOption_("osw_level");

    const String percolator_executable(getStringOption_("percolator_executable"));
    writeDebug_(String("Path to the percolator: ") + percolator_executable, 2);
    if (percolator_executable.empty())  //TODO? - TOPPBase::findExecutable after registerInputFile_("percolator_executable"... ???
    {
      writeLog_("No percolator executable specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    const String mzid_out(getStringOption_("mzid_out"));
    const String out(getStringOption_("out"));
    const String osw_out(getStringOption_("osw_out"));

    if (in_list.empty() && in_osw.empty())
    {
      writeLog_("Fatal error: no input file given (parameter 'in' or 'in_osw')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (mzid_out.empty() && out.empty() && osw_out.empty())
    {
      writeLog_("Fatal error: no output file given (parameter 'out' or 'mzid_out' or 'osw_out')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!in_osw.empty() && osw_out.empty())
    {
      writeLog_("Fatal error: OSW input requires OSW output.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!in_list.empty() && (out.empty() && mzid_out.empty()))
    {
      writeLog_("Fatal error: idXML/mzid input requires idXML/mzid output.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    bool peptide_level_fdrs = getFlag_("peptide-level-fdrs");
    bool protein_level_fdrs = getFlag_("protein-level-fdrs");  

    double ipf_max_peakgroup_pep = getDoubleOption_("ipf_max_peakgroup_pep");
    double ipf_max_transition_isotope_overlap = getDoubleOption_("ipf_max_transition_isotope_overlap");
    double ipf_min_transition_sn = getDoubleOption_("ipf_min_transition_sn");

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------

    string enz_str = getStringOption_("enzyme");
    
    // create temp directory to store percolator in file pin.tab temporarily
    String temp_directory_body = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory_body.toQString());
    }
    String txt_designator = File::getUniqueName();
    String pin_file(temp_directory_body + txt_designator + "_pin.tab");
    String pout_target_file(temp_directory_body + txt_designator + "_target_pout_psms.tab");
    String pout_decoy_file(temp_directory_body + txt_designator + "_decoy_pout_psms.tab");
    String pout_target_file_peptides(temp_directory_body + txt_designator + "_target_pout_peptides.tab");
    String pout_decoy_file_peptides(temp_directory_body + txt_designator + "_decoy_pout_peptides.tab");
    String pout_target_file_proteins(temp_directory_body + txt_designator + "_target_pout_proteins.tab");
    String pout_decoy_file_proteins(temp_directory_body + txt_designator + "_decoy_pout_proteins.tab");

    // prepare OSW I/O
    if (!in_osw.empty() && !osw_out.empty() && in_osw != osw_out)
    {
      // Copy input OSW to output OSW, because we want to retain all information
      remove(osw_out.c_str());
      if (!osw_out.empty())
      {
        std::ifstream  src(in_osw.c_str(), std::ios::binary);
        std::ofstream  dst(osw_out.c_str(), std::ios::binary);

        dst << src.rdbuf();
      }
    }

    // idXML or mzid input
    if (in_osw.empty())
    {
      //TODO introduce min/max charge to parameters for now take available range
      int max_charge = 0;
      int min_charge = 10;
      bool is_decoy = false;
      bool found_decoys = false;
      ExitCodes read_exit = readInputFiles_(in_list, all_peptide_ids, all_protein_ids, is_decoy, found_decoys, min_charge, max_charge);
      if (read_exit != EXECUTION_OK)
      {
        return read_exit;
      }
      
      if (!in_decoy.empty())
      {
        is_decoy = true;
        read_exit = readInputFiles_(in_decoy, all_peptide_ids, all_protein_ids, is_decoy, found_decoys, min_charge, max_charge);
        if (read_exit != EXECUTION_OK)
        {
          return read_exit;
        }
      }
      LOG_DEBUG << "Using min/max charges of " << min_charge << "/" << max_charge << endl;
      
      if (!found_decoys)
      {
        writeLog_("No decoys found, search results discrimination impossible. Aborting!");
        printUsage_();
        return INCOMPATIBLE_INPUT_DATA;
      }
      
      if (all_peptide_ids.empty())
      {
        writeLog_("No peptide hits found in input file. Aborting!");
        printUsage_();
        return INPUT_FILE_EMPTY;
      }
      
      if (all_protein_ids.empty())
      {
        writeLog_("No protein hits found in input file. Aborting!");
        printUsage_();
        return INPUT_FILE_EMPTY;
      }

      //-------------------------------------------------------------
      // prepare pin
      //-------------------------------------------------------------
      
      StringList feature_set;
      feature_set.push_back("SpecId");
      feature_set.push_back("Label");
      feature_set.push_back("ScanNr");
      feature_set.push_back("ExpMass");
      feature_set.push_back("CalcMass");
      feature_set.push_back("mass");
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
      
      ProteinIdentification::SearchParameters search_parameters = all_protein_ids.front().getSearchParameters();
      if (search_parameters.metaValueExists("extra_features"))
      {
        StringList extra_feature_set = ListUtils::create<String>(search_parameters.getMetaValue("extra_features").toString());
        feature_set.insert(feature_set.end(), extra_feature_set.begin(), extra_feature_set.end());
      }
      else if (getFlag_("generic-feature-set")) 
      {
        feature_set.push_back("score");
      } 
      else 
      {
        writeLog_("No search engine specific features found. Generate search engine specific features using PSMFeatureExtractor or set the -generic-features-set flag to override. Aborting!");
        printUsage_();
        return INCOMPATIBLE_INPUT_DATA;
      }
      
      feature_set.push_back("Peptide");
      feature_set.push_back("Proteins");
      
      LOG_DEBUG << "Writing percolator input file." << endl;
      TextFile txt;  
      txt.addLine(ListUtils::concatenate(feature_set, '\t'));
      preparePin_(all_peptide_ids, feature_set, enz_str, txt, min_charge, max_charge);
      txt.store(pin_file);
    }
    // OSW input
    else
    {
      LOG_DEBUG << "Writing percolator input file." << endl;
      TextFile txt;  
      std::stringstream pin_output;
      OSWFile().read(in_osw, osw_level, pin_output, ipf_max_peakgroup_pep, ipf_max_transition_isotope_overlap, ipf_min_transition_sn);
      txt << pin_output.str();
      txt.store(pin_file);
    }

    QStringList arguments;
    // Check all set parameters and get them into arguments StringList
    {    
      if (peptide_level_fdrs)
      { 
        arguments << "-r" << pout_target_file_peptides.toQString();
        arguments << "-B" << pout_decoy_file_peptides.toQString();
      }
      else
      {
        arguments << "-U";
      }
      arguments << "-m" << pout_target_file.toQString();
      arguments << "-M" << pout_decoy_file.toQString();
      
      if (protein_level_fdrs)
      {
        arguments << "-l" << pout_target_file_proteins.toQString();
        arguments << "-L" << pout_decoy_file_proteins.toQString();
        
        String fasta_file = getStringOption_("fasta");
        if (fasta_file.empty()) fasta_file = "auto";
        arguments << "-f" << fasta_file.toQString();
        
        String decoy_pattern = getStringOption_("decoy-pattern");
        if (decoy_pattern != "random") arguments << "-P" << decoy_pattern.toQString();
      }
      
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
      Int subset_max_train = getIntOption_("subset-max-train");
      if (subset_max_train > 0) arguments << "-N" << String(subset_max_train).toQString();
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

      arguments << pin_file.toQString();
    }
    writeLog_("Prepared percolator input.");

    //-------------------------------------------------------------
    // run percolator
    //-------------------------------------------------------------
    // Percolator execution with the executable and the arguments StringList
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
    map<String, PercolatorResult> pep_map;
    if (peptide_level_fdrs)
    {
      readPoutAsMap_(pout_target_file_peptides, pep_map);
      readPoutAsMap_(pout_decoy_file_peptides, pep_map);
    }
    else
    {
      readPoutAsMap_(pout_target_file, pep_map);
      readPoutAsMap_(pout_decoy_file, pep_map);
    }
    
    map<String, PercolatorProteinResult> protein_map;
    if (protein_level_fdrs)
    {
      readProteinPoutAsMap_(pout_target_file_proteins, protein_map);
      readProteinPoutAsMap_(pout_decoy_file_proteins, protein_map);
    }
    
    // As the percolator output file is not needed anymore, the temporary directory is going to be deleted
    if (this->debug_level_ < 5)
    {
      File::removeDirRecursively(temp_directory_body);
      LOG_WARN << "Removing temporary directory for Percolator in/output. Set debug level to >=5 to keep the temporary files." << endl;
    }
    else
    {
      LOG_WARN << "Keeping the temporary files at '" << temp_directory_body << "'. Set debug level to <5 to remove them." << endl;
    }

    // idXML or mzid input
    if (in_osw.empty())
    {
      // Add the percolator results to the peptide vector of the original input file
      //size_t c_debug = 0;
      size_t cnt = 0;
      String run_identifier = all_protein_ids.front().getIdentifier();
      for (vector<PeptideIdentification>::iterator it = all_peptide_ids.begin(); it != all_peptide_ids.end(); ++it)
      {
        it->setIdentifier(run_identifier);
        it->setScoreType("q-value");
        it->setHigherScoreBetter(false);
        
        String scan_identifier = getScanIdentifier_(it, all_peptide_ids.begin());
        
        //check each PeptideHit for compliance with one of the PercolatorResults (by sequence)
        for (vector<PeptideHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
        {
          String peptide_sequence = hit->getSequence().toString();
          String psm_identifier = scan_identifier + peptide_sequence;
          
          map<String, PercolatorResult>::iterator pr = pep_map.find(psm_identifier);
          if (pr != pep_map.end())
          {
            hit->setMetaValue("MS:1001492", pr->second.score);  // svm score
            hit->setMetaValue("MS:1001491", pr->second.qvalue);  // percolator q value
            hit->setMetaValue("MS:1001493", pr->second.posterior_error_prob);  // percolator pep
            hit->setScore(pr->second.qvalue);
            ++cnt;
          }
          else
          {
            hit->setScore(1.0); // set q-value to 1.0 if hit not found in results
          }
        }
      }
      //LOG_INFO << "No suitable PeptideIdentification for " << c_debug << " out of " << all_peptide_ids.size() << endl;
      LOG_INFO << "Suitable PeptideHits for " << cnt << " found." << endl;

      // TODO: There should only be 1 ProteinIdentification element in this vector, no need for a for loop
      for (vector<ProteinIdentification>::iterator it = all_protein_ids.begin(); it != all_protein_ids.end(); ++it)
      {      
        if (protein_level_fdrs)
        {
          //check each ProteinHit for compliance with one of the PercolatorProteinResults (by accession)
          for (vector<ProteinHit>::iterator hit = it->getHits().begin(); hit != it->getHits().end(); ++hit)
          {
            String protein_accession = hit->getAccession();        
            map<String, PercolatorProteinResult>::iterator pr = protein_map.find(protein_accession);
            if (pr != protein_map.end())
            {
              hit->setMetaValue("MS:1001491", pr->second.qvalue);  // percolator q value
              hit->setMetaValue("MS:1001493", pr->second.posterior_error_prob);  // percolator pep
              hit->setScore(pr->second.qvalue);
            }
            else
            {
              hit->setScore(1.0); // set q-value to 1.0 if hit not found in results
            }
          }
          it->setSearchEngine("Percolator");
          it->setScoreType("q-value");
          it->setHigherScoreBetter(false);
          it->sort();
        }
        
        //TODO add software percolator and PercolatorAdapter
        it->setMetaValue("percolator", "PercolatorAdapter");
        ProteinIdentification::SearchParameters search_parameters = it->getSearchParameters();
        
        search_parameters.setMetaValue("Percolator:peptide-level-fdrs", peptide_level_fdrs);
        search_parameters.setMetaValue("Percolator:protein-level-fdrs", protein_level_fdrs);
        search_parameters.setMetaValue("Percolator:generic-feature-set", getFlag_("generic-feature-set"));
        search_parameters.setMetaValue("Percolator:testFDR", getDoubleOption_("testFDR"));
        search_parameters.setMetaValue("Percolator:trainFDR", getDoubleOption_("trainFDR"));
        search_parameters.setMetaValue("Percolator:maxiter", getIntOption_("maxiter"));
        search_parameters.setMetaValue("Percolator:subset-max-train", getIntOption_("subset-max-train"));
        search_parameters.setMetaValue("Percolator:quick-validation", getFlag_("quick-validation"));
        search_parameters.setMetaValue("Percolator:weights", getStringOption_("weights"));
        search_parameters.setMetaValue("Percolator:init-weights", getStringOption_("init-weights"));
        search_parameters.setMetaValue("Percolator:default-direction", getStringOption_("default-direction"));
        search_parameters.setMetaValue("Percolator:cpos", getDoubleOption_("cpos"));
        search_parameters.setMetaValue("Percolator:cneg", getDoubleOption_("cneg"));
        search_parameters.setMetaValue("Percolator:unitnorm", getFlag_("unitnorm"));
        search_parameters.setMetaValue("Percolator:override", getFlag_("override"));
        search_parameters.setMetaValue("Percolator:seed", getIntOption_("seed"));
        search_parameters.setMetaValue("Percolator:doc", getIntOption_("doc"));
        search_parameters.setMetaValue("Percolator:klammer", getFlag_("klammer"));
        search_parameters.setMetaValue("Percolator:fasta", getStringOption_("fasta"));
        search_parameters.setMetaValue("Percolator:decoy-pattern", getStringOption_("decoy-pattern"));
        search_parameters.setMetaValue("Percolator:post-processing-tdc", getFlag_("post-processing-tdc"));
        
        it->setSearchParameters(search_parameters);
      }
      
      // Storing the PeptideHits with calculated q-value, pep and svm score
      if (!mzid_out.empty())
      {
        MzIdentMLFile().store(mzid_out.toQString().toStdString(), all_protein_ids, all_peptide_ids);
      }
      if (!out.empty())
      {
        IdXMLFile().store(out.toQString().toStdString(), all_protein_ids, all_peptide_ids);
      }
    }
    else
    {
      std::map< std::string, std::vector<double> > features;
      for (auto const &feat : pep_map)
      {
        features[feat.second.PSMId].push_back(feat.second.score);
        features[feat.second.PSMId].push_back(feat.second.qvalue);
        features[feat.second.PSMId].push_back(feat.second.posterior_error_prob);
      }
      OSWFile().write(osw_out, osw_level, features);
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
