// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/qfile.h>

#include <iostream>
#include <cmath>
#include <string>
#include <set>
//#include <typeinfo>

#include <boost/algorithm/clamp.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_PercolatorAdapter PercolatorAdapter

@brief PercolatorAdapter facilitates the input to, the call of and output integration of Percolator.
Percolator (http://percolator.ms/) is a tool to apply semi-supervised learning for peptide
identification from shotgun proteomics datasets.

@experimental This tool is work in progress and usage and input requirements might change.

<center>
  <table>
      <tr>
          <th ALIGN = "center"> pot. predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=2> &rarr; PercolatorAdapter &rarr;</td>
          <th ALIGN = "center"> pot. successor tools </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PSMFeatureExtractor </td>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
      </tr>
  </table>
</center>
<p>Percolator is search engine sensitive, i.e. it's input features vary,
depending on the search engine. Must be prepared beforehand. If you do not want
to use the specific features, use the generic_feature_set flag. Will incorporate
the score attribute of a PSM, so be sure, the score you want is set as main
score with @ref TOPP_IDScoreSwitcher . Be aware, that you might very well
experience a performance loss compared to the search engine specific features.
You can also perform protein inference with percolator when you activate the protein fdr parameter.
Additionally you need to set the enzyme setting.
We only read the q-value for protein groups since Percolator has a more elaborate FDR estimation.
For proteins we add q-value as main score and PEP as metavalue.
For PSMs you can choose the main score. Peptide level FDRs cannot be parsed and used yet.</p>

Multithreading: The thread parameter is passed to percolator.
Note: By default, a minimum of 3 threads is used (default of percolator) even if the number of threads
is set to e.g. 1 for backwards compatibility reasons. You can still force the usage of less than 3 threads
by setting the force flag.     

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
        std::size_t left_dot = row[4].find_first_of('.');
        std::size_t right_dot = row[4].find_last_of('.');
      
        OPENMS_PRECONDITION(left_dot < right_dot, "Peptide sequence encoding must have dot notation (e.g., A.PEPTIDER.C).")
 
        // retrieve pre and post AA, e.g., A and C in "A.PEPTIDE.C" or ".PEPTIDE."
        preAA = (left_dot == 0 || row[4][left_dot - 1] == '-') ? '[' : row[4][left_dot - 1];  // const char PeptideEvidence::N_TERMINAL_AA = '[';
        postAA = (right_dot + 1 < row[4].size() || row[4][right_dot + 1] == '-') ? ']' : row[4][right_dot + 1]; // const char PeptideEvidence::C_TERMINAL_AA = ']';

        // retrieve sequence between dots, e.g., PEPTIDE
        peptide = row[4].substr(left_dot + 1, (right_dot - 1) - (left_dot + 1) + 1);

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
          {
            return true;
          }
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
    static const bool force_openms_format(true);
        
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", !is_required);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML"));
    registerInputFileList_("in_decoy", "<files>", StringList(), "Input decoy file(s) in case of separate searches", !is_required);
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid,idXML"));
    registerInputFile_("in_osw", "<file>", "", "Input file in OSW format", !is_required);
    setValidFormats_("in_osw", ListUtils::create<String>("OSW"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,mzid,osw"));
    registerOutputFile_("out_pin", "<file>", "", "Write pin file (e.g., for debugging)", !is_required, is_advanced_option);
    setValidFormats_("out_pin", ListUtils::create<String>("tsv"), !force_openms_format);

    registerOutputFile_("out_pout_target", "<file>", "", "Write pout file (e.g., for debugging)", !is_required, is_advanced_option);
    setValidFormats_("out_pout_target", ListUtils::create<String>("tab"), !force_openms_format);
    registerOutputFile_("out_pout_decoy", "<file>", "", "Write pout file (e.g., for debugging)", !is_required, is_advanced_option);
    setValidFormats_("out_pout_decoy", ListUtils::create<String>("tab"), !force_openms_format);
    registerOutputFile_("out_pout_target_proteins", "<file>", "", "Write pout file (e.g., for debugging)", !is_required, is_advanced_option);
    setValidFormats_("out_pout_target_proteins", ListUtils::create<String>("tab"), !force_openms_format);
    registerOutputFile_("out_pout_decoy_proteins", "<file>", "", "Write pout file (e.g., for debugging)", !is_required, is_advanced_option);
    setValidFormats_("out_pout_decoy_proteins", ListUtils::create<String>("tab"), !force_openms_format);

    registerStringOption_("out_type", "<type>", "", "Output file type -- default: determined from file extension or content.", false);
    setValidStrings_("out_type", ListUtils::create<String>("mzid,idXML,osw"));
    String enzs = "no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin,trypsinp";
    registerStringOption_("enzyme", "<enzyme>", "trypsin", "Type of enzyme: "+enzs , !is_required);
    setValidStrings_("enzyme", ListUtils::create<String>(enzs));
    registerInputFile_("percolator_executable", "<executable>",
        // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "The Percolator executable. Provide a full or relative path, or make sure it can be found in your PATH environment.", is_required, !is_advanced_option, {"is_executable"}
    );
    registerFlag_("peptide_level_fdrs", "Calculate peptide-level FDRs instead of PSM-level FDRs.");
    registerFlag_("protein_level_fdrs", "Use the picked protein-level FDR to infer protein probabilities. Use the -fasta option and -decoy_pattern to set the Fasta file and decoy pattern.");
    
    registerStringOption_("osw_level", "<osw_level>", "ms2", "OSW: the data level selected for scoring.", !is_required);
    setValidStrings_("osw_level", StringList(OSWFile::names_of_oswlevel.begin(), OSWFile::names_of_oswlevel.end()));
    
    registerStringOption_("score_type", "<type>", "q-value", "Type of the peptide main score", false);
    setValidStrings_("score_type", ListUtils::create<String>("q-value,pep,svm"));

    //Advanced parameters
    registerFlag_("generic_feature_set", "Use only generic (i.e. not search engine specific) features. Generating search engine specific features for common search engines by PSMFeatureExtractor will typically boost the identification rate significantly.", is_advanced_option);
    registerIntOption_("subset_max_train", "<number>", 0, "Only train an SVM on a subset of <x> PSMs, and use the resulting score vector to evaluate the other PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal.", !is_required, is_advanced_option);
    registerDoubleOption_("cpos", "<value>", 0.0, "Cpos, penalty for mistakes made on positive examples. Set by cross validation if not specified.", !is_required, is_advanced_option);
    registerDoubleOption_("cneg", "<value>", 0.0, "Cneg, penalty for mistakes made on negative examples. Set by cross validation if not specified.", !is_required, is_advanced_option);
    registerDoubleOption_("testFDR", "<value>", 0.01, "False discovery rate threshold for evaluating best cross validation result and the reported end result.", !is_required, is_advanced_option);
    registerDoubleOption_("trainFDR", "<value>", 0.01, "False discovery rate threshold to define positive examples in training. Set to testFDR if 0.", !is_required, is_advanced_option);
    registerIntOption_("maxiter", "<number>", 10, "Maximal number of iterations", !is_required, is_advanced_option);
    registerIntOption_("nested_xval_bins", "<number>", 1, "Number of nested cross-validation bins in the 3 splits.", !is_required, is_advanced_option);
    registerFlag_("quick_validation", "Quicker execution by reduced internal cross-validation.", is_advanced_option);
    registerOutputFile_("weights", "<file>", "", "Output final weights to the given file", !is_required, is_advanced_option);
    setValidFormats_("weights", ListUtils::create<String>("tsv"), !force_openms_format);

    registerInputFile_("init_weights", "<file>", "", "Read initial weights to the given file", !is_required, is_advanced_option);
    setValidFormats_("init_weights", ListUtils::create<String>("tsv"), !force_openms_format);
    registerFlag_("static", "Use static model (requires init-weights parameter to be set)", is_advanced_option);

    registerStringOption_("default_direction", "<featurename>", "", "The most informative feature given as the feature name, can be negated to indicate that a lower value is better.", !is_required, is_advanced_option);
    registerIntOption_("verbose", "<level>", 2, "Set verbosity of output: 0=no processing info, 5=all.", !is_required, is_advanced_option);
    registerFlag_("unitnorm", "Use unit normalization [0-1] instead of standard deviation normalization", is_advanced_option);
    registerFlag_("test_each_iteration", "Measure performance on test set each iteration", is_advanced_option);
    registerFlag_("override", "Override error check and do not fall back on default score vector in case of suspect score vector", is_advanced_option);
    registerIntOption_("seed", "<value>", 1, "Setting seed of the random number generator.", !is_required, is_advanced_option);
    registerIntOption_("doc", "<value>", 0, "Include description of correct features", !is_required, is_advanced_option);
    registerFlag_("klammer", "Retention time features calculated as in Klammer et al. Only available if -doc is set", is_advanced_option);
    registerInputFile_("fasta", "<file>", "", "Provide the fasta file as the argument to this flag, which will be used for protein grouping based on an in-silico digest (only valid if option -protein_level_fdrs is active).", !is_required, is_advanced_option);
    setValidFormats_("fasta", ListUtils::create<String>("FASTA"));
    registerStringOption_("decoy_pattern", "<value>", "random", "Define the text pattern to identify the decoy proteins and/or PSMs, set this up if the label that identifies the decoys in the database is not the default (Only valid if option -protein_level_fdrs is active).", !is_required, is_advanced_option);
    registerFlag_("post_processing_tdc", "Use target-decoy competition to assign q-values and PEPs.", is_advanced_option);
    registerFlag_("train_best_positive", "Enforce that, for each spectrum, at most one PSM is included in the positive set during each training iteration. If the user only provides one PSM per spectrum, this filter will have no effect.", is_advanced_option);

    //OSW/IPF parameters
    registerDoubleOption_("ipf_max_peakgroup_pep", "<value>", 0.7, "OSW/IPF: Assess transitions only for candidate peak groups until maximum posterior error probability.", !is_required, is_advanced_option);
    registerDoubleOption_("ipf_max_transition_isotope_overlap", "<value>", 0.5, "OSW/IPF: Maximum isotope overlap to consider transitions in IPF.", !is_required, is_advanced_option);
    registerDoubleOption_("ipf_min_transition_sn", "<value>", 0, "OSW/IPF: Minimum log signal-to-noise level to consider transitions in IPF. Set -1 to disable this filter.", !is_required, is_advanced_option);
  }
  

  
  // Function adapted from Enzyme.h in Percolator converter
  // TODO: adapt to OpenMS enzymes. Use existing functionality in EnzymaticDigestion.
  bool isEnz_(const char& n, const char& c, string& enz)
  {
    if (enz == "trypsin")
    {
      return ((n == 'K' || n == 'R') && c != 'P') || n == '-' || c == '-';
    }
    else if (enz == "trypsinp")
    {
      return (n == 'K' || n == 'R') || n == '-' || c == '-';
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

  void readPoutAsMap_(const String& pout_file, std::map<String, PercolatorResult>& pep_map)
  {
    CsvFile csv_file(pout_file, '\t');
    StringList row;

    for (Size i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      PercolatorResult res(row);
      // note: Since we create our pin file in a way that the SpecID (=PSMId) is composed of filename + spectrum native id
      //  this will be passed through Percolator and we use it again to read it back in.
      String spec_ref = res.PSMId + res.peptide;
      writeDebug_("PSM identifier in pout file: " + spec_ref, 10);

      // retain only the best result in the unlikely case that a PSMId+peptide combination occurs multiple times
      if (pep_map.find(spec_ref) == pep_map.end())
      {
        pep_map.insert( map<String, PercolatorResult>::value_type ( spec_ref, res ) );
      }
    }
  }

  /// We only read the q-value for protein groups since Percolator has a more elaborate FDR estimation.
  /// For proteins we add q-value as main score and PEP as metavalue.
  void readProteinPoutAsMapAndAddGroups_(const String& pout_protein_file, std::map<String, PercolatorProteinResult>& protein_map, ProteinIdentification& protID_to_add_grps)
  {
    CsvFile csv_file(pout_protein_file, '\t');
    StringList row;
    std::vector<ProteinIdentification::ProteinGroup>& grps = protID_to_add_grps.getIndistinguishableProteins();

    for (Size i = 1; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, row);
      StringList protein_accessions;
      row[0].split(",", protein_accessions);
      double qvalue = row[2].toDouble();
      double posterior_error_prob = row[3].toDouble();
      for (const String& str : protein_accessions) 
      {
        protein_map.insert( map<String, PercolatorProteinResult>::value_type (str, PercolatorProteinResult(str, qvalue, posterior_error_prob ) ) );
      }

      ProteinIdentification::ProteinGroup grp;
      grp.probability = qvalue;
      grp.accessions = protein_accessions;
      grps.push_back(grp);
    }
  }
  
  ExitCodes readInputFiles_(const StringList& in_list, vector<PeptideIdentification>& all_peptide_ids, vector<ProteinIdentification>& all_protein_ids, bool isDecoy, bool& found_decoys, int& min_charge, int& max_charge)
  {
    for (StringList::const_iterator fit = in_list.begin(); fit != in_list.end(); ++fit)
    {
      String file_idx(distance(in_list.begin(), fit));
      vector<PeptideIdentification> peptide_ids;
      vector<ProteinIdentification> protein_ids;
      String in = *fit;
      FileTypes::Type in_type = FileHandler::getType(in);
      OPENMS_LOG_INFO << "Loading input file: " << in << endl;
      if (in_type == FileTypes::IDXML)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDXML});
      }
      else if (in_type == FileTypes::MZIDENTML)
      {
        OPENMS_LOG_WARN << "Converting from mzid: possible loss of information depending on target format." << endl;
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDXML});
      }
      //else catched by TOPPBase:registerInput being mandatory mzid or idxml
      if (protein_ids.empty())
      {
        throw Exception::ElementNotFound(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "File '" + in + "' has not ProteinIDRuns.");
      }
      else if (protein_ids.size() > 1)
      {
        throw Exception::InvalidValue(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "File '" + in + "' has more than one ProteinIDRun. This is currently not correctly handled."
            "Please use the merge_proteins_add_psms option if you used IDMerger. Alternatively, pass"
            " all original single-run idXML inputs as list to this tool.",
            "# runs: " + String(protein_ids.size()));
      }

      //being paranoid about the presence of target decoy denominations, which are crucial to the percolator process
      size_t index = 0;
      for (PeptideIdentification& pep_id : peptide_ids)
      {
        index++;
        if (in_list.size() > 1)
        {
          String scan_identifier = PercolatorInfile::getScanIdentifier(pep_id, index);
          scan_identifier = "file=" + file_idx + "," + scan_identifier;
          pep_id.setSpectrumReference( scan_identifier);
        }
        for (PeptideHit& hit : pep_id.getHits())
        {
          if (!hit.metaValueExists("target_decoy"))
          {
            if (isDecoy)
            {
              hit.setMetaValue("target_decoy", "decoy");
              found_decoys = true;
            }
            else
            {
              hit.setMetaValue("target_decoy", "target");
            }
          }
          else if (hit.getMetaValue("target_decoy").toString() == "decoy")
          {
            found_decoys = true;
          }
          
          if (hit.getCharge() > max_charge)
          {
            max_charge = hit.getCharge();
          }
          if (hit.getCharge() < min_charge)
          {
            min_charge = hit.getCharge();
          }

          // TODO: set min/max scores?
        }
      }
      
      //paranoia check if this comes from the same search engine! (only in the first proteinidentification of the first proteinidentifications vector vector)
      if (!all_protein_ids.empty()) 
      {
        if (protein_ids.front().getSearchEngine() != all_protein_ids.front().getSearchEngine())
        {
          writeLogError_("Input files are not all from the same search engine: " + protein_ids.front().getSearchEngine() + " and " + all_protein_ids.front().getSearchEngine() +
                         ". Use TOPP_PSMFeatureExtractor to merge results from different search engines if desired. Aborting!");
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
          writeLogError_("Input files do not have the same set of extra features from TOPP_PSMFeatureExtractor. Aborting!");
          return INCOMPATIBLE_INPUT_DATA;
        }
        
        if (protein_ids.front().getScoreType() != all_protein_ids.front().getScoreType())
        {
          OPENMS_LOG_WARN << "Warning: differing ScoreType between input files" << endl;
        }
        if (search_parameters.digestion_enzyme != all_search_parameters.digestion_enzyme)
        {
          OPENMS_LOG_WARN << "Warning: differing DigestionEnzyme between input files" << endl;
        }
        if (search_parameters.variable_modifications != all_search_parameters.variable_modifications)
        {
          OPENMS_LOG_WARN << "Warning: differing VarMods between input files" << endl;
        }
        if (search_parameters.fixed_modifications != all_search_parameters.fixed_modifications)
        {
          OPENMS_LOG_WARN << "Warning: differing FixMods between input files" << endl;
        }
        if (search_parameters.charges != all_search_parameters.charges)
        {
          OPENMS_LOG_WARN << "Warning: differing SearchCharges between input files" << endl;
        }
        if (search_parameters.fragment_mass_tolerance != all_search_parameters.fragment_mass_tolerance)
        {
          OPENMS_LOG_WARN << "Warning: differing FragTol between input files" << endl;
        }
        if (search_parameters.precursor_mass_tolerance != all_search_parameters.precursor_mass_tolerance)
        {
          OPENMS_LOG_WARN << "Warning: differing PrecTol between input files" << endl;
        }
      }
      OPENMS_LOG_INFO << "Merging peptide ids." << endl;
      all_peptide_ids.insert(all_peptide_ids.end(), peptide_ids.begin(), peptide_ids.end());
      OPENMS_LOG_INFO << "Merging protein ids." << endl;
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
    OPENMS_LOG_DEBUG << "Input file (of target?): " << ListUtils::concatenate(in_list, ",") << " & " << ListUtils::concatenate(in_decoy, ",") << " (decoy)" << endl;
    const String in_osw = getStringOption_("in_osw");
    const OSWFile::OSWLevel osw_level = (OSWFile::OSWLevel)Helpers::indexOf(OSWFile::names_of_oswlevel, getStringOption_("osw_level"));

    //output file names and types
    String out = getStringOption_("out");
    FileTypes::Type out_type = FileTypes::nameToType(getStringOption_("out_type"));

    if (out_type == FileTypes::UNKNOWN)
    {
      out_type = FileHandler::getTypeByFileName(out);
    }

    if (out_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Fatal error: Could not determine output file type!");
      return PARSE_ERROR;
    }

    const String percolator_executable(getStringOption_("percolator_executable"));
    
    if (in_list.empty() && in_osw.empty())
    {
      writeLogError_("Fatal error: no input file given (parameter 'in' or 'in_osw')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!in_list.empty() && !in_osw.empty())
    {
      writeLogError_("Fatal error: Provide either mzid/idXML or osw input files (parameter 'in' or 'in_osw')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (out.empty())
    {
      writeLogError_("Fatal error: no output file given (parameter 'out')");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!in_osw.empty() && out_type != FileTypes::OSW)
    {
      writeLogError_("Fatal error: OSW input requires OSW output.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }

    if (!in_list.empty() && out_type == FileTypes::OSW)
    {
      writeLogError_("Fatal error: idXML/mzid input requires idXML/mzid output.");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    bool peptide_level_fdrs = getFlag_("peptide_level_fdrs");
    bool protein_level_fdrs = getFlag_("protein_level_fdrs");  

    Int description_of_correct = getIntOption_("doc");

    double ipf_max_peakgroup_pep = getDoubleOption_("ipf_max_peakgroup_pep");
    double ipf_max_transition_isotope_overlap = getDoubleOption_("ipf_max_transition_isotope_overlap");
    double ipf_min_transition_sn = getDoubleOption_("ipf_min_transition_sn");

    //-------------------------------------------------------------
    // read input
    //-------------------------------------------------------------

    string enz_str = getStringOption_("enzyme");
    
    // create temp directory to store percolator in file pin.tab temporarily
    File::TempDir tmp_dir(debug_level_ >= 2);
    
    String txt_designator = File::getUniqueName();
    String pin_file;
    if (getStringOption_("out_pin").empty())
    {
      pin_file = tmp_dir.getPath() + txt_designator + "_pin.tab";
    }
    else
    {
      pin_file = getStringOption_("out_pin");
    }
    
    String pout_target_file(tmp_dir.getPath() + txt_designator + "_target_pout_psms.tab");
    String pout_decoy_file(tmp_dir.getPath() + txt_designator + "_decoy_pout_psms.tab");
    String pout_target_file_peptides(tmp_dir.getPath() + txt_designator + "_target_pout_peptides.tab");
    String pout_decoy_file_peptides(tmp_dir.getPath() + txt_designator + "_decoy_pout_peptides.tab");
    String pout_target_file_proteins(tmp_dir.getPath() + txt_designator + "_target_pout_proteins.tab");
    String pout_decoy_file_proteins(tmp_dir.getPath() + txt_designator + "_decoy_pout_proteins.tab");

    // prepare OSW I/O
    if (out_type == FileTypes::OSW && in_osw != out)
    {
      // Copy input OSW to output OSW, because we want to retain all information
      remove(out.c_str());
      if (!out.empty())
      {
        std::ifstream  src(in_osw.c_str(), std::ios::binary);
        std::ofstream  dst(out.c_str(), std::ios::binary);

        dst << src.rdbuf();
      }
    }

    // idXML or mzid input
    if (out_type != FileTypes::OSW)
    {
      //TODO introduce min/max charge to parameters. For now take available range
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
      OPENMS_LOG_DEBUG << "Using min/max charges of " << min_charge << "/" << max_charge << endl;
      
      if (!found_decoys)
      {
        writeLogError_("No decoys found, search results discrimination impossible. Aborting!");
        printUsage_();
        return INCOMPATIBLE_INPUT_DATA;
      }
      
      if (all_peptide_ids.empty())
      {
        writeLogError_("No peptide hits found in input file. Aborting!");
        printUsage_();
        return INPUT_FILE_EMPTY;
      }
      
      if (all_protein_ids.empty())
      {
        writeLogError_("No protein hits found in input file. Aborting!");
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
      if (description_of_correct != 0)
      {
        feature_set.push_back("retentiontime");
        feature_set.push_back("deltamass");
      }
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
      else if (getFlag_("generic_feature_set")) 
      {
        feature_set.push_back("score");
      } 
      else 
      {
        writeLogError_("No search engine specific features found. Generate search engine specific features using PSMFeatureExtractor or set the -generic-features-set flag to override. Aborting!");
        printUsage_();
        return INCOMPATIBLE_INPUT_DATA;
      }
      
      feature_set.push_back("Peptide");
      feature_set.push_back("Proteins");
      
      OPENMS_LOG_DEBUG << "Writing percolator input file." << endl;
      PercolatorInfile::store(pin_file, all_peptide_ids, feature_set, enz_str, min_charge, max_charge);
    }
    // OSW input
    else
    {
      OPENMS_LOG_DEBUG << "Writing percolator input file." << endl;
      std::ofstream pin_output(pin_file);
      OSWFile::readToPIN(in_osw, osw_level, pin_output, ipf_max_peakgroup_pep, ipf_max_transition_isotope_overlap, ipf_min_transition_sn);
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
        if (fasta_file.empty())
        {
          fasta_file = "auto";
        }
        arguments << "-f" << fasta_file.toQString();

        arguments << "-z" << String(enz_str).toQString();

        String decoy_pattern = getStringOption_("decoy_pattern");
        if (decoy_pattern != "random") arguments << "-P" << decoy_pattern.toQString();
      }
      
      int cv_threads = getIntOption_("threads"); // pass-through of OpenMS thread parameter

      if (cv_threads != 3) // default in percolator is 3
      {
        // If a lower than default value (3) is chosen the user needs to enforce it.
        // This ensures that existing workflows (which implicitly used 3 threads) don't slow down
        // if e.g. the OpenMS version and this adapter is updated.
        if (cv_threads > 3 || getFlag_("force"))
        { 
          arguments << "--num-threads" << String(cv_threads).toQString();
        }
      }
      
      double cpos = getDoubleOption_("cpos");
      double cneg = getDoubleOption_("cneg");
      if (cpos != 0.0)
      {
        arguments << "-p" << String(cpos).toQString();
      }
      if (cneg != 0.0)
      {
        arguments << "-n" << String(cneg).toQString();
      }
      double train_FDR = getDoubleOption_("trainFDR");
      double test_FDR = getDoubleOption_("testFDR");
      if (train_FDR != 0.01)
      {
        arguments << "-F" << String(train_FDR).toQString();
      }
      if (test_FDR != 0.01)
      {
        arguments << "-t" << String(test_FDR).toQString();
      }
      Int max_iter = getIntOption_("maxiter");
      if (max_iter != 10)
      {
        arguments << "-i" << String(max_iter).toQString();
      }
      Int subset_max_train = getIntOption_("subset_max_train");
      if (subset_max_train > 0)
      {
        arguments << "-N" << String(subset_max_train).toQString();
      }
      if (getFlag_("quick_validation"))
      {
        arguments << "-x";
      }
      if (getFlag_("post_processing_tdc"))
      {
        arguments << "-Y";
      }
      if (getFlag_("train_best_positive"))
      {
        arguments << "--train-best-positive";
      }
      if (getFlag_("static"))
      {
        arguments << "--static";
      }
      Int nested_xval_bins = getIntOption_("nested_xval_bins");
      if (nested_xval_bins > 1)
      {
        arguments << "--nested-xval-bins" << String(nested_xval_bins).toQString();
      }
      String weights_file = getStringOption_("weights");
      String init_weights_file = getStringOption_("init_weights");
      String default_search_direction = getStringOption_("default_direction");
      if (!weights_file.empty())
      {
        arguments << "-w" << weights_file.toQString();
      }
      if (!init_weights_file.empty())
      {
        arguments << "-W" << init_weights_file.toQString();
      }
      if (!default_search_direction.empty())
      {
        arguments << "-V" << default_search_direction.toQString();
      }
      Int verbose_level = getIntOption_("verbose");
      if (verbose_level != 2)
      {
        arguments << "-v" << String(verbose_level).toQString();
      }
      if (getFlag_("unitnorm"))
      {
        arguments << "-u";
      }
      if (getFlag_("test_each_iteration"))
      {
        arguments << "-R";
      }
      if (getFlag_("override"))
      {
        arguments << "-O";
      }
      Int seed = getIntOption_("seed");
      if (seed != 1)
      {
        arguments << "-S" << String(seed).toQString();
      }
      if (getFlag_("klammer"))
      {
        arguments << "-K";
      }
      if (description_of_correct != 0)
      {
        arguments << "-D" << String(description_of_correct).toQString();
      }
      arguments << pin_file.toQString();
    }
    writeLogInfo_("Prepared percolator input.");

    //-------------------------------------------------------------
    // run percolator
    //-------------------------------------------------------------
    // Percolator execution with the executable and the arguments StringList
    TOPPBase::ExitCodes exit_code = runExternalProcess_(percolator_executable.toQString(), arguments);
    if (exit_code != EXECUTION_OK)
    {
      return exit_code;
    }

    //-------------------------------------------------------------
    // reintegrate pout results
    //-------------------------------------------------------------
    // when percolator finished calculation, it stores the results -r option (with or without -U) or -m (which seems to be not working)
    //  WARNING: The -r option cannot be used in conjunction with -U: no peptide level statistics are calculated, redirecting PSM level statistics to provided file instead.
    map<String, PercolatorResult> pep_map;
    String pout_target = getStringOption_("out_pout_target");
    String pout_decoy = getStringOption_("out_pout_decoy");
    String pout_target_proteins = getStringOption_("out_pout_target_proteins");
    String pout_decoy_proteins = getStringOption_("out_pout_decoy_proteins");

    if (peptide_level_fdrs)
    {
      readPoutAsMap_(pout_target_file_peptides, pep_map);
      readPoutAsMap_(pout_decoy_file_peptides, pep_map);

      // copy file in tmp folder to output
      if (!pout_target.empty())
      {
        QFile::copy(pout_target_file_peptides.toQString(), pout_target.toQString());
      }
      if (!pout_decoy.empty())
      {
        QFile::copy(pout_decoy_file_peptides.toQString(), pout_decoy.toQString());
      }
    }
    else
    {
      readPoutAsMap_(pout_target_file, pep_map);
      readPoutAsMap_(pout_decoy_file, pep_map);

      // copy file in tmp folder to output
      if (!pout_target.empty())
      {
        QFile::copy(pout_target_file.toQString(), pout_target.toQString());
      }
      if (!pout_decoy.empty())
      {
        QFile::copy(pout_decoy_file.toQString(), pout_decoy.toQString());
      }
    }
    
    map<String, PercolatorProteinResult> protein_map;
    if (protein_level_fdrs)
    {
      readProteinPoutAsMapAndAddGroups_(pout_target_file_proteins, protein_map, all_protein_ids[0]);
      readProteinPoutAsMapAndAddGroups_(pout_decoy_file_proteins, protein_map, all_protein_ids[0] );

      // copy file in tmp folder to output filename
      if (!pout_target_proteins.empty())
      {
        QFile::copy(pout_target_file_proteins.toQString(), pout_target_proteins.toQString());
      }
      if (!pout_decoy_proteins.empty())
      {
        QFile::copy(pout_target_file_proteins.toQString(), pout_decoy_proteins.toQString());
      }
    }

    // idXML or mzid input
    if (in_osw.empty())
    {
      // Add the percolator results to the peptide vector of the original input file
      //size_t c_debug = 0;
      size_t cnt = 0;
      String run_identifier = all_protein_ids.front().getIdentifier();
      const String scoreType = getStringOption_("score_type");
      size_t index = 0;
      for (PeptideIdentification& pep_id : all_peptide_ids)
      {
        String old_score_type{pep_id.getScoreType()}; // copy because we modify the score type below
        index++;
        pep_id.setIdentifier(run_identifier);
        if (scoreType == "pep")
        {
          pep_id.setScoreType("Posterior Error Probability");
        }
        else
        {
          //TODO we should make a difference between peptide-level q-values and psm-level q-values!
          // I am just not changing it right now, because a lot of tools currently depend on
          // the score being exactly "q-value"
          pep_id.setScoreType(scoreType);
        }
        pep_id.setHigherScoreBetter(scoreType == "svm");
        
        String scan_identifier = PercolatorInfile::getScanIdentifier(pep_id, index);
        String file_identifier = pep_id.getMetaValue("file_origin", String());
        file_identifier += (String)pep_id.getMetaValue("id_merge_index", String());

        //check each PeptideHit for compliance with one of the PercolatorResults (by sequence)
        for (PeptideHit& hit : pep_id.getHits())
        {
          String peptide_sequence = hit.getSequence().toBracketString(false, true);
          String psm_identifier = file_identifier + scan_identifier + peptide_sequence;

          //Only for super debug
          writeDebug_("PSM identifier in PeptideHit: " + psm_identifier, 10);
 
          map<String, PercolatorResult>::iterator pr = pep_map.find(psm_identifier);
          if (pr != pep_map.end())
          {
            hit.setMetaValue(old_score_type, hit.getScore());  // old search engine "main" score as metavalue
            hit.setMetaValue("MS:1001492", pr->second.score);  // svm score
            hit.setMetaValue("MS:1001491", pr->second.qvalue);  // percolator q value
            hit.setMetaValue("MS:1001493", pr->second.posterior_error_prob);  // percolator pep

            if (scoreType == "q-value")
            {
              hit.setScore(pr->second.qvalue);
            }
            else if (scoreType == "pep")
            {
              hit.setScore(pr->second.posterior_error_prob);
            }
            else if (scoreType == "svm")
            {
              hit.setScore(pr->second.score);
            }

            ++cnt;
          }
          else
          {
            // If the input contains multiple PSMs per spectrum, Percolator only reports the top scoring PSM.
            // The remaining PSMs should be reported as not identified
            writeDebug_("PSM identifier " + psm_identifier + " not found in peptide map", 10);

            // Percolator's svm score is scaled such that 0.0 is the score at the chosen FDR threshold,
            // with positive scores representing PSMs under the FDR threshold (i.e. identified)
            // and negative scores PSMs above the FDR threshold (i.e. not identified);
            // -100.0 is typically more than low enough to represent a confidently non-identified PSM.
            hit.setMetaValue("MS:1001492", -100.0);  // svm score
            hit.setMetaValue("MS:1001491", 1.0);  // percolator q value
            hit.setMetaValue("MS:1001493", 1.0);  // percolator pep

            if (scoreType == "q-value" || scoreType == "pep")
            {
              hit.setScore(1.0); // set q-value or PEP to 1.0 if hit not found in results
            }
            else if (scoreType == "svm")
            {
              hit.setScore(-100.0); // set SVM score to -100.0 if hit not found in results
            }
          }
        }
      }

      if (!peptide_level_fdrs)
      {
      OPENMS_LOG_INFO << "PSM-level FDR: All PSMs are returned by percolator. Reannotating all PSMs in input data with percolator output." << endl;
      }
      else
      {
      OPENMS_LOG_INFO << "Peptide-level FDR: Only the best PSM per Peptide is returned by percolator. Reannotating the best PSM in input data with percolator output." << endl;
      }
      OPENMS_LOG_INFO << "Scores of all other PSMs will be set to 1.0." << endl;
      OPENMS_LOG_INFO << cnt << " suitable PeptideHits of " << all_peptide_ids.size() <<  " PSMs were reannotated." << endl;

      // TODO: There should only be 1 ProteinIdentification element in this vector, no need for a for loop
      for (ProteinIdentification& prot_id_run : all_protein_ids)
      {
        // it is not a real search engine but we set it so that we know that
        // scores were postprocessed
        prot_id_run.setSearchEngine("Percolator");
        prot_id_run.setSearchEngineVersion("3.05"); // TODO: read from percolator
        if (protein_level_fdrs)
        {
          //check each ProteinHit for compliance with one of the PercolatorProteinResults (by accession)
          for (ProteinHit& protein : prot_id_run.getHits())
          {
            String protein_accession = protein.getAccession();        
            map<String, PercolatorProteinResult>::iterator pr = protein_map.find(protein_accession);
            if (pr != protein_map.end())
            {
              protein.setMetaValue("MS:1001493", pr->second.posterior_error_prob);  // percolator pep
              protein.setScore(pr->second.qvalue);
              //remove to mark the protein as mapped. We can safely assume that every protein
              // only occurs once in Percolator output.
              protein_map.erase(pr);
            }
            else
            {
              protein.setScore(1.0); // set q-value to 1.0 if hit not found in results
              protein.setMetaValue("MS:1001493", 1.0);  // same for percolator pep
            }
          }
          if (protein_level_fdrs)
          {
            prot_id_run.setInferenceEngine("Percolator");
            prot_id_run.setInferenceEngineVersion("3.05");
          }
          prot_id_run.setScoreType("q-value");
          prot_id_run.setHigherScoreBetter(false);
          prot_id_run.sort();
        }
        
        if (!protein_map.empty())  //there remain unmapped proteins from Percolator
        {
          for (const auto& prot : protein_map)
          {
                  if (prot.second.posterior_error_prob < 1.0) //actually present according to Percolator
            {
                    OPENMS_LOG_WARN << "Warning: Protein " << prot.first << " reported by Percolator with non-zero probability was"
                "not present in the input idXML. Ignoring to keep consistency of the PeptideIndexer settings..";
            }
          }
          // filter groups that might contain these unmapped proteins so we do not get errors while writing our output.
          IDFilter::updateProteinGroups(all_protein_ids[0].getIndistinguishableProteins(), all_protein_ids[0].getHits());
        }

        //TODO add software percolator and PercolatorAdapter
        prot_id_run.setMetaValue("percolator", "PercolatorAdapter");
        ProteinIdentification::SearchParameters search_parameters = prot_id_run.getSearchParameters();
        
        search_parameters.setMetaValue("Percolator:peptide_level_fdrs", peptide_level_fdrs);
        search_parameters.setMetaValue("Percolator:protein_level_fdrs", protein_level_fdrs);
        search_parameters.setMetaValue("Percolator:generic_feature_set", getFlag_("generic_feature_set"));
        search_parameters.setMetaValue("Percolator:testFDR", getDoubleOption_("testFDR"));
        search_parameters.setMetaValue("Percolator:trainFDR", getDoubleOption_("trainFDR"));
        search_parameters.setMetaValue("Percolator:maxiter", getIntOption_("maxiter"));
        search_parameters.setMetaValue("Percolator:subset_max_train", getIntOption_("subset_max_train"));
        search_parameters.setMetaValue("Percolator:quick_validation", getFlag_("quick_validation"));
        search_parameters.setMetaValue("Percolator:static", getFlag_("static"));
        search_parameters.setMetaValue("Percolator:weights", getStringOption_("weights"));
        search_parameters.setMetaValue("Percolator:init_weights", getStringOption_("init_weights"));
        search_parameters.setMetaValue("Percolator:default_direction", getStringOption_("default_direction"));
        search_parameters.setMetaValue("Percolator:cpos", getDoubleOption_("cpos"));
        search_parameters.setMetaValue("Percolator:cneg", getDoubleOption_("cneg"));
        search_parameters.setMetaValue("Percolator:unitnorm", getFlag_("unitnorm"));
        search_parameters.setMetaValue("Percolator:override", getFlag_("override"));
        search_parameters.setMetaValue("Percolator:seed", getIntOption_("seed"));
        search_parameters.setMetaValue("Percolator:doc", getIntOption_("doc"));
        search_parameters.setMetaValue("Percolator:klammer", getFlag_("klammer"));
        search_parameters.setMetaValue("Percolator:fasta", getStringOption_("fasta"));
        search_parameters.setMetaValue("Percolator:decoy_pattern", getStringOption_("decoy_pattern"));
        search_parameters.setMetaValue("Percolator:post_processing_tdc", getFlag_("post_processing_tdc"));
        search_parameters.setMetaValue("Percolator:train_best_positive", getFlag_("train_best_positive"));
        
        prot_id_run.setSearchParameters(search_parameters);
      }
      // Storing the PeptideHits with calculated q-value, pep and svm score
      FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {FileTypes::IDXML, FileTypes::MZIDENTML});
    }
    else
    {
      std::map< std::string, OSWFile::PercolatorFeature > features;
      for (auto const &feat : pep_map)
      {
        features.emplace(std::piecewise_construct,
                         std::forward_as_tuple(feat.second.PSMId),
                         std::forward_as_tuple(feat.second.score, feat.second.qvalue, feat.second.posterior_error_prob));
      }
      OSWFile::writeFromPercolator(out, osw_level, features);
    }

    writeLogInfo_("PercolatorAdapter finished successfully!");
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  PercolatorAdapter tool;

  return tool.main(argc, argv);
}

/// @endcond
