// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Andreas Simon, Mathias Walzer, Matthew The $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorFeatureSetHelper.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/PercolatorInfile.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <boost/regex.hpp>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>


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

  /// Strip PIN-format meta values that PercolatorInfile::stampPinFeaturesOnHits
  /// leaves on each hit during the in-process feature-matrix build. Without
  /// this, the in-process output idXML carries 17+ internal columns (SpecId,
  /// ScanNr, Label, ExpMass, CalcMass, mass, peplen, charge<N>, enzN/C/Int,
  /// dm, absdm, deltamass, retentiontime, score, Peptide, Proteins) that
  /// downstream tools shouldn't see. The subprocess path doesn't need this:
  /// it writes those fields to a temporary PIN file and never stamps them
  /// onto hits in the first place.
  static void stripPinFeatureMetaValues_(
    PeptideIdentificationList& pep_ids,
    int min_charge, int max_charge)
  {
    static const std::vector<String> KEYS_TO_STRIP = {
      "SpecId", "ScanNr", "Label", "ExpMass", "CalcMass",
      "mass", "peplen", "deltamass", "retentiontime", "score",
      "enzN", "enzC", "enzInt", "dm", "absdm",
      "Peptide", "Proteins"
    };
    for (auto& pid : pep_ids.getData())
    {
      for (auto& hit : pid.getHits())
      {
        for (const String& k : KEYS_TO_STRIP) hit.removeMetaValue(k);
        for (int c = min_charge; c <= max_charge; ++c)
        {
          hit.removeMetaValue("charge" + String(c));
        }
      }
    }
  }

  /// Stamp "this run was post-processed by Percolator via PercolatorAdapter"
  /// metadata on every ProteinIdentification: search_engine, version, the
  /// "percolator" marker UserParam, and the 23 Percolator:* SearchParameter
  /// UserParams recording the adapter's CLI configuration. Both backends
  /// (in-process and subprocess) must produce equivalent metadata so that
  /// downstream tools see the same provenance regardless of which path ran.
  ///
  /// @param prot_ids       all ProteinIdentifications from the input
  /// @param peptide_level_fdrs / @param protein_level_fdrs  CLI flag values
  /// @param version_string  e.g. "3.08" or "in-process eb157f7"
  /// @param protein_map     populated only when protein_level_fdrs=true on
  ///                        the subprocess path (provides per-protein qvalue/
  ///                        PEP from Percolator's protein-level results to
  ///                        be stamped back onto each ProteinHit). nullptr
  ///                        on the in-process path skips that block.
  void stampPercolatorAdapterMetadata_(
    vector<ProteinIdentification>& prot_ids,
    bool peptide_level_fdrs,
    bool protein_level_fdrs,
    const String& version_string,
    map<String, PercolatorProteinResult>* protein_map = nullptr) const
  {
    // Guard: it would be a silent metadata lie to stamp
    // `Percolator:protein_level_fdrs=true` without actually running the
    // protein-level FDR mapping block (which only happens when protein_map
    // is non-null). Catch any future caller that mismatches the two.
    assert(!protein_level_fdrs || protein_map != nullptr);

    for (ProteinIdentification& prot_id_run : prot_ids)
    {
      // It is not a real search engine but we set it so downstream tools
      // can see that scores were post-processed.
      prot_id_run.setSearchEngine("Percolator");
      prot_id_run.setSearchEngineVersion(version_string);

      if (protein_level_fdrs && protein_map)
      {
        // Map Percolator protein-level qvalue/PEP back onto each ProteinHit
        // by accession.
        for (ProteinHit& protein : prot_id_run.getHits())
        {
          String protein_accession = protein.getAccession();
          map<String, PercolatorProteinResult>::iterator pr =
            protein_map->find(protein_accession);
          if (pr != protein_map->end())
          {
            protein.setMetaValue("MS:1001493", pr->second.posterior_error_prob);
            protein.setScore(pr->second.qvalue);
            // Mark mapped — safe because each protein appears once in the
            // Percolator output.
            protein_map->erase(pr);
          }
          else
          {
            protein.setScore(1.0);
            protein.setMetaValue("MS:1001493", 1.0);
          }
        }
        prot_id_run.setInferenceEngine("Percolator");
        prot_id_run.setInferenceEngineVersion(version_string);
        prot_id_run.setScoreType("q-value");
        prot_id_run.setHigherScoreBetter(false);
        prot_id_run.sort();

        if (!protein_map->empty())
        {
          for (const auto& prot : *protein_map)
          {
            if (prot.second.posterior_error_prob < 1.0)
            {
              OPENMS_LOG_WARN << "Warning: Protein " << prot.first
                << " reported by Percolator with non-zero probability was"
                << " not present in the input idXML. Ignoring to keep "
                << "consistency of the PeptideIndexer settings.";
            }
          }
          IDFilter::updateProteinGroups(
            prot_ids[0].getIndistinguishableProteins(), prot_ids[0].getHits());
        }
      }

      prot_id_run.setMetaValue("percolator", "PercolatorAdapter");
      ProteinIdentification::SearchParameters sp = prot_id_run.getSearchParameters();
      sp.setMetaValue("Percolator:peptide_level_fdrs",   peptide_level_fdrs);
      sp.setMetaValue("Percolator:protein_level_fdrs",   protein_level_fdrs);
      sp.setMetaValue("Percolator:generic_feature_set",  getFlag_("generic_feature_set"));
      sp.setMetaValue("Percolator:testFDR",              getDoubleOption_("testFDR"));
      sp.setMetaValue("Percolator:trainFDR",             getDoubleOption_("trainFDR"));
      sp.setMetaValue("Percolator:maxiter",              getIntOption_("maxiter"));
      sp.setMetaValue("Percolator:subset_max_train",     getIntOption_("subset_max_train"));
      sp.setMetaValue("Percolator:quick_validation",     getFlag_("quick_validation"));
      sp.setMetaValue("Percolator:static",               getFlag_("static"));
      sp.setMetaValue("Percolator:weights",              getStringOption_("weights"));
      sp.setMetaValue("Percolator:init_weights",         getStringOption_("init_weights"));
      sp.setMetaValue("Percolator:default_direction",    getStringOption_("default_direction"));
      sp.setMetaValue("Percolator:cpos",                 getDoubleOption_("cpos"));
      sp.setMetaValue("Percolator:cneg",                 getDoubleOption_("cneg"));
      sp.setMetaValue("Percolator:unitnorm",             getFlag_("unitnorm"));
      sp.setMetaValue("Percolator:override",             getFlag_("override"));
      sp.setMetaValue("Percolator:seed",                 getIntOption_("seed"));
      sp.setMetaValue("Percolator:doc",                  getIntOption_("doc"));
      sp.setMetaValue("Percolator:klammer",              getFlag_("klammer"));
      sp.setMetaValue("Percolator:fasta",                getStringOption_("fasta"));
      sp.setMetaValue("Percolator:decoy_pattern",        getStringOption_("decoy_pattern"));
      sp.setMetaValue("Percolator:post_processing_tdc",  getFlag_("post_processing_tdc"));
      sp.setMetaValue("Percolator:train_best_positive",  getFlag_("train_best_positive"));
      prot_id_run.setSearchParameters(sp);
    }
  }

  /// Apply the user-facing post-filters (-score:fdr, -best_per_spectrum_only)
  /// to in-memory identifications before they are written.  Used by both the
  /// subprocess path and the in-process path so the two backends produce
  /// equivalent output.
  void applyPostFilters_(PeptideIdentificationList& all_peptide_ids)
  {
    const double score_fdr = getDoubleOption_("score:fdr");
    if (score_fdr < 1.0)
    {
      // -score:fdr is always interpreted as a q-value cutoff, independent of -score_type.
      // The 3-arg overload uses IDScoreSwitcherAlgorithm: if the main score already is the
      // q-value it filters on it directly; otherwise it locates the q-value in the hit's
      // meta values (PercolatorAdapter stamps it as MS:1001491 in both backends, see the
      // svm/PEP code paths) and filters on that. This keeps the user's chosen main score
      // (e.g. svm or PEP) for downstream tools while still applying the FDR cutoff.
      IDFilter::filterHitsByScore(all_peptide_ids.getData(), score_fdr,
                                  IDScoreSwitcherAlgorithm::ScoreType::QVAL);
      IDFilter::removeEmptyIdentifications(all_peptide_ids);
      OPENMS_LOG_INFO << "Applied score:fdr cutoff " << score_fdr
                      << "; remaining peptide identifications: " << all_peptide_ids.size() << std::endl;
      if (all_peptide_ids.empty())
      {
        OPENMS_LOG_WARN << "score:fdr cutoff " << score_fdr << " dropped all PSMs. "
                        << "Output will contain no peptide identifications." << std::endl;
      }
    }

    if (getFlag_("best_per_spectrum_only"))
    {
      IDFilter::keepBestPeptideHits(all_peptide_ids);
    }
  }

  void registerOptionsAndFlags_() override
  {
    static const bool is_required(true);
    static const bool is_advanced_option(true);
    static const bool force_openms_format(true);
        
    registerInputFileList_("in", "<files>", StringList(), "Input file(s)", !is_required);
    setValidFormats_("in", ListUtils::create<String>("mzid,idXML,idparquet"));
    registerInputFileList_("in_decoy", "<files>", StringList(), "Input decoy file(s) in case of separate searches", !is_required);
    setValidFormats_("in_decoy", ListUtils::create<String>("mzid,idXML,idparquet"));
    registerInputFile_("in_osw", "<file>", "", "Input file in OSW format", !is_required);
    setValidFormats_("in_osw", ListUtils::create<String>("OSW"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,mzid,osw,idparquet"));
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
    setValidStrings_("out_type", ListUtils::create<String>("mzid,idXML,osw,idparquet"));
    String enzs = "no_enzyme,elastase,pepsin,proteinasek,thermolysin,chymotrypsin,lys-n,lys-c,arg-c,asp-n,glu-c,trypsin,trypsinp";
    registerStringOption_("enzyme", "<enzyme>", "trypsin", "Type of enzyme: "+enzs , !is_required);
    setValidStrings_("enzyme", ListUtils::create<String>(enzs));
    registerStringOption_("use_subprocess", "<choice>", "false",
        "Run the external 'percolator' binary instead of the in-process OpenMS::Percolator library. "
        "The in-process backend covers the idXML/mzid + PSM-level FDR path; "
        "OSW input, protein-level FDR, and peptide-level FDR still require the subprocess.", false);
    setValidStrings_("use_subprocess", {"true","false"});

    registerInputFile_("percolator_executable", "<executable>",
        // choose the default value according to the platform where it will be executed
        #ifdef OPENMS_WINDOWSPLATFORM
                       "percolator.exe",
        #else
                       "percolator",
        #endif
                       "The Percolator executable. Required only when -use_subprocess=true; "
                       "the in-process backend doesn't need it.",
                       !is_required, !is_advanced_option, {"is_executable"}
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

    //Post-filter parameters
    registerTOPPSubsection_("score", "Post-filter parameters applied to Percolator output");
    registerDoubleOption_("score:fdr", "<value>",
      1.0,
      "FDR cutoff applied to the Percolator q-value before writing output. "
      "PSMs with q-value > cutoff are dropped. 1.0 disables the filter.",
      false, false);
    setMinFloat_("score:fdr", 0.0);
    setMaxFloat_("score:fdr", 1.0);

    registerFlag_("best_per_spectrum_only",
      "After applying score:fdr, retain only the best-scoring PSM per spectrum (default: keep all hits, matching legacy PercolatorAdapter behaviour).",
      false);
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
  
  ExitCodes readInputFiles_(const StringList& in_list, PeptideIdentificationList& all_peptide_ids, vector<ProteinIdentification>& all_protein_ids, bool isDecoy, bool& found_decoys, int& min_charge, int& max_charge)
  {
    for (StringList::const_iterator fit = in_list.begin(); fit != in_list.end(); ++fit)
    {
      String file_idx(distance(in_list.begin(), fit));
      PeptideIdentificationList peptide_ids;
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
      else if (in_type == FileTypes::IDPARQUET)
      {
        FileHandler().loadIdentifications(in, protein_ids, peptide_ids, {FileTypes::IDPARQUET});
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
          else if (hit.isDecoy())
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
    PeptideIdentificationList all_peptide_ids;
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
      
      StringList feature_set = PercolatorInfile::getStandardFeatureSet(min_charge, max_charge);
      if (description_of_correct != 0)
      {
        feature_set.push_back("retentiontime");
        feature_set.push_back("deltamass");
      }

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

      // In-process path: handle the simple idXML/mzid + PSM-level FDR case
      // without spawning the external 'percolator' binary. Falls through to
      // the subprocess path for anything the in-process wrapper can't do
      // (protein-level FDR, peptide-level FDR, init_weights, etc.).
      const bool use_subprocess = (getStringOption_("use_subprocess") == "true");
      const bool in_process_ok = !use_subprocess
                                 && !protein_level_fdrs
                                 && !peptide_level_fdrs
                                 && description_of_correct == 0
                                 && getStringOption_("init_weights").empty();

      if (in_process_ok)
      {
        // Stamp PIN meta values on all hits — this mirrors what the subprocess
        // path does at .pin write time. After this, every hit carries the
        // full PIN feature set (CalcMass, ExpMass, mass, peplen, chargeN,
        // enzN/enzC/enzInt, dm/absdm, score, etc.) plus any search-engine-
        // specific extra_features already present. Training on this set
        // matches the subprocess path's feature vectors row-for-row.
        const auto skipped = PercolatorInfile::stampPinFeaturesOnHits(
          all_peptide_ids, enz_str, min_charge, max_charge);
        OPENMS_LOG_INFO << "Stamped PIN feature meta values; "
                        << skipped.size() << " PSMs skipped (no evidence or unknown TD)."
                        << std::endl;

        // After stamping, feature_set entries SpecId/Label/Peptide/Proteins/ScanNr
        // are either non-numeric (string) or bookkeeping metadata, not training
        // features. ExpMass is a mass value Percolator uses for sort hashing,
        // not a training feature either (it goes into RescoreInput.exp_masses).
        // Filter feature_set to the numeric columns that actually discriminate.
        const std::set<String> pin_metadata_not_feature {
          "SpecId", "Label", "ScanNr", "ExpMass", "Peptide", "Proteins"
        };
        StringList numeric_features;
        for (const String& f : feature_set)
        {
          if (pin_metadata_not_feature.count(f)) continue;
          numeric_features.push_back(f);
        }
        if (all_peptide_ids.empty() || all_peptide_ids.front().getHits().empty()
            || numeric_features.empty())
        {
          writeLogError_("No usable PSMs/features for in-process path; "
                         "use -use_subprocess true.");
          return INCOMPATIBLE_INPUT_DATA;
        }
        OPENMS_LOG_INFO << "Rescoring " << all_peptide_ids.size()
                        << " PSMs in-process with " << numeric_features.size()
                        << " features (skipping external percolator binary)."
                        << std::endl;

        // Build a RescoreInput manually with PIN-compatible fields so the
        // in-process path mirrors the PSM sort order the subprocess path
        // would produce. We do the derivation inline here so rows stay
        // aligned when some hits get skipped (missing target_decoy, missing
        // feature) — fillPINCompatibleFields doesn't know about those skips.
        RescoreInput ri;
        ri.feature_names = numeric_features;
        std::vector<std::pair<size_t, size_t>> hit_locs;

        // After stamping, every kept hit carries ScanNr / ExpMass / CalcMass
        // / Label meta values with PIN-equivalent derivations. Read them back
        // out to build RescoreInput, and skip the same hits the stamp skipped.
        for (size_t i = 0; i < all_peptide_ids.size(); ++i)
        {
          const auto& pid = all_peptide_ids[i];
          for (size_t j = 0; j < pid.getHits().size(); ++j)
          {
            if (skipped.count({i, j})) continue;
            const PeptideHit& hit = pid.getHits()[j];

            // Build feature row from stamped meta values.
            std::vector<double> row;
            row.reserve(numeric_features.size());
            bool ok = true;
            for (const String& f : numeric_features)
            {
              if (!hit.metaValueExists(f)) { ok = false; break; }
              row.push_back(static_cast<double>(hit.getMetaValue(f)));
            }
            if (!ok) continue;

            ri.features.push_back(std::move(row));
            ri.is_decoy.push_back(hit.isDecoy());
            ri.scan_numbers.push_back(
              static_cast<int>(hit.getMetaValue("ScanNr")));
            // Derive stable specFileNr from ScanNr (good enough for sort order
            // since all hits with same scan share the same spec file).
            ri.spec_file_numbers.push_back(0);
            ri.exp_masses.push_back(
              static_cast<double>(hit.getMetaValue("ExpMass")));
            ri.calc_masses.push_back(
              static_cast<double>(hit.getMetaValue("CalcMass")));
            hit_locs.emplace_back(i, j);
          }
        }

        Percolator perc;
        Param pp = perc.getDefaults();
        if (getDoubleOption_("cpos") > 0.0)   pp.setValue("c_pos",    getDoubleOption_("cpos"));
        if (getDoubleOption_("cneg") > 0.0)   pp.setValue("c_neg",    getDoubleOption_("cneg"));
        pp.setValue("test_fdr",       getDoubleOption_("testFDR"));
        pp.setValue("train_fdr",      getDoubleOption_("trainFDR"));
        pp.setValue("num_iterations", getIntOption_("maxiter"));
        pp.setValue("seed",           getIntOption_("seed"));
        pp.setValue("normalizer",     getFlag_("unitnorm") ? "uni" : "stdv");
        pp.setValue("train_best_positive",
                    getFlag_("train_best_positive") ? "true" : "false");
        pp.setValue("post_processing_tdc",
                    getFlag_("post_processing_tdc") ? "true" : "false");
        pp.setValue("nested_xval_bins", getIntOption_("nested_xval_bins"));
        pp.setValue("subset_max_train", getIntOption_("subset_max_train"));
        pp.setValue("initial_direction", getStringOption_("default_direction"));
        pp.setValue("pep_method", "nonparametric");  // match external percolator binary's PEP algorithm
        {
          const int in_process_threads = std::min(std::max(getIntOption_("threads"), 1), 3);
          pp.setValue("num_threads", in_process_threads);  // mirror subprocess --num-threads
        }
        perc.setParameters(pp);
        // Call the low-level API so the PIN-compatible fields on `ri` are
        // actually used (the high-level rescore(peptide_ids, …) ignores them).
        RescoreOutput ro = perc.rescore(ri);

        // Stamp percolator_* meta values back onto each hit from ro.
        for (size_t row = 0; row < ro.scores.size(); ++row)
        {
          const auto [pid_i, hit_i] = hit_locs[row];
          PeptideHit& hit = all_peptide_ids[pid_i].getHits()[hit_i];
          hit.setMetaValue("percolator_score",   ro.scores[row]);
          hit.setMetaValue("percolator_q_value", ro.q_values[row]);
          hit.setMetaValue("percolator_pep",     ro.peps[row]);
        }

        // Transfer percolator_* meta values into the canonical score fields
        // expected by downstream idXML/mzid consumers. Mirror the subprocess
        // path's identifier normalization + meta-value stamps so the idXML
        // writer's strict peptide↔protein identifier cross-check passes.
        const String score_type = getStringOption_("score_type");
        const String run_identifier = all_protein_ids.front().getIdentifier();
        for (auto& pid : all_peptide_ids.getData())
        {
          const String old_score_type = pid.getScoreType();
          pid.setIdentifier(run_identifier);  // align with the (single) IdentificationRun
          if (score_type == "pep")
          {
            pid.setScoreType("Posterior Error Probability");
            pid.setHigherScoreBetter(false);
          }
          else if (score_type == "svm")
          {
            pid.setScoreType("svm");
            pid.setHigherScoreBetter(true);
          }
          else // "q-value" (default)
          {
            pid.setScoreType("q-value");
            pid.setHigherScoreBetter(false);
          }

          for (auto& hit : pid.getHits())
          {
            if (!hit.metaValueExists("percolator_score")) continue;
            const double svm  = hit.getMetaValue("percolator_score");
            const double qval = hit.getMetaValue("percolator_q_value");
            const double pep  = hit.getMetaValue("percolator_pep");

            // Mirror subprocess path's PSI-MS CV meta values
            hit.setMetaValue(old_score_type, hit.getScore());  // preserve original
            hit.setMetaValue("MS:1001492", svm);    // percolator:score
            hit.setMetaValue("MS:1001491", qval);   // percolator:PEP / q-value
            hit.setMetaValue("MS:1001493", pep);    // PEP

            if (score_type == "q-value")      hit.setScore(qval);
            else if (score_type == "pep")     hit.setScore(pep);
            else                              hit.setScore(svm);
          }
        }

        // Clean up the PIN-format meta values that stampPinFeaturesOnHits
        // left on each hit, then stamp PercolatorAdapter provenance metadata
        // (search engine identity + 23 SearchParameter UserParams) so the
        // in-process output matches the historical subprocess output's
        // metadata contract that downstream tools rely on.
        stripPinFeatureMetaValues_(all_peptide_ids, min_charge, max_charge);
        stampPercolatorAdapterMetadata_(all_protein_ids,
          peptide_level_fdrs, protein_level_fdrs,
          /*version_string=*/"3.08-vendored");

        // Apply the same post-filters the subprocess path uses, so the two
        // backends produce equivalent output for -score:fdr / -best_per_spectrum_only.
        applyPostFilters_(all_peptide_ids);

        // Write output and return — skipping the pin / subprocess / pout path.
        FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {out_type});
        return EXECUTION_OK;
      }

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

    std::vector<String> arguments;
    // Check all set parameters and get them into arguments StringList
    {
      if (peptide_level_fdrs)
      {
        arguments.push_back("-r"); arguments.push_back(pout_target_file_peptides);
        arguments.push_back("-B"); arguments.push_back(pout_decoy_file_peptides);
      }
      else
      {
        arguments.push_back("-U");
      }
      arguments.push_back("-m"); arguments.push_back(pout_target_file);
      arguments.push_back("-M"); arguments.push_back(pout_decoy_file);

      if (protein_level_fdrs)
      {
        arguments.push_back("-l"); arguments.push_back(pout_target_file_proteins);
        arguments.push_back("-L"); arguments.push_back(pout_decoy_file_proteins);

        String fasta_file = getStringOption_("fasta");
        if (fasta_file.empty())
        {
          fasta_file = "auto";
        }
        arguments.push_back("-f"); arguments.push_back(fasta_file);

        arguments.push_back("-z"); arguments.push_back(String(enz_str));

        String decoy_pattern = getStringOption_("decoy_pattern");
        if (decoy_pattern != "random") { arguments.push_back("-P"); arguments.push_back(decoy_pattern); }
      }

      int cv_threads = getIntOption_("threads"); // pass-through of OpenMS thread parameter

      if (cv_threads != 3) // default in percolator is 3
      {
        // If a lower than default value (3) is chosen the user needs to enforce it.
        // This ensures that existing workflows (which implicitly used 3 threads) don't slow down
        // if e.g. the OpenMS version and this adapter is updated.
        if (cv_threads > 3 || getFlag_("force"))
        {
          arguments.push_back("--num-threads"); arguments.push_back(String(cv_threads));
        }
      }

      double cpos = getDoubleOption_("cpos");
      double cneg = getDoubleOption_("cneg");
      if (cpos != 0.0)
      {
        arguments.push_back("-p"); arguments.push_back(String(cpos));
      }
      if (cneg != 0.0)
      {
        arguments.push_back("-n"); arguments.push_back(String(cneg));
      }
      double train_FDR = getDoubleOption_("trainFDR");
      double test_FDR = getDoubleOption_("testFDR");
      if (train_FDR != 0.01)
      {
        arguments.push_back("-F"); arguments.push_back(String(train_FDR));
      }
      if (test_FDR != 0.01)
      {
        arguments.push_back("-t"); arguments.push_back(String(test_FDR));
      }
      Int max_iter = getIntOption_("maxiter");
      if (max_iter != 10)
      {
        arguments.push_back("-i"); arguments.push_back(String(max_iter));
      }
      Int subset_max_train = getIntOption_("subset_max_train");
      if (subset_max_train > 0)
      {
        arguments.push_back("-N"); arguments.push_back(String(subset_max_train));
      }
      if (getFlag_("quick_validation"))
      {
        arguments.push_back("-x");
      }
      if (getFlag_("post_processing_tdc"))
      {
        arguments.push_back("-Y");
      }
      if (getFlag_("train_best_positive"))
      {
        arguments.push_back("--train-best-positive");
      }
      if (getFlag_("static"))
      {
        arguments.push_back("--static");
      }
      Int nested_xval_bins = getIntOption_("nested_xval_bins");
      if (nested_xval_bins > 1)
      {
        arguments.push_back("--nested-xval-bins"); arguments.push_back(String(nested_xval_bins));
      }
      String weights_file = getStringOption_("weights");
      String init_weights_file = getStringOption_("init_weights");
      String default_search_direction = getStringOption_("default_direction");
      if (!weights_file.empty())
      {
        arguments.push_back("-w"); arguments.push_back(weights_file);
      }
      if (!init_weights_file.empty())
      {
        arguments.push_back("-W"); arguments.push_back(init_weights_file);
      }
      if (!default_search_direction.empty())
      {
        arguments.push_back("-V"); arguments.push_back(default_search_direction);
      }
      Int verbose_level = getIntOption_("verbose");
      if (verbose_level != 2)
      {
        arguments.push_back("-v"); arguments.push_back(String(verbose_level));
      }
      if (getFlag_("unitnorm"))
      {
        arguments.push_back("-u");
      }
      if (getFlag_("test_each_iteration"))
      {
        arguments.push_back("-R");
      }
      if (getFlag_("override"))
      {
        arguments.push_back("-O");
      }
      Int seed = getIntOption_("seed");
      if (seed != 1)
      {
        arguments.push_back("-S"); arguments.push_back(String(seed));
      }
      if (getFlag_("klammer"))
      {
        arguments.push_back("-K");
      }
      if (description_of_correct != 0)
      {
        arguments.push_back("-D"); arguments.push_back(String(description_of_correct));
      }
      arguments.push_back(pin_file);
    }
    writeLogInfo_("Prepared percolator input.");

    //-------------------------------------------------------------
    // run percolator
    //-------------------------------------------------------------
    // Percolator execution with the executable and the arguments StringList
    TOPPBase::ExitCodes exit_code = runExternalProcess_(percolator_executable, arguments);
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
        File::copy(pout_target_file_peptides, pout_target);
      }
      if (!pout_decoy.empty())
      {
        File::copy(pout_decoy_file_peptides, pout_decoy);
      }
    }
    else
    {
      readPoutAsMap_(pout_target_file, pep_map);
      readPoutAsMap_(pout_decoy_file, pep_map);

      // copy file in tmp folder to output
      if (!pout_target.empty())
      {
        File::copy(pout_target_file, pout_target);
      }
      if (!pout_decoy.empty())
      {
        File::copy(pout_decoy_file, pout_decoy);
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
        File::copy(pout_target_file_proteins, pout_target_proteins);
      }
      if (!pout_decoy_proteins.empty())
      {
        File::copy(pout_decoy_file_proteins, pout_decoy_proteins);
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
      stampPercolatorAdapterMetadata_(all_protein_ids,
        peptide_level_fdrs, protein_level_fdrs,
        /*version_string=*/"3.07",  // TODO: read from percolator binary --version
        protein_level_fdrs ? &protein_map : nullptr);

      applyPostFilters_(all_peptide_ids);

      // Storing the PeptideHits with calculated q-value, pep and svm score
      FileHandler().storeIdentifications(out, all_protein_ids, all_peptide_ids, {FileTypes::IDXML, FileTypes::MZIDENTML, FileTypes::IDPARQUET});
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
