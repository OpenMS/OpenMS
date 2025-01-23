// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/ConsensusMapMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <vector>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_Epifany Epifany

@brief EPIFANY - Efficient protein inference for any peptide-protein network is a Bayesian
protein inference engine. It uses PSM (posterior) probabilities from Percolator, OpenMS' IDPosteriorErrorProbability
or similar tools to calculate posterior probabilities for proteins and protein groups.

@experimental This tool is work in progress and usage and input requirements might change.

<center>
  <table>
      <tr>
          <th ALIGN = "center"> pot. predecessor tools </td>
          <td VALIGN="middle" ROWSPAN=2> &rarr; Epifany &rarr;</td>
          <th ALIGN = "center"> pot. successor tools </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PercolatorAdapter </td>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
      </tr>
      <tr>
          <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDPosteriorErrorProbability </td>
      </tr>
  </table>
</center>
<p>It is a protein inference engine based on a Bayesian network. Currently the same model like
Fido is used with the main parameters alpha (pep_emission), beta (pep_spurious_emission) and gamma (prot_prior).
If not specified,
these parameters are trained based on their classification performance and calibration via a grid search
by simply running with several possible combinations and evaluating. Unless you see very extreme output
probabilities (e.g. many close to 1.0) or you know good parameters (e.g. from an earlier run),
grid search is recommended although slower. The tool will merge multiple idXML files (union of proteins
and concatenation of PSMs) when given more than one. It assumes one search engine run per input file but
might work on more. Proteins need to be indexed by OpenMS's PeptideIndexer but this is usually done before
Percolator/IDPEP since target/decoy associations are needed there already. Make sure that the input PSM
probabilities are not too extreme already (garbage in - garbage out). After merging the input probabilities
are preprocessed with a low posterior probability cutoff to neglect very unreliable matches. Then
the probabilities are aggregated with the maximum per peptide and the graph is built and split into
connected components. When compiled with the OpenMP
flag (default enabled in the release binaries) the tool is multi-threaded which can
be activated at runtime by the threads parameter. Note that peak memory requirements
may rise significantly when processing multiple components of the graph at the same time.
</p>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_Epifany.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_Epifany.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class Epifany :
public TOPPBase
{
public:
  Epifany() :
  TOPPBase("Epifany", "Runs a Bayesian protein inference.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    //TODO support separate runs
    registerInputFileList_("in", "<file>", StringList(), "Input: identification results");
    setValidFormats_("in", {"idXML","consensusXML"});
    registerInputFile_("exp_design", "<file>", "", "(Currently unused) Input: experimental design", false);
    setValidFormats_("exp_design", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output: identification results with scored/grouped proteins");
    setValidFormats_("out", {"idXML","consensusXML"});
    registerStringOption_("out_type", "<file>", "", "Output type: auto detected by file extension but can be overwritten here.", false);
    setValidStrings_("out_type", {"idXML","consensusXML"});

    registerStringOption_("protein_fdr",
                          "<option>",
                          "false",
                          "Additionally calculate the target-decoy FDR on protein-level based on the posteriors", false, false);
    setValidStrings_("protein_fdr", {"true","false"});

    registerStringOption_("conservative_fdr",
                          "<option>",
                          "true",
                          "Use (D+1)/(T) instead of (D+1)/(T+D) for reporting protein FDRs.", false, true);
    setValidStrings_("conservative_fdr", {"true","false"});

    registerStringOption_("picked_fdr",
                          "<option>",
                          "true",
                          "Use picked protein FDRs.", false, true);
    setValidStrings_("picked_fdr", {"true","false"});
    registerStringOption_("picked_decoy_string",
                          "<decoy_string>",
                          "",
                          "If using picked protein FDRs, which decoy string was used? Leave blank for auto-detection.", false, true);
    registerStringOption_("picked_decoy_prefix",
                          "<option>",
                          "prefix",
                          "If using picked protein FDRs, was the decoy string a prefix or suffix? Ignored during auto-detection.", false, true);
    setValidStrings_("picked_decoy_prefix", {"prefix","suffix"});

    registerStringOption_("greedy_group_resolution",
                       "<option>",
                       "none",
                       "Post-process inference output with greedy resolution of shared peptides based on the parent protein probabilities."
                       " Also adds the resolved ambiguity groups to output.", false, false);
    setValidStrings_("greedy_group_resolution", {"none","remove_associations_only","remove_proteins_wo_evidence"});

    registerDoubleOption_("min_psms_extreme_probability",
        "<float>",
        0.0,
        "Set PSMs with probability lower than this to this minimum probability.", false, true);
    registerDoubleOption_("max_psms_extreme_probability",
        "<float>",
        1.0,
        "Set PSMs with probability higher than this to this maximum probability.", false, false);

    addEmptyLine_();

    registerSubsection_("algorithm", "Parameters for the Algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return BayesianProteinInferenceAlgorithm().getParameters();
  }

  pair<double,double> checkExtremePSMScores_(vector<PeptideIdentification>& mergedpeps)
  {
    double minscore = 2.;
    double maxscore = -1.;
    //convert all scores to PPs
    for (auto& pep_id : mergedpeps)
    {
      for (auto& pep_hit : pep_id.getHits())
      {
        double newScore = pep_hit.getScore();
        if (newScore > 0)
        {
          minscore = std::min(minscore, newScore);
        }
        if (newScore < 1.)
        {
          maxscore = std::max(maxscore, newScore);
        }
      }
    }
    return {minscore, maxscore};
  }

  void convertPSMScores_(vector<PeptideIdentification>& mergedpeps)
  {
    //convert all scores to PPs
    for (auto& pep_id : mergedpeps)
    {
      String score_l = pep_id.getScoreType();
      score_l = score_l.toLower();
      if (score_l == "pep" || score_l == "posterior error probability")
      {
        for (auto& pep_hit : pep_id.getHits())
        {
          double newScore = 1. - pep_hit.getScore();
          pep_hit.setScore(newScore);
        }
        pep_id.setScoreType("Posterior Probability");
        pep_id.setHigherScoreBetter(true);
      }
      else
      {
        if (score_l != "Posterior Probability")
        {
          throw OpenMS::Exception::InvalidParameter(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Epifany needs Posterior (Error) Probabilities in the Peptide Hits. Use Percolator with PEP score"
              "or run IDPosteriorErrorProbability first.");
        }
      }
    }
  }

  void removeExtremeValues_(vector<PeptideIdentification>& mergedpeps, double minscore, double maxscore)
  {
    //convert all scores to PPs
    for (auto& pep_id : mergedpeps)
    {
      for (auto& pep_hit : pep_id.getHits())
      {
        double score = pep_hit.getScore();
        pep_hit.setScore(std::min(std::max(score,minscore),maxscore));
      }
    }
  }

  ExitCodes main_(int, const char**) override
  {
    // get parameters specific for the algorithm underneath
    Param epifany_param = getParam_().copy("algorithm:", true);
    bool greedy_group_resolution = getStringOption_("greedy_group_resolution") != "none";
    bool remove_prots_wo_evidence = getStringOption_("greedy_group_resolution") == "remove_proteins_wo_evidence";

    //writeDebug_("Parameters passed to Epifany", epifany_param, 3);
    StringList files = getStringList_("in");
    if (files.empty())
    {
      OPENMS_LOG_ERROR << "No files given.\n";
    }

    FileTypes::Type in_type = FileHandler::getType(files[0]);
    String exp_des = getStringOption_("exp_design");

    StopWatch sw;
    sw.start();

    String out_file = getStringOption_("out");
    String out_type = getStringOption_("out_type");

    if (!files.empty() && (in_type == FileTypes::CONSENSUSXML))
    {
      if (FileHandler::getTypeByFileName(out_file) != FileTypes::CONSENSUSXML &&
          FileTypes::nameToType(out_type) != FileTypes::CONSENSUSXML)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Running on consensusXML requires output as consensusXML. Please change the "
                                  "output type.\n";
      }
      OPENMS_LOG_INFO << "Loading input..." << std::endl;

      if (files.size() > 1)
      {
        OPENMS_LOG_FATAL_ERROR << "Error: Multiple inputs only supported for idXML\n";
      }
      ConsensusMapMergerAlgorithm cmerge;
      ConsensusMap cmap;
      FileHandler().loadConsensusFeatures(files[0], cmap, {FileTypes::CONSENSUSXML});
      std::optional<const ExperimentalDesign> edopt = maybeGetExpDesign_(exp_des);
      if (!exp_des.empty())
      {
        cmerge.mergeProteinsAcrossFractionsAndReplicates(cmap, edopt.value());
      }
      else
      {
        cmerge.mergeAllIDRuns(cmap);
      }

      OPENMS_LOG_INFO << "Loading took " << sw.toString() << std::endl;
      sw.reset();

      BayesianProteinInferenceAlgorithm bpi1(getIntOption_("debug"));
      bpi1.setParameters(epifany_param);
      bpi1.inferPosteriorProbabilities(cmap, greedy_group_resolution, edopt);
      OPENMS_LOG_INFO << "Inference total took " << sw.toString() << std::endl;
      sw.stop();

      if (remove_prots_wo_evidence)
      {
        OPENMS_LOG_INFO << "Postprocessing: Removing proteins without associated evidence..." << std::endl;
        IDFilter::removeUnreferencedProteins(cmap, true);
        for (auto& run : cmap.getProteinIdentifications())
        {
          IDFilter::updateProteinGroups(run.getIndistinguishableProteins(), run.getHits());
        }
      }

      for (auto& run : cmap.getProteinIdentifications())
      {
        std::sort(run.getHits().begin(), run.getHits().end(),
                  [](const ProteinHit& f,
                     const ProteinHit& g)
                  {return f.getAccession() < g.getAccession();});
        //sort for output because they might have been added in a different order
        std::sort(
            run.getIndistinguishableProteins().begin(),
            run.getIndistinguishableProteins().end(),
            [](const ProteinIdentification::ProteinGroup& f,
               const ProteinIdentification::ProteinGroup& g)
            {return f.accessions < g.accessions;});
      }

      bool calc_protFDR = getStringOption_("protein_fdr") == "true";
      if (calc_protFDR)
      {
        OPENMS_LOG_INFO << "Calculating target-decoy q-values..." << std::endl;
        FalseDiscoveryRate fdr;
        Param fdrparam = fdr.getParameters();
        fdrparam.setValue("conservative", getStringOption_("conservative_fdr"));
        fdrparam.setValue("add_decoy_proteins","true");
        fdr.setParameters(fdrparam);
        for (auto& run : cmap.getProteinIdentifications())
        {
          if (getStringOption_("picked_fdr") == "true")
          {
            fdr.applyPickedProteinFDR(run, getStringOption_("picked_decoy_string"), getStringOption_("picked_decoy_prefix") == "prefix");
          }
          else
          {
            fdr.applyBasic(run, true);
          }
        }
      }

      FileHandler().storeConsensusFeatures(out_file, cmap, {FileTypes::CONSENSUSXML});
    }
    else // ----------------------------   IdXML   -------------------------------------
    {
      IDMergerAlgorithm merger{};
      OPENMS_LOG_INFO << "Loading input..." << std::endl;
      vector<ProteinIdentification> mergedprots{1};
      vector<PeptideIdentification> mergedpeps;
      if (files.size() > 1)
      {
        for (String& file : files)
        {
          vector<ProteinIdentification> prots;
          vector<PeptideIdentification> peps;
          FileHandler().loadIdentifications(file, prots, peps, {FileTypes::IDXML});
          //TODO merger does not support groups yet, so clear them here right away.
          // Not so easy to implement at first sight. Merge groups whenever one protein overlaps?
          prots[0].getIndistinguishableProteins().clear();
          prots[0].getProteinGroups().clear();
          merger.insertRuns(prots, peps);
        }
        merger.returnResultsAndClear(mergedprots[0], mergedpeps);
      }
      else
      {
        FileHandler().loadIdentifications(files[0], mergedprots, mergedpeps, {FileTypes::IDXML});
        //TODO For now we delete because we want to add new groups here.
        // Think about:
        // 1) keeping the groups and allow them to be used as a prior grouping (e.g. gene based)
        // 2) keeping the groups and store them in a separate group object to output both and compare.
        mergedprots[0].getIndistinguishableProteins().clear();
        mergedprots[0].getProteinGroups().clear();
      }

      // Currently this is needed because otherwise there might be proteins with a previous score
      // that get evaluated during FDR without a new posterior being set. (since components of size 1 are skipped)
      // Alternative would be to reset scores but this does not work well if you wanna work with i.e. user priors
      // However, this is done additionally in the Inference class after filtering, so maybe not necessary.

      IDFilter::removeUnreferencedProteins(mergedprots, mergedpeps);

      OPENMS_LOG_INFO << "Loading took " << sw.toString() << std::endl;
      sw.reset();

      //Check if score types are valid.
      try
      {
        IDScoreSwitcherAlgorithm switcher;
        Size c = 0;
        switcher.switchToGeneralScoreType(mergedpeps, IDScoreSwitcherAlgorithm::ScoreType::PEP, c);
      }
      catch(Exception::MissingInformation&)
      {
        OPENMS_LOG_FATAL_ERROR <<
              "Epifany expects a Posterior Error Probability score in all Peptide IDs." << endl;
             return ExitCodes::INCOMPATIBLE_INPUT_DATA;
      }

      BayesianProteinInferenceAlgorithm bpi1(getIntOption_("debug"));
      bpi1.setParameters(epifany_param);
      bpi1.inferPosteriorProbabilities(mergedprots, mergedpeps, greedy_group_resolution);
      OPENMS_LOG_INFO << "Inference total took " << sw.toString() << std::endl;
      sw.stop();

      if (remove_prots_wo_evidence)
      {
        OPENMS_LOG_INFO << "Postprocessing: Removing proteins without associated evidence..." << std::endl;
        IDFilter::removeUnreferencedProteins(mergedprots, mergedpeps);
        IDFilter::updateProteinGroups(mergedprots[0].getIndistinguishableProteins(), mergedprots[0].getHits());
        IDFilter::updateProteinGroups(mergedprots[0].getProteinGroups(), mergedprots[0].getHits());
      }

      bool calc_protFDR = getStringOption_("protein_fdr") == "true";
      if (calc_protFDR)
      {
        OPENMS_LOG_INFO << "Calculating target-decoy q-values..." << std::endl;
        FalseDiscoveryRate fdr;
        Param fdrparam = fdr.getParameters();
        fdrparam.setValue("conservative", getStringOption_("conservative_fdr"));
        fdrparam.setValue("add_decoy_proteins","true");
        fdr.setParameters(fdrparam);
        if (getStringOption_("picked_fdr") == "true")
        {
          fdr.applyPickedProteinFDR(mergedprots[0], getStringOption_("picked_decoy_string"), getStringOption_("picked_decoy_prefix") == "prefix");
        }
        else
        {
          fdr.applyBasic(mergedprots[0], true);
        }
      }

      OPENMS_LOG_INFO << "Writing inference run as first ProteinIDRun with " <<
                      mergedprots[0].getHits().size() << " proteins in " <<
                      mergedprots[0].getIndistinguishableProteins().size() <<
                      " indist. groups." << std::endl;

      //sort for output because they might have been added in a different order
      std::sort(
          mergedprots[0].getIndistinguishableProteins().begin(),
          mergedprots[0].getIndistinguishableProteins().end(),
          [](const ProteinIdentification::ProteinGroup& f,
             const ProteinIdentification::ProteinGroup& g)
          {return f.accessions < g.accessions;});

      FileHandler().storeIdentifications(out_file, mergedprots, mergedpeps, {FileTypes::IDXML});
    }
    return ExitCodes::EXECUTION_OK;
  }

  // Some thoughts about how to leverage info from different runs.
  //Fractions: Always merge (not much to leverage, maybe agreement at borders)
  // - Think about only allowing one/the best PSM per peptidoform across fractions
  //Replicates: Use matching ID and quant, also a always merge
  //Samples: In theory they could yield different proteins/pep-protein-associations
  // 3 options:
  // - don't merge: -> don't leverage peptide quant profiles (or use them repeatedly -> same as second opt.?)
  // - merge and assume same proteins: -> You can use the current graph and weigh the associations
  //   based on deviation in profiles
  // - merge and don't assume same proteins: -> We need an extended graph, that has multiple versions
  //   of the proteins for every sample

  static std::optional<const ExperimentalDesign> maybeGetExpDesign_(const String& filename)
  {
    if (filename.empty()) return std::nullopt;
    return std::optional<const ExperimentalDesign>(ExperimentalDesignFile::load(filename, false));
  }
};



int main(int argc, const char** argv)
{
  Epifany tool;

  return tool.main(argc, argv);
}

/// @endcond
