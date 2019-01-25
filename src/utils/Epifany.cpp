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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <vector>

using namespace OpenMS;
using namespace std;

class TOPPEpifany :
public TOPPBase
{
public:
  TOPPEpifany() :
  TOPPBase("EPIFANY", "Runs a Bayesian protein inference.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file>", StringList(), "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML,consensusXML,fasta"));
    registerInputFile_("exp_design", "<file>", "", "Input: experimental design", false);
    setValidFormats_("exp_design", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output: identification results with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML,consensusXML"));
    registerFlag_("separate_runs", "Process multiple protein identification runs in the input separately,"
                                   " don't merge them. Merging results in loss of descriptive information"
                                   " of the single protein identification runs.", false);
    registerStringOption_("protein_fdr",
                          "<option>",
                          "false",
                          "Additionally calculate the target-decoy FDR on protein-level based on the posteriors", false, false);
    setValidStrings_("protein_fdr", {"true","false"});

    registerStringOption_("greedy_group_resolution",
                       "<option>",
                       "none",
                       "Post-process inference output with greedy resolution of shared peptides based on the parent protein probabilities. Also adds the resolved ambiguity groups to output.", false, true);
    setValidStrings_("greedy_group_resolution", {"none","remove_associations_only","remove_proteins_wo_evidence"});
    registerDoubleOption_("evidence_cutoff",
                          "<option>",
                          0.001,
                          "After greedy resolution, before filtering proteins evidence-based, remove very low-scoring PSMs under this posterior probability cutoff.", false, true);
    addEmptyLine_();

    registerSubsection_("algorithm", "Parameters for the Algorithm section");
  }

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    return BayesianProteinInferenceAlgorithm().getParameters();
  }

  ExitCodes main_(int, const char**) override
  {
    // get parameters specific for the feature finder
    Param epifany_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to Epifany", epifany_param, 3);

    StringList files = getStringList_("in");
    IdXMLFile idXMLf;
    IDMergerAlgorithm merger{};
    StopWatch sw;
    sw.start();
    vector<ProteinIdentification> mergedprots{1};
    vector<PeptideIdentification> mergedpeps;

    if (files.size() > 1)
    {
      for (String& file : files)
      {
        vector<ProteinIdentification> prots;
        vector<PeptideIdentification> peps;
        idXMLf.load(file,prots,peps);
        merger.insertRun(prots,peps);
      }

      merger.returnResultsAndClear(mergedprots[0], mergedpeps);
    }
    else
    {
      idXMLf.load(files[0],mergedprots,mergedpeps);
    }

    IDFilter::filterBestPerPeptide(mergedpeps, true, true, 1);

    IDFilter::filterEmptyPeptideIDs(mergedpeps);

    //TODO get from data or set as parameter
    double min_nonnull_obs_probability = 0.0001;
    double pre_filter = 0.05;
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
          if (newScore <= min_nonnull_obs_probability)
          {
            newScore = min_nonnull_obs_probability;
          }
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
        for (auto& pep_hit : pep_id.getHits())
        {
          double newScore = pep_hit.getScore();
          if (newScore <= min_nonnull_obs_probability)
          {
            newScore = min_nonnull_obs_probability;
          }
          pep_hit.setScore(newScore);
        }
      }
    }
    if (pre_filter > 0.0)
    {
      IDFilter::filterHitsByScore(mergedpeps, pre_filter);
      IDFilter::filterEmptyPeptideIDs(mergedpeps);
    }

    // Currently this is needed because otherwise there might be proteins with a previous score
    // that get evaluated during FDR without a new posterior being set.
    // Alternative would be to reset scores but this does not work well if you wanna work with i.e. user priors
    IDFilter::removeUnreferencedProteins(mergedprots, mergedpeps);

    LOG_INFO << "Loading and merging took " << sw.toString() << std::endl;
    sw.reset();

    BayesianProteinInferenceAlgorithm bpi1;
    bpi1.setParameters(epifany_param);
    bpi1.inferPosteriorProbabilities(mergedprots, mergedpeps);
    LOG_INFO << "Inference total took " << sw.toString() << std::endl;
    sw.stop();

    bool greedy_group_resolution = getStringOption_("greedy_group_resolution") != "none";
    bool remove_prots_wo_evidence = getStringOption_("greedy_group_resolution") == "remove_proteins_wo_evidence";

    if (greedy_group_resolution)
    {
      LOG_INFO << "Postprocessing: Removing associations from spectrum to all but the best protein group..." << std::endl;
      //TODO add group resolution to the IDBoostGraph class so we do not
      // unnecessarily build a second (old) data structure
      mergedprots[0].fillIndistinguishableGroupsWithSingletons();
      PeptideProteinResolution ppr;
      ppr.buildGraph(mergedprots[0], mergedpeps);
      ppr.resolveGraph(mergedprots[0], mergedpeps);

      //PeptideProteinResolution::resolve(mergedprots[0], mergedpeps, true, false);
    }
    if (remove_prots_wo_evidence)
    {
      LOG_INFO << "Postprocessing: Removing proteins without associated evidence..." << std::endl;
      //TODO filtering as parameter? Do it before (remember that before the inference PSMs often hold PEPs not PPs)
      // Should be bigger than min_nonnull_obs_probability!!
      double cutoff = getDoubleOption_("evidence_cutoff");
      if (cutoff <= min_nonnull_obs_probability)
      {
        LOG_WARN << "Post-greedy, pre-filter cutoff smaller than minimum PSM probability. Maybe not useful." << std::endl;
      }
      IDFilter::filterHitsByScore(mergedpeps, cutoff);
      IDFilter::removeUnreferencedProteins(mergedprots, mergedpeps);
      IDFilter::updateProteinGroups(mergedprots[0].getIndistinguishableProteins(), mergedprots[0].getHits());
      IDFilter::updateProteinGroups(mergedprots[0].getProteinGroups(), mergedprots[0].getHits());
    }

    bool calc_protFDR = getStringOption_("protein_fdr") == "true";
    if (calc_protFDR)
    {
      LOG_INFO << "Calculating target-decoy q-values..." << std::endl;
      FalseDiscoveryRate fdr;
      fdr.applyBasic(mergedprots[0], true);
    }

    LOG_INFO << "Writing inference run as first ProteinIDRun with " <<
             mergedprots[0].getHits().size() << " proteins in " <<
             mergedprots[0].getIndistinguishableProteins().size() <<
             " indist. groups." << std::endl;
    idXMLf.store(getStringOption_("out"),mergedprots,mergedpeps);
    return ExitCodes::EXECUTION_OK;

    /* That was test code for the old merger
    ProteinIdentification pi1{};
    ProteinIdentification pi2{};
    ProteinIdentification pi3{};
    ProteinIdentification pi4{};
    pi1.setPrimaryMSRunPath({"a.mzml", "b.mzml", "c.mzml"});
    pi2.setPrimaryMSRunPath({"d.mzml", "c.mzml", "e.mzml"});
    pi3.setPrimaryMSRunPath({"g.mzml", "f.mzml", "e.mzml"});
    pi4.setPrimaryMSRunPath({"h.mzml", "e.mzml"});
    vector<ProteinIdentification> pis{pi1,pi2,pi3,pi4};
    vector<PeptideIdentification> peps1{};
    IDMergerAlgorithm merger1;
    vector<vector<PeptideIdentification>> newpeps1{};
    merger1.mergeIDRunsAndSplitPeptides(pis, peps1, {{0,0},{1,0},{2,1},{3,1}}, newpeps1);

    return ExitCodes::EXECUTION_OK;


    //TODO remove until ExitCode and allow a (merged) ConsensusXML as input
    vector<PeptideIdentification> peps;
    vector<ProteinIdentification> prots;
    ConsensusMap cmap;
    ConsensusXMLFile cXML;
    cXML.load(getStringOption_("in"), cmap);
    IDMergerAlgorithm merger;
    ExperimentalDesignFile expFile;
    ExperimentalDesign expDesign = expFile.load(getStringOption_("exp_design"), false);
    merger.mergeProteinsAcrossFractionsAndReplicates(cmap, expDesign);
    cXML.store(getStringOption_("out"), cmap);

    return ExitCodes::EXECUTION_OK;
     */

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

    vector<PeptideIdentification> peps;
    vector<ProteinIdentification> prots;
    IdXMLFile idXML;
    idXML.load(getStringOption_("in"), prots, peps);
    //TODO filter unmatched proteins and peptides before!
    //TODO check t+d annotations

    BayesianProteinInferenceAlgorithm bpi;
    bpi.inferPosteriorProbabilities(prots, peps);

    idXML.store(getStringOption_("out"),prots, peps);

    return ExitCodes::EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPEpifany tool;

  return tool.main(argc, argv);
}
