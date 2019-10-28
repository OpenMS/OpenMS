// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
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


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_Epifany Epifany

  @brief EPIFANY - Efficient protein inference for any peptide-protein network is a Bayesian
  protein inference engine. It uses PSM (posterior) probabilities from Percolator, OpenMS' IDPosteriorErrorProbability
  or similar tools to calculate posterior probabilities for proteins and protein groups.

  @experimental This tool is work in progress and usage and input requirements might change.

  <center>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ Epifany \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
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
  @verbinclude UTILS_Epifany.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_Epifany.html

*/

class Epifany :
public TOPPBase
{
public:
  Epifany() :
  TOPPBase("Epifany", "Runs a Bayesian protein inference.", false)
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    //TODO support separate runs
    registerInputFileList_("in", "<file>", StringList(), "Input: identification results");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    //registerInputFile_("exp_design", "<file>", "", "(Currently unused) Input: experimental design", false);
    //setValidFormats_("exp_design", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Output: identification results with scored/grouped proteins");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

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
    for (auto &pep_id : mergedpeps)
    {
      for (auto &pep_hit : pep_id.getHits())
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
    for (auto &pep_id : mergedpeps)
    {
      String score_l = pep_id.getScoreType();
      score_l = score_l.toLower();
      if (score_l == "pep" || score_l == "posterior error probability")
      {
        for (auto &pep_hit : pep_id.getHits())
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
    for (auto &pep_id : mergedpeps)
    {
      for (auto &pep_hit : pep_id.getHits())
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
    writeDebug_("Parameters passed to Epifany", epifany_param, 3);

    StringList files = getStringList_("in");
    IdXMLFile idXMLf;
    IDMergerAlgorithm merger{};
    StopWatch sw;
    sw.start();
    vector<ProteinIdentification> mergedprots{1};
    vector<PeptideIdentification> mergedpeps;

    OPENMS_LOG_INFO << "Loading input..." << std::endl;
    if (files.size() > 1)
    {
      for (String& file : files)
      {
        //TODO this only works for idXML
        vector<ProteinIdentification> prots;
        vector<PeptideIdentification> peps;
        idXMLf.load(file, prots, peps);
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
      idXMLf.load(files[0], mergedprots, mergedpeps);
      //TODO For now we delete because we want to add new groups here.
      // Think about:
      // 1) keeping the groups and allow them to be used as a prior grouping (e.g. gene based)
      // 2) keeping the groups and store them in a separate group object to output both and compare.
      mergedprots[0].getIndistinguishableProteins().clear();
      mergedprots[0].getProteinGroups().clear();
    }

    // Currently this is needed because otherwise there might be proteins with a previous score
    // that get evaluated during FDR without a new posterior being set.
    // Alternative would be to reset scores but this does not work well if you wanna work with i.e. user priors
    IDFilter::removeUnreferencedProteins(mergedprots, mergedpeps);

    OPENMS_LOG_INFO << "Loading took " << sw.toString() << std::endl;
    sw.reset();

    BayesianProteinInferenceAlgorithm bpi1(getIntOption_("debug"));
    bpi1.setParameters(epifany_param);
    bpi1.inferPosteriorProbabilities(mergedprots, mergedpeps);
    OPENMS_LOG_INFO << "Inference total took " << sw.toString() << std::endl;
    sw.stop();

    // Let's always add all the proteins to the protein group section, easier in postprocessing.
    // PeptideProteinResolution needs it anyway.
    //TODO check if still needed after adding the addSingleton option to the IDGraph function
    mergedprots[0].fillIndistinguishableGroupsWithSingletons();

    bool greedy_group_resolution = getStringOption_("greedy_group_resolution") != "none";
    bool remove_prots_wo_evidence = getStringOption_("greedy_group_resolution") == "remove_proteins_wo_evidence";

    if (greedy_group_resolution)
    {
      OPENMS_LOG_INFO << "Postprocessing: Removing associations from spectrum via best PSM to all but the best protein group..." << std::endl;
      //TODO add group resolution to the IDBoostGraph class so we do not
      // unnecessarily build a second (old) data structure

      PeptideProteinResolution ppr;
      ppr.buildGraph(mergedprots[0], mergedpeps);
      ppr.resolveGraph(mergedprots[0], mergedpeps);

      //PeptideProteinResolution::resolve(mergedprots[0], mergedpeps, true, false);
    }
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
      fdr.applyBasic(mergedprots[0], true);
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


    idXMLf.store(getStringOption_("out"),mergedprots,mergedpeps);
    return ExitCodes::EXECUTION_OK;


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
  }

};

int main(int argc, const char** argv)
{
  Epifany tool;

  return tool.main(argc, argv);
}
