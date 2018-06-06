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
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <set>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

using namespace std;

namespace OpenMS
{
  /// A functor that specifies what to do on a connected component (IDBoostGraph::FilteredGraph)
  class BayesianProteinInferenceAlgorithm::FilteredGraphInferenceFunctor :
      public std::function<void(IDBoostGraph::FilteredGraph&)>
  {
  public:
    const Param& param_;

    explicit FilteredGraphInferenceFunctor(const Param& param):
        param_(param)
    {}

    void operator() (IDBoostGraph::FilteredGraph& fg) {
      //------------------ Now actual inference ------------------- //
      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types or skip nodes without neighbors inside the if instead
      //TODO do quick bruteforce calculation if the cc is really small
      if (boost::num_vertices(fg) >= 3)
      {
        MessagePasserFactory<unsigned long> mpf (param_.getValue("model_parameters:pep_emission"),
                                                 param_.getValue("model_parameters:pep_spurious_emission"),
                                                 param_.getValue("model_parameters:prot_prior"),
                                                 1.0); // the p used for marginalization: 1 = sum product, inf = max product
        BetheInferenceGraphBuilder<unsigned long> bigb;

        boost::filtered_graph<IDBoostGraph::Graph, boost::function<bool(IDBoostGraph::edge_t)>, boost::function<bool(IDBoostGraph::vertex_t)> >::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end (usually at least proteins)
        // Maybe later peptides (e.g. for an iterative procedure)
        vector<vector<unsigned long>> posteriorVars;

        // direct neighbors are proteins on the "left" side and peptides on the "right" side
        // TODO use directed graph? But finding "non-strong" connected components in a directed graph is not
        // out of the box supported by boost
        std::vector<IDBoostGraph::vertex_t> in{};
        //std::vector<IDBoostGraph::vertex_t> out{};

        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::FilteredGraph::adjacency_iterator nbIt, nbIt_end;
          boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

          in.clear();
          //out.clear();

          for (; nbIt != nbIt_end; ++nbIt)
          {
            if (fg[*nbIt].which() < fg[*ui].which())
            {
              in.push_back(*nbIt);
            }
            /*else
            {
              out.push_back(*nbIt);
            }*/
          }

          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.

          if (fg[*ui].which() == 3) // pep
          {
            bigb.insert_dependency(mpf.createSumEvidenceFactor(boost::get<PeptideHit*>(fg[*ui])->getPeptideEvidences().size(), in[0], *ui));
            bigb.insert_dependency(mpf.createPeptideEvidenceFactor(*ui, boost::get<PeptideHit*>(fg[*ui])->getScore()));
          }
          else if (fg[*ui].which() == 2) // pep group
          {
            bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(in, *ui));
          }
          else if (fg[*ui].which() == 1) // prot group
          {
            bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(in, *ui));
          }
          else if (fg[*ui].which() == 0) // prot
          {
            //TODO allow an already present prior probability here
            bigb.insert_dependency(mpf.createProteinFactor(*ui));
            posteriorVars.push_back({*ui});
          }

        }

        // create factor graph for Bayesian network
        InferenceGraph<unsigned long> ig = bigb.to_graph();

        //TODO parametrize the type of scheduler.
        PriorityScheduler<unsigned long> scheduler(param_.getValue("loopy_belief_propagation:dampening_lambda"),
                                                   param_.getValue("loopy_belief_propagation:convergence_threshold"),
                                                   param_.getValue("loopy_belief_propagation:max_nr_iterations"));
        scheduler.add_ab_initio_edges(ig);

        BeliefPropagationInferenceEngine<unsigned long> bpie(scheduler, ig);
        auto posteriorFactors = bpie.estimate_posteriors(posteriorVars);

        //TODO play around with adjusted peptide posteriors that you could also easily request here.
        for (auto const& posteriorFactor : posteriorFactors)
        {
          double posterior = 0.0;
          IDBoostGraph::SetPosteriorVisitor pv;
          unsigned long nodeId = posteriorFactor.ordered_variables()[0];
          const PMF& pmf = posteriorFactor.pmf();
          // If Index 1 is in the range of this result PMFFactor it is non-zero
          if (1 >= pmf.first_support()[0] && 1 <= pmf.last_support()[0]) {
            posterior = pmf.table()[1 - pmf.first_support()[0]];
          }
          auto bound_visitor = std::bind(pv, std::placeholders::_1, posterior);
          boost::apply_visitor(bound_visitor, fg[nodeId]);
        }
      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  };

  struct BayesianProteinInferenceAlgorithm::GridSearchEvaluator
  {
    Param& param_;
    IDBoostGraph& ibg_;
    const ProteinIdentification& prots_;

    explicit GridSearchEvaluator(Param& param, IDBoostGraph& ibg, const ProteinIdentification& prots):
        param_(param),
        ibg_(ibg),
        prots_(prots)
    {}

    double operator() (double alpha, double beta, double gamma)
    {
      std::cout << "Evaluating: " << alpha << " " << beta << " " << gamma << std::endl;
      param_.setValue("model_parameters:prot_prior", gamma);
      param_.setValue("model_parameters:pep_emission", alpha);
      param_.setValue("model_parameters:pep_spurious_emission", beta);
      ibg_.applyFunctorOnCCs(FilteredGraphInferenceFunctor(const_cast<const Param&>(param_)));
      FalseDiscoveryRate fdr;
      return fdr.applyEvaluateProteinIDs(prots_);
    }
  };


  BayesianProteinInferenceAlgorithm::BayesianProteinInferenceAlgorithm() :
      DefaultParamHandler("BayesianProteinInference"),
      ProgressLogger()
  {
    // set default parameter values

    /* More parameter TODOs:
     * - grid search settings: e.g. fine, coarse, prob. threshold, lower convergence crit.
     * - use own groups (and regularize)
     * - use own priors
     * - multiple runs
     * - what to do about multiple charge states or modded peptides
     * - use add. pep. infos (rt, ms1dev)
     * - add dependencies on peptides in same feature and psms to same peptide (so that there is competition)
     * - ...
     */

/* TODO not yet implemented
 * defaults_.setValue("keep_threshold",
                       "false",
                       "Keep only proteins and protein groups with estimated probability higher than this threshold");
    defaults_.setValue("greedy_group_resolution",
                       "false",
                       "Post-process inference output with greedy resolution of shared peptides based on the parent protein probabilities. Also adds the resolved ambiguity groups to output.");
    defaults_.setValue("combine_indist_groups",
                       "false",
                       "Combine indistinguishable protein groups beforehand to only perform inference on them (probability for the whole group = is ANY of them present).");*/
    defaults_.setValue("annotate_groups_only",
                       "false",
                       "Skips complex inference completely and just annotates indistinguishable groups.");
    defaults_.setValue("all_PSMs",
                       "false",
                       "Consider all PSMs of each peptide, instead of only the best one.");
    defaults_.addSection("model_parameters","Model parameters for the Bayesian network");
    defaults_.setValue("model_parameters:prot_prior",
                       0.9,
                       "Protein prior probability ('gamma' parameter).");
    defaults_.setMinFloat("model_parameters:prot_prior", 0.0);
    defaults_.setMaxFloat("model_parameters:prot_prior", 1.0);
    defaults_.setValue("model_parameters:pep_emission",
                       0.1,
                       "Peptide emission probability ('alpha' parameter)");
    defaults_.setMinFloat("model_parameters:pep_emission", 0.0);
    defaults_.setMaxFloat("model_parameters:pep_emission", 1.0);
    defaults_.setValue("model_parameters:pep_spurious_emission",
                       0.001,
                       "Spurious peptide identification probability ('beta' parameter). Usually much smaller than emission from proteins");
    defaults_.setMinFloat("model_parameters:pep_spurious_emission", 0.0);
    defaults_.setMaxFloat("model_parameters:pep_spurious_emission", 1.0);

    defaults_.addSection("loopy_belief_propagation","Settings for the loopy belief propagation algorithm.");
    defaults_.setValue("loopy_belief_propagation:scheduling_type",
                       "priority",
                       "How to pick the next message:"
                           " priority = based on difference to last message (higher = more important)."
                           " fifo = first in first out."
                           " random_spanning_tree = message passing follows a random spanning tree in each iteration");
    defaults_.setValidStrings("loopy_belief_propagation:scheduling_type", {"priority","fifo","random_spanning_tree"});

    //TODO not yet implemented
/*    defaults_.setValue("loopy_belief_propagation:message_difference",
                       "MSE",
                       "How to calculate the difference of distributions in updated messages.");
    defaults_.setValidStrings("loopy_belief_propagation:message_difference", {"MSE"});*/
    defaults_.setValue("loopy_belief_propagation:convergence_threshold",
                       1e-5,
                       "Under which threshold is a message considered to be converged.");
    defaults_.setValue("loopy_belief_propagation:dampening_lambda",
                       1e-3,
                       "How strongly should messages be updated in each step. "
                           "0 = new message overwrites old completely (no dampening),"
                           "1 = old message stays (no convergence, don't do that)"
                           "In-between it will be a convex combination of both. Prevents oscillations but hinders convergence");
    defaults_.setValue("loopy_belief_propagation:max_nr_iterations",
                       1ul<<32,
                       "If not all messages converge, how many iterations should be done at max?");

    defaults_.addSection("param_optimize","Settings for the parameter optimization.");
    defaults_.setValue("param_optimize:aucweight",
                       0.2,
                       "How important is AUC vs calibration of the posteriors?"
                       " 0 = maximize calibration only,"
                       " 1 = maximize AUC only,"
                       " between = convex combination.");
    defaults_.setMinFloat("param_optimize:aucweight", 0.0);
    defaults_.setMaxFloat("param_optimize:aucweight", 1.0);


    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs)
  {
    // get enzyme settings from peptideID
    const DigestionEnzymeProtein enzyme = proteinIDs[0].getSearchParameters().digestion_enzyme;
    Size missed_cleavages = proteinIDs[0].getSearchParameters().missed_cleavages;
    EnzymaticDigestion ed{};
    ed.setEnzyme(&enzyme);
    ed.setMissedCleavages(missed_cleavages);

    std::vector<StringView> tempDigests{};
    // if not annotated, assign max nr of digests
    for (auto& protein : proteinIDs[0].getHits())
    {
      // check for existing max nr peptides metavalue annotation
      if (!protein.metaValueExists("missingTheorDigests"))
      {
        if(!protein.getSequence().empty())
        {
          tempDigests.clear();
          //TODO check which peptide lengths we should support. Parameter?
          Size nrDiscarded = ed.digestUnmodified(protein.getSequence(), tempDigests);
          //TODO add the discarded digestions products, too?
          protein.setMetaValue("missingTheorDigests", tempDigests.size());
        }
        else
        {
          //TODO Exception
          std::cerr << "Protein sequence not annotated" << std::endl;
        }
      }
    }

    //TODO would be better if we set this after inference but only here we currently have
    // non-const access.
    proteinIDs[0].setScoreType("Posterior Probability");
    proteinIDs[0].setHigherScoreBetter(true);

    // init empty graph
    IDBoostGraph ibg(proteinIDs[0], peptideIDs);
    ibg.buildGraph(param_.getValue("all_PSMs").toBool());
    ibg.computeConnectedComponents();
    ibg.annotateIndistinguishableGroups();

    //TODO how to perform group inference
    // Three options:
    // -collapse proteins to groups beforehand and run inference
    // -use the automatically created indist. groups and report their posterior
    // -calculate prior from proteins for the group beforehand and remove proteins from network (saves computation
    //  because messages are not passed from prots to groups anymore.


    //TODO Use gold search that goes deeper into the grid where it finds the best value.
    //We have to do it on a whole dataset basis though (all CCs). -> I have to refactor to actually store as much
    //as possible (it would be cool to store the inference graph but this is probably not possible bc that is why
    //I split up in CCs.
    // OR I could save the outputs! One value for every protein, per parameter set.

    vector<double> gamma_search{0.5};
    vector<double> beta_search{0.001};
    vector<double> alpha_search{0.1, 0.3, 0.5, 0.7, 0.9};
    //Percolator settings
    //vector<double> alpha_search{0.008, 0.032, 0.128};

    GridSearch<double,double,double> gs{alpha_search, beta_search, gamma_search};

    std::array<size_t, 3> bestParams{};
    //TODO run grid search on reduced graph?
    //TODO if not, think about storing results temporary and only keep the best in the end
    gs.evaluate(GridSearchEvaluator(param_, ibg, proteinIDs[0]), -1.0, bestParams);

    std::cout << "Best params found at " << bestParams[0] << "," << bestParams[1] << "," << bestParams[2] << std::endl;

    //TODO write graphfile?
    //TODO let user modify Grid for GridSearch and/or provide some more default settings
  }
}
