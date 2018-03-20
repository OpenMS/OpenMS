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

using namespace std;

namespace OpenMS
{

  /// A functor that specifies what to do on a connected component (IDBoostGraph::FilteredGraph)
  class BayesianProteinInferenceAlgorithm::FilteredGraphInferenceFunctorNoGroups :
  public std::function<void(IDBoostGraph::FilteredGraph&)>
  {
  public:
    const Param& param_;
    vector<ProteinIdentification::ProteinGroup>& indistGroups_;

    explicit FilteredGraphInferenceFunctorNoGroups(const Param& param, vector<ProteinIdentification::ProteinGroup>& indistGroups):
        param_(param),
        indistGroups_(indistGroups)
    {}

    void operator() (IDBoostGraph::FilteredGraph& fg) {
      //------------------ Now actual inference ------------------- //
      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types or skip nodes without neighbors inside the if instead
      //TODO do quick bruteforce calculation if the cc is really small
      if (boost::num_edges(fg) >= 1)
      {
        //TODO allow estimation via gold/grid search
        MessagePasserFactory<unsigned long> mpf (param_.getValue("model_parameters:prot_prior"),
                                                 param_.getValue("model_parameters:pep_emission"),
                                                 param_.getValue("model_parameters:pep_spurious_emission"),
                                                 1.0); // the p used for marginalization: 1 = sum product, inf = max product
        BetheInferenceGraphBuilder<unsigned long> bigb;

        // Cluster peptides with same parents
        //TODO this could be sped up by a good hashing function for sets of uints and using unordered_map
        map< set<IDBoostGraph::vertex_t>, set<IDBoostGraph::vertex_t> > pepClusters; //maps the parent (protein) set to peptides that have the same
        boost::filtered_graph<IDBoostGraph::Graph, boost::function<bool(IDBoostGraph::edge_t)>, boost::function<bool(IDBoostGraph::vertex_t)> >::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end (usually at least proteins)
        // Maybe later peptides (e.g. for an iterative procedure)
        vector<vector<unsigned long>> posteriorVars;

        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = fg[*ui];
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          if (curr_idObj.which() == 0) //it's a peptide
          {
            //TODO assert that there is at least one protein mapping to this peptide! Eg. Require IDFilter removeUnmatched before.
            //Or just check rigorously here.
            set<IDBoostGraph::vertex_t> parents;
            IDBoostGraph::FilteredGraph::adjacency_iterator adjIt, adjIt_end;
            boost::tie(adjIt, adjIt_end) = boost::adjacent_vertices(*ui, fg);
            for (; adjIt != adjIt_end; ++adjIt)
            {
              if (fg[*adjIt].which() == 1) //if there are only two types (pep,prot) this check for prot is actually unnecessary
              {
                parents.insert(*adjIt);
              }
            }

            auto clusterIt = pepClusters.find(parents);
            if (clusterIt != pepClusters.end())
            {
              clusterIt->second.insert(*ui);
            }
            else
            {
              pepClusters[parents] = set<IDBoostGraph::vertex_t>({*ui});
            }
          }
          else if (curr_idObj.which() == 1)
          {
            //TODO allow an already present prior probability here
            bigb.insert_dependency(mpf.createProteinFactor(*ui));
            posteriorVars.push_back({*ui});
          }
        }

        bool annotateIndistGroupsOnly = param_.getValue("annotate_groups_only").toBool();
        if (!annotateIndistGroupsOnly)
        {
          int count = 1;
          //TODO all setpair.first entries are actually indistinguishable prot groups.
          //Add them to the group object (with param?)
          //Allow user to use the whole group as one entity (the old Fido group option).
          for (auto const& setpair : pepClusters)
          {
            #ifdef INFERENCE_DEBUG
            for (auto const& j : setpair.first)
            std::cout << j << ",";

          std::cout << ": ";
          for (auto const& j : setpair.second)
            std::cout << j << ",";

          std::cout << std::endl;
            #endif
            unsigned long label = boost::num_vertices(fg) + (count++);
            bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(setpair.first, label));

            for (auto const& j : setpair.second) // foreach peptide
            {
              //TODO assert that the peptide score is of type PEP! Probably enough in the beginning.
              bigb.insert_dependency(mpf.createSumEvidenceFactor(setpair.first.size(), label, j));
              IDBoostGraph::IDPointer p = fg[j];
              bigb.insert_dependency(mpf.createPeptideEvidenceFactor(j, boost::get<PeptideHit*>(p)->getScore()));
            }

          }
          InferenceGraph<unsigned long> ig = bigb.to_graph();

          //TODO parametrize the type of scheduler.
          PriorityScheduler<unsigned long> scheduler(param_.getValue("loopy_belief_propagation:dampening_lambda"),
                                                     param_.getValue("loopy_belief_propagation:convergence_threshold"),
                                                     param_.getValue("loopy_belief_propagation:max_nr_iterations"));
          scheduler.add_ab_initio_edges(ig);

          BeliefPropagationInferenceEngine<unsigned long> bpie(scheduler, ig);
          auto posteriorFactors = bpie.estimate_posteriors(posteriorVars);

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

        // Afterwards go through all pepClusters, i.e. indist. groups and add them
        for (auto const& setpair : pepClusters)
        {
          ProteinIdentification::ProteinGroup pg;
          double groupProb = 0.0;
          for (auto const& proteinVID : setpair.first)
          {
            // Usually the probabilities in a group should be the same. But there could be slight inexactnesses due to
            // message passing order and especially user priors.
            ProteinHit* proteinPtr = boost::get<ProteinHit*>(fg[proteinVID]);
            double currScore = proteinPtr->getScore();
            if (currScore > 0.001)
            {
              pg.accessions.push_back(proteinPtr->getAccession());
              if (currScore > groupProb) //use max. prob. for the group
              {
                groupProb = currScore;
              }
            }
          }
          if (!pg.accessions.empty())
          {
            pg.probability = groupProb;
            indistGroups_.push_back(pg);
          }

        }

      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
      }
    }
  };

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
        MessagePasserFactory<unsigned long> mpf (param_.getValue("model_parameters:prot_prior"),
                                                 param_.getValue("model_parameters:pep_emission"),
                                                 param_.getValue("model_parameters:pep_spurious_emission"),
                                                 1.0); // the p used for marginalization: 1 = sum product, inf = max product
        BetheInferenceGraphBuilder<unsigned long> bigb;

        boost::filtered_graph<IDBoostGraph::Graph, boost::function<bool(IDBoostGraph::edge_t)>, boost::function<bool(IDBoostGraph::vertex_t)> >::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end (usually at least proteins)
        // Maybe later peptides (e.g. for an iterative procedure)
        vector<vector<unsigned long>> posteriorVars;

        //int count = 1; // count for the (hidden) sumfactors
        for (; ui != ui_end; ++ui)
        {
          IDBoostGraph::IDPointer curr_idObj = fg[*ui];
          //TODO introduce an enum for the types to make it more clear.
          //Or use the static_visitor pattern: You have to pass the vertex with its neighbors as a second arg though.
          //TODO we could actually filter for current connected component AND protein group type and then work our way
          //to the PSMs from there in the original graph.
          if (curr_idObj.which() == 2) //it's a protein group
          {
            vector<IDBoostGraph::vertex_t> indistProts;
            vector<IDBoostGraph::vertex_t> sharedPeps;

            // direct neighbors of indist. groups are proteins
            IDBoostGraph::FilteredGraph::adjacency_iterator protVIt, protVIt_end;
            boost::tie(protVIt, protVIt_end) = boost::adjacent_vertices(*ui, fg);

            // neighbors of the first protein in the group are the shared peptides or the group again.
            IDBoostGraph::FilteredGraph::adjacency_iterator pepVIt, pepVIt_end;
            boost::tie(pepVIt, pepVIt_end) = boost::adjacent_vertices(*protVIt, fg);

            for (; protVIt != protVIt_end; ++protVIt)
            {
              //TODO allow an already present prior probability here
              bigb.insert_dependency(mpf.createProteinFactor(*protVIt));
              posteriorVars.push_back({*protVIt});
              indistProts.push_back(*protVIt);
            }

            bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(indistProts, *ui));

            for (int j = 0; pepVIt != pepVIt_end; ++pepVIt, ++j)
            {
              //TODO if we use a directed graph, this can be avoided. Might actually be the better choice.
              if (fg[*pepVIt].which() == 0)
              {
                sharedPeps.emplace_back(*pepVIt);
                bigb.insert_dependency(mpf.createSumEvidenceFactor(indistProts.size(), *ui, *pepVIt));
              }
            }

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
     * - use add. pep. infos (rt, ms1dev)
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

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs)
  {
    // init empty graph
    IDBoostGraph ibg(proteinIDs[0], peptideIDs);
    ibg.buildGraph(param_.getValue("all_PSMs").toBool());
    ibg.computeConnectedComponents();
    ibg.annotateIndistinguishableGroups();

    //TODO if we only perform group inference I think we should use another functor.
    // It is too different from the normal type of inference, since we kind of create a new "protein-like" entity
    // for every indist group

    //TODO perform parameter search here with Statistics::rocN(50), a lambda for convex combination and the difference
    //between target-decoy and posterior FDR (to be implemented)
    //Use gold search that goes deeper into the grid where it finds the best value.
    //We have to do it on a whole dataset basis though (all CCs). -> I have to refactor to actually store as much
    //as possible (it would be cool to store the inference graph but this is probably not possible bc that is why
    //I split up in CCs.
    // OR I could save the outputs! One value for every protein, per parameter set.

    vector<double> gamma_search{0.5};
    vector<double> beta_search{0.001};
    vector<double> alpha_search{0.008, 0.032, 0.128};

    GridSearch<double,double,double> gs{alpha_search, beta_search, gamma_search};

    std::array<size_t, 3> bestParams{};
    gs.evaluate(GridSearchEvaluator(param_, ibg, proteinIDs[0]), -1.0, bestParams);

    // what to do, with which params and where to write additional output (here groups) not accessible via the graphs.
    FilteredGraphInferenceFunctor f{param_};
    // apply functor
    ibg.applyFunctorOnCCs(f);
    //TODO write graphfile?
    //TODO let user modify Grid for GridSearch and/or provide some more default settings
  }
}
