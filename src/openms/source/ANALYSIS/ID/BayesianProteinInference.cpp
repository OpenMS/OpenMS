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
#include <OpenMS/ANALYSIS/ID/BayesianProteinInference.h>

#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <set>

using namespace std;

namespace OpenMS
{
  class BayesianProteinInference::FilteredGraphInferenceFunctor :
  public std::function<void(IDBoostGraph::FilteredGraph&)>
  {
  public:
    FilteredGraphInferenceFunctor() = default;
    void operator() (IDBoostGraph::FilteredGraph& fg) {
      //------------------ Now actual inference ------------------- //
      // Skip cc without peptide or protein
      //TODO this currently does not work because we do not filter edges I think
      //TODO introduce edge types or skip single nodes instead
      if (boost::num_edges(fg) >= 1)
      {
        //TODO make them parameters and/or estimate with gold/grid search
        MessagePasserFactory<unsigned long> mpf (0.9,0.01,0.5,1.0);
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
            //TODO assert that there is at least one protein mapping to this peptide! Require IDFilter removeUnmatched before.
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

        int count = 1;
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
            //TODO assert that the peptide score is of type PEP!
            bigb.insert_dependency(mpf.createSumEvidenceFactor(setpair.first.size(), label, j));
            IDBoostGraph::IDPointer p = fg[j];
            bigb.insert_dependency(mpf.createPeptideEvidenceFactor(j, boost::get<PeptideHit*>(p)->getScore()));
          }

        }
        InferenceGraph<unsigned long> ig = bigb.to_graph();
        //TODO we should parametrize the type of scheduler and its params.
        PriorityScheduler<unsigned long> scheduler(0.001, 1e-8, 1ul<<32);
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

  BayesianProteinInference::BayesianProteinInference(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs) :
      DefaultParamHandler("BayesianProteinInference"),
      ProgressLogger()
  {
    // set default parameter values
    defaults_.setValue("keep_zero_group",
                       false,
                       "Keep the group of proteins with estimated probability of zero, which is otherwise removed (it may be very large)");
    defaults_.setValue("greedy_group_resolution",
                       false,
                       "Post-process Fido output with greedy resolution of shared peptides based on the protein probabilities. Also adds the resolved ambiguity groups to output.");
    defaults_.setValue("all_PSMs",
                       false,
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
    defaults_.setValue("loopy_belief_propagation:message_difference",
                       "MSE",
                       "How to calculate the difference of distributions in updated messages.");
    defaults_.setValidStrings("loopy_belief_propagation:message_difference", {"MSE"});
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

    IDBoostGraph ibg;
    FilteredGraphInferenceFunctor f;
    ibg.applyFunctorOnCCs(proteinIDs[0], peptideIDs, f);
    std::string file = "";
  }
}
