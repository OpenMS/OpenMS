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
#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/openms_package_version.h>

#include <set>


using namespace std;

namespace OpenMS
{
  /// Only works if ProteinGroup nodes are present, which is the case when used in this class
  //TODO private
  class BayesianProteinInferenceAlgorithm::AnnotateIndistGroupsFunctor :
      public std::function<void(IDBoostGraph::Graph&)>
  {
  private:
    ProteinIdentification& prots;
  public:
    explicit AnnotateIndistGroupsFunctor(ProteinIdentification& proteinIDToAnnotateGroups):
    prots(proteinIDToAnnotateGroups)
    {}

    void operator() (IDBoostGraph::Graph& fg) {
      // this skips CCs with just peps or prots. We only add edges between different types.
      // and if there were no edges, it would not be a CC.
      if (boost::num_vertices(fg) >= 2)
      {
        IDBoostGraph::Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        for (; ui != ui_end; ++ui)
        {
          if (fg[*ui].which() == 1) //prot group
          {
            ProteinIdentification::ProteinGroup pg{};
            pg.probability = (double) boost::get<IDBoostGraph::ProteinGroup>(fg[*ui]); //init
            IDBoostGraph::Graph::adjacency_iterator nbIt, nbIt_end;
            boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

            ProteinHit *proteinPtr = nullptr;
            for (; nbIt != nbIt_end; ++nbIt)
            {
              if (fg[*nbIt].which() == 0)
              {
                proteinPtr = boost::get<ProteinHit*>(fg[*nbIt]);
                pg.accessions.push_back(proteinPtr->getAccession());
              }
            }
            prots.getIndistinguishableProteins().push_back(pg);
          }
        }
      }
    }
  };

  /// A functor that specifies what to do on a connected component (IDBoostGraph::FilteredGraph)
  class BayesianProteinInferenceAlgorithm::GraphInferenceFunctor :
      public std::function<void(IDBoostGraph::Graph&)>
  {
  public:
    //TODO think about restructuring params (we do not need every param from the BPI class here.
    const Param& param_;
    unsigned int debug_lvl_;
    unsigned long cnt_;

    explicit GraphInferenceFunctor(const Param& param, unsigned int debug_lvl):
        param_(param),
        debug_lvl_(debug_lvl),
        cnt_(0)
    {}

    unsigned long operator() (IDBoostGraph::Graph& fg) {
      //TODO do quick bruteforce calculation if the cc is really small?
      cnt_++;
      // this skips CCs with just peps or prots. We only add edges between different types.
      // and if there were no edges, it would not be a CC.
      if (boost::num_vertices(fg) >= 2)
      {
        bool graph_mp_ownership_acquired = false;
        bool update_PSM_probabilities = param_.getValue("update_PSM_probabilities").toBool();
        bool annotate_group_posterior = param_.getValue("annotate_group_probabilities").toBool();
        bool user_defined_priors = param_.getValue("user_defined_priors").toBool();
        bool regularize = param_.getValue("model_parameters:regularize").toBool();
        double pnorm = param_.getValue("loopy_belief_propagation:p_norm_inference");
        if (pnorm <= 0)
        {
          pnorm = std::numeric_limits<double>::infinity();
        }

        MessagePasserFactory<IDBoostGraph::vertex_t> mpf (param_.getValue("model_parameters:pep_emission"),
                                                 param_.getValue("model_parameters:pep_spurious_emission"),
                                                 param_.getValue("model_parameters:prot_prior"),
                                                 pnorm,
                                                 param_.getValue("model_parameters:pep_prior")); // the p used for marginalization: 1 = sum product, inf = max product
        BetheInferenceGraphBuilder<IDBoostGraph::vertex_t> bigb;

        IDBoostGraph::Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end
        vector<vector<IDBoostGraph::vertex_t>> posteriorVars;

        // direct neighbors are proteins on the "left" side and peptides on the "right" side
        // TODO Can be sped up using directed graph. Needs some restructuring in IDBoostGraph class first tho.
        std::vector<IDBoostGraph::vertex_t> in{};
        //std::vector<IDBoostGraph::vertex_t> out{};

        //TODO the try section could in theory be slimmed down a little bit. First use of insertDependency maybe.
        // check performance impact.
        try
        {
          for (; ui != ui_end; ++ui)
          {
            IDBoostGraph::Graph::adjacency_iterator nbIt, nbIt_end;
            boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

            in.clear();
            //out.clear(); // we dont need out edges currently

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

            if (fg[*ui].which() == 6) // pep hit = psm
            {
              if (regularize)
              {
                bigb.insert_dependency(mpf.createRegularizingSumEvidenceFactor(boost::get<PeptideHit *>(fg[*ui])
                                                                                   ->getPeptideEvidences().size(), in[0], *ui));
              }
              else
              {
                bigb.insert_dependency(mpf.createSumEvidenceFactor(boost::get<PeptideHit *>(fg[*ui])
                                                                                   ->getPeptideEvidences().size(), in[0], *ui));
              }

              bigb.insert_dependency(mpf.createPeptideEvidenceFactor(*ui,
                                                                     boost::get<PeptideHit *>(fg[*ui])->getScore()));
              if (update_PSM_probabilities)
              {
                posteriorVars.push_back({*ui});
              }
            }
            else if (fg[*ui].which() == 2) // pep group
            {
              bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(in, *ui));
            }
            else if (fg[*ui].which() == 1) // prot group
            {
              bigb.insert_dependency(mpf.createPeptideProbabilisticAdderFactor(in, *ui));
              if (annotate_group_posterior)
              {
                posteriorVars.push_back({*ui});
              }
            }
            else if (fg[*ui].which() == 0) // prot
            {
              //TODO modify createProteinFactor to start with a modified prior based on the number of missing
              // peptides (later tweak to include conditional prob. for that peptide
              if (user_defined_priors)
              {
                bigb.insert_dependency(mpf.createProteinFactor(*ui,
                                                               (double) boost::get<ProteinHit *>(fg[*ui])
                                                                   ->getMetaValue("Prior")));
              }
              else
              {
                bigb.insert_dependency(mpf.createProteinFactor(*ui));
              }
              posteriorVars.push_back({*ui});
            }
          }

          // create factor graph for Bayesian network
          InferenceGraph < IDBoostGraph::vertex_t > ig = bigb.to_graph();
          graph_mp_ownership_acquired = true;

          unsigned long maxMessages = param_
              .getValue("loopy_belief_propagation:max_nr_iterations");
          double initDampeningLambda = param_
              .getValue("loopy_belief_propagation:dampening_lambda");
          double initConvergenceThreshold = param_.getValue(
              "loopy_belief_propagation:convergence_threshold");
          unsigned long nrEdges = boost::num_edges(fg);

          //TODO parametrize the type of scheduler.
          PriorityScheduler<IDBoostGraph::vertex_t> scheduler(initDampeningLambda,
                                                              initConvergenceThreshold,
                                                              maxMessages);
          scheduler.add_ab_initio_edges(ig);

          BeliefPropagationInferenceEngine<IDBoostGraph::vertex_t> bpie(scheduler, ig);

          auto posteriorFactors = bpie.estimate_posteriors_in_steps(posteriorVars,
              {
                  std::make_tuple(std::max<unsigned long>(10000ul, nrEdges*nrEdges*2ul), initDampeningLambda, initConvergenceThreshold),
                  std::make_tuple(nrEdges*nrEdges, std::min(0.5,initDampeningLambda*10), std::min(0.01,initConvergenceThreshold*10)),
                  std::make_tuple(nrEdges*nrEdges/2ul, std::min(0.5,initDampeningLambda*100), std::min(0.01,initConvergenceThreshold*100))
              });

          // TODO move the writing of statistics from IDBoostGraph here and write more stats
          //  like nr messages and failure/success
          unsigned long nrMessagesNeeded = bpie.getNrMessagesPassed();

          for (auto const &posteriorFactor : posteriorFactors)
          {
            double posterior = 1.0;
            IDBoostGraph::SetPosteriorVisitor pv;
            IDBoostGraph::vertex_t nodeId = posteriorFactor.ordered_variables()[0];
            const PMF &pmf = posteriorFactor.pmf();
            // If Index 0 is in the range of this result PMFFactor is probability is non-zero
            // and the prob of presence is 1-P(p=0). Important in multi-value factors like protein groups.
            if (0 >= pmf.first_support()[0] && 0 <= pmf.last_support()[0])
            {
              posterior = 1. - pmf.table()[0ul];
            }
            auto bound_visitor = std::bind(pv, std::placeholders::_1, posterior);
            boost::apply_visitor(bound_visitor, fg[nodeId]);
          }
          //TODO we could write out/save the posteriors here,
          // so we can easily read them later for the best params of the grid search
          return nrMessagesNeeded;
        }
        catch (const std::runtime_error& /*e*/)
        {
          //TODO print failing component and implement the following options
          // 1) Leave posteriors (e.g. if Percolator was ran before. Make sure they are PPs not PEPs)
          // 2) set posteriors to priors
          // 3) try another type of inference on that connected component. Different scheduler,
          //    different extreme probabilities or maybe best: trivial aggregation-based inference.
          // 4) Cancelling this and all other threads/ the loop and call this set of parameters invalid

          //For now we just warn and continue with the rest of the iterations. Might still be a valid run.

          // Graph builder needs to build otherwise it leaks memory.
          if (!graph_mp_ownership_acquired) bigb.to_graph();

          if (debug_lvl_ > 2)
          {
            std::ofstream ofs;
            ofs.open ("failed_cc_a"+ String(param_.getValue("model_parameters:pep_emission")) +
                "_b" + String(param_.getValue("model_parameters:pep_spurious_emission")) + "_g" +
                String(param_.getValue("model_parameters:prot_prior")) + "_c" +
                String(param_.getValue("model_parameters:pep_prior")) + "_p" + String(pnorm) + "_"
                + String(cnt_) + ".graphviz"
                , std::ofstream::out | std::ofstream::app);
            IDBoostGraph::printGraph(ofs, fg);
          }
          LOG_WARN << "Warning: Loopy belief propagation encountered a problem in a connected component. Skipping"
                      " inference there." << std::endl;
          return 0;
        }
      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
        return 0;
      }
    }
  };

  /// A functor that specifies what to do on a connected component (IDBoostGraph::FilteredGraph)
  class BayesianProteinInferenceAlgorithm::ExtendedGraphInferenceFunctor :
      public std::function<void(IDBoostGraph::Graph&)>
  {
  public:
    const Param& param_;

    explicit ExtendedGraphInferenceFunctor(const Param& param):
        param_(param)
    {}

    unsigned long operator() (IDBoostGraph::Graph& fg) {
      //TODO do quick bruteforce calculation if the cc is really small

      double pnorm = param_.getValue("loopy_belief_propagation:p_norm_inference");
      if (pnorm <= 0)
      {
        pnorm = std::numeric_limits<double>::infinity();
      }

      // this skips CCs with just peps or prots. We only add edges between different types.
      // and if there were no edges, it would not be a CC.
      if (boost::num_vertices(fg) >= 2)
      {
        MessagePasserFactory<IDBoostGraph::vertex_t> mpf (param_.getValue("model_parameters:pep_emission"),
                                                 param_.getValue("model_parameters:pep_spurious_emission"),
                                                 param_.getValue("model_parameters:prot_prior"),
                                                 pnorm,
                                                 param_.getValue("model_parameters:pep_prior")); // the p used for marginalization: 1 = sum product, inf = max product

        BetheInferenceGraphBuilder<IDBoostGraph::vertex_t> bigb;

        IDBoostGraph::Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end (usually at least proteins)
        // Maybe later peptides (e.g. for an iterative procedure)
        vector<vector<IDBoostGraph::vertex_t>> posteriorVars;

        // direct neighbors are proteins on the "left" side and peptides on the "right" side
        // TODO can be sped up using directed graph. Requires some restructuring first.
        std::vector<IDBoostGraph::vertex_t> in{};
        //std::vector<IDBoostGraph::vertex_t> out{};

        //TODO the try section could in theory be slimmed down a little bit. First use of insertDependency maybe.
        // check performance impact.
        try
        {
          for (; ui != ui_end; ++ui)
          {
            IDBoostGraph::Graph::adjacency_iterator nbIt, nbIt_end;
            boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

            in.clear();
            //out.clear(); // we dont need out edges currently

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

            if (fg[*ui].which() == 6) // pep hit = psm
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
              //TODO modify createProteinFactor to start with a modified prior based on the number of missing
              // peptides (later tweak to include conditional prob. for that peptide
              bigb.insert_dependency(mpf.createProteinFactor(*ui));
              posteriorVars.push_back({*ui});
            }
          }

          // create factor graph for Bayesian network
          InferenceGraph<IDBoostGraph::vertex_t> ig = bigb.to_graph();

          //TODO parametrize the type of scheduler.
          PriorityScheduler<IDBoostGraph::vertex_t> scheduler(param_.getValue("loopy_belief_propagation:dampening_lambda"),
                                                     param_.getValue("loopy_belief_propagation:convergence_threshold"),
                                                     param_.getValue("loopy_belief_propagation:max_nr_iterations"));
          scheduler.add_ab_initio_edges(ig);

          BeliefPropagationInferenceEngine<IDBoostGraph::vertex_t> bpie(scheduler, ig);

          auto posteriorFactors = bpie.estimate_posteriors(posteriorVars);

          //TODO you could also save the indices of the peptides here and request + update their posteriors, too.
          for (auto const &posteriorFactor : posteriorFactors)
          {
            double posterior = 0.0;
            IDBoostGraph::SetPosteriorVisitor pv;
            unsigned long nodeId = posteriorFactor.ordered_variables()[0];
            const PMF &pmf = posteriorFactor.pmf();
            // If Index 1 is in the range of this result PMFFactor it is non-zero
            if (1 >= pmf.first_support()[0] && 1 <= pmf.last_support()[0])
            {
              posterior = pmf.table()[1 - pmf.first_support()[0]];
            }
            auto bound_visitor = std::bind(pv, std::placeholders::_1, posterior);
            boost::apply_visitor(bound_visitor, fg[nodeId]);
          }
          //TODO update to use actual nr of messages
          return 1;
        }
        catch (const std::runtime_error& /*e*/)
        {
          //TODO print failing component
          // set posteriors to priors or try another type of inference?
          // Think about cancelling all other threads/ the loop
          //For now we just warn and continue with the rest of the iterations. Might still be a valid run.

          // Graph builder needs to build otherwise it leaks memory.
          bigb.to_graph();
          LOG_WARN << "Warning: Loopy belief propagation encountered a problem in a connected component. Skipping"
                      "inference there." << std::endl;
          return 0;
        }
        //TODO we could write out the posteriors here, so we can easily read them for the best params of the grid search
      }
      else
      {
        std::cout << "Skipped cc with only one type (proteins or peptides)" << std::endl;
        return 0;
      }
    }
  };

  struct BayesianProteinInferenceAlgorithm::GridSearchEvaluator
  {
    Param& param_;
    IDBoostGraph& ibg_;
    const ProteinIdentification& prots_;
    const unsigned int debug_lvl_;

    explicit GridSearchEvaluator(Param& param, IDBoostGraph& ibg, const ProteinIdentification& prots, unsigned int debug_lvl):
        param_(param),
        ibg_(ibg),
        prots_(prots),
        debug_lvl_(debug_lvl)
    {}

    double operator() (double alpha, double beta, double gamma)
    {
      std::cout << "Evaluating: " << alpha << " " << beta << " " << gamma << std::endl;
      param_.setValue("model_parameters:prot_prior", gamma);
      param_.setValue("model_parameters:pep_emission", alpha);
      param_.setValue("model_parameters:pep_spurious_emission", beta);
      ibg_.applyFunctorOnCCs(GraphInferenceFunctor(const_cast<const Param&>(param_), debug_lvl_));
      FalseDiscoveryRate fdr;
      Param fdrparam = fdr.getParameters();
      fdrparam.setValue("conservative",param_.getValue("param_optimize:conservative_fdr"));
      fdr.setParameters(fdrparam);
      return fdr.applyEvaluateProteinIDs(prots_, 1.0, 50, static_cast<double>(param_.getValue("param_optimize:aucweight")));
    }
  };


  BayesianProteinInferenceAlgorithm::BayesianProteinInferenceAlgorithm(unsigned int debug_lvl) :
      DefaultParamHandler("BayesianProteinInferenceAlgorithm"),
      ProgressLogger(),
      debug_lvl_(debug_lvl)
  {
    // set default parameter values

    /* More parameter TODOs:
     * - grid search settings: e.g. fine, coarse, prob. threshold, lower convergence crit., own lists
     * - use own groups (and regularize)
     * - multiple runs
     * - what to do about multiple charge states or modded peptides
     * - use add. pep. infos (rt, ms1dev)
     * - add dependencies on peptides in same feature and psms to same peptide (so that there is competition)
     * - option to write graphfile?
     */
    /*
    defaults_.setValue("combine_indist_groups",
                       "false",
                       "Combine indistinguishable protein groups beforehand to only perform inference on them (probability for the whole group = is ANY of them present).");*/

    defaults_.setValue("top_PSMs",
                       1,
                       "Consider only top X PSMs per spectrum. 0 considers all.");
    defaults_.setMinInt("top_PSMs", 0);

    defaults_.setValue("update_PSM_probabilities",
                       "true",
                       "(Experimental:) Update PSM probabilities with their posteriors under consideration of the protein probabilities.");
    defaults_.setValidStrings("update_PSM_probabilities", {"true","false"});

    defaults_.setValue("user_defined_priors",
                       "false",
                       "(Experimental:) Uses the current protein scores as user-defined priors.");
    defaults_.setValidStrings("user_defined_priors", {"true","false"});

    defaults_.setValue("annotate_group_probabilities",
                       "true",
                       "Annotates group probabilities for indistinguishable protein groups (indistinguishable by "
                       "experimentally observed PSMs).");
    defaults_.setValidStrings("annotate_group_probabilities", {"true","false"});

    defaults_.addSection("model_parameters","Model parameters for the Bayesian network");

    defaults_.setValue("model_parameters:prot_prior",
                       -1.,
                       "Protein prior probability ('gamma' parameter). Negative values enable grid search for this param.");
    defaults_.setMinFloat("model_parameters:prot_prior", -1.0);
    defaults_.setMaxFloat("model_parameters:prot_prior", 1.0);

    defaults_.setValue("model_parameters:pep_emission",
                       -1.,
                       "Peptide emission probability ('alpha' parameter). Negative values enable grid search for this param.");
    defaults_.setMinFloat("model_parameters:pep_emission", -1.0);
    defaults_.setMaxFloat("model_parameters:pep_emission", 1.0);

    defaults_.setValue("model_parameters:pep_spurious_emission",
                       -1.,
                       "Spurious peptide identification probability ('beta' parameter)."
                       " Usually much smaller than emission from proteins. "
                       "Negative values enable grid search for this param.");
    defaults_.setMinFloat("model_parameters:pep_spurious_emission", -1.0);
    defaults_.setMaxFloat("model_parameters:pep_spurious_emission", 1.0);

    defaults_.setValue("model_parameters:pep_prior",
                       0.1,
                       "Peptide prior probability (experimental, should be covered by combinations of the other params).");
    defaults_.setMinFloat("model_parameters:pep_prior", 0.0);
    defaults_.setMaxFloat("model_parameters:pep_prior", 1.0);

    defaults_.setValue("model_parameters:regularize",
                       "false",
                       "Regularize the number of proteins that produce a peptide together (experimental, should be activated when using higher p-norms).");
    defaults_.setValidStrings("model_parameters:regularize",{"true","false"});

    defaults_.addSection("loopy_belief_propagation","Settings for the loopy belief propagation algorithm.");

    defaults_.setValue("loopy_belief_propagation:scheduling_type",
                       "priority",
                       "(Not used yet) How to pick the next message:"
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
                       "Initial threshold under which MSE difference a message is considered to be converged.");
    defaults_.setMinFloat("loopy_belief_propagation:convergence_threshold", 1e-7);
    defaults_.setMaxFloat("loopy_belief_propagation:convergence_threshold", 1.0);

    defaults_.setValue("loopy_belief_propagation:dampening_lambda",
                       1e-3,
                       "Initial value for how strongly should messages be updated in each step. "
                           "0 = new message overwrites old completely (no dampening),"
                           "1 = old message stays (no convergence, don't do that)"
                           "In-between it will be a convex combination of both. Prevents oscillations but hinders convergence.");
    defaults_.setMinFloat("loopy_belief_propagation:dampening_lambda", 1e-7);
    defaults_.setMaxFloat("loopy_belief_propagation:dampening_lambda", 0.5);

    defaults_.setValue("loopy_belief_propagation:max_nr_iterations",
                       1ul<<31,
                       "(Unused, autodetermined) If not all messages converge, how many iterations should be done at max?");
    //I think restricting does not work because it only works for type Int (= int)
    //defaults_.setMinInt("loopy_belief_propagation:max_nr_iterations", 10);

    defaults_.setValue("loopy_belief_propagation:p_norm_inference",
                       1.0,
                       "P-norm used for marginalization of multidimensional factors. "
                       "1 == sum-product inference (all configurations vote equally) (default),"
                       "<= 0 == infinity = max-product inference (only best configurations propagate)"
                       "The higher the value the more important high probability configurations get."
                       );

    defaults_.addSection("param_optimize","Settings for the parameter optimization.");
    defaults_.setValue("param_optimize:aucweight",
                       0.2,
                       "How important is AUC vs calibration of the posteriors?"
                       " 0 = maximize calibration only,"
                       " 1 = maximize AUC only,"
                       " between = convex combination.");
    defaults_.setMinFloat("param_optimize:aucweight", 0.0);
    defaults_.setMaxFloat("param_optimize:aucweight", 1.0);

    defaults_.setValue("param_optimize:conservative_fdr",
                          "true",
                          "Use (D+1)/(T) instead of (D+1)/(T+D) for parameter estimation.");
    defaults_.setValidStrings("param_optimize:conservative_fdr", {"true","false"});


    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }


/* TODO under construction
  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(
      ConsensusMap cmap,
      const ExperimentalDesign& expDesign,
      const String& db,
      const ProteaseDigestion pd,
      ProteinIdentification& proteinIds
      )
  {
    proteinIds.getHits().clear();

    FASTAContainer<TFI_File> fastaDB(db);
    //assume that fasta is given for now. (you could use protein ids in the input idXMLs instead and merge)
    //run peptideindexer (with integrated missing peptide counting or after?)
    PeptideIndexing pi{};
    Param piparams = pi.getDefaults();
    piparams.setValue("annotate_nr_theoretical_peptides", true);
    piparams.setValue("enzyme:name", pd.getEnzymeName());
    //TODO expose enzyme specificity. Default Full
    //piparams.setValue("enzyme:specificity", full)
    pi.setParameters(piparams);
    std::vector<ProteinIdentification> singletonProteinId{proteinIds};

    //Next line: PepIdxr needs to accept ConsensusMap and ExpDesign and annotate one ProtID per sample (frac group/rep combo)
    //pi.run<FASTAContainer<TFI_File>>(fastaDB, singletonProteinId, pepIdConcatReplicates);


    //TODO would be better if we set this after inference but only here we currently have
    // non-const access.
    proteinIds.setScoreType("Posterior Probability");
    proteinIds.setHigherScoreBetter(true);

    // init empty graph
    IDBoostGraph ibg(proteinIds, pepIdConcatReplicates);
    ibg.buildGraph(param_.getValue("nr_PSMs").toBool());
    ibg.computeConnectedComponents();
    ibg.clusterIndistProteinsAndPeptides();

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
    gs.evaluate(GridSearchEvaluator(param_, ibg, proteinIds), -1.0, bestParams);

    std::cout << "Best params found at " << bestParams[0] << "," << bestParams[1] << "," << bestParams[2] << std::endl;

    //TODO write graphfile?
    //TODO let user modify Grid for GridSearch and/or provide some more default settings
  }*/

/*    void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs, OpenMS::ExperimentalDesign expDesign)
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
      if (!protein.metaValueExists("maxNrTheoreticalDigests"))
      {
        if(!protein.getSequence().empty())
        {
          tempDigests.clear();
          //TODO check which peptide lengths we should support. Parameter?
          ed.digestUnmodified(protein.getSequence(), tempDigests);
          //TODO add the discarded digestions products, too?
          protein.setMetaValue("maxNrTheoreticalDigests", tempDigests.size());
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
    ibg.buildGraph(param_.getValue("nr_PSMs").toInt());
    ibg.computeConnectedComponents();
    ibg.clusterIndistProteinsAndPeptides();

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
  }*/


  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs)
  {

    //TODO The following is a sketch to think about how to include missing peptides
    // Requirement: Datastructures for peptides first
    // Options:
    // - Require annotation from peptideindexer
    // - Require sequence from peptideindexer
    // - Require fasta file and annotate here
    /*
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
          //Size nrDiscarded =
          ed.digestUnmodified(protein.getSequence(), tempDigests);
          //TODO add the discarded digestions products, too?
          protein.setMetaValue("missingTheorDigests", tempDigests.size());
        }
        else
        {
          //TODO Exception
          std::cerr << "Protein sequence not annotated" << std::endl;
        }
      }
    }*/

    //TODO would be better if we set this after inference but only here we currently have
    // non-const access.
    //TODO if you do peptide rescoring, you could/should actually set the search_engine parameter
    //TODO add the most important settings of Epifany in SearchParameter metavalues, like PercolatorAdapter.
    proteinIDs[0].setScoreType("Posterior Probability");
    proteinIDs[0].setInferenceEngine("Epifany");
    proteinIDs[0].setInferenceEngineVersion(OPENMS_PACKAGE_VERSION);
    proteinIDs[0].setHigherScoreBetter(true);

    // init empty graph
    IDBoostGraph ibg(proteinIDs[0], peptideIDs);

    FalseDiscoveryRate pepFDR;
    Param p = pepFDR.getParameters();

    // On second thought, I think it is best to always use the best PSM only
    // inference might change the ranking.
    // OLDTODO adapt, so it matches top_PSMs parameter (which is not bool).
    // or just do a filter according to top_PSMs in the beginning and then
    // use all here.
    p.setValue("use_all_hits", "false");
    pepFDR.setParameters(p);
    LOG_INFO << "Peptide FDR AUC before protein inference: " << pepFDR.rocN(peptideIDs, 0) << std::endl;

    bool user_defined_priors = param_.getValue("user_defined_priors").toBool();
    if (user_defined_priors)
    {
      // Save current protein score into a metaValue
      for (auto& prot_hit : proteinIDs[0].getHits())
      {
        prot_hit.setMetaValue("Prior", prot_hit.getScore());
      }
    }

    bool oldway = true;
    if (oldway)
    {
      ibg.buildGraph(param_.getValue("top_PSMs"));
      ibg.computeConnectedComponents();
      ibg.clusterIndistProteinsAndPeptides();

      //Note: How to perform group inference
      // Three options:
      // -(implemented) use the automatically created indist. groups and report their posterior
      // - collapse proteins to groups beforehand and run inference
      // - calculate prior from proteins for the group beforehand and remove proteins from network (saves computation
      //  because messages are not passed from prots to groups anymore.


      //TODO Use gold search that goes deeper into the grid where it finds the best value.
      //We have to do it on a whole dataset basis though (all CCs). -> I have to refactor to actually store as much
      //as possible (it would be cool to store the inference graph but this is probably not possible bc that is why
      //I split up in CCs.
      // OR you could save the outputs! One value for every protein, per parameter set.

      // Do not expand gamma_search when user_defined_priors is on. Would be unused.
      double alpha = param_.getValue("model_parameters:pep_emission");
      double beta = param_.getValue("model_parameters:pep_spurious_emission");
      double gamma = param_.getValue("model_parameters:prot_prior");
      vector<double> gamma_search;
      vector<double> beta_search;
      vector<double> alpha_search;
      if (gamma > 1.0 || gamma < 0.0)
      {
        gamma_search = {0.2, 0.5, 0.7};
      }
      else
      {
        gamma_search = {gamma};
      }
      if (beta > 1.0 || beta < 0.0)
      {
        beta_search = {0.001};
      }
      else
      {
        beta_search = {beta};
      }
      if (alpha > 1.0 || alpha < 0.0)
      {
        alpha_search = {0.1, 0.25, 0.5, 0.65, 0.8};
      }
      else
      {
        alpha_search = {alpha};
      }

      GridSearch<double,double,double> gs{alpha_search, beta_search, gamma_search};

      std::array<size_t, 3> bestParams{{0, 0, 0}};

      //Save initial settings and deactivate certain features to save time during grid search and to not
      // interfere with later runs.
      // TODO We could think about optimizing PSM FDR as another goal though.
      bool update_PSM_probabilities = param_.getValue("update_PSM_probabilities").toBool();
      param_.setValue("update_PSM_probabilities","false");

      bool annotate_group_posteriors = param_.getValue("annotate_group_probabilities").toBool();
      param_.setValue("annotate_group_probabilities","false");

      //TODO run grid search on reduced graph? Then make sure, untouched protein/peps do not affect evaluation results.
      //TODO if not, think about storing results temporary (file? mem?) and only keep the best in the end
      //TODO think about running grid search on the small CCs only (maybe it's enough)
      if (gs.getNrCombos() > 1)
      {
        std::cout << "Testing " << gs.getNrCombos() << " param combinations." << std::endl;
        /*double res =*/ gs.evaluate(GridSearchEvaluator(param_, ibg, proteinIDs[0], debug_lvl_), -1.0, bestParams);
      }
      else
      {
        std::cout << "Only one combination specified: Skipping grid search." << std::endl;
      }

      double bestGamma = gamma_search[bestParams[2]];
      double bestBeta = beta_search[bestParams[1]];
      double bestAlpha = alpha_search[bestParams[0]];
      std::cout << "Best params found at a=" << bestAlpha << ", b=" << bestBeta << ", g=" << bestGamma << std::endl;
      std::cout << "Running with best parameters:" << std::endl;
      param_.setValue("model_parameters:prot_prior", bestGamma);
      param_.setValue("model_parameters:pep_emission", bestAlpha);
      param_.setValue("model_parameters:pep_spurious_emission", bestBeta);
      // Reset original values for those two options
      param_.setValue("update_PSM_probabilities", update_PSM_probabilities ? "true" : "false");
      param_.setValue("annotate_group_probabilities", annotate_group_posteriors ? "true" : "false");
      ibg.applyFunctorOnCCs(GraphInferenceFunctor(const_cast<const Param&>(param_), debug_lvl_));


      LOG_INFO << "Peptide FDR AUC after protein inference: " << pepFDR.rocN(peptideIDs, 0) << std::endl;
      ibg.applyFunctorOnCCsST(AnnotateIndistGroupsFunctor(proteinIDs[0]));

      // rename score_type in PepIDs? I think not. Posterior Probability is still fine. You can
      // get the type from search_engine = Epifany + setting = on.
      // TODO But we could set the "search engine". E.g. Percolator sets itself, too and is just a rescoring.
      //if (update_PSM_probabilities)
      //{}
      //TODO set all unused (= not top) PSMs to 0 or remove! Currently not so bad because FDR also can take just the best.

    }
    else
    {
      //TODO create run info parameter or even a different tool/class.
      ibg.buildGraphWithRunInfo(param_.getValue("top_PSMs"));
      ibg.computeConnectedComponents();
      ibg.clusterIndistProteinsAndPeptidesAndExtendGraph();

      vector<double> gamma_search{0.5};
      vector<double> beta_search{0.001};
      vector<double> alpha_search{0.1, 0.3, 0.5, 0.7, 0.9};

      GridSearch<double,double,double> gs{alpha_search, beta_search, gamma_search};

      std::array<size_t, 3> bestParams{{0, 0, 0}};
      //TODO run grid search on reduced graph?
      //TODO if not, think about storing results temporary (file? mem?) and only keep the best in the end
      gs.evaluate(GridSearchEvaluator(param_, ibg, proteinIDs[0], debug_lvl_), -1.0, bestParams);

      std::cout << "Best params found at " << bestParams[0] << "," << bestParams[1] << "," << bestParams[2] << std::endl;
      double bestGamma = gamma_search[bestParams[2]];
      double bestBeta = beta_search[bestParams[1]];
      double bestAlpha = alpha_search[bestParams[0]];
      std::cout << "Running with best parameters again." << std::endl;
      param_.setValue("model_parameters:prot_prior", bestGamma);
      param_.setValue("model_parameters:pep_emission", bestAlpha);
      param_.setValue("model_parameters:pep_spurious_emission", bestBeta);
      ibg.applyFunctorOnCCs(ExtendedGraphInferenceFunctor(const_cast<const Param&>(param_)));
      ibg.applyFunctorOnCCsST(AnnotateIndistGroupsFunctor(proteinIDs[0]));
    }
  }
}
