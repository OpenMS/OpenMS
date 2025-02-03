// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/BayesianProteinInferenceAlgorithm.h>

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/ANALYSIS/ID/IDScoreGetterSetter.h>
#include <OpenMS/ANALYSIS/ID/IDBoostGraph.h>
#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <set>

using namespace std;
using namespace OpenMS::Internal;

namespace OpenMS
{

  /// A functor that specifies what to do on a connected component (IDBoostGraph::FilteredGraph)
  class BayesianProteinInferenceAlgorithm::GraphInferenceFunctor
      //: public std::function<unsigned long(IDBoostGraph::Graph&)>
  {
  public:
    //TODO think about restructuring the passed params (we do not need every param from the BPI class here.
    const Param& param_;
    unsigned int debug_lvl_;
    unsigned long cnt_;

    explicit GraphInferenceFunctor(const Param& param, unsigned int debug_lvl):
        param_(param),
        debug_lvl_(debug_lvl),
        cnt_(0)
    {}

    unsigned long operator() (IDBoostGraph::Graph& fg, unsigned int idx) {
      //TODO do quick brute-force calculation if the cc is really small?

      // this skips CCs with just peps or prots. We only add edges between different types.
      // and if there were no edges, it would not be a CC.
      if (boost::num_vertices(fg) >= 2)
      {
        unsigned long nrEdges = boost::num_edges(fg);

        // avoid critical sections if not needed
        if (debug_lvl_ > 1)
        {
          // we do not need information about file and line so use LOG_INFO instead
          OPENMS_LOG_INFO << "Running cc " << String(idx) << "...\n";
          OPENMS_LOG_INFO << "CC " << String(idx) << " has " << String(nrEdges) << " edges.\n";
        }

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
        evergreen::BetheInferenceGraphBuilder<IDBoostGraph::vertex_t> bigb;

        IDBoostGraph::Graph::vertex_iterator ui, ui_end;
        boost::tie(ui,ui_end) = boost::vertices(fg);

        // Store the IDs of the nodes for which you want the posteriors in the end
        vector<vector<IDBoostGraph::vertex_t>> posteriorVars;

        // direct neighbors are proteins on the "left" side and peptides on the "right" side
        // TODO Can be sped up using directed graph. Needs some restructuring in IDBoostGraph class first tho.
        vector<IDBoostGraph::vertex_t> in{};
        //std::vector<IDBoostGraph::vertex_t> out{};

        //TODO the try section could in theory be slimmed down a little bit. Start at first use of insertDependency maybe.
        // check performance impact.
        try
        {
          for (; ui != ui_end; ++ui)
          {
            IDBoostGraph::Graph::adjacency_iterator nbIt, nbIt_end;
            boost::tie(nbIt, nbIt_end) = boost::adjacent_vertices(*ui, fg);

            in.clear();
            //out.clear(); // we don't need out edges currently

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
          evergreen::InferenceGraph <IDBoostGraph::vertex_t> ig = bigb.to_graph();
          graph_mp_ownership_acquired = true;

          unsigned long maxMessages = param_
              .getValue("loopy_belief_propagation:max_nr_iterations");
          double initDampeningLambda = param_
              .getValue("loopy_belief_propagation:dampening_lambda");
          double initConvergenceThreshold = param_.getValue(
              "loopy_belief_propagation:convergence_threshold");
          std::string scheduler_type = param_.getValue(
              "loopy_belief_propagation:scheduling_type");

          std::unique_ptr<evergreen::Scheduler<IDBoostGraph::vertex_t>> scheduler;
          if (scheduler_type == "priority")
          {
            scheduler = std::unique_ptr<evergreen::Scheduler<IDBoostGraph::vertex_t>>(
                new evergreen::PriorityScheduler<IDBoostGraph::vertex_t>(
                 initDampeningLambda,
                 initConvergenceThreshold,
                 maxMessages));
          }
          else if (scheduler_type == "subtree")
          {
            scheduler = std::unique_ptr<evergreen::Scheduler<IDBoostGraph::vertex_t>>(
                new evergreen::RandomSubtreeScheduler<IDBoostGraph::vertex_t>(
                    initDampeningLambda,
                    initConvergenceThreshold,
                    maxMessages));
          }
          else if (scheduler_type == "fifo")
          {
            scheduler = std::unique_ptr<evergreen::Scheduler<IDBoostGraph::vertex_t>>(
                new evergreen::FIFOScheduler<IDBoostGraph::vertex_t>(
                  initDampeningLambda,
                  initConvergenceThreshold,
                  maxMessages));
          }
          else
          {
            scheduler = std::unique_ptr<evergreen::Scheduler<IDBoostGraph::vertex_t>>(
                new evergreen::PriorityScheduler<IDBoostGraph::vertex_t>(
                    initDampeningLambda,
                    initConvergenceThreshold,
                    maxMessages));
          }
          scheduler->add_ab_initio_edges(ig);

          evergreen::BeliefPropagationInferenceEngine<IDBoostGraph::vertex_t> bpie(*scheduler, ig);

          vector<evergreen::LabeledPMF<IDBoostGraph::vertex_t>> posteriorFactors;
          unsigned long nrEdgesSq = nrEdges*nrEdges;
          if (maxMessages < nrEdgesSq * 3ul)
          {
            posteriorFactors = bpie.estimate_posteriors_in_steps(posteriorVars,
            {
                std::make_tuple(maxMessages, initDampeningLambda, initConvergenceThreshold)});
          }
          else
          {
            posteriorFactors = bpie.estimate_posteriors_in_steps(posteriorVars,
            {
                std::make_tuple(std::max<unsigned long>(10000ul, nrEdgesSq*2ul), initDampeningLambda, initConvergenceThreshold),
                std::make_tuple(nrEdgesSq, std::min(0.49,initDampeningLambda*10), std::min(0.01,initConvergenceThreshold*10)),
                std::make_tuple(nrEdgesSq/2ul, std::min(0.49,initDampeningLambda*100), std::min(0.01,initConvergenceThreshold*100))
            });
          }

          // TODO move the writing of statistics from IDBoostGraph here and write more stats
          //  like nr messages and failure/success
          unsigned long nrMessagesNeeded = bpie.getNrMessagesPassed();

          for (auto const &posteriorFactor : posteriorFactors)
          {
            double posterior = 1.0;
            IDBoostGraph::SetPosteriorVisitor pv;
            IDBoostGraph::vertex_t nodeId = posteriorFactor.ordered_variables()[0];
            const evergreen::PMF &pmf = posteriorFactor.pmf();
            // If Index 0 is in the range of this result PMFFactor its probability is non-zero
            // and the prob of presence is 1-P(p=0). Important in multi-value factors like protein groups.
            if (0 >= pmf.first_support()[0] && 0 <= pmf.last_support()[0])
            {
              posterior = 1. - pmf.table()[0ul];
            }
            auto bound_visitor = std::bind(pv, std::placeholders::_1, posterior);
            boost::apply_visitor(bound_visitor, fg[nodeId]);
          }

          // avoid critical sections if not needed
          if (debug_lvl_ > 1)
          {
            // we do not need information about file and line so use LOG_INFO instead
            OPENMS_LOG_INFO << "Finished cc " << String(idx) << "after " << String(nrMessagesNeeded) << " messages\n";
          }

          //TODO we could write out/save the posteriors here,
          // so we can easily read them later for the best params of the grid search
          return nrMessagesNeeded;
        }
        catch (const std::runtime_error& /*e*/)
        {
          //TODO print failing component and implement the following options
          // 1) Leave posteriors (e.g. if Percolator was ran before. Make sure they are PPs not PEPs)
          // 2) set posteriors to priors (implicitly done right now)
          // 3) try another type of inference on that connected component. Different scheduler,
          //    different extreme probabilities or maybe best: trivial aggregation-based inference.
          // 4) Cancelling this and all other threads/ the loop and call this set of parameters invalid

          //For now we just warn and continue with the rest of the iterations. Might still be a valid run.

          // Graph builder needs to build otherwise it leaks memory.
          if (!graph_mp_ownership_acquired) bigb.to_graph();

          if (debug_lvl_ > 2)
          {
            std::ofstream ofs;
            ofs.open (std::string("failed_cc_a") + 
                param_.getValue("model_parameters:pep_emission").toString() + "_b" + 
                param_.getValue("model_parameters:pep_spurious_emission").toString() + "_g" +
                param_.getValue("model_parameters:prot_prior").toString() + "_c" +
                param_.getValue("model_parameters:pep_prior").toString() + "_p" + String(pnorm) + "_"
                + String(idx) + ".dot"
                , std::ofstream::out);
            IDBoostGraph::printGraph(ofs, fg);
            //TODO print graph with peptide probabilities to see which evidences cause problems with which params
          }
          OPENMS_LOG_WARN << "Warning: Loopy belief propagation encountered a problem in a connected component. Skipping"
                      " inference there." << std::endl;
          return 0;
        }
      }
      else
      {
        OPENMS_LOG_WARN << "Skipped cc with only one type (proteins or peptides)" << std::endl;
        return 0;
      }
    }
  };

  /// A functor that specifies what to do on a connected component with additional layers (i.e. implicitly extended
  /// graph. @TODO static type checking
  class BayesianProteinInferenceAlgorithm::ExtendedGraphInferenceFunctor
     //: public std::function<unsigned long(IDBoostGraph::Graph&)>
  {
  public:
    const Param& param_;

    explicit ExtendedGraphInferenceFunctor(const Param& param):
        param_(param)
    {}

    unsigned long operator() (IDBoostGraph::Graph& fg, unsigned int /*idx*/) {
      //TODO do quick brute-force calculation if the cc is really small
      //TODO make use of idx

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

        evergreen::BetheInferenceGraphBuilder<IDBoostGraph::vertex_t> bigb;

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
            //out.clear(); // we don't need out edges currently

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
          evergreen::InferenceGraph<IDBoostGraph::vertex_t> ig = bigb.to_graph();

          //TODO parametrize the type of scheduler.
          evergreen::PriorityScheduler<IDBoostGraph::vertex_t> scheduler(param_.getValue("loopy_belief_propagation:dampening_lambda"),
                                                     param_.getValue("loopy_belief_propagation:convergence_threshold"),
                                                     param_.getValue("loopy_belief_propagation:max_nr_iterations"));
          scheduler.add_ab_initio_edges(ig);

          evergreen::BeliefPropagationInferenceEngine<IDBoostGraph::vertex_t> bpie(scheduler, ig);

          auto posteriorFactors = bpie.estimate_posteriors(posteriorVars);

          //TODO you could also save the indices of the peptides here and request + update their posteriors, too.
          for (auto const &posteriorFactor : posteriorFactors)
          {
            double posterior = 0.0;
            IDBoostGraph::SetPosteriorVisitor pv;
            unsigned long nodeId = posteriorFactor.ordered_variables()[0];
            const evergreen::PMF &pmf = posteriorFactor.pmf();
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
          std::cout << "Warning: Loopy belief propagation encountered a problem in a connected component. Skipping"
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
    const unsigned int debug_lvl_;

    explicit GridSearchEvaluator(Param& param, IDBoostGraph& ibg, unsigned int debug_lvl):
        param_(param),
        ibg_(ibg),
        debug_lvl_(debug_lvl)
    {}

    double operator() (double alpha, double beta, double gamma)
    {
      OPENMS_LOG_INFO << "Evaluating: " << alpha << " " << beta << " " << gamma << std::endl;
      if (beta - alpha >= 0.3 && alpha + beta <= 1.0)
      {
        OPENMS_LOG_INFO << "Skipping improbable parameter combination.. " << std::endl;
        return 0.;
      }
      param_.setValue("model_parameters:prot_prior", gamma);
      param_.setValue("model_parameters:pep_emission", alpha);
      param_.setValue("model_parameters:pep_spurious_emission", beta);
      GraphInferenceFunctor gif {param_, debug_lvl_};
      ibg_.applyFunctorOnCCs(gif);

      FalseDiscoveryRate fdr;
      Param fdrparam = fdr.getParameters();
      fdrparam.setValue("conservative",param_.getValue("param_optimize:conservative_fdr"));
      fdrparam.setValue("add_decoy_proteins","true");
      fdr.setParameters(fdrparam);

      if (param_.getValue("annotate_group_probabilities").toBool())
      {
        ScoreToTgtDecLabelPairs scores_and_tgt_fraction{};
        ibg_.getProteinGroupScoresAndTgtFraction(scores_and_tgt_fraction);
        return fdr.applyEvaluateProteinIDs(scores_and_tgt_fraction, 1.0, 100, static_cast<double>(param_.getValue("param_optimize:aucweight")));
      }
      // The following line assumes that ALL proteins in the ID structure are used in the graph
      return fdr.applyEvaluateProteinIDs(ibg_.getProteinIDs(), 1.0, 100, static_cast<double>(param_.getValue("param_optimize:aucweight")));
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

    defaults_.setValue("psm_probability_cutoff",
                          0.001,
                          "Remove PSMs with probabilities less than this cutoff");
    defaults_.setMinFloat("psm_probability_cutoff", 0.0);
    defaults_.setMaxFloat("psm_probability_cutoff", 1.0);

    defaults_.setValue("top_PSMs",
                       1,
                       "Consider only top X PSMs per spectrum. 0 considers all.");
    defaults_.setMinInt("top_PSMs", 0);

    defaults_.setValue("keep_best_PSM_only",
                       "true",
                       "Epifany uses the best PSM per peptide for inference. Discard the rest (true) or keep"
                       "e.g. for quantification/reporting?");
    defaults_.setValidStrings("keep_best_PSM_only", {"true","false"});

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

    defaults_.setValue("use_ids_outside_features",
                       "false",
                       "(Only consensusXML) Also use IDs without associated features for inference?");
    defaults_.setValidStrings("use_ids_outside_features", {"true","false"});

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

    defaults_.setValue("model_parameters:extended_model",
                       "false",
                       "Uses information from different peptidoforms also across runs"
                       " (automatically activated if an experimental design is given!)");
    defaults_.setValidStrings("model_parameters:extended_model", {"true","false"});

    defaults_.addSection("loopy_belief_propagation","Settings for the loopy belief propagation algorithm.");

    defaults_.setValue("loopy_belief_propagation:scheduling_type",
                       "priority",
                       "(Not used yet) How to pick the next message:"
                           " priority = based on difference to last message (higher = more important)."
                           " fifo = first in first out."
                           " subtree = message passing follows a random spanning tree in each iteration");
    defaults_.setValidStrings("loopy_belief_propagation:scheduling_type", {"priority","fifo","subtree"});

    //TODO not yet implemented
/*    defaults_.setValue("loopy_belief_propagation:message_difference",
                       "MSE",
                       "How to calculate the difference of distributions in updated messages.");
    defaults_.setValidStrings("loopy_belief_propagation:message_difference", {"MSE"});*/
    defaults_.setValue("loopy_belief_propagation:convergence_threshold",
                       1e-5,
                       "Initial threshold under which MSE difference a message is considered to be converged.");
    defaults_.setMinFloat("loopy_belief_propagation:convergence_threshold", 1e-9);
    defaults_.setMaxFloat("loopy_belief_propagation:convergence_threshold", 1.0);

    defaults_.setValue("loopy_belief_propagation:dampening_lambda",
                       1e-3,
                       "Initial value for how strongly should messages be updated in each step. "
                           "0 = new message overwrites old completely (no dampening; only recommended for trees),"
                           "0.5 = equal contribution of old and new message (stay below that),"
                           "In-between it will be a convex combination of both. Prevents oscillations but hinders convergence.");
    defaults_.setMinFloat("loopy_belief_propagation:dampening_lambda", 0.0);
    defaults_.setMaxFloat("loopy_belief_propagation:dampening_lambda", 0.49999);

    defaults_.setValue("loopy_belief_propagation:max_nr_iterations",
                       (1ul<<31)-1,
                       "(Usually auto-determined by estimated but you can set a hard limit here)."
                       " If not all messages converge, how many iterations should be done at max per connected component?");
    //I think restricting does not work because it only works for type Int (= int), not unsigned long
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
                       0.3,
                       "How important is target decoy AUC vs calibration of the posteriors?"
                       " 0 = maximize calibration only,"
                       " 1 = maximize AUC only,"
                       " between = convex combination.");
    defaults_.setMinFloat("param_optimize:aucweight", 0.0);
    defaults_.setMaxFloat("param_optimize:aucweight", 1.0);

    defaults_.setValue("param_optimize:conservative_fdr",
                          "true",
                          "Use (D+1)/(T) instead of (D+1)/(T+D) for parameter estimation.");
    defaults_.setValidStrings("param_optimize:conservative_fdr", {"true","false"});

    defaults_.setValue("param_optimize:regularized_fdr",
                       "true",
                       "Use a regularized FDR for proteins without unique peptides.");
    defaults_.setValidStrings("param_optimize:regularized_fdr", {"true","false"});


    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
   }

  void BayesianProteinInferenceAlgorithm::updateMembers_()
  {
    //Note: the following lambda function can be changed, e.g. when we want to do a extremum removal etc. beforehand
    /*
    double min_nonnull_obs_probability = getDoubleOption_("min_psms_extreme_probability");
    double max_nonone_obs_probability = getDoubleOption_("max_psms_extreme_probability");
    // Currently unused
    bool datadependent_extrema_removal = false;
    if (datadependent_extrema_removal)
    {
      pair<double,double> minmax = checkExtremePSMScores_(mergedpeps);
      min_nonnull_obs_probability = minmax.first;
      max_nonone_obs_probability = minmax.second;
    }

    if (min_nonnull_obs_probability > 0.0 || max_nonone_obs_probability < 1.0 )
    {
      removeExtremeValues_(mergedpeps, min_nonnull_obs_probability, max_nonone_obs_probability);
    }
    */
    //TODO also convert potential PEPs to PPs in ProteinHits? In case you want to use them as priors or
    // emergency posteriors?
    //TODO test performance of getting the probability cutoff every time vs capture free lambda
    double probability_cutoff = param_.getValue("psm_probability_cutoff");
    checkConvertAndFilterPepHits_ = [probability_cutoff](PeptideIdentification &pep_id/*, const String& run_id*/)
    {
      //if (pep_id.getIdentifier() == run_id)
      //{
      String score_l = pep_id.getScoreType();
      score_l = score_l.toLower();
      if (score_l == "pep" || score_l == "posterior error probability" || score_l == "ms:1001493")
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
        if (score_l != "posterior probability")
        {
          throw OpenMS::Exception::InvalidParameter(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Epifany needs Posterior (Error) Probabilities in the Peptide Hits. Use Percolator with PEP score"
              " or run IDPosteriorErrorProbability first.");
        }
      }
      //TODO remove hits "on-the-go"?
      IDFilter::removeMatchingItems(pep_id.getHits(),
                                    [&probability_cutoff](PeptideHit &hit)
                                    { return hit.getScore() < probability_cutoff; });
      //}
    };
  }

  void BayesianProteinInferenceAlgorithm::setScoreTypeAndSettings_(ProteinIdentification& proteinIDs)
  {
    proteinIDs.setScoreType("Posterior Probability");
    proteinIDs.setInferenceEngine("Epifany");
    proteinIDs.setInferenceEngineVersion(VersionInfo::getVersion());
    proteinIDs.setHigherScoreBetter(true);
  }

  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(
      ConsensusMap& cmap,
      bool greedy_group_resolution, // TODO probably better to add it as a Param
      std::optional<const ExperimentalDesign> exp_des)
  {
    IDScoreSwitcherAlgorithm switcher;
    Size counter(0);
    try
    {
      switcher.switchToGeneralScoreType(cmap, IDScoreSwitcherAlgorithm::ScoreType::PEP, counter);
    }
    catch (OpenMS::Exception::MissingInformation& /*e*/)
    {
      throw OpenMS::Exception::MissingInformation(
          __FILE__,
          __LINE__,
          OPENMS_PRETTY_FUNCTION,
          "Epifany needs Posterior Error Probabilities in the Peptide Hits. Use Percolator with PEP score"
          " or run IDPosteriorErrorProbability first.");
    }

    //TODO filtering needs to account for run info if we allow running on a subset.
    cmap.applyFunctionOnPeptideIDs(checkConvertAndFilterPepHits_);
    //TODO BIG filter empty PeptideIDs afterwards

    bool keep_all_psms = param_.getValue("keep_best_PSM_only").toString() == "false";
    bool user_defined_priors = param_.getValue("user_defined_priors").toBool();
    bool use_unannotated_ids = param_.getValue("use_ids_outside_features").toBool();
    bool use_run_info = param_.getValue("model_parameters:extended_model").toBool();
    Size nr_top_psms = static_cast<Size>(param_.getValue("top_PSMs"));

    FalseDiscoveryRate pepFDR;
    Param p = pepFDR.getParameters();

    // I think it is best to always use the best PSM only for comparing PSM FDR before-after
    // since inference might change the ranking.
    p.setValue("use_all_hits", "false");
    pepFDR.setParameters(p);

    //TODO allow use unassigned everywhere
    //TODO actually if we just want to use replicate information, we can still filter for best per run,
    // but the extended model is currently coupled to multiple charge and mod states (which would be removed)
    if (!use_run_info)
    {
      if (!keep_all_psms)
      {
        IDFilter::keepBestPerPeptidePerRun(cmap, true,
                                           true, static_cast<unsigned int>(nr_top_psms));
        IDFilter::removeEmptyIdentifications(cmap);
      }
      else
      {
        IDFilter::annotateBestPerPeptidePerRun(cmap, true,
                                               true, static_cast<unsigned int>(nr_top_psms));
      }
    }
    // filter in both cases since the PSM filtering by score is always active
    IDFilter::removeUnreferencedProteins(cmap, true);

    // extract proteins that are "theoretically" unreferenced, since
    // unassigned PSMs might not be considered in inference (depending on param).
    // NOTE: this would in theory not be necessary if we calculate the FDR based on
    // the graph only (because then only the used proteins are in the graph).
    // But FDR on the ProteinID data structure should be faster.
    std::map<String, vector<ProteinHit>> unassigned{};
    if (!use_unannotated_ids)
    {
      unassigned = IDFilter::extractUnassignedProteins(cmap);
    }

    vector<ProteinIdentification>& proteinIDs = cmap.getProteinIdentifications();
    if (proteinIDs.size() == 1) // could be merged with the general case, but we can save the runid lookup here.
    {
      resetProteinScores_(proteinIDs[0], user_defined_priors);

      // TODO try to calc AUC partial only (e.g. up to 5% FDR)
      if (!keep_all_psms)
        OPENMS_LOG_INFO << "Peptide FDR AUC before protein inference: " << pepFDR.rocN(cmap, 0) << std::endl;

      setScoreTypeAndSettings_(proteinIDs[0]);
      IDBoostGraph ibg(proteinIDs[0], cmap, nr_top_psms, use_run_info, use_unannotated_ids, keep_all_psms, exp_des);
      inferPosteriorProbabilities_(ibg);
      if (greedy_group_resolution) ibg.resolveGraphPeptideCentric(true);

      if (!keep_all_psms)
        OPENMS_LOG_INFO << "Peptide FDR AUC after protein inference: " << pepFDR.rocN(cmap, 0) << std::endl;

      if (!use_unannotated_ids)
      {
        auto& unassigned_for_run = unassigned.at(proteinIDs[0].getIdentifier());
        for (auto& h : unassigned_for_run) h.setScore(0.);
        proteinIDs[0].getHits().reserve(proteinIDs[0].getHits().size() + unassigned_for_run.size());
        std::move(std::begin(unassigned_for_run), std::end(unassigned_for_run), std::back_inserter(proteinIDs[0].getHits()));
        unassigned_for_run.clear();
      }

      proteinIDs[0].fillIndistinguishableGroupsWithSingletons();
    }
    else if (cmap.getProteinIdentifications().size() > 1)
    {
      for (auto& proteinID : cmap.getProteinIdentifications())
      {
        //TODO use if case for unassigned
        //TODO Log the currently processed run
        resetProteinScores_(proteinID, user_defined_priors);

        //TODO try to calc AUC partial only (e.g. up to 5% FDR)
        if (!keep_all_psms)
          OPENMS_LOG_INFO << "Peptide FDR AUC before protein inference: " << pepFDR.rocN(cmap, 0, proteinID.getIdentifier()) << std::endl;

        setScoreTypeAndSettings_(proteinID);
        IDBoostGraph ibg(proteinID, cmap, nr_top_psms, use_run_info, use_unannotated_ids, keep_all_psms);
        inferPosteriorProbabilities_(ibg);
        if (greedy_group_resolution) ibg.resolveGraphPeptideCentric(true);

        if (!keep_all_psms)
          OPENMS_LOG_INFO << "Peptide FDR AUC after protein inference: " << pepFDR.rocN(cmap, 0, proteinID.getIdentifier()) << std::endl;

        if (!use_unannotated_ids)
        {
          auto& unassigned_for_run = unassigned.at(proteinIDs[0].getIdentifier());
          for (auto& h : unassigned_for_run) h.setScore(0.);
          proteinID.getHits().reserve(proteinID.getHits().size() + unassigned_for_run.size());
          std::move(std::begin(unassigned_for_run), std::end(unassigned_for_run), std::back_inserter(proteinID.getHits()));
          unassigned_for_run.clear();
        }

        proteinID.fillIndistinguishableGroupsWithSingletons();
      }
    }
  }

  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities_(
      IDBoostGraph& ibg)
  {
    bool use_run_info = param_.getValue("model_parameters:extended_model").toBool();

    ibg.computeConnectedComponents();
    ibg.clusterIndistProteinsAndPeptides();

    vector<double> gamma_search;
    vector<double> beta_search;
    vector<double> alpha_search;
    GridSearch<double,double,double> gs = initGridSearchFromParams_(alpha_search, beta_search, gamma_search);

    std::array<size_t, 3> bestParams{{0, 0, 0}};

    //Save initial settings and deactivate certain features to save time during grid search and to not
    // interfere with later runs.
    // TODO We could think about optimizing PSM FDR as another goal though.
    bool update_PSM_probabilities = param_.getValue("update_PSM_probabilities").toBool();
    param_.setValue("update_PSM_probabilities","false");

    bool annotate_group_posteriors = param_.getValue("annotate_group_probabilities").toBool();
    // during grid search we evaluate on single protein-level
    param_.setValue("annotate_group_probabilities","false");

    //TODO run grid search on reduced graph? Then make sure, untouched protein/peps do not affect evaluation results.
    //TODO if not, think about storing results temporary (file? mem?) and only keep the best in the end
    //TODO think about running grid search on the small CCs only (maybe it's enough)
    if (gs.getNrCombos() > 1)
    {
     OPENMS_LOG_INFO << "Testing " << gs.getNrCombos() << " param combinations." << std::endl;
      /*double res =*/ gs.evaluate(GridSearchEvaluator(param_, ibg, debug_lvl_), -1.0, bestParams);
    }
    else
    {
     OPENMS_LOG_INFO << "Only one combination specified: Skipping grid search." << std::endl;
    }

    double bestGamma = gamma_search[bestParams[2]];
    double bestBeta = beta_search[bestParams[1]];
    double bestAlpha = alpha_search[bestParams[0]];
    OPENMS_LOG_INFO << "Best params found at a=" << bestAlpha << ", b=" << bestBeta << ", g=" << bestGamma << std::endl;
    OPENMS_LOG_INFO << "Running with best parameters:" << std::endl;
    param_.setValue("model_parameters:prot_prior", bestGamma);
    param_.setValue("model_parameters:pep_emission", bestAlpha);
    param_.setValue("model_parameters:pep_spurious_emission", bestBeta);
    // Reset original values for those two options
    param_.setValue("update_PSM_probabilities", update_PSM_probabilities ? "true" : "false");
    param_.setValue("annotate_group_probabilities", annotate_group_posteriors ? "true" : "false");

    if (!use_run_info)
    {
      GraphInferenceFunctor gif {param_, debug_lvl_};
      ibg.applyFunctorOnCCs(gif);
    }
    else
    {
      //TODO under construction
      ExtendedGraphInferenceFunctor gif {param_};
      ibg.applyFunctorOnCCs(gif);
    }

    //uses the existing protein group nodes in the graph
    ibg.annotateIndistProteins(true);
  }

  GridSearch<double,double,double> BayesianProteinInferenceAlgorithm::initGridSearchFromParams_(
      vector<double>& alpha_search,
      vector<double>& beta_search,
      vector<double>& gamma_search
      )
  {
    // Do not expand gamma_search when user_defined_priors is on. Would be unused.
    double alpha = param_.getValue("model_parameters:pep_emission");
    double beta = param_.getValue("model_parameters:pep_spurious_emission");
    double gamma = param_.getValue("model_parameters:prot_prior");

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
      beta_search = {0.01, 0.2, 0.4};
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

    return GridSearch<double,double,double>{alpha_search, beta_search, gamma_search};
  }

  void BayesianProteinInferenceAlgorithm::inferPosteriorProbabilities(
      std::vector<ProteinIdentification>& proteinIDs,
      std::vector<PeptideIdentification>& peptideIDs,
      bool greedy_group_resolution,
      std::optional<const ExperimentalDesign> exp_des)
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

    //TODO actually loop over all proteinID runs.
    if (proteinIDs.size() > 1)
    {
      OPENMS_LOG_WARN << "Warning: more than one protein identification run provided for inference. Only "
                         "the first will be processed for now." << std::endl;
    }

    // groups will be reannotated
    proteinIDs[0].getIndistinguishableProteins().clear();

    bool use_run_info = param_.getValue("model_parameters:extended_model").toBool();

    //TODO BIG filtering needs to account for run info if only a subset is to be processed!
    std::for_each(peptideIDs.begin(), peptideIDs.end(), checkConvertAndFilterPepHits_);
    IDFilter::removeEmptyIdentifications(peptideIDs);

    Size nr_top_psms = static_cast<Size>(param_.getValue("top_PSMs"));
    bool keep_all_psms = param_.getValue("keep_best_PSM_only").toString() == "false";

    //TODO actually if we just want to use replicate information, we can still filter for best per run,
    // but the extended model is currently coupled to multiple charge and mod states (which would be removed)
    if (!use_run_info)
    {
      if (!keep_all_psms)
      {
        IDFilter::keepBestPerPeptidePerRun(proteinIDs, peptideIDs, true,
                                           true, static_cast<unsigned int>(nr_top_psms));
        IDFilter::removeEmptyIdentifications(peptideIDs);
      }
      else
      {
        IDFilter::annotateBestPerPeptidePerRun(proteinIDs, peptideIDs, true,
                                           true, static_cast<unsigned int>(nr_top_psms));
      }
    }

    IDFilter::removeUnreferencedProteins(proteinIDs, peptideIDs);

    FalseDiscoveryRate pepFDR;
    Param p = pepFDR.getParameters();

    // I think it is best to always use the best PSM only for comparing PSM FDR before-after
    // since inference might change the ranking.
    p.setValue("use_all_hits", "false");
    pepFDR.setParameters(p);

    bool user_defined_priors = param_.getValue("user_defined_priors").toBool();
    resetProteinScores_(proteinIDs[0], user_defined_priors);

    if (!keep_all_psms)
      OPENMS_LOG_INFO << "Peptide FDR AUC before protein inference: " << pepFDR.rocN(peptideIDs, 0, proteinIDs[0].getIdentifier()) << std::endl;

    setScoreTypeAndSettings_(proteinIDs[0]);
    IDBoostGraph ibg(proteinIDs[0], peptideIDs, nr_top_psms, use_run_info, keep_all_psms, exp_des);
    inferPosteriorProbabilities_(ibg);
    if (greedy_group_resolution) ibg.resolveGraphPeptideCentric(true);
    proteinIDs[0].fillIndistinguishableGroupsWithSingletons();

    if (!keep_all_psms)
      OPENMS_LOG_INFO << "Peptide FDR AUC after protein inference: " << pepFDR.rocN(peptideIDs, 0, proteinIDs[0].getIdentifier()) << std::endl;
  }

  void BayesianProteinInferenceAlgorithm::resetProteinScores_(ProteinIdentification& protein_id, bool keep_old_as_prior)
  {
    // Save current scores as priors if requested
    if (keep_old_as_prior)
    {
      // Save current protein score into a metaValue
      for (auto& prot_hit : protein_id.getHits())
      {
        prot_hit.setMetaValue("Prior", prot_hit.getScore());
        prot_hit.setScore(0.);
      }
    }
    else
    {
      // Save current protein score into a metaValue
      for (auto& prot_hit : protein_id.getHits())
      {
        prot_hit.setScore(0.);
      }
    }
  }

}
