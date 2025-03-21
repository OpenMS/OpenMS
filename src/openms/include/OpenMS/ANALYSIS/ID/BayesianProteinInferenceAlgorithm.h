// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
#pragma once

//#define INFERENCE_BENCH

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ML/GRIDSEARCH/GridSearch.h>

#include <vector>
#include <functional>
#include <optional>

namespace OpenMS
{
  class ConsensusMap;
  namespace Internal
  {
    class IDBoostGraph;
  }
  class PeptideIdentification;
  class ProteinIdentification;

  /**
   * @brief Performs a Bayesian protein inference on Protein/Peptide identifications or ConsensusMap (experimental).
   * - Filters for best n PSMs per spectrum.
   * - Calculates and filters for best peptide per spectrum.
   * - Builds a k-partite graph from the structures.
   * - Finds and splits into connected components by DFS
   * - Extends the graph by adding layers from indist. protein groups, peptides with the same parents and optionally
   *   some additional layers (peptide sequence, charge, replicate -> extended model = experimental)
   * - Builds a factor graph representation of a Bayesian network using the Evergreen library
   *   See model param section. It is based on the Fido noisy-OR model with an option for
   *   regularizing the number of proteins per peptide.
   * - Performs loopy belief propagation on the graph and queries protein, protein group and/or peptide posteriors
   *   See loopy_belief_propagation param section.
   * - Learns best parameters via grid search if the parameters were not given in the param section.
   * - Writes posteriors to peptides and/or proteins and adds indistinguishable protein groups to the underlying
   *   data structures.
   * - Can make use of OpenMP to parallelize over connected components.
   */
  class OPENMS_DLLAPI BayesianProteinInferenceAlgorithm :
      public DefaultParamHandler,
      public ProgressLogger
  {
  public:
    /// Constructor @todo is there a better way to pass the debug level from TOPPBase?
    explicit BayesianProteinInferenceAlgorithm(unsigned int debug_lvl = 0);

    /// Destructor
    ~BayesianProteinInferenceAlgorithm() override = default;

    void updateMembers_() override;

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components
    class GraphInferenceFunctor;

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components. This can make use of additional layers. @todo use static type checking
    /// by using two different Graph types
    class ExtendedGraphInferenceFunctor;

    /// A function object to pass into the GridSearch class
    struct GridSearchEvaluator;


    /**
     * @brief Perform inference. Filter, build graph, run the private inferPosteriorProbabilities_ function.
     *  Writes its results into protein and (optionally also) peptide hits (as new score).
     *  Optionally adds indistinguishable protein groups with separate scores, too.
     *  Output scores are always posterior probabilities. Input can be posterior or error probabilities.
     *  See Param object defaults_ within the BayesianProteinInferenceAlgorithm for more settings.
     *  Currently only takes first proteinID run and all peptides (irrespective of getIdentifier()).
     * @param proteinIDs Input/output proteins
     * @param peptideIDs Input/output peptides
     * @param greedy_group_resolution Do greedy group resolution? Remove all but best association for "razor" peptides.
     * @param exp_des Experimental design can be used to create an extended graph with replicate information. (experimental)
     * 
     * @todo loop over all runs
     * 
     */
    void inferPosteriorProbabilities(
        std::vector<ProteinIdentification>& proteinIDs,
        std::vector<PeptideIdentification>& peptideIDs,
        bool greedy_group_resolution,
        std::optional<const ExperimentalDesign> exp_des = std::optional<const ExperimentalDesign>());

    /**
     * @brief Perform inference. Filter, build graph, run the private inferPosteriorProbabilities_ function.
     *  Writes its results into protein and (optionally also) peptide hits (as new score).
     *  Optionally adds indistinguishable protein groups with separate scores, too.
     *  Output scores are always posterior probabilities. Input can be posterior or error probabilities.
     *  See Param object defaults_ within the BayesianProteinInferenceAlgorithm for more settings.
     *  Currently only takes first proteinID run and all peptides (irrespective of getIdentifier()).
     * @param cmap Features with input/output peptides and proteins (from getProteinIdentifications)
     * @param greedy_group_resolution Do greedy group resolution? Remove all but best association for "razor" peptides.
     * @param exp_des Experimental design can be used to create an extended graph with replicate information. (experimental)
     */
    void inferPosteriorProbabilities(
        ConsensusMap& cmap,
        bool greedy_group_resolution,
        std::optional<const ExperimentalDesign> exp_des = std::optional<const ExperimentalDesign>());

  private:

    /// after a graph was built, use this method to perform inference and write results to the structures
    /// with which the graph was built
    void inferPosteriorProbabilities_(Internal::IDBoostGraph& ibg);

    /// read Param object and set the grid
    GridSearch<double,double,double> initGridSearchFromParams_(
        std::vector<double>& alpha_search,
        std::vector<double>& beta_search,
        std::vector<double>& gamma_search
        );

    /// set score type and settings for every ProteinID run processed
    void setScoreTypeAndSettings_(ProteinIdentification& proteinIDs);

    /// reset all protein scores to 0.0, save old ones as Prior MetaValue if requested
    // TODO double-check if -1 is maybe the better option
    //  to distinguish between "untouched/unused/unreferenced" (e.g. if somehow
    //  not removed/filtered) and an inferred probability of 0.0. But it might give
    //  problems in FDR algorithms if not ignored/removed correctly
    void resetProteinScores_(ProteinIdentification& protein_id, bool keep_old_as_prior);

    /// function initialized based on the algorithm parameters that is used to filter PeptideHits
    /// @todo extend to allow filtering only for the current run
    std::function<void(PeptideIdentification&/*, const String& run_id*/)> checkConvertAndFilterPepHits_;

    unsigned int debug_lvl_;

    #ifdef INFERENCE_BENCH
    std::vector<std::pair<double,Size>> debug_times_;
    #endif

  };
}
