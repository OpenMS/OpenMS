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
#pragma once

//#define INFERENCE_BENCH

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/GridSearch.h>

#include <vector>
#include <functional>
#include <boost/optional.hpp>

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

    /// Perform inference. Filter, build graph, run the private inferPosteriorProbabilities_ function.
    /// Writes its results into protein and (optionally also) peptide hits (as new score).
    /// Optionally adds indistinguishable protein groups with separate scores, too. See Param object of class.
    /// Currently only takes first proteinID run and all peptides.
    /// Experimental design can be used to create an extended graph with replicate information. (experimental)
    /// @TODO loop over all runs
    void inferPosteriorProbabilities(
        std::vector<ProteinIdentification>& proteinIDs,
        std::vector<PeptideIdentification>& peptideIDs,
        boost::optional<const ExperimentalDesign> exp_des = boost::optional<const ExperimentalDesign>());

    /// Perform inference. Filter, build graph, run the private inferPosteriorProbabilities_ function.
    /// Writes its results into protein and (optionally also) peptide hits (as new score).
    /// Optionally adds indistinguishable protein groups with separate scores, too. See Param object of class.
    /// Loops over all runs in the ConsensusMaps' protein IDs. (experimental)
    void inferPosteriorProbabilities(
        ConsensusMap& cmap,
        boost::optional<const ExperimentalDesign> exp_des = boost::optional<const ExperimentalDesign>());

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

    /// function initialized based on the algorithm parameters that is used to filter PeptideHits
    /// @todo extend to allow filtering only for the current run
    std::function<void(PeptideIdentification&/*, const String& run_id*/)> checkConvertAndFilterPepHits_;

    unsigned int debug_lvl_;

    #ifdef INFERENCE_BENCH
    std::vector<std::pair<double,Size>> debug_times_;
    #endif

  };
}
