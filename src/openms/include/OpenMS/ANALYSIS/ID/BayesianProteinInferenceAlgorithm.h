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
#ifndef OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H
#define OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H

//#define INFERENCE_BENCH

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/GridSearch.h>
#include <vector>
#include <functional>
#include <boost/optional.hpp>

namespace OpenMS
{
  class ConsensusMap;
  class IDBoostGraph;
  class PeptideIdentification;
  class ProteinIdentification;
  class ExperimentalDesign;

  class OPENMS_DLLAPI BayesianProteinInferenceAlgorithm :
      public DefaultParamHandler,
      public ProgressLogger
  {
  public:
    /// Constructor
    explicit BayesianProteinInferenceAlgorithm(unsigned int debug_lvl = 0);

    /// Destructor
    ~BayesianProteinInferenceAlgorithm() override = default;

    void updateMembers_() override;

    //Note: How to perform group inference
    // Three options:
    // - (as implemented) use the automatically created indist. groups and report their posterior
    // - (can be done additionally) collapse proteins to groups beforehand and run inference
    // - (if no single protein scores wanted at all) calculate prior from proteins for the group
    //  beforehand and remove proteins from network (saves computation
    //  because messages are not passed from prots to groups anymore.

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components
    class GraphInferenceFunctor;

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components
    class ExtendedGraphInferenceFunctor;

    /// A function object to pass into the GridSearch;
    struct GridSearchEvaluator;

    /// Perform inference. Writes its results into protein and (optionally) peptide hits (new score).
    /// Optionally adds indistinguishable protein groups with separate scores, too.
    /// Currently only takes first proteinID run.
    /// TODO loop over all runs
    void inferPosteriorProbabilities(
        std::vector<ProteinIdentification>& proteinIDs,
        std::vector<PeptideIdentification>& peptideIDs,
        boost::optional<const ExperimentalDesign&> = boost::optional<const ExperimentalDesign&>());

    void inferPosteriorProbabilities(
        ConsensusMap& cmap,
        boost::optional<const ExperimentalDesign&> = boost::optional<const ExperimentalDesign&>());

  private:

    void inferPosteriorProbabilities_(IDBoostGraph& ibg);

    GridSearch<double,double,double> initGridSearchFromParams_(
        std::vector<double>& alpha_search,
        std::vector<double>& beta_search,
        std::vector<double>& gamma_search
        );

    void setScoreTypeAndSettings_(ProteinIdentification& proteinIDs);

    std::function<void(PeptideIdentification&/*, const String& run_id*/)> checkConvertAndFilterPepHits_;

    /// The grid search object initialized with a default grid
    GridSearch<double,double,double> grid{{0.008,0.032,0.128},{0.001},{0.5}};

    unsigned int debug_lvl_;

    #ifdef INFERENCE_BENCH
    std::vector<std::pair<double,Size>> debug_times_;
    #endif

  };
}
#endif // OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H
