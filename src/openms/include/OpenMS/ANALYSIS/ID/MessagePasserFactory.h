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

#include <cmath>
#include <Evergreen/evergreen.hpp>

typedef unsigned long int uiint;

namespace OpenMS
{
  namespace Internal
{
  /// Produces MessagePassers (nodes in a factor graph = bayesian network) for use with Evergreen library,
  /// based on a parameterization of the Protein-Peptide Bayesian network.
  /// Those MessagePassers can be tables or convolution trees. Labels are used to associate the variables they are
  /// working on. They can be integers (for speed) or strings (for readability/debug)
  template <typename Label>
  class MessagePasserFactory {
  private:
    //const int minInputsPAF = 3; //@todo could be used to decide when brute force is better.

    /// the model parameters
    double alpha_, beta_, gamma_, p_, pepPrior_;

    /// Likelihoods for the charge states given presence of the peptide sequence (@todo could be calculated from IDPEP
    /// if we do per charge state fitting) or empirically estimated from the input PSMs
    std::map<int, double> chgLLhoods = {{1, 0.7}, {2, 0.9}, {3, 0.7}, {4, 0.5}, {5, 0.5}};

    /// to fill the noisy-OR table for a peptide given parent proteins
    /// TODO precompute for like a hundred parent proteins
    /// TODO introduce special case for alpha or beta = 1. The log formula does not work otherwise.
    inline double notConditionalGivenSum(unsigned long summ) {
      // use log for better precision
      return std::pow(2., log2(1. - beta_) + summ * log2(1. - alpha_));
      //return std::pow((1.0 - alpha_), summ) * (1.0 - beta_); // standard way
    }

  public:
    /// Protein Factor initialized with model prior (missing peps are experimental)
    evergreen::TableDependency<Label> createProteinFactor(Label id, int nrMissingPeps = 0);
    /// Protein Factor initialized with user prior (missing peps are experimental)
    evergreen::TableDependency<Label> createProteinFactor(Label id, double prior, int nrMissingPeps = 0);

    /// Peptide Factor initialized with:
    /// @param prob peptide evidence probability
    evergreen::TableDependency<Label> createPeptideEvidenceFactor(Label id, double prob);

    /// Conditional probability table of peptide given number of parent proteins, based on model params.
    /// Additionally regularizes on the amount of parent proteins.
    /// @param nrParents (maximum) number of parent proteins
    evergreen::TableDependency<Label> createRegularizingSumEvidenceFactor(size_t nrParents, Label nId, Label pepId);

    /// Conditional probability table of peptide given number of parent proteins, based on model params.
    /// @param nrParents (maximum) number of parent proteins
    evergreen::TableDependency<Label> createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId);

    //For extended model. @todo currently unused
    evergreen::TableDependency<Label> createSumFactor(size_t nrParents, Label nId);
    evergreen::TableDependency<Label> createReplicateFactor(Label seqId, Label repId);
    evergreen::TableDependency<Label> createChargeFactor(Label repId, Label chargeId, int chg);

    /// To sum up distributions for the number of parent proteins of a peptide with convolution trees
    evergreen::AdditiveDependency<Label> createPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId);
    /// To sum up distributions for the number of parent proteins of a peptide with convolution trees
    evergreen::AdditiveDependency<Label> createPeptideProbabilisticAdderFactor(const std::vector<Label> & parentProteinIDs, Label nId);
    /// To sum up distributions for the number of parent proteins of a peptide brute-force
    evergreen::PseudoAdditiveDependency<Label> createBFPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId, const std::vector<evergreen::TableDependency <Label> > & deps);

    /**
     * @brief Constructor
     * @param alpha Peptide emission probability
     * @param beta Spurious peptide emission probability
     * @param gamma Protein prior
     * @param p Marginalization norm
     * @param pepPrior Peptide prior (defines at which evidence probability, additional evidence is beneficial)
     */
    MessagePasserFactory<Label>(double alpha, double beta, double gamma, double p, double pep_prior);



    /// Works on a vector of protein indices (potentially not consecutive)
    // TODO we could recollect the protIDs from the union of parents.
    void fillVectorsOfMessagePassers(const std::vector<Label> & protIDs,
                                     const std::vector<std::vector<Label>> & parentsOfPeps,
                                     const std::vector<double> & pepEvidences,
                                     evergreen::InferenceGraphBuilder<Label> & igb);

    //void fillVectorsOfMessagePassersBruteForce(const std::vector<Label> & protIDs,
    //                                 const std::vector<std::vector<Label>> & parentsOfPeps,
    //                                 const std::vector<double> & pepEvidences,
    //                                 evergreen::InferenceGraphBuilder<Label> & igb);
  };

  //IMPLEMENTATIONS:

  template <typename L>
  MessagePasserFactory<L>::MessagePasserFactory(double alpha, double beta, double gamma, double p, double pep_prior) {
    assert(0. <= alpha && alpha <= 1.);
    assert(0. <= beta && beta <= 1.);
    assert(0. <= gamma && gamma <= 1.);
    //Note: smaller than 1 might be possible but is untested right now.
    assert(p >= 1.);
    assert(0. < pep_prior && pep_prior < 1.);
    alpha_ = alpha;
    beta_ = beta;
    gamma_ = gamma;
    p_ = p;
    pepPrior_ = pep_prior;
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createProteinFactor(L id, int nrMissingPeps) {
    double prior = gamma_;
    if (nrMissingPeps > 0)
    {
      double powFactor = std::pow(1.0 - alpha_, -nrMissingPeps);
      prior = -prior/(prior * powFactor - prior - powFactor);
    }
    double table[] = {1.0 - prior, prior};
    evergreen::LabeledPMF<L> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createProteinFactor(L id, double prior, int nrMissingPeps) {
    if (nrMissingPeps > 0)
    {
      double powFactor = std::pow(1.0 - alpha_, -nrMissingPeps);
      prior = -prior/(prior * powFactor - prior - powFactor);
    }
    double table[] = {1.0 - prior, prior};
    evergreen::LabeledPMF<L> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createPeptideEvidenceFactor(L id, double prob) {
    double table[] = {(1 - prob) * (1 - pepPrior_), prob * pepPrior_};
    evergreen::LabeledPMF<L> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<L>(lpmf,p_);
  }


  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createSumEvidenceFactor(size_t nrParents, L nId, L pepId) {
    evergreen::Tensor<double> table({static_cast<unsigned long>(nrParents + 1) , 2});
    for (unsigned long i=0; i <= nrParents; ++i) {
      double notConditional = notConditionalGivenSum(i);
      unsigned long indexArr[2] = {i,0ul};
      table[indexArr] = notConditional;
      unsigned long indexArr2[2] = {i,1ul};
      table[indexArr2] = 1.0 - notConditional;
    }
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<L> lpmf({nId, pepId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createRegularizingSumEvidenceFactor(size_t nrParents, L nId, L pepId) {
    evergreen::Tensor<double> table({static_cast<unsigned long>(nrParents + 1) , 2});
    unsigned long z[2]{0ul,0ul};
    unsigned long z1[2]{0ul,1ul};
    table[z] = 1. - beta_;
    table[z1] = beta_;
    for (unsigned long i=1; i <= nrParents; ++i) {
      double notConditional = notConditionalGivenSum(i);
      unsigned long indexArr[2] = {i,0ul};
      table[indexArr] = notConditional / i;
      unsigned long indexArr2[2] = {i,1ul};
      table[indexArr2] = (1.0 - notConditional) / i;
    }
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<L> lpmf({nId, pepId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createSumFactor(size_t nrParents, L nId) {
    evergreen::Tensor<double> table({nrParents+1});
    for (unsigned long i=0; i <= nrParents; ++i) {
      table[i] = 1.0/(nrParents+1.);
    }
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<L> lpmf({nId}, evergreen::PMF({0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createReplicateFactor(L seqId, L repId) {
    using arr = unsigned long[2];
    evergreen::Tensor<double> table({2,2});
    table[arr{0,0}] = 0.999;
    table[arr{0,1}] = 0.001;
    table[arr{1,0}] = 0.1;
    table[arr{1,1}] = 0.9;
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<L> lpmf({seqId,repId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::TableDependency<L> MessagePasserFactory<L>::createChargeFactor(L repId, L chgId, int chg) {
    double chgPrior = chgLLhoods[chg];
    using arr = unsigned long[2];
    evergreen::Tensor<double> table({2,2});
    table[arr{0,0}] = 0.999;
    table[arr{0,1}] = 0.001;
    table[arr{1,0}] = 0.1;
    table[arr{1,1}] = chgPrior;
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<L> lpmf({repId,chgId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<L>(lpmf,p_);
  }

  template <typename L>
  evergreen::AdditiveDependency<L> MessagePasserFactory<L>::createPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId) {
    std::vector<std::vector<L>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
    return evergreen::AdditiveDependency<L>(parents, {nId}, p_);
  }

  template <typename L>
  evergreen::AdditiveDependency<L> MessagePasserFactory<L>::createPeptideProbabilisticAdderFactor(const std::vector<L> & parentProteinIDs, L nId) {
    std::vector<std::vector<L>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
    return evergreen::AdditiveDependency<L>(parents, {nId}, p_);
  }

  template <typename L>
  evergreen::PseudoAdditiveDependency<L> MessagePasserFactory<L>::createBFPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId, const std::vector<evergreen::TableDependency<L>> & deps) {
    std::vector<std::vector<L>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
    return evergreen::PseudoAdditiveDependency<L>(parents, {nId}, deps, p_);
  }

  /// Works on a vector of protein indices (potentially not consecutive)
  // TODO we could recollect the protIDs from the union of parents.
  template <typename L>
  void MessagePasserFactory<L>::fillVectorsOfMessagePassers(const std::vector<L> & protIDs,
                                                            const std::vector<std::vector<L>> & parentsOfPeps,
                                                            const std::vector<double> & pepEvidences,
                                                            evergreen::InferenceGraphBuilder<L> & igb)
  {
    //TODO asserts could be loosened
    assert(parentsOfPeps.size() == pepEvidences.size());
    for (const std::vector<uiint>& parents : parentsOfPeps)
      for (L parent : parents)
        assert(std::find(protIDs.begin(), protIDs.end(), parent) != protIDs.end());

    for (uiint pid : protIDs)
      igb.insert_dependency(createProteinFactor(pid));

    for (uiint j = 0; j < parentsOfPeps.size(); j++)
    {
      igb.insert_dependency(createPeptideEvidenceFactor(j,pepEvidences[j]));
      igb.insert_dependency(createSumEvidenceFactor(parentsOfPeps[j],j,j));
      igb.insert_dependency(createPeptideProbabilisticAdderFactor(parentsOfPeps[j],j));
    }
  }

  /* unused but working
  template <typename L>
  void MessagePasserFactory<L>::fillVectorsOfMessagePassersBruteForce(const std::vector<L> & protIDs,
                                                                      const std::vector<std::vector<L>> & parentsOfPeps,
                                                                      const std::vector<double> & pepEvidences,
                                                                      evergreen::InferenceGraphBuilder<L> & igb)
  {
    assert(parentsOfPeps.size() == pepEvidences.size());
    for (std::vector<uiint> parents : parentsOfPeps)
      for (uiint parent : parents)
        assert(std::find(protIDs.begin(), protIDs.end(), parent) != protIDs.end());

    for (uiint pid : protIDs)
      igb.insert_dependency(createProteinFactor(pid));

    for (uiint j = 0; j < parentsOfPeps.size(); j++)
    {
      std::vector<evergreen::TableDependency<std::string> > deps;
      auto pepdep = createSumEvidenceFactor(parentsOfPeps[j],j,j);
      auto sumdep = createSumFactor(parentsOfPeps[j],j);
      igb.insert_dependency(createPeptideEvidenceFactor(j,pepEvidences[j]));
      igb.insert_dependency(pepdep);
      deps.push_back(sumdep);
      for (const auto& parent : parentsOfPeps[j]) {
        deps.push_back(createProteinFactor(parent));
      }

      //igb.insert_dependency(createEmptyPeptideProbabilisticAdderFactor(parentsOfPeps[j],j));
      igb.insert_dependency(createBFPeptideProbabilisticAdderFactor(parentsOfPeps[j],j,deps));
    }
  }
   */

  } // namespace Internal
} // namespace OpenMS
