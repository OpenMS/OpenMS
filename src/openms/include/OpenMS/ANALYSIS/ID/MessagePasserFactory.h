// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
    /// TODO pre-compute for like a hundred parent proteins
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
    /// @param id ID for the LabeledPMF
    /// @param prob peptide evidence probability
    evergreen::TableDependency<Label> createPeptideEvidenceFactor(Label id, double prob);

    /// Conditional probability table of peptide given number of parent proteins, based on model params.
    /// Additionally regularizes on the amount of parent proteins.
    /// @param nr_parents (maximum) number of parent proteins
    /// @param id ID for the LabeledPMF    
    /// @param pep_id ID for the LabeledPMF    
    evergreen::TableDependency<Label> createRegularizingSumEvidenceFactor(size_t nr_parents, Label id, Label pep_id);

    /// Conditional probability table of peptide given number of parent proteins, based on model params.
    /// @param nr_parents (maximum) number of parent proteins
    /// @param id ID for the LabeledPMF    
    /// @param pep_id ID for the LabeledPMF    
    evergreen::TableDependency<Label> createSumEvidenceFactor(size_t nr_parents, Label id, Label pep_id);

    //For extended model. @todo currently unused
    evergreen::TableDependency<Label> createSumFactor(size_t nr_parents, Label nId);
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
     * @param pep_prior Peptide prior (defines at which evidence probability, additional evidence is beneficial)
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

  template <typename Label>
  MessagePasserFactory<Label>::MessagePasserFactory(double alpha, double beta, double gamma, double p, double pep_prior) {
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

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createProteinFactor(Label id, int nrMissingPeps) {
    double prior = gamma_;
    if (nrMissingPeps > 0)
    {
      double powFactor = std::pow(1.0 - alpha_, -nrMissingPeps);
      prior = -prior/(prior * powFactor - prior - powFactor);
    }
    double table[] = {1.0 - prior, prior};
    evergreen::LabeledPMF<Label> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createProteinFactor(Label id, double prior, int nrMissingPeps) {
    if (nrMissingPeps > 0)
    {
      double powFactor = std::pow(1.0 - alpha_, -nrMissingPeps);
      prior = -prior/(prior * powFactor - prior - powFactor);
    }
    double table[] = {1.0 - prior, prior};
    evergreen::LabeledPMF<Label> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createPeptideEvidenceFactor(Label id, double prob) {
    double table[] = {(1 - prob) * (1 - pepPrior_), prob * pepPrior_};
    evergreen::LabeledPMF<Label> lpmf({id}, evergreen::PMF({0L}, evergreen::Tensor<double>::from_array(table)));
    return evergreen::TableDependency<Label>(lpmf,p_);
  }


  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId) {
    evergreen::Tensor<double> table({static_cast<unsigned long>(nrParents + 1) , 2});
    for (unsigned long i=0; i <= nrParents; ++i) {
      double notConditional = notConditionalGivenSum(i);
      unsigned long indexArr[2] = {i,0ul};
      table[indexArr] = notConditional;
      unsigned long indexArr2[2] = {i,1ul};
      table[indexArr2] = 1.0 - notConditional;
    }
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<Label> lpmf({nId, pepId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createRegularizingSumEvidenceFactor(size_t nrParents, Label nId, Label pepId) {
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
    evergreen::LabeledPMF<Label> lpmf({nId, pepId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createSumFactor(size_t nrParents, Label nId) {
    evergreen::Tensor<double> table({nrParents+1});
    for (unsigned long i=0; i <= nrParents; ++i) {
      table[i] = 1.0/(nrParents+1.);
    }
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<Label> lpmf({nId}, evergreen::PMF({0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createReplicateFactor(Label seqId, Label repId) {
    using arr = unsigned long[2];
    evergreen::Tensor<double> table({2,2});
    table[arr{0,0}] = 0.999;
    table[arr{0,1}] = 0.001;
    table[arr{1,0}] = 0.1;
    table[arr{1,1}] = 0.9;
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<Label> lpmf({seqId,repId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::TableDependency<Label> MessagePasserFactory<Label>::createChargeFactor(Label repId, Label chgId, int chg) {
    double chgPrior = chgLLhoods[chg];
    using arr = unsigned long[2];
    evergreen::Tensor<double> table({2,2});
    table[arr{0,0}] = 0.999;
    table[arr{0,1}] = 0.001;
    table[arr{1,0}] = 0.1;
    table[arr{1,1}] = chgPrior;
    //std::cout << table << std::endl;
    evergreen::LabeledPMF<Label> lpmf({repId,chgId}, evergreen::PMF({0L,0L}, table));
    //std::cout << lpmf << std::endl;
    return evergreen::TableDependency<Label>(lpmf,p_);
  }

  template <typename Label>
  evergreen::AdditiveDependency<Label> MessagePasserFactory<Label>::createPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId) {
    std::vector<std::vector<Label>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const Label& l){return std::vector<Label>{l};});
    return evergreen::AdditiveDependency<Label>(parents, {nId}, p_);
  }

  template <typename Label>
  evergreen::AdditiveDependency<Label> MessagePasserFactory<Label>::createPeptideProbabilisticAdderFactor(const std::vector<Label> & parentProteinIDs, Label nId) {
    std::vector<std::vector<Label>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const Label& l){return std::vector<Label>{l};});
    return evergreen::AdditiveDependency<Label>(parents, {nId}, p_);
  }

  template <typename Label>
  evergreen::PseudoAdditiveDependency<Label> MessagePasserFactory<Label>::createBFPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId, const std::vector<evergreen::TableDependency<Label>> & deps) {
    std::vector<std::vector<Label>> parents;
    std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const Label& l){return std::vector<Label>{l};});
    return evergreen::PseudoAdditiveDependency<Label>(parents, {nId}, deps, p_);
  }

  /// Works on a vector of protein indices (potentially not consecutive)
  // TODO we could recollect the protIDs from the union of parents.
  template <typename Label>
  void MessagePasserFactory<Label>::fillVectorsOfMessagePassers(const std::vector<Label>& protIDs,
                                                            const std::vector<std::vector<Label>>& parentsOfPeps,
                                                            const std::vector<double>& pepEvidences,
                                                            evergreen::InferenceGraphBuilder<Label>& igb)
  {
    //TODO asserts could be loosened
    assert(parentsOfPeps.size() == pepEvidences.size());
    for (const std::vector<uiint>& parents : parentsOfPeps)
      for (Label parent : parents)
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
  template <typename Label>
  void MessagePasserFactory<Label>::fillVectorsOfMessagePassersBruteForce(const std::vector<Label> & protIDs,
                                                                      const std::vector<std::vector<Label>> & parentsOfPeps,
                                                                      const std::vector<double> & pepEvidences,
                                                                      evergreen::InferenceGraphBuilder<Label> & igb)
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
