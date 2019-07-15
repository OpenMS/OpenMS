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

#ifndef OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP
#define OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP

#include <cmath>
#include <Evergreen/evergreen.hpp>
#include <Utility/inference_utilities.hpp>

typedef unsigned long int uiint;

template <typename Label>
class MessagePasserFactory {
private:
    const int minInputsPAF = 3;
    double alpha, beta, gamma, p, pepPrior;
    Label offset;
    std::map<int, double> chgPriors = {{1,0.7},{2,0.9},{3,0.7},{4,0.5},{5,0.5}};

    inline double notConditionalGivenSum(unsigned long summ) {
        // use log for better precision
        return pow(2., log2(1. - beta) + summ * log2(1. - alpha));
        //return std::pow((1.0 - alpha), summ) * (1.0 - beta);
    }

public:
    TableDependency<Label> createProteinFactor(Label id, int nrMissingPeps = 0);
    TableDependency<Label> createProteinFactor(Label id, double prior, int nrMissingPeps = 0);

    TableDependency<Label> createPeptideEvidenceFactor(Label id, double prob);

    TableDependency<Label> createRegularizingSumEvidenceFactor(size_t nrParents, Label nId, Label pepId);
    TableDependency<Label> createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId);

    TableDependency<Label> createSumFactor(size_t nrParents, Label nId);

    TableDependency<Label> createReplicateFactor(Label seqID, Label repID);
    TableDependency<Label> createChargeFactor(Label repID, Label chargeID, int chg);

    AdditiveDependency<Label> createPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId);
    AdditiveDependency<Label> createPeptideProbabilisticAdderFactor(const std::vector<Label> & parentProteinIDs, Label nId);

    PseudoAdditiveDependency<Label> createBFPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId, const std::vector<TableDependency <Label> > & deps);

    MessagePasserFactory<Label>(double alpha, double beta, double gamma, double p, double pepPrior);



    /// Works on a vector of protein indices (potentially not consecutive)
    // TODO we could recollect the protIDs from the union of parents.
    void fillVectorsOfMessagePassers(const std::vector<Label> & protIDs,
                                     const std::vector<std::vector<Label>> & parentsOfPeps,
                                     const std::vector<double> & pepEvidences,
                                     InferenceGraphBuilder<Label> & igb);

    //void fillVectorsOfMessagePassersBruteForce(const std::vector<Label> & protIDs,
    //                                 const std::vector<std::vector<Label>> & parentsOfPeps,
    //                                 const std::vector<double> & pepEvidences,
    //                                 InferenceGraphBuilder<Label> & igb);

    //const std::vector<std::set<Label>> getPosteriorVariables(const std::vector<uiint> & protIDs);
    //const std::vector<std::vector<Label>> getPosteriorVariablesVectors(const std::vector<uiint> & protIDs);
    //const std::vector<std::set<Label>> getPosteriorVariables(uiint rangeProtIDs);
};

//IMPLEMENTATIONS:

template <typename L>
MessagePasserFactory<L>::MessagePasserFactory(double alpha_, double beta_, double gamma_, double p_, double pep_prior_) {
  assert(0. <= alpha_ && alpha_ <= 1.);
  assert(0. <= beta_ && beta_ <= 1.);
  assert(0. <= gamma_ && gamma_ <= 1.);
  //Note: smaller than 1 might be possible but is untested right now.
  assert(p_ >= 1.);
  assert(0. < pep_prior_ && pep_prior_ < 1.);
  alpha = alpha_;
  beta = beta_;
  gamma = gamma_;
  p = p_;
  pepPrior = pep_prior_;
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createProteinFactor(L id, int nrMissingPeps) {
  double prior = gamma;
  if (nrMissingPeps > 0)
  {
    double powFactor = std::pow(1.0 - alpha, -nrMissingPeps);
    prior = -prior/(prior * powFactor - prior - powFactor);
  }
  double table[] = {1.0 - prior, prior};
  LabeledPMF<L> lpmf({id}, PMF({0L}, Tensor<double>::from_array(table)));
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createProteinFactor(L id, double prior, int nrMissingPeps) {
  if (nrMissingPeps > 0)
  {
    double powFactor = std::pow(1.0 - alpha, -nrMissingPeps);
    prior = -prior/(prior * powFactor - prior - powFactor);
  }
  double table[] = {1.0 - prior, prior};
  LabeledPMF<L> lpmf({id}, PMF({0L}, Tensor<double>::from_array(table)));
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createPeptideEvidenceFactor(L id, double prob) {
  double table[] = {(1 - prob) * (1 - pepPrior), prob * pepPrior};
  LabeledPMF<L> lpmf({id}, PMF({0L}, Tensor<double>::from_array(table)));
  return TableDependency<L>(lpmf,p);
}


template <typename L>
TableDependency<L> MessagePasserFactory<L>::createSumEvidenceFactor(size_t nrParents, L nId, L pepId) {
  Tensor<double> table({static_cast<unsigned long>(nrParents + 1) , 2});
  for (unsigned long i=0; i <= nrParents; ++i) {
    double notConditional = notConditionalGivenSum(i);
    unsigned long indexArr[2] = {i,0ul};
    table[indexArr] = notConditional;
    unsigned long indexArr2[2] = {i,1ul};
    table[indexArr2] = 1.0 - notConditional;
  }
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({nId, pepId}, PMF({0L,0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createRegularizingSumEvidenceFactor(size_t nrParents, L nId, L pepId) {
  Tensor<double> table({static_cast<unsigned long>(nrParents + 1) , 2});
  unsigned long z[2]{0ul,0ul};
  unsigned long z1[2]{0ul,1ul};
  table[z] = 1. - beta;
  table[z1] = beta;
  for (unsigned long i=1; i <= nrParents; ++i) {
    double notConditional = notConditionalGivenSum(i);
    unsigned long indexArr[2] = {i,0ul};
    table[indexArr] = notConditional / i;
    unsigned long indexArr2[2] = {i,1ul};
    table[indexArr2] = (1.0 - notConditional) / i;
  }
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({nId, pepId}, PMF({0L,0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createSumFactor(size_t nrParents, L nId) {
  Tensor<double> table({nrParents+1});
  for (unsigned long i=0; i <= nrParents; ++i) {
    table[i] = 1.0/(nrParents+1);
  }
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({nId}, PMF({0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createReplicateFactor(L seqId, L repId) {
  using arr = unsigned long[2];
  Tensor<double> table({2,2});
  table[arr{0,0}] = 0.999;
  table[arr{0,1}] = 0.001;
  table[arr{1,0}] = 0.1;
  table[arr{1,1}] = 0.9;
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({seqId,repId}, PMF({0L,0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createChargeFactor(L repId, L chgId, int chg) {
  double chgPrior = chgPriors[chg];
  using arr = unsigned long[2];
  Tensor<double> table({2,2});
  table[arr{0,0}] = 0.999;
  table[arr{0,1}] = 0.001;
  table[arr{1,0}] = 0.1;
  table[arr{1,1}] = chgPrior;
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({repId,chgId}, PMF({0L,0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
AdditiveDependency<L> MessagePasserFactory<L>::createPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId) {
  std::vector<std::vector<L>> parents;
  std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
  return AdditiveDependency<L>(parents, {nId}, p);
}

template <typename L>
AdditiveDependency<L> MessagePasserFactory<L>::createPeptideProbabilisticAdderFactor(const std::vector<L> & parentProteinIDs, L nId) {
  std::vector<std::vector<L>> parents;
  std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
  return AdditiveDependency<L>(parents, {nId}, p);
}

template <typename L>
PseudoAdditiveDependency<L> MessagePasserFactory<L>::createBFPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId, const std::vector<TableDependency<L>> & deps) {
  std::vector<std::vector<L>> parents;
  std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
  return PseudoAdditiveDependency<L>(parents, {nId}, deps, p);
}

/// Works on a vector of protein indices (potentially not consecutive)
// TODO we could recollect the protIDs from the union of parents.
template <typename L>
void MessagePasserFactory<L>::fillVectorsOfMessagePassers(const std::vector<L> & protIDs,
                                                          const std::vector<std::vector<L>> & parentsOfPeps,
                                                          const std::vector<double> & pepEvidences,
                                                          InferenceGraphBuilder<L> & igb)
{
  //TODO asserts could be loosened
  assert(parentsOfPeps.size() == pepEvidences.size());
  for (std::vector<uiint> parents : parentsOfPeps)
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
                                                                    InferenceGraphBuilder<L> & igb)
{
  assert(parentsOfPeps.size() == pepEvidences.size());
  for (std::vector<uiint> parents : parentsOfPeps)
    for (uiint parent : parents)
      assert(std::find(protIDs.begin(), protIDs.end(), parent) != protIDs.end());

  for (uiint pid : protIDs)
    igb.insert_dependency(createProteinFactor(pid));

  for (uiint j = 0; j < parentsOfPeps.size(); j++)
  {
    std::vector<TableDependency<std::string> > deps;
    auto pepdep = createSumEvidenceFactor(parentsOfPeps[j],j,j);
    auto sumdep = createSumFactor(parentsOfPeps[j],j);
    igb.insert_dependency(createPeptideEvidenceFactor(j,pepEvidences[j]));
    igb.insert_dependency(pepdep);
    deps.push_back(sumdep);
    for (auto parent : parentsOfPeps[j]) {
      deps.push_back(createProteinFactor(parent));
    }

    //igb.insert_dependency(createEmptyPeptideProbabilisticAdderFactor(parentsOfPeps[j],j));
    igb.insert_dependency(createBFPeptideProbabilisticAdderFactor(parentsOfPeps[j],j,deps));
  }
}
 */

/* Not needed anymore. We use indices directly now.
template <typename L>
const std::vector<std::set<L>> MessagePasserFactory<L>::getPosteriorVariables(const std::vector<L> & protIDs){
    std::vector<std::set<L>> varSets{};
    for (L protID : protIDs){
        std::set<L> varSet{"Pr" + std::to_string(protID)};
        varSets.push_back(varSet);
    }
    return varSets;
}

template <typename L>
const std::vector<std::vector<std::string>> MessagePasserFactory<L>::getPosteriorVariablesVectors(const std::vector<uiint> & protIDs){
  std::vector<std::vector<std::string>> varVecs{};
  for (uiint protID : protIDs){
    std::vector<std::string> varVec{"Pr" + std::to_string(protID)};
    varVecs.push_back(varVec);
  }
  return varVecs;
}

template <typename L>
const std::vector<std::set<std::string>> MessagePasserFactory<L>::getPosteriorVariables(uiint rangeProtIDs){
    std::vector<std::set<std::string>> varSets{};
    for (uiint i=0; i < rangeProtIDs; ++i){
        std::set<std::string> varSet{"Pr" + std::to_string(i)};
        varSets.push_back(varSet);
    }
    return varSets;
}*/

#endif //OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP
