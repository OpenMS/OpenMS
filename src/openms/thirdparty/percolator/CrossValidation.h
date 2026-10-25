/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef CROSSVALIDATION_H_
#define CROSSVALIDATION_H_


#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "DataSet.h"
#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "Globals.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "ssl.h"
#include "newton_svm.h"

namespace OpenMS { namespace Internal { namespace Percolator {

struct CandidateCposCfrac {
  double cpos;
  double cfrac;
  unsigned int set;
  int nestedSet;
  vector<double> ww;
  int tp;
};

class CrossValidation {
 public:
  CrossValidation(bool quickValidation,
                  bool reportPerformanceEachIteration,
                  double testFdr,
                  double selectionFdr,
                  double initialSelectionFdr,
                  double selectedCpos,
                  double selectedCneg,
                  unsigned int niter,
                  bool usePi0,
                  unsigned int nestedXvalBins,
                  bool trainBestPositive,
                  unsigned int numThreads,
                  bool skipNormalizeScores,
                  double decoyFractionTraining,
                  unsigned int numFolds);

  ~CrossValidation();

  int preIterationSetup(Scores& fullset,
                        SanityCheck* pCheck,
                        const Normalizer* pNorm,
                        FeatureMemoryPool& featurePool);

  void train(const Normalizer* pNorm);

  // OpenMS (experimental options of the training; the defaults are Percolator's training): TRON's
  // stopping tolerance (default EPSILON), warm starts (each (Cpos, Cneg) pair starts TRON from its
  // weights of the previous iteration instead of 0), one Cneg candidate (3 x the target/decoy
  // ratio) instead of three. Set before preIterationSetup. A tolerance that is not finite and
  // positive keeps EPSILON; warm starts apply without nested cross-validation only (its nested
  // sets are drawn anew in each iteration).
  void setTrainingOptions(double tronEpsilon, bool warmStart, bool singleCneg) {
    tronEpsilon_ = tronEpsilon > 0.0 && tronEpsilon < std::numeric_limits<double>::infinity()
                       ? tronEpsilon : EPSILON;
    warmStart_ = warmStart;
    singleCneg_ = singleCneg;
  }

  void postIterationProcessing(Scores& fullset, SanityCheck* pCheck);

  void printAllWeights(ostream& weightStream, const Normalizer* pNorm);

  void getAvgWeights(std::vector<double>& weights, const Normalizer* pNorm);

  // OpenMS (out-of-fold capped training): the folds of the PSMs are given by foldOfPsm (0 to
  // numFolds - 1, the same for all PSMs of a spectrum) instead of being drawn; set before
  // preIterationSetup. Nested cross-validation sets are still drawn.
  void setFixedFolds(std::function<unsigned int(const PSMDescription&)> foldOfPsm) {
    fixedFold_ = std::move(foldOfPsm);
  }
  // OpenMS: the weights of each fold in the raw feature space (bias last), as getAvgWeights
  // un-normalizes them, without averaging
  void getRawFoldWeights(std::vector<std::vector<double> >& weights,
                         const Normalizer* pNorm) const;
  // OpenMS: the PSMs of the training set of fold set (after preIterationSetup)
  void getTrainingPsms(std::size_t set,
                       std::vector<const PSMDescription*>& psms) const;

  void inline setSelectedCpos(double cpos) { selectedCpos_ = cpos; }
  double inline getSelectedCpos() { return selectedCpos_; }
  void inline setSelectedCneg(double cneg) { selectedCneg_ = cneg; }
  double inline getSelectedCneg() { return selectedCneg_; }

  void inline setSelectionFdr(double fdr) { selectionFdr_ = fdr; }
  double inline getSelectionFdr() { return selectionFdr_; }
  void inline setTestFdr(double fdr) { testFdr_ = fdr; }

  void inline setNiter(unsigned int n) { niter_ = n; }
  unsigned int inline getNiter() { return niter_; }
  void inline setQuickValidation(bool on) { quickValidation_ = on; }
  void inline setReportPerformanceEachIteration(bool on) {
    reportPerformanceEachIteration_ = on;
  }
  // OpenMS: the solver of the SVMs (SvmSolver, newton_svm.h; default kSolverTron)
  void inline setSolver(int solver) { solver_ = solver; }

 protected:
  std::vector<AlgIn*> svmInputs_;
  std::vector<std::vector<double> > weights_;  // svm weights for each fold
  std::vector<CandidateCposCfrac>
      classWeightsPerFold_;  // cpos, cneg pairs to train for each nested CV
                             // fold

  bool quickValidation_;
  bool usePi0_;
  bool reportPerformanceEachIteration_;

  unsigned int numThreads_;

  double testFdr_;       // fdr used for cross validation performance measuring
  double selectionFdr_;  // fdr used for determining positive training set
  double initialSelectionFdr_;  // fdr used for determining positive training
                                // set in first iteration
  double selectedCpos_;  // soft margin parameter for positive training set
  double selectedCneg_;  // soft margin parameter for negative training set

  unsigned int niter_;
  unsigned int nestedXvalBins_;

  bool trainBestPositive_;
  bool skipNormalizeScores_;

  double decoyFractionTraining_;

  const static double requiredIncreaseOver2Iterations_;

  unsigned int numFolds_;  // number of folds for cross validation
  // OpenMS: see setTrainingOptions
  double tronEpsilon_ = EPSILON;
  bool warmStart_ = false;
  bool singleCneg_ = false;
  std::vector<Scores> trainScores_, testScores_;
  // OpenMS: whether the scores of trainScores_ are those of weights_ (calcScores, at the end
  // of mergeCpCnPairs) and have not changed since; see doStep
  bool trainScoresOfWeights_ = false;
  // OpenMS: see setFixedFolds (empty: the folds are drawn)
  std::function<unsigned int(const PSMDescription&)> fixedFold_;
  // OpenMS: the solver of the SVMs (SvmSolver)
  int solver_ = kSolverTron;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  void initializeGridSearch(double targetDecoySizeRatio);
  // OpenMS: returns whether the weights' ranking does not depend on the costs (trainSvm)
  bool trainCpCnPair(CandidateCposCfrac& cpCnFold,
                     options& pOptions,
                     AlgIn* svmInput);

  // OpenMS: solves the SVM of svmInput with the solver (solver_), from the weights in w
  // (Newton, LDA; TRON starts from 0); returns whether LDA gave the weights
  bool trainSvm(const AlgIn& svmInput, options& pOptions, double cpos, double cneg,
                vector_double& w, vector_double& outputs);

  int mergeCpCnPairs(double selectionFdr,
                     options& pOptions,
                     const std::vector<int>& pairTruePoses,
                     const vector<double>& cpos_vec,
                     const vector<double>& cfrac_vec);
  int doStep(const Normalizer* pNorm, double selectionFdr);
  // OpenMS: see CrossValidation.cpp
  void setupFoldInputs(int set, double selectionFdr, std::vector<Scores>& nestedTestScores,
                       const Scores** testSets);
  void trainFoldChains();

  void printSetWeights(ostream& weightStream, unsigned int set);
  void printRawSetWeights(ostream& weightStream,
                          unsigned int set,
                          const Normalizer* pNorm);

  void printAllWeightsColumns(ostream& weightStream);
  void printAllRawWeightsColumns(ostream& weightStream,
                                 const Normalizer* pNorm);
  void printAllWeightsColumns(
      std::vector<std::vector<double> > weightMatrix,
      ostream& weightStream);
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*CROSSVALIDATION_H_*/
