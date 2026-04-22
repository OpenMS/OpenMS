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

#include <iomanip>
#include <iostream>
#include <vector>

#include "DataSet.h"
#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "Globals.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Scores.h"
#include "ssl.h"

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

  void postIterationProcessing(Scores& fullset, SanityCheck* pCheck);

  void printAllWeights(ostream& weightStream, const Normalizer* pNorm);

  void getAvgWeights(std::vector<double>& weights, const Normalizer* pNorm);

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
  std::vector<Scores> trainScores_, testScores_;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  void initializeGridSearch(double targetDecoySizeRatio);
  void trainCpCnPair(CandidateCposCfrac& cpCnFold,
                     options& pOptions,
                     AlgIn* svmInput);

  int mergeCpCnPairs(double selectionFdr,
                     options& pOptions,
                     std::vector<std::vector<Scores> >& nestedTestScoresVec,
                     const vector<double>& cpos_vec,
                     const vector<double>& cfrac_vec);
  int doStep(const Normalizer* pNorm, double selectionFdr);

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

#endif /*CROSSVALIDATION_H_*/
