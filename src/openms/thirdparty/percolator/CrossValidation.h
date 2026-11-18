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

#include <cstdint>
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
#include "newton_svm.h"  // OpenMS: SvmSolver

namespace OpenMS { namespace Internal { namespace Percolator {


struct CandidateCposCfrac {
  double cpos;
  double cfrac;
  unsigned int set;
  int nestedSet;
  vector<double> ww;
  int tp;  // validation count; 0 for the only trained pair of a fold
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
                  double decoyFractionTraining,
                  unsigned int numFolds);

  int preIterationSetup(Scores& fullset,
                        SanityCheck& check,
                        const Normalizer* pNorm);

  void train(const Normalizer* pNorm);

  void postIterationProcessing(Scores& fullset, SanityCheck& check);

  void printAllWeights(ostream& weightStream, const Normalizer* pNorm);

  void getAvgWeights(std::vector<double>& weights, const Normalizer* pNorm);
  // the weights of each fold for unnormalized features, as getAvgWeights
  void getRawFoldWeights(std::vector<std::vector<double> >& weights,
                         const Normalizer* pNorm);
  // (input line, fold) of every PSM in the test sets, sorted by line, as
  // recorded by postIterationProcessing when setRecordTestFolds(true) was
  // called; the test sets themselves are emptied by the merge. Empty if
  // the direction check failed, after a reset and with --override, since
  // all PSMs are then scored with the averaged weights. Moves the list out.
  void getTestFolds(
      std::vector<std::pair<unsigned int, unsigned int> >& lineFolds);

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
  // train as Percolator 3.09 did, see Caller's --classic-training
  void inline setClassicTraining(bool on) { classicTraining_ = on; }
  // OpenMS: the solver of the SVMs (SvmSolver, newton_svm.h; default kSolverTron, which ssl.h provides in place of
  // upstream's L2_SVM_MFN)
  void inline setSolver(int solver) { solver_ = solver; }
  // OpenMS (options of the training; the defaults are Percolator's): TRON's stopping tolerance (default EPSILON),
  // warm starts with classic training (each (Cpos, Cneg) pair starts TRON from its weights of the previous
  // iteration instead of 0, without nested cross-validation, whose nested sets are drawn anew in each iteration;
  // the training of setClassicTraining(false) starts every solve from the previous weights anyway), one Cneg
  // candidate (3 x the target/decoy ratio) instead of three. Set before preIterationSetup. A tolerance that is not
  // finite and positive keeps EPSILON.
  // OpenMS (out-of-fold capped training): the folds of the PSMs are given by foldOfPsm (0 to numFolds - 1, the same
  // for all PSMs of a spectrum) instead of being drawn; set before preIterationSetup. Nested cross-validation sets
  // are still drawn (inside a fold's training set, from a generator state per fold and iteration), and each fold's
  // Cneg candidates come from its training set's target/decoy ratio.
  void setFixedFolds(std::function<unsigned int(const PSMDescription&)> foldOfPsm) {
    fixedFold_ = std::move(foldOfPsm);
  }
  // OpenMS (learning curve): train only the folds marked active (empty: all). An inactive fold keeps the weights of
  // its initial direction, and its model is not to be used; the counts of the training then cover the active folds.
  void setActiveFolds(std::vector<char> active) { activeFolds_ = std::move(active); }
  // OpenMS: the PSMs of the training set of fold set, between preIterationSetup and postIterationProcessing (which
  // releases the training sets)
  void getTrainingPsms(std::size_t set, std::vector<const PSMDescription*>& psms) const;
  void setTrainingOptions(double tronEpsilon, bool warmStart, bool singleCneg) {
    tronEpsilon_ = tronEpsilon > 0.0 && tronEpsilon < std::numeric_limits<double>::infinity()
                       ? tronEpsilon : EPSILON;
    warmStart_ = warmStart;
    singleCneg_ = singleCneg;
  }
  // train the SVM on the feature rows in place instead of on copies (see
  // AlgIn), which needs rows that end with the intercept 1.0, as those of a
  // FeatureMemoryPool for which canTrainInPlace(numFeatures) holds; set
  // before preIterationSetup
  void inline setSvmRowsInPlace(bool on) { svmRowsInPlace_ = on; }
  // The number of folds whose training sets are held as copies of the
  // feature rows at once (see AlgIn): the SVMs of these folds are trained
  // in parallel, then the buffers are reused for the next folds. The
  // results do not depend on it. 0, the default, holds all folds at once
  // for up to kMaxFoldsAllAtOnce folds and one fold with more; a larger
  // number than there are folds means all. No effect with
  // setSvmRowsInPlace, which makes no copies. Set before preIterationSetup.
  void inline setFoldsAtOnce(unsigned int n) { foldsAtOnce_ = n; }
  static constexpr unsigned int kMaxFoldsAllAtOnce = 3u;
  // the number of folds held at once for the setting foldsAtOnce and
  // numFolds folds
  static unsigned int foldsHeldAtOnce(unsigned int foldsAtOnce,
                                      unsigned int numFolds);
  // With fewer folds held at once than there are, the passes in which each
  // training set serves a single solve (after the grid search is fixed, with
  // a single candidate pair, and the solve on a whole fold of nested cross
  // validation) would run their solves one after the other. With this
  // setting they train all folds in parallel on the rows in place instead,
  // which needs rows as for setSvmRowsInPlace.
  void inline setSingleSolvesInPlace(bool on) { singleSolvesInPlace_ = on; }
  // let postIterationProcessing list the input line and fold of every test
  // PSM for getTestFolds, before the merge empties the test sets (-N); after
  // a failed direction check under test, a reset or --override, it then
  // scores every test bin with the first fold's vector, since -N scores all
  // PSMs with the averaged weights
  void inline setRecordTestFolds(bool on) { recordTestFolds_ = on; }
  // How the scores of the test bins are put on a common scale before they
  // are merged, see MergePlan and Caller's --merge-normalization. The plan's
  // fdr is set here, to the selection FDR of the training.
  void setMergePlan(const MergePlan& plan) {
    mergePlan_ = plan;
    mergePlan_.fdr = selectionFdr_;
    usedMergeNormalization_ = mergePlan_.usedScheme(false);
  }
  const MergePlan& mergePlan() const { return mergePlan_; }
  // What the direction check of postIterationProcessing decided, and the
  // scheme it then merged with (MergePlan::usedScheme): the plan's scheme,
  // or its secondChanceScheme when the second chance kept the trained SVMs;
  // always test when the scores are not rescaled (--static)
  SanityCheck::DirectionOutcome directionOutcome() const {
    return directionOutcome_;
  }
  MergeNormalization usedMergeNormalization() const {
    return usedMergeNormalization_;
  }
  // The scheme for numFolds folds when none is chosen: test up to three
  // folds, the default, whose test bins hold a third of the PSMs or more;
  // train-decoys with more folds, since the q = fdr anchor of test is noisy
  // in smaller test bins (see Caller's --num-folds).
  static MergeNormalization defaultMergeNormalization(unsigned int numFolds);

 protected:
  // the training sets of the folds held at once, one per nested fold
  std::vector<AlgIn> svmInputs_;
  // the training sets of all folds as row pointers, see
  // setSingleSolvesInPlace
  std::vector<AlgIn> svmRowLists_;
  std::vector<std::vector<double> > weights_;  // svm weights for each fold
  std::vector<CandidateCposCfrac>
      classWeightsPerFold_;  // cpos, cneg pairs to train for each nested CV
                             // fold

  bool quickValidation_;
  bool classicTraining_ = false;
  // OpenMS: see setFixedFolds (empty: the folds are drawn)
  std::function<unsigned int(const PSMDescription&)> fixedFold_;
  // OpenMS: see setActiveFolds
  std::vector<char> activeFolds_;
  // OpenMS: with fixed folds, the folds whose training has stopped (see train) and the targets with q < test FDR in
  // the training set of each fold after the last step
  std::vector<char> stoppedFolds_;
  std::vector<int> foldPositives_;
  bool isActiveFold(unsigned int set) const {
    return (activeFolds_.empty() || activeFolds_[set] != 0) &&
           (stoppedFolds_.empty() || stoppedFolds_[set] == 0);
  }
  // OpenMS: with fixed folds, the generator state of the nested sets of each fold and iteration derives from this
  // state (at preIterationSetup) and the iteration (see doStep)
  uint64_t nestedSeedBase_ = 0;
  // OpenMS: per fold, with fixed folds (see initializeGridSearch); empty otherwise
  std::vector<std::vector<double> > candidatesCfracPerFold_;
  // OpenMS: see setSolver and setTrainingOptions
  int solver_ = kSolverTron;
  double tronEpsilon_ = EPSILON;
  bool warmStart_ = false;
  bool singleCneg_ = false;
  bool svmRowsInPlace_ = false;
  unsigned int foldsAtOnce_ = 0u;
  bool singleSolvesInPlace_ = false;
  bool recordTestFolds_ = false;
  std::vector<std::pair<unsigned int, unsigned int> > testFolds_;
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
  MergePlan mergePlan_;
  SanityCheck::DirectionOutcome directionOutcome_ =
      SanityCheck::DirectionOutcome::KEPT;
  MergeNormalization usedMergeNormalization_ = MergeNormalization::TEST;

  double decoyFractionTraining_;

  const static double requiredIncreaseOver2Iterations_;
  const static double requiredIncreaseOver2IterationsDefault_;
  const static unsigned int gridSearchIterations_;

  // after gridSearchIterations_ only each fold's selected (cpos, cfrac)
  // pair is trained
  unsigned int iteration_ = 0;
  bool gridFixed_ = false;
  std::vector<std::pair<double, double> > selectedPairs_;
  bool isSelectedPair(const CandidateCposCfrac& pair) const {
    return !gridFixed_ || (pair.cpos == selectedPairs_[pair.set].first &&
                           pair.cfrac == selectedPairs_[pair.set].second);
  }

  unsigned int numFolds_;  // number of folds for cross validation
  std::vector<Scores> trainScores_, testScores_;
  // whether trainScores_ are scored and sorted with the current weights_
  bool trainScoresCurrent_ = false;
  std::vector<double> candidatesCpos_, candidatesCfrac_;

  // buffers of the SVM solves of one thread, reused for its pairs
  struct SvmTrainingBuffers {
    SolverWorkspace solver;
    std::vector<double> outputs;
  };

  // the number of folds whose training sets are held at once
  unsigned int foldsAtOnce() const {
    return foldsHeldAtOnce(foldsAtOnce_, numFolds_);
  }
  // The training sets to fill for a pass over the folds and the number of
  // folds per group of that pass: the set of nested fold b of the i-th fold
  // of a group is element i * nestedXvalBins_ + b. singleSolves: each
  // training set serves one solve.
  std::vector<AlgIn>& svmInputsOfPass(bool singleSolves,
                                      unsigned int& foldsPerPass);
  void initializeGridSearch(double targetDecoySizeRatio);
  // OpenMS: whether train() runs the folds as chains of steps (trainFoldChains): with OpenMP tasks and more than one
  // thread, without nested cross-validation, -R or verbosity above 2, with all folds active and their training sets
  // held at once
  bool trainsFoldChains() const;
  // OpenMS: the iterations of train(), each fold as a chain of tasks that waits only for its own (Cpos, Cneg) pairs,
  // with the results and reports of the iterations of train()
  void trainFoldChains();
  // OpenMS: group, the TronPairGroup of svmInput (or nullptr)
  void trainCpCnPair(CandidateCposCfrac& cpCnFold,
                     const options& pOptions,
                     const AlgIn& svmInput,
                     SvmTrainingBuffers& buffers,
                     TronPairGroup* group);
  // OpenMS: solves the SVM of svmInput with the solver of setSolver, from the weights in w; the TRON solves of the
  // pairs of one training set at a time share its copy of the set in group (ssl.h)
  void trainSvm(const AlgIn& svmInput, const options& pOptions, double cpos,
                double cneg, std::vector<double>& w,
                std::vector<double>& outputs,
                TronPairGroup* group = nullptr) const;

  int mergeCpCnPairs(double selectionFdr,
                     const options& pOptions,
                     std::vector<std::vector<Scores> >& nestedTestScoresVec,
                     const vector<double>& cpos_vec,
                     const vector<double>& cfrac_vec);
  int doStep(const Normalizer* pNorm, double selectionFdr);
  // throws, or warns with --no-terminate, if a weight of a fold is not finite
  void checkWeightsFinite() const;
  // lists the input line and fold of every test PSM in testFolds_
  void recordTestFolds();

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
