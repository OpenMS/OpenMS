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

#include "CrossValidation.h"
#include "PosteriorEstimator.h"
#include "Timer.h"

#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "OpenMPThreads.h"  // OpenMS: maxThreads()
#include "OpenMPTasks.h"    // OpenMS: PERCOLATOR_OMP_TASKS, PERCOLATOR_TASK_PRAGMA (trainFoldChains)

#include <atomic>
#include <memory>
#include <mutex>

namespace OpenMS { namespace Internal { namespace Percolator {

// checks cross validation convergence in case of quickValidation_
const double CrossValidation::requiredIncreaseOver2Iterations_ = 0.01;
// checks convergence in normal (not classic) training
const double CrossValidation::requiredIncreaseOver2IterationsDefault_ = 0.001;
// number of iterations in which Cpos and Cneg are selected; afterwards each
// fold keeps its pair (not in classic training)
const unsigned int CrossValidation::gridSearchIterations_ = 3u;

// OpenMS: Newton's tolerance (as TRON's, the gradient norm relative to that at 0) in the iterations before the last
// one of classic training without quick validation, which solves to the tolerance of the options: the earlier
// iterations only choose the positives of the next. (A training that can stop early does not know its last iteration.)
constexpr double kNewtonEarlyEpsilon = 1e-4;

// the number of threads a parallel region can have, and the number of the
// calling thread in its region
// OpenMS: maxThreads() (OpenMPThreads.h), omp_get_max_threads() at most the cap of the calling thread; the regions
// below get it as num_threads, so that the thread numbers stay below it
static std::size_t maxNumThreads() {
  return static_cast<std::size_t>(maxThreads());
}

static std::size_t threadNum() {
#ifdef _OPENMP
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0u;
#endif
}

CrossValidation::CrossValidation(bool quickValidation,
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
                                 unsigned int numFolds)
    : quickValidation_(quickValidation),
      usePi0_(usePi0),
      reportPerformanceEachIteration_(reportPerformanceEachIteration),
      numThreads_(numThreads),
      testFdr_(testFdr),
      selectionFdr_(selectionFdr),
      initialSelectionFdr_(initialSelectionFdr),
      selectedCpos_(selectedCpos),
      selectedCneg_(selectedCneg),
      niter_(niter),
      nestedXvalBins_(nestedXvalBins),
      trainBestPositive_(trainBestPositive),
      decoyFractionTraining_(decoyFractionTraining),
      numFolds_(numFolds) {
  // initialize selectionFdr_ and initialSelectionFdr_ if they have not been set
  // explicitly by the user.
  if (selectionFdr_ <= 0.0) {
    selectionFdr_ = testFdr_;
    if (initialSelectionFdr_ <= 0.0) {
      initialSelectionFdr_ = testFdr_;
    }
  }
  mergePlan_.fdr = selectionFdr_;
}

/**
 * Sets up the SVM classifier:
 * - divide dataset into training and test sets for each fold
 * - set parameters (fdr, soft margin)
 *
 * Side effects: updates the following member variables:
 * - `svmInputs_`
 * - `trainScores_`, `testScores_` (in createXvalSetsBySpectrum)
 * - `weights_`
 * - `candidatesCpos_`, `candidatesCfrac_` and `classWeightsPerFold_` (in
 * `initializeGridSearch`)
 * - empties `fullset`, whose PSMs are then held by the folds
 *
 * @param w_ vector of SVM weights
 * @return number of positives for initial setup
 */
int CrossValidation::preIterationSetup(Scores& fullset,
                                       SanityCheck& check,
                                       const Normalizer* pNorm) {
  assert(nestedXvalBins_ >= 1u);
  // The SVM inputs are made by svmInputsOfPass when they are first needed.
  if (VERB > 1 && !svmRowsInPlace_ && foldsAtOnce() < numFolds_) {
    cerr << "Training the SVMs on copies of the feature rows of "
         << foldsAtOnce() << " of the " << numFolds_ << " folds at a time"
         << (singleSolvesInPlace_
                 ? ", and on the rows in place where a fold trains one SVM"
                 : "")
         << " (--folds-in-memory)." << endl;
  }

  unsigned int decoysPerTarget = 1u;  // TODO: make this a command line argument
  double targetDecoySizeRatio = 1.0;
  const size_t numPsms = fullset.size();
  {
    StageTimer stageTimer("cross-validation split");
    nestedSeedBase_ = PseudoRandom::getSeed();  // OpenMS: see doStep (fixed folds)
    fullset.createXvalSetsBySpectrum(trainScores_, testScores_, numFolds_,
                                     decoyFractionTraining_, decoysPerTarget,
                                     &fixedFold_);
    // the folds hold copies of all PSMs; postIterationProcessing merges them
    // back into fullset
    targetDecoySizeRatio = fullset.getTargetDecoySizeRatio();
    fullset.reset();
  }
  // A fold without PSMs cannot be scored or merged; this happens when the
  // input has about as many spectra as folds, or fewer.
  for (unsigned int set = 0; set < numFolds_; ++set) {
    if (testScores_[set].size() == 0u) {
      ostringstream oss;
      oss << "Error: cross-validation fold " << set + 1 << " of " << numFolds_
          << " received no PSMs, the input has too few spectra for "
          << numFolds_ << " folds. Use fewer folds (--num-folds).\n";
      Globals::failUnlessNoTerminate(oss.str(),
                                     "No-terminate flag set: ignoring error.");
      break;
    }
  }
  // With very few spectra all of them can fall into the test set of one
  // split, whose SVM then has nothing to train on. With two or more folds
  // the other test sets are then empty, so this is only reached with the
  // no-terminate flag; one fold trains on all PSMs.
  for (unsigned int set = 0; set < numFolds_; ++set) {
    if (trainScores_[set].size() == 0 && numPsms > 0 && VERB > 0) {
      cerr << "Warning: cross validation split " << set + 1 << " of "
           << numFolds_ << " has no PSMs to train on, as all "
           << numPsms << " PSMs fall into its test set. The input "
           << "has too few spectra for " << numFolds_ << "-fold cross "
           << "validation." << endl;
    }
  }

  // initialize weights vector for all folds
  weights_ = vector<vector<double> >(
      numFolds_, vector<double>(FeatureNames::getNumFeatures() + 1));
  StageTimer stageTimer("initial direction");
  // the merged counts of the direction check merge the test bins as
  // postIterationProcessing does
  check.setMergePlan(mergePlan_);
  int numPositive =
      check.getInitDirection(testScores_, trainScores_, pNorm, weights_,
                               testFdr_, initialSelectionFdr_);

  initializeGridSearch(targetDecoySizeRatio);
  return numPositive;
}

/**
 * Fill member variables with grid search parameters.
 *
 * Side effects: updates the candidatesCpos_, candidatesCfrac_ and
 * classWeightsPerFold_ member variables.
 *
 * @param targetDecoySizeRatio ratio of target to decoy PSMs
 */
void CrossValidation::initializeGridSearch(double targetDecoySizeRatio) {
  if (selectedCpos_ > 0) {
    candidatesCpos_.push_back(selectedCpos_);
  } else {
    candidatesCpos_.push_back(10);
    candidatesCpos_.push_back(1);
    candidatesCpos_.push_back(0.1);
    if (VERB > 0) {
      cerr << "Selecting Cpos by cross-validation." << endl;
    }
  }
  if (selectedCpos_ > 0 && selectedCneg_ > 0) {
    candidatesCfrac_.push_back(selectedCneg_ / selectedCpos_);
  } else {
    // OpenMS: with singleCneg_, the middle candidate only
    if (!singleCneg_) candidatesCfrac_.push_back(1.0 * targetDecoySizeRatio);
    candidatesCfrac_.push_back(3.0 * targetDecoySizeRatio);
    if (!singleCneg_) candidatesCfrac_.push_back(10.0 * targetDecoySizeRatio);
    if (VERB > 0) {
      cerr << "Selecting Cneg by cross-validation." << endl;
    }
    // OpenMS: with fixed folds (out-of-fold capped training), each fold's candidates from the target/decoy ratio of
    // its own training set (the ratio of fullset counts the fold's own labels)
    if (fixedFold_) {
      for (unsigned int set = 0; set < numFolds_; ++set) {
        const double ratio = trainScores_[set].getTargetDecoySizeRatio();
        std::vector<double> cfracs;
        if (!singleCneg_) cfracs.push_back(1.0 * ratio);
        cfracs.push_back(3.0 * ratio);
        if (!singleCneg_) cfracs.push_back(10.0 * ratio);
        candidatesCfracPerFold_.push_back(cfracs);
      }
    }
  }

  // Quick validation trains only Cpos = Cneg = 1. The candidate lists are
  // replaced as well, since the selection after nested cross validation
  // loops over them and would otherwise pick a pair that was never trained.
  if (quickValidation_) {
    candidatesCpos_.assign(1, 1.0);
    candidatesCfrac_.assign(1, 1.0);
    candidatesCfracPerFold_.clear();  // OpenMS
  }

  // Form cpos, cneg pairs per nested CV fold
  CandidateCposCfrac cpCnFold;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    for (unsigned int nestedSet = 0; nestedSet < nestedXvalBins_; nestedSet++) {
      for (double cpos : candidatesCpos_) {
        // OpenMS: the fold's own Cneg candidates with fixed folds
        for (double cfrac : candidatesCfracPerFold_.empty() ? candidatesCfrac_ : candidatesCfracPerFold_[set]) {
          cpCnFold.cpos = cpos;
          cpCnFold.cfrac = cfrac;
          cpCnFold.set = set;
          cpCnFold.nestedSet = nestedSet;
          cpCnFold.tp = 0;
          cpCnFold.ww.clear();
          for (int i = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
               i--;) {
            cpCnFold.ww.push_back(0);
          }
          classWeightsPerFold_.push_back(cpCnFold);
        }
      }
    }
  }
}

/**
 * Train the SVM using several cross validation iterations
 * @param pNorm Normalization object
 */
void CrossValidation::train(const Normalizer* pNorm) {
  // OpenMS: the solves of the training reuse the buffers of TRON's workspaces (ssl.h), also the large ones
  TronWorkspaceHold tronWorkspaceHold;
  if (VERB > 0) {
    cerr << "---Training with Cpos";
    if (selectedCpos_ > 0) {
      cerr << "=" << selectedCpos_;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", Cneg";
    if (selectedCneg_ > 0) {
      cerr << "=" << selectedCneg_;
    } else {
      cerr << " selected by cross validation";
    }
    cerr << ", initial_fdr=" << initialSelectionFdr_;
    cerr << ", fdr=" << selectionFdr_ << endl;
  }

  // iterate
  int foundPositivesOldOld = 0, foundPositivesOld = 0, foundPositives = 0;
  // OpenMS: with fixed folds (out-of-fold capped training), each fold stops by the rule below on its own count, the
  // targets with q < test FDR in its training set: the training sets of the other folds hold its PSMs, so their
  // sum would let the fold's held-out labels decide when its model stops training. A stopped fold keeps its weights.
  const bool stopPerFold = static_cast<bool>(fixedFold_);
  std::vector<int> foldPositivesOld(numFolds_, 0), foldPositivesOldOld(numFolds_, 0);
  stoppedFolds_.assign(numFolds_, 0);
  foldPositives_.assign(numFolds_, 0);
  // OpenMS: the folds as chains of tasks (trainFoldChains), with the results and reports of this loop
  const bool foldChains = trainsFoldChains();
  if (foldChains) {
    trainFoldChains();
  }
  for (unsigned int i = 0; !foldChains && i < niter_; i++) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << ":\t";
    }

    double selectionFdr = selectionFdr_;
    if (i == 0u) {
      selectionFdr = initialSelectionFdr_;
    }
    foundPositives = doStep(pNorm, selectionFdr);

    if (reportPerformanceEachIteration_) {
      int foundTestPositives = 0;
      for (size_t set = 0; set < numFolds_; ++set) {
        foundTestPositives +=
            testScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
      }
      if (VERB > 1) {
        std::cerr << "Found " << foundTestPositives << " test set PSMs with q<"
                  << testFdr_ << endl;
      }
    } else if (VERB > 1) {
      cerr << "Estimated " << foundPositives << " PSMs with q<" << testFdr_
           << endl;
    }

    if (VERB > 2) {
      printAllWeightsColumns(cerr);
    }
    if (VERB > 3) {
      printAllRawWeightsColumns(cerr, pNorm);
    }
    // stop when the estimated positives grew by less than a fraction over the
    // last two iterations; classic training runs all iterations unless
    // quick validation is used
    const double requiredIncrease = quickValidation_
                                        ? requiredIncreaseOver2Iterations_
                                        : requiredIncreaseOver2IterationsDefault_;
    if (stopPerFold) {
      bool allStopped = true;
      for (unsigned int set = 0; set < numFolds_; ++set) {
        if (!isActiveFold(set)) continue;
        const int positives = foldPositives_[set];
        if (positives > 0 && foldPositivesOldOld[set] > 0 &&
            (quickValidation_ || !classicTraining_) &&
            (static_cast<double>(positives - foldPositivesOldOld[set]) <=
             foldPositivesOldOld[set] * requiredIncrease)) {
          stoppedFolds_[set] = 1;
          if (VERB > 1) {
            std::cerr << "Split " << set + 1 << " has converged (" << positives
                      << " vs " << foldPositivesOldOld[set] << ")" << std::endl;
          }
        } else {
          allStopped = false;
        }
        foldPositivesOldOld[set] = foldPositivesOld[set];
        foldPositivesOld[set] = positives;
      }
      if (allStopped) break;
      continue;
    }
    if (foundPositives > 0 && foundPositivesOldOld > 0 &&
        (quickValidation_ || !classicTraining_) &&
        (static_cast<double>(foundPositives - foundPositivesOldOld) <=
         foundPositivesOldOld * requiredIncrease)) {
      if (VERB > 0) {
        std::cerr << "Performance increase over the last two iterations "
                  << "indicate that the algorithm has converged\n"
                  << "(" << foundPositives << " vs " << foundPositivesOldOld
                  << ")" << std::endl;
      }
      break;
    }
    foundPositivesOldOld = foundPositivesOld;
    foundPositivesOld = foundPositives;
  }
  if (VERB == 2) {
    printAllWeightsColumns(cerr);
  }
  if (VERB == 3) {
    printAllRawWeightsColumns(cerr, pNorm);
  }
  checkWeightsFinite();
  foundPositives = 0;
  for (size_t set = 0; set < numFolds_; ++set) {
    foundPositives +=
        testScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
  }
  if (VERB > 0) {
    std::cerr << "Found " << foundPositives << " test set PSMs with q<"
              << testFdr_ << "." << std::endl;
  }
}

void CrossValidation::checkWeightsFinite() const {
  for (size_t set = 0; set < weights_.size(); ++set) {
    for (size_t ix = 0; ix < weights_[set].size(); ++ix) {
      if (std::isfinite(weights_[set][ix])) continue;
      ostringstream oss;
      oss << "ERROR: The SVM weight of "
          << (ix < FeatureNames::getNumFeatures()
                  ? "feature " + DataSet::getFeatureNames().describeFeature(ix)
                  : std::string("the bias term m0"))
          << " in cross validation split " << set + 1 << " is "
          << weights_[set][ix] << " after training." << std::endl;
      Globals::failUnlessNoTerminate(
          oss.str(), "No-terminate flag set: ignoring the error.",
          "Terminating.\n");
      return;
    }
  }
}

unsigned int CrossValidation::foldsHeldAtOnce(unsigned int foldsAtOnce,
                                              unsigned int numFolds) {
  if (foldsAtOnce == 0u) {
    return numFolds > kMaxFoldsAllAtOnce ? 1u : numFolds;
  }
  return (std::min)(foldsAtOnce, numFolds);
}

std::vector<AlgIn>& CrossValidation::svmInputsOfPass(
    bool singleSolves,
    unsigned int& foldsPerPass) {
  const int n = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
  // row pointers need no copies, so all folds are held at once
  foldsPerPass = svmRowsInPlace_ ? numFolds_ : foldsAtOnce();
  assert(foldsPerPass >= 1u);
  if (singleSolvesInPlace_ && singleSolves && foldsPerPass < numFolds_) {
    // Solves that would run one after the other on the copies of a fold at
    // a time: all folds in parallel on the rows in place.
    foldsPerPass = numFolds_;
    if (gridFixed_) {
      std::vector<AlgIn>().swap(svmInputs_);  // no further solves on copies
    }
    if (svmRowLists_.empty()) {
      svmRowLists_.assign(numFolds_ * nestedXvalBins_, AlgIn(n, true));
    }
    return svmRowLists_;
  }
  const std::size_t numInputs = foldsPerPass * nestedXvalBins_;
  if (svmInputs_.size() != numInputs) {
    svmInputs_.assign(numInputs, AlgIn(n, svmRowsInPlace_));
  }
  return svmInputs_;
}

/**
 * Executes a cross validation step
 * @param w_ list of the bins' normal vectors (in linear algebra sense) of the
 *        hyperplane from SVM
 * @return Estimation of number of true positives
 */
int CrossValidation::doStep(const Normalizer* /* pNorm */, double selectionFdr) {
  // Setup
  options pOptions;
  pOptions.lambda = 1.0;
  pOptions.lambda_u = 1.0;
  pOptions.epsilon = tronEpsilon_;  // OpenMS: setTrainingOptions (default EPSILON)
  pOptions.cgitermax = CGITERMAX;
  pOptions.mfnitermax = MFNITERMAX;
  pOptions.directSolve = !classicTraining_;
  // OpenMS: the training uses the weights of the solves only (tron() then skips the outputs and the convergence check)
  pOptions.weightsOnly = true;
  // OpenMS: Newton solves the iterations before the last one of classic training to a looser tolerance (not with
  // quick validation: its training can stop after any iteration)
  if (solver_ == kSolverNewton && classicTraining_ && !quickValidation_ && iteration_ + 1u < niter_) {
    pOptions.epsilon = kNewtonEarlyEpsilon;
  }
  int estTruePos = 0;

  // for determining an appropriate positive training set, the decoys+1 in the
  // FDR estimates is too restrictive for small datasets
  bool skipDecoysPlusOne = true;
  {
    StageTimer stageTimer("rescore training sets");
    for (std::size_t set = 0; set < numFolds_; ++set) {
      if (trainScoresCurrent_) {
        // already scored and sorted with weights_ at the end of the last step
        trainScores_[set].calcQvals(selectionFdr, skipDecoysPlusOne);
      } else {
        trainScores_[set].calcScoresAndQvals(weights_[set], selectionFdr,
                                             skipDecoysPlusOne);
      }
    }
  }

  // Below implements the series of speedups detailed in the following:
  // ////////////////////////////////
  // Speeding Up Percolator
  // John T. Halloran, Hantian Zhang, Kaan Kara, Cédric Renggli, Matthew The, Ce
  // Zhang, David M. Rocke, Lukas Käll, and William Stafford Noble
  //   Journal of Proteome Research 2019 18 (9), 3353-3359
  // ////////////////////////////////
  // Note that the implementation further improves on the speedups in the paper
  // by:
  //   -implementing a single threadpool using OMP for SVM training per each
  //   cpos,cneg pair per nested CV fold -has a much smaller memory footprint by
  //   fixing memory leaks in L2_SVM_MFN and more efficient validation of the
  //   learned SVM parameters
  // ////

  // Create SVM input data for parallelization: the training sets of
  // foldsPerPass folds at a time, whose pairs are then trained in parallel.
  // The folds are built in the same order and every pair is trained on the
  // same rows for any number of folds per pass.
  std::vector<std::vector<Scores> > nestedTestScoresVec(numFolds_);
  const unsigned int pairsPerFold =
      static_cast<unsigned int>(classWeightsPerFold_.size() / numFolds_);
  // whether every training set serves a single solve, as it does once the
  // grid is fixed or with a single candidate pair
  const bool singleSolves = gridFixed_ || pairsPerFold == nestedXvalBins_;
  unsigned int foldsPerPass = 0u;
  std::vector<AlgIn>& svmInputs = svmInputsOfPass(singleSolves, foldsPerPass);
  {
    // the buffers of the threads are released before the grid evaluation
    std::vector<SvmTrainingBuffers> buffers(maxNumThreads());
    for (unsigned int firstSet = 0; firstSet < numFolds_;
         firstSet += foldsPerPass) {
      const unsigned int endSet =
          (std::min)(numFolds_, firstSet + foldsPerPass);
      {
        StageTimer stageTimer("build SVM input");
        for (unsigned int set = firstSet; set < endSet; ++set) {
          // Without nested cross validation, the training set is also used
          // for validation (sub-optimal cross validation).
          std::vector<Scores> nestedTrainScores, nestedTestScores;
          if (!isActiveFold(set)) continue;  // OpenMS: setActiveFolds
          if (nestedXvalBins_ > 1) {
            // OpenMS: with fixed folds, each fold draws its nested sets from a generator state of its own (from the
            // state at preIterationSetup, the iteration and the fold), so the draws of the other folds, whose
            // training sets hold this fold's PSMs, do not reach it
            if (fixedFold_) {
              uint64_t x = nestedSeedBase_ ^ (static_cast<uint64_t>(iteration_) << 8) ^ static_cast<uint64_t>(set);
              x += 0x9E3779B97F4A7C15ULL;
              x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
              x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
              x ^= x >> 31;
              PseudoRandom::setSeed(static_cast<unsigned long>(1u + x % (PseudoRandom::kRandMax - 1u)));
            }
            trainScores_[set].createXvalSetsBySpectrum(
                nestedTrainScores, nestedTestScores, nestedXvalBins_);
          }
          nestedTestScoresVec[set] = std::move(nestedTestScores);
          // Set SVM input data for L2-SVM-MFN
          for (std::size_t nestedFold = 0; nestedFold < nestedXvalBins_;
               ++nestedFold) {
            AlgIn& svmInput =
                svmInputs[(set - firstSet) * nestedXvalBins_ + nestedFold];
            Scores& nestedTrainSet = nestedTrainScores.empty()
                                         ? trainScores_[set]
                                         : nestedTrainScores[nestedFold];
            nestedTrainSet.reserveTrainingSet(svmInput, selectionFdr);
            nestedTrainSet.generateNegativeTrainingSet(svmInput, 1.0);
            nestedTrainSet.generatePositiveTrainingSet(
                svmInput, selectionFdr, 1.0, trainBestPositive_);
            // with nested cross validation, the training set of the whole
            // fold is reported when it is built, in mergeCpCnPairs
            if (VERB > 2) {
              cerr << "Split " << set + 1;
              if (nestedXvalBins_ > 1u) {
                cerr << " nested fold " << nestedFold;
              }
              cerr << ": Training with " << svmInput.positives
                   << " positives and " << svmInput.negatives << " negatives"
                   << std::endl;
            }
          }
        }
      }

      StageTimer stageTimer("SVM training");
      // OpenMS: the TRON solves of the pairs of a training set, which run at the same time, share one copy of the
      // set in the solver's layout (ssl.h), instead of one per solve and thread
      std::vector<std::unique_ptr<TronPairGroup> > tronGroups;
      if (!singleSolves) {
        tronGroups.resize(static_cast<std::size_t>(endSet - firstSet) * nestedXvalBins_);
        for (std::size_t k = 0; k < tronGroups.size(); ++k) {
          tronGroups[k] = std::make_unique<TronPairGroup>(svmInputs[k]);
        }
      }
      // the pairs are ordered by fold
      const int firstPair = static_cast<int>(firstSet * pairsPerFold);
      const int endPair = static_cast<int>(endSet * pairsPerFold);
      LoopErrors pairErrors(static_cast<std::size_t>(endPair - firstPair));  // OpenMS
#pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(maxNumThreads()))
      for (int pairIdx = firstPair; pairIdx < endPair; pairIdx++) {
        try {  // OpenMS: LoopErrors
        CandidateCposCfrac* cpCnFold = &classWeightsPerFold_[pairIdx];
        if (!isSelectedPair(*cpCnFold) || !isActiveFold(cpCnFold->set)) {  // OpenMS: setActiveFolds
          continue;
        }
        const std::size_t input =
            (cpCnFold->set - firstSet) * nestedXvalBins_ +
            static_cast<unsigned int>(cpCnFold->nestedSet);
        trainCpCnPair(*cpCnFold, pOptions, svmInputs[input], buffers[threadNum()],
                      tronGroups.empty() ? nullptr : tronGroups[input].get());
        } catch (...) {
          pairErrors.capture(static_cast<std::size_t>(pairIdx - firstPair));
        }
      }
      pairErrors.rethrowFirst();
    }
  }

  StageTimer stageTimer("grid evaluation");
  estTruePos = mergeCpCnPairs(selectionFdr, pOptions, nestedTestScoresVec,
                              candidatesCpos_, candidatesCfrac_);
  return estTruePos;
}

bool CrossValidation::trainsFoldChains() const {
#if PERCOLATOR_OMP_TASKS
  return maxNumThreads() > 1u && numFolds_ > 1u && niter_ > 0u && nestedXvalBins_ == 1u &&
         !reportPerformanceEachIteration_ && VERB <= 2 && activeFolds_.empty() &&
         (svmRowsInPlace_ || foldsAtOnce() == numFolds_);
#else
  return false;
#endif
}

/**
 * OpenMS: the iterations of train() as one chain of tasks per fold. A step of a fold is what doStep and
 * mergeCpCnPairs do for it, in the same order: its training set scored with its weights (rescored in the first
 * step), its SVM input, the solves of its (Cpos, Cneg) pairs (one task each, with the count of the pair's positives
 * right after its solve), the choice of the best pair and the count of the fold's positives with its weights. The
 * last pair of a step to finish runs the rest of it and starts the fold's next step. A fold thus waits only for its
 * own pairs, not for the other folds at each stage, and the threads run whatever task is ready. Only the stop rule of
 * the training on the sum of the folds (new training without fixed folds) waits for all folds of an iteration: the
 * last fold to finish it decides and starts the next steps. With fixed folds each fold stops on its own count, and
 * classic training runs all iterations, so the chains never wait for each other.
 *
 * The results are those of the loop of train(): a step of a fold reads and writes only the fold's own data, and its
 * computations do not depend on the threads (the scores of a sorted set, its order and its q-values are the same on
 * one thread). The reports of the iterations are printed afterwards as the loop prints them. An error ends its fold's
 * chain; the first error (by iteration, then stage, then fold) is rethrown.
 */
void CrossValidation::trainFoldChains() {
  const unsigned int numPairsPerFold =
      static_cast<unsigned int>(classWeightsPerFold_.size() / numFolds_);
  const bool stopPerFold = static_cast<bool>(fixedFold_);
  const bool stopsEarly = quickValidation_ || !classicTraining_;
  // the stop rule on the sum of the folds (train)
  const bool coupled = !stopPerFold && stopsEarly;
  const double requiredIncrease = quickValidation_ ? requiredIncreaseOver2Iterations_
                                                   : requiredIncreaseOver2IterationsDefault_;
  const bool canFixGrid = !classicTraining_ && classWeightsPerFold_.size() > numFolds_ * nestedXvalBins_;
  const bool skipDecoysPlusOne = true;  // as doStep and mergeCpCnPairs
  const int numValues = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
  if (svmInputs_.size() != numFolds_) {
    svmInputs_.assign(numFolds_, AlgIn(numValues, svmRowsInPlace_));
  }
  const std::size_t teamThreads = maxNumThreads();

  struct Fold {
    unsigned int iteration = 0;  // the step that runs
    options stepOptions;
    bool needCounts = false;
    bool fixed = false;  // the grid of the fold is fixed (after gridSearchIterations_ steps)
    std::pair<double, double> selected{1.0, 1.0};
    std::unique_ptr<TronPairGroup> group;
    Scores::FeatureRows rows;
    std::vector<int> pairPositives;
    std::vector<std::exception_ptr> pairErrors;
    std::atomic<unsigned int> pendingPairs{0u};
    // per iteration: the positives of the fold (-1: no step), the best pair
    std::vector<int> positives;
    std::vector<std::pair<double, double> > best;
    int positivesOld = 0, positivesOldOld = 0;
    int stoppedAt = -1;  // the iteration after which the fold stopped (fixed folds)
    std::exception_ptr error;
    unsigned int errorIteration = 0;
    int errorStage = 0;
  };
  std::vector<Fold> folds(numFolds_);
  for (Fold& fold : folds) {
    fold.positives.assign(niter_, -1);
    fold.best.assign(niter_, std::make_pair(1.0, 1.0));
  }
  std::vector<SvmTrainingBuffers> svmBuffers(teamThreads);
  std::vector<Scores::CountBuffers> countBuffers(teamThreads);
  // the stop rule on the sum of the folds
  std::mutex decisionMutex;
  std::vector<unsigned int> foldsDone(niter_, 0u);
  std::vector<int> found(niter_, 0);
  int foundOld = 0, foundOldOld = 0;
  int stoppedAfter = -1;  // the iteration after which the rule stopped the training

  auto failed = [&](Fold& fold, int stage) {
    fold.error = std::current_exception();
    fold.errorIteration = fold.iteration;
    fold.errorStage = stage;
  };
  std::function<void(unsigned int)> startStep;
  std::function<void(unsigned int)> finishStep;

  // the second part of a step: after the solves of the fold's pairs
  finishStep = [&](unsigned int set) {
    const TaskLoopScope taskLoops;  // its scoring of the training set as a taskloop
    Fold& fold = folds[set];
    const unsigned int i = fold.iteration;
    int stage = 3;
    try {
      fold.group.reset();
      for (const std::exception_ptr& error : fold.pairErrors) {
        if (error) {
          stage = 2;
          std::rethrow_exception(error);
        }
      }
      // mergeCpCnPairs without nested cross validation, for this fold
      int bestTruePoses = 0;
      double bestCpos = 1.0, bestCfrac = 1.0;
      for (unsigned int k = 0; k < numPairsPerFold; ++k) {
        CandidateCposCfrac& cpCnFold = classWeightsPerFold_[set * numPairsPerFold + k];
        const int tp = fold.pairPositives[k];
        cpCnFold.tp = tp;
        if (tp >= bestTruePoses) {
          bestTruePoses = tp;
          weights_[set] = cpCnFold.ww;
          bestCpos = cpCnFold.cpos;
          bestCfrac = cpCnFold.cfrac;
        }
      }
      fold.best[i] = std::make_pair(bestCpos, bestCfrac);
      if (canFixGrid && !fold.fixed && i + 1u >= gridSearchIterations_) {
        fold.fixed = true;
        fold.selected = fold.best[i];
      }
      stage = 4;
      const int positives = trainScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
      foldPositives_[set] = positives;
      fold.positives[i] = positives;
    } catch (...) {
      failed(fold, stage);
      if (!coupled) return;
    }
    if (!coupled) {
      bool stop = i + 1u >= niter_;
      if (stopPerFold) {
        const int positives = fold.positives[i];
        if (positives > 0 && fold.positivesOldOld > 0 && stopsEarly &&
            (static_cast<double>(positives - fold.positivesOldOld) <=
             fold.positivesOldOld * requiredIncrease)) {
          stoppedFolds_[set] = 1;
          fold.stoppedAt = static_cast<int>(i);
          stop = true;
        }
        fold.positivesOldOld = fold.positivesOld;
        fold.positivesOld = positives;
      }
      if (!stop) {
        fold.iteration = i + 1u;
        PERCOLATOR_TASK_PRAGMA(omp task default(shared) firstprivate(set))
        startStep(set);
      }
      return;
    }
    // the stop rule of train() on the sum of the folds, by the last fold to finish the iteration
    bool next = false;
    {
      std::lock_guard<std::mutex> lock(decisionMutex);
      if (++foldsDone[i] < numFolds_) return;
      bool anyError = false;
      double positives = 0.0;
      for (const Fold& other : folds) {
        anyError = anyError || static_cast<bool>(other.error);
        positives += other.positives[i];
      }
      if (anyError) return;
      // each PSM is in numFolds_ - 1 training sets (doStep)
      found[i] = static_cast<int>(positives / std::max(1u, numFolds_ - 1u));
      if (found[i] > 0 && foundOldOld > 0 &&
          (static_cast<double>(found[i] - foundOldOld) <= foundOldOld * requiredIncrease)) {
        stoppedAfter = static_cast<int>(i);
      } else {
        foundOldOld = foundOld;
        foundOld = found[i];
        next = i + 1u < niter_;
      }
    }
    if (next) {
      for (unsigned int other = 0; other < numFolds_; ++other) {
        folds[other].iteration = i + 1u;
        PERCOLATOR_TASK_PRAGMA(omp task default(shared) firstprivate(other))
        startStep(other);
      }
    }
  };

  // the first part of a step: the training set, the SVM input and the solves of the pairs (as tasks)
  startStep = [&](unsigned int set) {
    const TaskLoopScope taskLoops;  // its scoring and the copies of the training rows as taskloops
    Fold& fold = folds[set];
    const unsigned int i = fold.iteration;
    int stage = 0;
    std::vector<unsigned int> pairs;
    try {
      const double selectionFdr = i == 0u ? initialSelectionFdr_ : selectionFdr_;
      // doStep's options for the step
      options& pOptions = fold.stepOptions;
      pOptions.lambda = 1.0;
      pOptions.lambda_u = 1.0;
      pOptions.epsilon = tronEpsilon_;
      pOptions.cgitermax = CGITERMAX;
      pOptions.mfnitermax = MFNITERMAX;
      pOptions.directSolve = !classicTraining_;
      pOptions.weightsOnly = true;
      if (solver_ == kSolverNewton && classicTraining_ && !quickValidation_ && i + 1u < niter_) {
        pOptions.epsilon = kNewtonEarlyEpsilon;
      }
      if (i > 0u || trainScoresCurrent_) {
        // scored and sorted with weights_ at the end of the fold's last step
        trainScores_[set].calcQvals(selectionFdr, skipDecoysPlusOne);
      } else {
        trainScores_[set].calcScoresAndQvals(weights_[set], selectionFdr, skipDecoysPlusOne);
      }
      stage = 1;
      AlgIn& svmInput = svmInputs_[set];
      trainScores_[set].reserveTrainingSet(svmInput, selectionFdr);
      trainScores_[set].generateNegativeTrainingSet(svmInput, 1.0);
      trainScores_[set].generatePositiveTrainingSet(svmInput, selectionFdr, 1.0, trainBestPositive_);
      fold.needCounts = !fold.fixed && numPairsPerFold > nestedXvalBins_;
      if (fold.needCounts) {
        trainScores_[set].getFeatureRows(fold.rows);
      }
      const bool singleSolves = fold.fixed || numPairsPerFold == nestedXvalBins_;
      if (!singleSolves) {
        fold.group = std::make_unique<TronPairGroup>(svmInput);
      }
      stage = 2;
      // the pairs of the step (isSelectedPair); a pair that is not trained counts -1 (mergeCpCnPairs)
      fold.pairPositives.assign(numPairsPerFold, -1);
      fold.pairErrors.assign(numPairsPerFold, nullptr);
      for (unsigned int k = 0; k < numPairsPerFold; ++k) {
        const CandidateCposCfrac& cpCnFold = classWeightsPerFold_[set * numPairsPerFold + k];
        if (!fold.fixed || (cpCnFold.cpos == fold.selected.first && cpCnFold.cfrac == fold.selected.second)) {
          pairs.push_back(k);
        }
      }
    } catch (...) {
      failed(fold, stage);
      fold.group.reset();
      if (coupled) {
        std::lock_guard<std::mutex> lock(decisionMutex);
        ++foldsDone[i];
      }
      return;
    }
    if (pairs.empty()) {  // (not reached: a fixed grid keeps one of the fold's pairs)
      finishStep(set);
      return;
    }
    fold.pendingPairs.store(static_cast<unsigned int>(pairs.size()), std::memory_order_relaxed);
    for (unsigned int k : pairs) {
      PERCOLATOR_TASK_PRAGMA(omp task default(shared) firstprivate(set, k))
      {
        Fold& pairFold = folds[set];
        try {
          CandidateCposCfrac& cpCnFold = classWeightsPerFold_[set * numPairsPerFold + k];
          const std::size_t thread = threadNum();
          trainCpCnPair(cpCnFold, pairFold.stepOptions, svmInputs_[set], svmBuffers[thread],
                        pairFold.group.get());
          pairFold.pairPositives[k] =
              pairFold.needCounts
                  ? trainScores_[set].countTargetsBelowFdr(pairFold.rows, cpCnFold.ww, testFdr_,
                                                           skipDecoysPlusOne, countBuffers[thread])
                  : 0;  // the only trained pair of its fold
        } catch (...) {
          pairFold.pairErrors[k] = std::current_exception();
        }
        if (pairFold.pendingPairs.fetch_sub(1u, std::memory_order_acq_rel) == 1u) {
          finishStep(set);
        }
      }
    }
  };

#pragma omp parallel num_threads(static_cast<int>(teamThreads))
#pragma omp single
  {
    for (unsigned int set = 0; set < numFolds_; ++set) {
      PERCOLATOR_TASK_PRAGMA(omp task default(shared) firstprivate(set))
      startStep(set);
    }
  }

  // the steps that ran: all folds for the iterations up to the last one, and the state of the loop after them
  unsigned int iterations = 0;
  for (const Fold& fold : folds) {
    for (unsigned int i = 0; i < niter_; ++i) {
      if (fold.positives[i] >= 0) iterations = std::max(iterations, i + 1u);
    }
  }
  iteration_ = iterations;
  if (canFixGrid && iterations >= gridSearchIterations_) {
    gridFixed_ = true;
    selectedPairs_.resize(numFolds_);
    for (unsigned int set = 0; set < numFolds_; ++set) {
      selectedPairs_[set] = folds[set].best[gridSearchIterations_ - 1u];
    }
  }
  trainScoresCurrent_ = true;
  // a fold that stopped before the last iteration had its q-values recalculated in the iterations after (doStep)
  for (unsigned int set = 0; set < numFolds_; ++set) {
    if (folds[set].stoppedAt >= 0 && static_cast<unsigned int>(folds[set].stoppedAt) + 1u < iterations) {
      trainScores_[set].calcQvals(selectionFdr_, skipDecoysPlusOne);
    }
  }

  // the reports of the iterations, as train() prints them
  for (unsigned int i = 0; i < iterations; ++i) {
    if (VERB > 1) {
      cerr << "Iteration " << i + 1 << ":\t";
    }
    if (canFixGrid && i + 1u == gridSearchIterations_ && VERB > 1) {
      for (unsigned int set = 0; set < numFolds_; ++set) {
        const std::pair<double, double>& best = folds[set].best[i];
        std::cerr << "Split " << set + 1 << ": using Cpos=" << best.first
                  << ", Cneg=" << best.first * best.second
                  << " in the remaining iterations" << std::endl;
      }
    }
    double positives = 0.0;
    for (const Fold& fold : folds) {
      if (fold.positives[i] >= 0) positives += fold.positives[i];
    }
    const int foundPositives = static_cast<int>(positives / std::max(1u, numFolds_ - 1u));
    if (VERB > 1) {
      cerr << "Estimated " << foundPositives << " PSMs with q<" << testFdr_ << endl;
    }
    if (stopPerFold) {
      for (unsigned int set = 0; set < numFolds_; ++set) {
        if (folds[set].stoppedAt == static_cast<int>(i) && VERB > 1) {
          const int before = i >= 2u ? folds[set].positives[i - 2u] : 0;
          std::cerr << "Split " << set + 1 << " has converged (" << folds[set].positives[i]
                    << " vs " << before << ")" << std::endl;
        }
      }
    } else if (stoppedAfter == static_cast<int>(i) && VERB > 0) {
      std::cerr << "Performance increase over the last two iterations "
                << "indicate that the algorithm has converged\n"
                << "(" << found[i] << " vs " << (i >= 2u ? found[i - 2u] : 0) << ")" << std::endl;
    }
  }

  // the first error by iteration, stage and fold
  const Fold* first = nullptr;
  for (const Fold& fold : folds) {
    if (fold.error && (!first || fold.errorIteration < first->errorIteration ||
                       (fold.errorIteration == first->errorIteration && fold.errorStage < first->errorStage))) {
      first = &fold;
    }
  }
  if (first) {
    std::rethrow_exception(first->error);
  }
}

/**
 * Train SVM over a single (cpos, cneg) pair
 * @param cpCnFold contains cpos, cneg pair and SVM learned weights
 * @param pOptions options for the SVM algorithm
 * @param svmInput training data for this particular nested CV fold
 * @param buffers buffers of the calling thread
 */
void CrossValidation::trainCpCnPair(CandidateCposCfrac& cpCnFold,
                                    const options& pOptions,
                                    const AlgIn& svmInput,
                                    SvmTrainingBuffers& buffers,
                                    TronPairGroup* group) {
  double cpos = cpCnFold.cpos;
  double cfrac = cpCnFold.cfrac;

  if (VERB > 3)
    cerr << "- cross-validation with Cpos=" << cpos << ", Cneg=" << cfrac * cpos
         << endl;
  // The SVM algorithm (see ssl.cpp) starts from the weights this pair had
  // after the previous iteration (zero in the first one), which are close
  // to the new solution, or from zero with classic training.
  // OpenMS: TRON (ssl.h) starts from the weights alone, so the outputs X w of upstream's solver are not computed.
  // With classic training, TRON starts from the previous weights as well with warmStart_ (and no nested
  // cross-validation); Newton and LDA (setSolver) always do.
  std::vector<double>& weights = cpCnFold.ww;
  if (classicTraining_ && solver_ == kSolverTron && !(warmStart_ && nestedXvalBins_ == 1u)) {
    std::fill(weights.begin(), weights.end(), 0.0);
  }
  trainSvm(svmInput, pOptions, cpos, cfrac * cpos, weights, buffers.outputs, group);
}

// OpenMS: the solver of setSolver. Newton (newton_svm.h) and LDA start from the weights in w; LDA needs two rows of
// each class, Newton at most kNewtonMaxColumns columns: otherwise, and if Newton fails or does not converge, TRON
// solves the SVM from 0. TRON (kSolverTron) starts from w as given.
void CrossValidation::trainSvm(const AlgIn& svmInput, const options& pOptions,
                               double cpos, double cneg, std::vector<double>& w,
                               std::vector<double>& outputs, TronPairGroup* group) const {
  if (solver_ == kSolverLda && ldaSvm(svmInput, pOptions, cpos, cneg, w.data())) {
    return;
  }
  if (solver_ != kSolverTron) {
    if (newtonSvm(svmInput, pOptions, cpos, cneg, pOptions.epsilon, w.data(), nullptr) == 1) {
      return;
    }
    std::fill(w.begin(), w.end(), 0.0);  // the fallback: TRON from 0
  }
  L2_SVM_MFN(svmInput, pOptions, w, outputs, cpos, cneg, group);
}

/**
 * Validate and merge weights learned per cpos,cneg pairs per nested CV fold per
 * CV fold
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
 * @param nestedTestScoresVec 2D vector, test sets per nested CV fold per CV
 * fold
 */
int CrossValidation::mergeCpCnPairs(
    double selectionFdr,
    const options& pOptions,
    vector<vector<Scores> >& nestedTestScoresVec,
    const vector<double>& cposCandidates,
    const vector<double>& cfracCandidates) {
  // for determining the number of positives, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true;

  vector<int> bestTruePoses(numFolds_, 0);
  vector<double> bestCposes(numFolds_, 1);
  vector<double> bestCfracs(numFolds_, 1);
  
  // Validate learned parameters per (cpos,cneg) pair per nested CV fold
  // Note: this cannot be done in trainCpCnPair without setting a critical
  // pragma, due to the
  //       scoring calculation in calcScores.
  unsigned int numCpCnPairsPerSet =
      static_cast<unsigned int>(classWeightsPerFold_.size() / numFolds_);

  // The validation counts of all pairs are independent, so they are computed
  // in parallel. Without nested cross validation the training sets are the
  // validation sets.
  auto testSet = [&](unsigned int set,
                     unsigned int nestedSet) -> const Scores& {
    const std::vector<Scores>& nestedTestScores = nestedTestScoresVec[set];
    return nestedTestScores.empty() ? trainScores_[set]
                                    : nestedTestScores[nestedSet];
  };
  // When a fold trains only one pair per nested fold, as it does once the
  // grid is fixed or with a single candidate pair, the selection below
  // picks that pair for any count >= 0, with and without nested cross
  // validation. The counts are then not computed, and the pair gets 0.
  const bool needCounts = !gridFixed_ && numCpCnPairsPerSet > nestedXvalBins_;
  const int numTestSets =
      needCounts ? static_cast<int>(numFolds_ * nestedXvalBins_) : 0;
  std::vector<Scores::FeatureRows> testRows(numTestSets);
  LoopErrors rowErrors(static_cast<std::size_t>(numTestSets));  // OpenMS
#pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(maxNumThreads()))
  for (int idx = 0; idx < numTestSets; ++idx) {
    try {  // OpenMS: LoopErrors
    if (!isActiveFold(idx / nestedXvalBins_)) continue;  // OpenMS: setActiveFolds
    testSet(idx / nestedXvalBins_, idx % nestedXvalBins_)
        .getFeatureRows(testRows[idx]);
    } catch (...) {
      rowErrors.capture(static_cast<std::size_t>(idx));
    }
  }
  rowErrors.rethrowFirst();

  std::vector<int> truePositives(classWeightsPerFold_.size());
  LoopErrors countErrors(classWeightsPerFold_.size());  // OpenMS
#pragma omp parallel num_threads(static_cast<int>(maxNumThreads()))
  {
    Scores::CountBuffers buffers;  // one per thread, for all its pairs
#pragma omp for schedule(dynamic, 1)
    for (int pairIdx = 0;
         pairIdx < static_cast<int>(classWeightsPerFold_.size()); ++pairIdx) {
      try {  // OpenMS: LoopErrors
      const CandidateCposCfrac& cpCnPair = classWeightsPerFold_[pairIdx];
      const unsigned int nestedSet =
          static_cast<unsigned int>(cpCnPair.nestedSet);
      if (!isSelectedPair(cpCnPair) || !isActiveFold(cpCnPair.set)) {  // OpenMS: setActiveFolds
        truePositives[pairIdx] = -1;  // not trained in this iteration
        continue;
      }
      if (!needCounts) {
        truePositives[pairIdx] = 0;  // the only trained pair of its fold
        continue;
      }
      truePositives[pairIdx] =
          testSet(cpCnPair.set, nestedSet)
              .countTargetsBelowFdr(
                  testRows[cpCnPair.set * nestedXvalBins_ + nestedSet],
                  cpCnPair.ww, testFdr_, skipDecoysPlusOne, buffers);
      } catch (...) {
        countErrors.capture(static_cast<std::size_t>(pairIdx));
      }
    }
  }
  countErrors.rethrowFirst();

  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    if (!isActiveFold(static_cast<unsigned int>(set))) continue;  // OpenMS: keeps its weights
    unsigned int a = set * numCpCnPairsPerSet;
    unsigned int b = (set + 1) * numCpCnPairsPerSet;
    int tp = 0;
    std::vector<CandidateCposCfrac>::iterator itCpCnPair;
    std::map<std::pair<double, double>, int> intermediateResults;
    for (itCpCnPair = classWeightsPerFold_.begin() + a;
         itCpCnPair < classWeightsPerFold_.begin() + b; itCpCnPair++) {
      tp = truePositives[itCpCnPair - classWeightsPerFold_.begin()];
      intermediateResults[std::pair<double,double>{itCpCnPair->cpos,
                                                   itCpCnPair->cfrac}] += tp;
      itCpCnPair->tp = tp;
      if (nestedXvalBins_ <= 1) {
        if (tp >= bestTruePoses[set]) {
          bestTruePoses[set] = tp;
          weights_[set] = itCpCnPair->ww;
          bestCposes[set] = itCpCnPair->cpos;
          bestCfracs[set] = itCpCnPair->cfrac;
        }
      }
    }
    if (nestedXvalBins_ > 1) {  // Check nestedXvalBins, which collapse
                                // (accumulate) tp estimated for each CV bin
      // Now check which achieved best performance among cpos, cneg pairs
      std::vector<double>::const_iterator itCpos = cposCandidates.begin();
      for (; itCpos != cposCandidates.end(); ++itCpos) {
        double cpos = *itCpos;
        // OpenMS: the fold's own Cneg candidates with fixed folds (initializeGridSearch)
        const vector<double>& setCfracCandidates =
            candidatesCfracPerFold_.empty() ? cfracCandidates : candidatesCfracPerFold_[set];
        std::vector<double>::const_iterator itCfrac = setCfracCandidates.begin();
        for (; itCfrac != setCfracCandidates.end(); ++itCfrac) {
          double cfrac = *itCfrac;
          tp = intermediateResults[std::pair<double,double>{cpos, cfrac}];
          if (tp >= bestTruePoses[set]) {
            bestTruePoses[set] = tp;
            bestCposes[set] = cpos;
            bestCfracs[set] = cfrac;
          }
        }
      }
    }
  }

  if (nestedXvalBins_ > 1) {
    // one solve per fold on its whole training set, for foldsPerPass folds
    // at a time as in doStep
    unsigned int foldsPerPass = 0u;
    std::vector<AlgIn>& svmInputs = svmInputsOfPass(true, foldsPerPass);
    for (unsigned int firstSet = 0; firstSet < numFolds_;
         firstSet += foldsPerPass) {
      const unsigned int endSet =
          (std::min)(numFolds_, firstSet + foldsPerPass);
      LoopErrors refitErrors(numFolds_);  // OpenMS
#pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(maxNumThreads()))
      for (int set = static_cast<int>(firstSet);
           set < static_cast<int>(endSet); ++set) {
        try {  // OpenMS: LoopErrors
        if (!isActiveFold(static_cast<unsigned int>(set))) continue;  // OpenMS: setActiveFolds
        AlgIn& svmInput = svmInputs[
            (static_cast<unsigned int>(set) - firstSet) * nestedXvalBins_];
        trainScores_[set].reserveTrainingSet(svmInput, selectionFdr);
        trainScores_[set].generateNegativeTrainingSet(svmInput, 1.0);
        trainScores_[set].generatePositiveTrainingSet(
            svmInput, selectionFdr, 1.0, trainBestPositive_);

        // The SVM algorithm (see ssl.cpp) starts from the fold's weights of
        // the previous iteration, which are close to the new solution, as
        // the grid solves do; in the first iteration, which holds the
        // initial direction, and with classic training it starts from zero.
        // OpenMS: TRON (ssl.h) starts from the weights alone (no outputs X w); the solver of setSolver
        std::vector<double> outputs;
        if (classicTraining_ || iteration_ == 0u) {
          std::fill(weights_[set].begin(), weights_[set].end(), 0.0);
        }
        trainSvm(svmInput, pOptions, bestCposes[set],
                 bestCposes[set] * bestCfracs[set], weights_[set], outputs);
        } catch (...) {
          refitErrors.capture(static_cast<std::size_t>(set));
        }
      }
      refitErrors.rethrowFirst();
      if (VERB > 2) {
        for (unsigned int set = firstSet; set < endSet; ++set) {
          const AlgIn& svmInput =
              svmInputs[(set - firstSet) * nestedXvalBins_];
          cerr << "Split " << set + 1 << ": Training with "
               << svmInput.positives << " positives and "
               << svmInput.negatives << " negatives" << std::endl;
        }
      }
    }
  }

  // After the first iterations the selected pairs hardly change, and their
  // number of positives differs from the other pairs' only by noise. From
  // then on only the selected pair of each fold is trained.
  ++iteration_;
  if (!classicTraining_ && !gridFixed_ && iteration_ >= gridSearchIterations_ &&
      classWeightsPerFold_.size() > numFolds_ * nestedXvalBins_) {
    gridFixed_ = true;
    selectedPairs_.resize(numFolds_);
    for (unsigned int set = 0; set < numFolds_; ++set) {
      selectedPairs_[set] = std::make_pair(bestCposes[set], bestCfracs[set]);
      if (VERB > 1) {
        std::cerr << "Split " << set + 1 << ": using Cpos=" << bestCposes[set]
                  << ", Cneg=" << bestCposes[set] * bestCfracs[set]
                  << " in the remaining iterations" << std::endl;
      }
    }
  }

  double bestTruePos = 0;
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    if (!isActiveFold(static_cast<unsigned int>(set))) continue;  // OpenMS: the counts of the active folds
    const int positives = trainScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
    if (static_cast<std::size_t>(set) < foldPositives_.size()) {
      foldPositives_[static_cast<std::size_t>(set)] = positives;  // OpenMS: see train
    }
    bestTruePos += positives;
  }
  trainScoresCurrent_ = true;
  // each PSM is in numFolds_ - 1 training sets, or in the only one
  return static_cast<int>(bestTruePos / std::max(1u, numFolds_ - 1));
}

/**
 * Scores the test sets with the trained weights and merges them into
 * fullset. Releases the training data, so train() cannot be called again.
 */
void CrossValidation::postIterationProcessing(Scores& fullset,
                                              SanityCheck& check) {
  // The training is done. Its data are released before the merge, which
  // copies all PSMs once more and would otherwise set the peak memory.
  std::vector<AlgIn>().swap(svmInputs_);
  std::vector<AlgIn>().swap(svmRowLists_);
  if (decoyFractionTraining_ < 1.0) {
    // with RESET 1 decoy: decoysPerTarget = 1 and decoyFractionTraining = 0.5,
    // then factor = 1.5 and testNullTargetWinProb = 2/3 (and decoyFactor = 2)
    //
    // with RESET 2 decoys: decoysPerTarget = 2 and decoyFractionTraining = 0.5,
    // then factor = 2 and testNullTargetWinProb = 1/2 (and decoyFactor = 1)
    unsigned int decoysPerTarget =
        1u;  // TODO: make this a command line argument
    double factor =
        (decoysPerTarget + 1) - decoyFractionTraining_ * decoysPerTarget;
    double testNullTargetWinProb = 1.0 / factor;

    for (std::size_t set = 0; set < numFolds_; ++set) {
      Scores newTestScores(false);
      for (ScoreHolder& sh : testScores_[set]) {
        if (sh.label == LabelType::DECOY)
          continue;

        if (sh.label == LabelType::PSEUDO_TARGET) {
          sh.label = LabelType::DECOY;
        }
        newTestScores.addScoreHolder(sh);
      }
      newTestScores.setNullTargetWinProb(testNullTargetWinProb);
      newTestScores.recalculateSizes();
      testScores_[set] = std::move(newTestScores);
    }
    fullset.setNullTargetWinProb(testNullTargetWinProb);
  }

  directionOutcome_ =
      check.validateDirection(testScores_, trainScores_, weights_);
  const bool directionValid =
      directionOutcome_ != SanityCheck::DirectionOutcome::FAILED;
  // When the second chance of the direction check kept the trained SVMs,
  // whose test bins were too small for the per-bin count, the bins are
  // merged as it counted them, with the training-set decoys, unless the
  // caller chose another scheme for this case, whose list the second chance
  // has then checked as well (MergePlan::usedScheme).
  usedMergeNormalization_ = mergePlan_.usedScheme(
      directionOutcome_ ==
      SanityCheck::DirectionOutcome::KEPT_BY_SECOND_CHANCE);
  if (!directionValid) {
    // Each test bin is scored with its own fold's vector: after a reset the
    // vector the fallback chose for that fold (under --init-weights the given
    // vector or the fold's default direction, whichever gives more targets
    // on the fold's training set; otherwise the fold's default vector, and
    // with --static the given vector in every fold); with --override the
    // trained vector, since the other folds' vectors were trained on the
    // bin. The weights written with -w then describe the scores, but for
    // the exception that follows.
    //
    // Except with -N under test, after a reset and with --override: -N then
    // scores all PSMs with the averaged weights (see below), so the merged
    // scores of the bins are not reported, and the merge rescales each
    // fold's vector by the anchors of a bin before they are averaged. A bin
    // with no target at q < fdr is anchored at its best score, so with its
    // own fold's vector a real score with a long tail got a smaller factor
    // than another fold's noise feature, which then dominated the average.
    // All bins are scored with the first fold's vector instead, as before
    // the folds kept their own vectors, so that the folds' vectors get
    // nearly the same factor and count alike in the average.
    const bool averagedUnderTest =
        recordTestFolds_ &&
        usedMergeNormalization_ == MergeNormalization::TEST;
    for (std::size_t set = 0; set < numFolds_; ++set) {
      testScores_[set].calcScoresAndQvals(
          averagedUnderTest ? weights_[0] : weights_[set], selectionFdr_);
    }
  }

  // With the TRAIN and TRAIN_DECOYS schemes, each test bin is rescaled by
  // anchors taken from its fold's training set, scored with the same final
  // weights as the test bin. The training set holds (numFolds_ - 1) times
  // as many PSMs as the test bin, or the same PSMs with one fold, so the
  // anchors exist and are less noisy for small test bins.
  const std::vector<ScoreScale> binScales = Scores::calcBinScales(
      trainScores_, weights_, usedMergeNormalization_, mergePlan_.fdr,
      mergePlan_.skipNormalizeScores);
  // the default direction of validateDirection and the anchors above were
  // the last use of the training sets; merge empties the test sets
  std::vector<Scores>().swap(trainScores_);
  // After a reset the weights of each fold are a fallback picked per fold,
  // e.g. the best single feature of its training set, which can be a noise
  // feature in one fold and the real score in another; all PSMs are then
  // scored with the averaged weights, as without the test folds. So they
  // are after a failed check with --override, whose trained vectors the
  // check did not accept.
  if (recordTestFolds_ && directionValid) recordTestFolds();
  fullset.merge(testScores_, mergePlan_.fdr, mergePlan_.skipNormalizeScores,
                weights_, usedMergeNormalization_, binScales);
}

MergeNormalization CrossValidation::defaultMergeNormalization(
    unsigned int numFolds) {
  return numFolds <= 3u ? MergeNormalization::TEST
                        : MergeNormalization::TRAIN_DECOYS;
}

/**
 * Prints the weights of the normalized vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printSetWeights(ostream& weightStream, unsigned int set) {
  weightStream << DataSet::getFeatureNames().getFeatureNames() << "\tm0"
               << std::endl;
  weightStream.precision(3);
  weightStream << weights_[set][0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << weights_[set][ix];
  }
  weightStream << endl;
}

/**
 * Prints the weights of the raw vector to a stream
 * @param weightStream stream where the weights are written to
 * @param w_ normal vector
 */
void CrossValidation::printRawSetWeights(ostream& weightStream,
                                         unsigned int set,
                                         const Normalizer* pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);
  pNorm->unnormalizeweight(weights_[set], ww);
  weightStream << ww[0];
  for (unsigned int ix = 1; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    weightStream << "\t" << fixed << setprecision(4) << ww[ix];
  }
  weightStream << endl;
}

// used to save weights for reuse as weight input file
void CrossValidation::printAllWeights(ostream& weightStream,
                                      const Normalizer* pNorm) {
  weightStream << "# This file contains the weights from each cross validation "
                  "bin from percolator training"
               << endl;
  weightStream << "# First line is the feature names, followed by normalized "
                  "weights, and the raw weights of bin 1"
               << endl;
  weightStream << "# This is repeated for the other bins" << endl;
  for (unsigned int ix = 0; ix < numFolds_; ++ix) {
    printSetWeights(weightStream, ix);
    printRawSetWeights(weightStream, ix, pNorm);
  }
}

void CrossValidation::getAvgWeights(std::vector<double>& weights,
                                    const Normalizer* pNorm) {
  vector<double> ww(FeatureNames::getNumFeatures() + 1);

  weights.resize(FeatureNames::getNumFeatures() + 1);
  for (size_t set = 0; set < weights_.size(); ++set) {
    pNorm->unnormalizeweight(weights_[set], ww);
    for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
      weights[ix] += ww[ix] / static_cast<double>(weights_.size());
    }
  }
}

void CrossValidation::getRawFoldWeights(
    std::vector<std::vector<double> >& weights,
    const Normalizer* pNorm) {
  weights.assign(weights_.size(),
                 vector<double>(FeatureNames::getNumFeatures() + 1));
  for (size_t set = 0; set < weights_.size(); ++set) {
    pNorm->unnormalizeweight(weights_[set], weights[set]);
  }
}

void CrossValidation::getTrainingPsms(
    std::size_t set, std::vector<const PSMDescription*>& psms) const {
  psms.clear();
  if (set >= trainScores_.size()) return;
  psms.reserve(trainScores_[set].size());
  for (const ScoreHolder& sh : trainScores_[set]) {
    psms.push_back(sh.pPSM);
  }
}

void CrossValidation::getTestFolds(
    std::vector<std::pair<unsigned int, unsigned int> >& lineFolds) {
  lineFolds.clear();
  lineFolds.swap(testFolds_);
}

void CrossValidation::recordTestFolds() {
  testFolds_.clear();
  size_t numPsms = 0;
  for (const Scores& testSet : testScores_) {
    numPsms += testSet.size();
  }
  testFolds_.reserve(numPsms);
  for (unsigned int set = 0; set < testScores_.size(); ++set) {
    for (const ScoreHolder& sh : testScores_[set]) {
      testFolds_.emplace_back(sh.pPSM->getLineNr(), set);
    }
  }
  std::sort(testFolds_.begin(), testFolds_.end());
}

void CrossValidation::printAllWeightsColumns(
    std::vector<std::vector<double> > weightMatrix,
    ostream& outputStream) {
  // write to intermediate stream to prevent the fixed precision from sticking
  ostringstream weightStream;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    weightStream << " Split" << set + 1 << '\t';  // right-align with weights
  }
  weightStream << "FeatureName" << std::endl;
  size_t numRows = FeatureNames::getNumFeatures() + 1;
  for (unsigned int ix = 0; ix < numRows; ix++) {
    for (unsigned int set = 0; set < numFolds_; ++set) {
      // align positive and negative weights
      if (weightMatrix[set][ix] >= 0)
        weightStream << " ";
      weightStream << fixed << setprecision(4) << weightMatrix[set][ix];
      weightStream << "\t";
    }
    if (ix == numRows - 1) {
      weightStream << "m0";
    } else {
      weightStream << DataSet::getFeatureNames().getFeatureName(ix);
    }
    weightStream << std::endl;
  }
  outputStream << weightStream.str();
}

void CrossValidation::printAllWeightsColumns(ostream& weightStream) {
  std::cerr << "Learned normalized SVM weights for the " << numFolds_
            << " cross-validation splits:" << std::endl;
  printAllWeightsColumns(weights_, weightStream);
}

void CrossValidation::printAllRawWeightsColumns(ostream& weightStream,
                                                const Normalizer* pNorm) {
  std::cerr << "Learned raw SVM weights for the " << numFolds_
            << " cross-validation splits:" << std::endl;
  std::vector<std::vector<double> > ww = weights_;
  for (size_t set = 0; set < weights_.size(); ++set) {
    pNorm->unnormalizeweight(weights_[set], ww[set]);
  }
  printAllWeightsColumns(ww, weightStream);
}

}}}  // namespace OpenMS::Internal::Percolator
