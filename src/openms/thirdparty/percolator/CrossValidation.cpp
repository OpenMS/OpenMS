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

#include <exception>

#ifdef _OPENMP
#include <omp.h>
#include <algorithm>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// checks cross validation convergence in case of quickValidation_
const double CrossValidation::requiredIncreaseOver2Iterations_ = 0.01;

// OpenMS: the cross-validation loops run on at most numThreads threads, and never on
// more than OpenMP allows here (a TOPP tool's -threads). Upstream calls
// omp_set_num_threads() instead, which would change the count for the whole process.
// Without an OpenMP runtime (e.g. -fopenmp-simd, which still parses the clauses) 1.
// A loop over the folds (or pairs) uses no more threads than it has iterations. The loops
// whose iterations sort (Scores, see SortLikeStd.h) run them as tasks of a parallel region on
// more threads instead (numTeamThreads): the threads without an iteration of their own run
// the tasks of the sorts.
[[maybe_unused]] static int numLoopThreads(unsigned int numThreads,
                                           std::size_t iterations = std::size_t(-1)) {
#ifdef _OPENMP
  const std::size_t threads = std::min<std::size_t>(
      std::min(std::max(numThreads, 1u),
               static_cast<unsigned int>(omp_get_max_threads())),
      std::max<std::size_t>(iterations, 1));
  return static_cast<int>(threads);
#else
  (void)numThreads;
  (void)iterations;
  return 1;
#endif
}

// OpenMS: the parallel regions whose iterations run as tasks (with the tasks of their sorts,
// see SortLikeStd.h) use at most kMaxTeamThreads threads: more threads would mostly wait for
// tasks, polling the team's task queue.
constexpr unsigned int kMaxTeamThreads = 12;
[[maybe_unused]] static int numTeamThreads(unsigned int numThreads) {
  return numLoopThreads(std::min(numThreads, kMaxTeamThreads));
}

// OpenMS: rethrows the first error, in iteration order, caught in the iterations of a
// parallel loop (an exception must not leave an OpenMP region).
static void rethrowFirstError(const std::vector<std::exception_ptr>& errors) {
  for (const std::exception_ptr& error : errors) {
    if (error) std::rethrow_exception(error);
  }
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
                                 bool skipNormalizeScores,
                                 double decoyFractionTraining,
                                 unsigned int numFolds)
    : quickValidation_(quickValidation),
      usePi0_(usePi0),
      reportPerformanceEachIteration_(reportPerformanceEachIteration),
      testFdr_(testFdr),
      selectionFdr_(selectionFdr),
      initialSelectionFdr_(initialSelectionFdr),
      selectedCpos_(selectedCpos),
      selectedCneg_(selectedCneg),
      niter_(niter),
      nestedXvalBins_(nestedXvalBins),
      trainBestPositive_(trainBestPositive),
      numThreads_(numThreads),
      skipNormalizeScores_(skipNormalizeScores),
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
}

CrossValidation::~CrossValidation() {
  for (unsigned int set = 0; set < numFolds_ * nestedXvalBins_; ++set) {
    if (svmInputs_[set]) {
      delete svmInputs_[set];
    }
    svmInputs_[set] = NULL;
  }
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
 *
 * @param w_ vector of SVM weights
 * @return number of positives for initial setup
 */
int CrossValidation::preIterationSetup(Scores& fullset,
                                       SanityCheck* pCheck,
                                       const Normalizer* pNorm,
                                       FeatureMemoryPool& featurePool) {
  assert(nestedXvalBins_ >= 1u);
  // One input set per (nested) cross validation bin, to be reused multiple
  // times
  for (unsigned int set = 0; set < numFolds_ * nestedXvalBins_; ++set) {
    svmInputs_.push_back(new AlgIn(
        fullset.size(), static_cast<int>(FeatureNames::getNumFeatures()) + 1));
    assert(svmInputs_.back());
  }

  unsigned int decoysPerTarget = 1u;  // TODO: make this a command line argument
  // OpenMS: in a parallel region, whose threads run the tasks of the sort (see SortLikeStd.h)
  std::exception_ptr xvalError;
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
  try {
    fullset.createXvalSetsBySpectrum(trainScores_, testScores_, numFolds_,
                                     featurePool, decoyFractionTraining_,
                                     decoysPerTarget);
  } catch (...) {
    xvalError = std::current_exception();
  }
  if (xvalError) std::rethrow_exception(xvalError);

  // initialize weights vector for all folds
  weights_ = vector<vector<double> >(
      numFolds_, vector<double>(FeatureNames::getNumFeatures() + 1));
  pCheck->setNumThreads(static_cast<unsigned int>(numTeamThreads(numThreads_)));
  int numPositive =
      pCheck->getInitDirection(testScores_, trainScores_, pNorm, weights_,
                               testFdr_, initialSelectionFdr_);

  initializeGridSearch(fullset.getTargetDecoySizeRatio());

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
    candidatesCfrac_.push_back(1.0 * targetDecoySizeRatio);
    candidatesCfrac_.push_back(3.0 * targetDecoySizeRatio);
    candidatesCfrac_.push_back(10.0 * targetDecoySizeRatio);
    if (VERB > 0) {
      cerr << "Selecting Cneg by cross-validation." << endl;
    }
  }

  // Form cpos, cneg pairs per nested CV fold
  CandidateCposCfrac cpCnFold;
  for (unsigned int set = 0; set < numFolds_; ++set) {
    for (int nestedSet = 0; nestedSet < nestedXvalBins_; nestedSet++) {
      std::vector<double> candidatesCpos = candidatesCpos_;
      std::vector<double> candidatesCfrac = candidatesCfrac_;
      if (quickValidation_) {  // in quick validation,
        candidatesCpos = std::vector<double>(1, 1);
        candidatesCfrac = std::vector<double>(1, 1);
      }
      for (double cpos : candidatesCpos) {
        for (double cfrac : candidatesCfrac) {
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
  for (unsigned int i = 0; i < niter_; i++) {
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
    if (foundPositives > 0 && foundPositivesOldOld > 0 && quickValidation_ &&
        (static_cast<double>(foundPositives - foundPositivesOldOld) <=
         foundPositivesOldOld * requiredIncreaseOver2Iterations_)) {
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
  // OpenMS: the folds in parallel; the counts are added up in fold order
  vector<int> testPositives(numFolds_);
  vector<std::exception_ptr> testErrors(numFolds_);
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
#pragma omp taskloop grainsize(1)
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    try {
      testPositives[set] =
          testScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
    } catch (...) {
      testErrors[set] = std::current_exception();
    }
  }
  rethrowFirstError(testErrors);
  foundPositives = 0;
  for (size_t set = 0; set < numFolds_; ++set) {
    foundPositives += testPositives[set];
  }
  if (VERB > 0) {
    std::cerr << "Found " << foundPositives << " test set PSMs with q<"
              << testFdr_ << "." << std::endl;
  }
}

/**
 * Executes a cross validation step
 * @param w_ list of the bins' normal vectors (in linear algebra sense) of the
 *        hyperplane from SVM
 * @return Estimation of number of true positives
 */
int CrossValidation::doStep(const Normalizer* pNorm, double selectionFdr) {
  // Setup
  options pOptions;
  pOptions.lambda = 1.0;
  pOptions.lambda_u = 1.0;
  pOptions.epsilon = EPSILON;
  pOptions.cgitermax = CGITERMAX;
  pOptions.mfnitermax = MFNITERMAX;
  int estTruePos = 0;

  // for determining an appropriate positive training set, the decoys+1 in the
  // FDR estimates is too restrictive for small datasets
  bool skipDecoysPlusOne = true;
  // OpenMS: the folds in parallel; each scores and sorts only its own training set
  vector<std::exception_ptr> trainErrors(numFolds_);
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
#pragma omp taskloop grainsize(1)
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    try {
      trainScores_[set].calcScoresAndQvals(weights_[set], selectionFdr,
                                           skipDecoysPlusOne);
    } catch (...) {
      trainErrors[set] = std::current_exception();
    }
  }
  rethrowFirstError(trainErrors);

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

  // Create SVM input data for parallelization
  // OpenMS: the folds in parallel (each fills its own SVM inputs from its own training set),
  // unless nested cross-validation draws the nested sets at random
  std::vector<AlgIn*> svmInputsVec(numFolds_ * nestedXvalBins_);
  std::vector<std::vector<Scores> > nestedTestScoresVec(numFolds_);
  vector<std::exception_ptr> setupErrors(numFolds_);
  const auto setupFold = [&](int set) {
    try {
      std::vector<Scores> nestedTrainScores(nestedXvalBins_, usePi0_),
          nestedTestScores(nestedXvalBins_, usePi0_);
      if (nestedXvalBins_ > 1) {
        FeatureMemoryPool featurePool;
        trainScores_[set].createXvalSetsBySpectrum(
            nestedTrainScores, nestedTestScores, nestedXvalBins_, featurePool);
      } else {
        // sub-optimal cross validation
        nestedTrainScores[0] = trainScores_[set];
        nestedTestScores[0] = trainScores_[set];
      }
      nestedTestScoresVec[set] = nestedTestScores;
      // Set SVM input data for L2-SVM-MFN
      for (std::size_t nestedFold = 0; nestedFold < nestedXvalBins_;
           ++nestedFold) {
        AlgIn* svmInput = svmInputs_[set * nestedXvalBins_ + nestedFold];
        if ((VERB > 2) && (nestedFold == 0)) {
          cerr << "Split " << set + 1 << ": Training with " << svmInput->positives
               << " positives and " << svmInput->negatives << " negatives"
               << std::endl;
        }
        nestedTrainScores[nestedFold].generateNegativeTrainingSet(*svmInput, 1.0);
        nestedTrainScores[nestedFold].generatePositiveTrainingSet(
            *svmInput, selectionFdr, 1.0, trainBestPositive_);
        if ((VERB > 2) && (nestedXvalBins_ > 1u)) {
          cerr << "Split " << set + 1 << " nested fold " << nestedFold
               << ": Training with " << svmInput->positives << " positives and "
               << svmInput->negatives << " negatives" << std::endl;
        }
        svmInputsVec[set * nestedXvalBins_ + nestedFold] = svmInput;
      }
    } catch (...) {
      setupErrors[set] = std::current_exception();
    }
  };
  if (nestedXvalBins_ > 1) {
    // nested cross-validation draws the nested sets at random: the folds in their order
    for (int set = 0; set < static_cast<int>(numFolds_); ++set) setupFold(set);
  } else {
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
#pragma omp taskloop grainsize(1)
    for (int set = 0; set < static_cast<int>(numFolds_); ++set) setupFold(set);
  }
  rethrowFirstError(setupErrors);

  // OpenMS: the (cpos, cneg) pairs as tasks of one parallel region (see numTeamThreads). Each
  // task then validates its pair: the positives of the pair's weights on the pair's test set,
  // counted without changing the set (Scores::countPositives), so that the pairs of a fold can
  // do so at the same time (the set is not used afterwards). The results do not depend on the
  // number of threads.
  // for determining the number of positives, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  const bool skipDecoysPlusOneValidation = true;
  vector<int> pairTruePoses(classWeightsPerFold_.size());
  vector<std::exception_ptr> pairErrors(classWeightsPerFold_.size());
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
  for (int pairIdx = 0; pairIdx < static_cast<int>(classWeightsPerFold_.size()); pairIdx++) {
#pragma omp task firstprivate(pairIdx) shared(pairErrors, pairTruePoses, svmInputsVec, nestedTestScoresVec, pOptions)
    {
      try {
        CandidateCposCfrac* cpCnFold = &classWeightsPerFold_[pairIdx];
        AlgIn* svmInput =
            svmInputsVec[cpCnFold->set * nestedXvalBins_ +
                         static_cast<unsigned int>(cpCnFold->nestedSet)];
        trainCpCnPair(*cpCnFold, pOptions, svmInput);
        pairTruePoses[pairIdx] =
            nestedTestScoresVec[cpCnFold->set][static_cast<std::size_t>(cpCnFold->nestedSet)]
                .countPositives(cpCnFold->ww, testFdr_, skipDecoysPlusOneValidation);
      } catch (...) {
        pairErrors[pairIdx] = std::current_exception();
      }
    }
  }
  rethrowFirstError(pairErrors);

  estTruePos = mergeCpCnPairs(selectionFdr, pOptions, pairTruePoses,
                              candidatesCpos_, candidatesCfrac_);
  return estTruePos;
}

/**
 * Train SVM over a single (cpos, cneg) pair
 * @param cpCnFold contains cpos, cneg pair and SVM learned weights
 * @param pOptions options for the SVM algorithm
 * @param svmInput training data for this particular nested CV fold
 */
void CrossValidation::trainCpCnPair(CandidateCposCfrac& cpCnFold,
                                    options& pOptions,
                                    AlgIn* svmInput) {
  vector_double pWeights;
  pWeights.d = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
  pWeights.vec = new double[pWeights.d];

  double cpos = cpCnFold.cpos;
  double cfrac = cpCnFold.cfrac;

  // Create storage vector for SVM algorithm
  vector_double Outputs;
  size_t numInputs =
      static_cast<std::size_t>(svmInput->positives + svmInput->negatives);
  Outputs.vec = new double[numInputs];
  Outputs.d = static_cast<int>(numInputs);

  if (VERB > 3)
    cerr << "- cross-validation with Cpos=" << cpos << ", Cneg=" << cfrac * cpos
         << endl;
  for (int ix = 0; ix < pWeights.d; ix++) {
    pWeights.vec[ix] = 0;
  }
  for (int ix = 0; ix < Outputs.d; ix++) {
    Outputs.vec[ix] = 0;
  }

  // Call SVM algorithm (see ssl.cpp)
  L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs, cpos, cfrac * cpos);

  for (std::size_t i = FeatureNames::getNumFeatures() + 1; i--;) {
    cpCnFold.ww[i] = pWeights.vec[i];
  }
}

/**
 * Validate and merge weights learned per cpos,cneg pairs per nested CV fold per
 * CV fold
 * @param pWeights results vector from the SVM algorithm
 * @param pOptions options for the SVM algorithm
 * @param pairTruePoses positives per (cpos,cneg) pair per nested CV fold per
 * CV fold, on its test set (OpenMS: computed by doStep)
 */
int CrossValidation::mergeCpCnPairs(
    double selectionFdr,
    options& pOptions,
    const vector<int>& pairTruePoses,
    const vector<double>& cposCandidates,
    const vector<double>& cfracCandidates) {
  vector<int> bestTruePoses(numFolds_, 0);
  vector<double> bestCposes(numFolds_, 1);
  vector<double> bestCfracs(numFolds_, 1);
  
  // Validate learned parameters per (cpos,cneg) pair per nested CV fold
  // OpenMS: the positives of the pairs are validated in doStep, in the pairs' tasks;
  // the pairs are merged here in their order, as before.
  unsigned int numCpCnPairsPerSet =
      static_cast<unsigned int>(classWeightsPerFold_.size() / numFolds_);
  vector<std::exception_ptr> mergeErrors(numFolds_);
#pragma omp parallel for schedule(dynamic, 1) ordered num_threads(numLoopThreads(std::min(numThreads_, kMaxTeamThreads), numFolds_))
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    try {
      unsigned int a = set * numCpCnPairsPerSet;
      unsigned int b = (set + 1) * numCpCnPairsPerSet;
      int tp = 0;
      std::vector<CandidateCposCfrac>::iterator itCpCnPair;
      std::map<std::pair<double, double>, int> intermediateResults;
      for (itCpCnPair = classWeightsPerFold_.begin() + a;
           itCpCnPair < classWeightsPerFold_.begin() + b; itCpCnPair++) {
        tp = pairTruePoses[static_cast<std::size_t>(itCpCnPair - classWeightsPerFold_.begin())];
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
          std::vector<double>::const_iterator itCfrac = cfracCandidates.begin();
          for (; itCfrac != cfracCandidates.end(); ++itCfrac) {
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
    } catch (...) {
      mergeErrors[set] = std::current_exception();
    }
  }
  rethrowFirstError(mergeErrors);

  if (nestedXvalBins_ > 1) {
    vector<std::exception_ptr> refitErrors(numFolds_);
#pragma omp parallel for schedule(dynamic, 1) ordered num_threads(numLoopThreads(std::min(numThreads_, kMaxTeamThreads), numFolds_))
    for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
      try {
        vector_double pWeights;
        pWeights.d = static_cast<int>(FeatureNames::getNumFeatures()) + 1;
        pWeights.vec = new double[pWeights.d];

        AlgIn* svmInput = svmInputs_[set * nestedXvalBins_];
        trainScores_[set].generateNegativeTrainingSet(*svmInput, 1.0);
        trainScores_[set].generatePositiveTrainingSet(*svmInput, selectionFdr,
                                                      1.0, trainBestPositive_);

        // Create storage vector for SVM algorithm
        vector_double Outputs;
        size_t numInputs =
            static_cast<std::size_t>(svmInput->positives + svmInput->negatives);
        Outputs.vec = new double[numInputs];
        Outputs.d = static_cast<int>(numInputs);

        for (int ix = 0; ix < pWeights.d; ix++) {
          pWeights.vec[ix] = 0;
        }
        for (int ix = 0; ix < Outputs.d; ix++) {
          Outputs.vec[ix] = 0;
        }
        // Call SVM algorithm (see ssl.cpp)
        L2_SVM_MFN(*svmInput, pOptions, pWeights, Outputs, bestCposes[set],
                   bestCposes[set] * bestCfracs[set]);

        for (std::size_t i = FeatureNames::getNumFeatures() + 1; i--;) {
          weights_[set][i] = pWeights.vec[i];
        }
      } catch (...) {
        refitErrors[set] = std::current_exception();
      }
    }
    rethrowFirstError(refitErrors);
  }

  // OpenMS: the folds in parallel; the counts are added up in fold order
  vector<int> truePoses(numFolds_);
  vector<std::exception_ptr> foldErrors(numFolds_);
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
#pragma omp taskloop grainsize(1)
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    try {
      truePoses[set] =
          trainScores_[set].calcScoresAndQvals(weights_[set], testFdr_);
    } catch (...) {
      foldErrors[set] = std::current_exception();
    }
  }
  rethrowFirstError(foldErrors);
  double bestTruePos = 0;
  for (int set = 0; set < static_cast<int>(numFolds_); ++set) {
    bestTruePos += truePoses[set];
  }
  return static_cast<int>(bestTruePos / (numFolds_ - 1));
}

void CrossValidation::postIterationProcessing(Scores& fullset,
                                              SanityCheck* pCheck) {
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

  if (!pCheck->validateDirection(weights_)) {
    for (std::size_t set = 0; set < numFolds_; ++set) {
      testScores_[set].calcScoresAndQvals(weights_[0], selectionFdr_);
    }
  }
  // OpenMS: in a parallel region, whose threads run the tasks of the sorts (see SortLikeStd.h)
  std::exception_ptr mergeError;
#pragma omp parallel num_threads(numTeamThreads(numThreads_))
#pragma omp single
  try {
    fullset.merge(testScores_, selectionFdr_, skipNormalizeScores_, weights_);
  } catch (...) {
    mergeError = std::current_exception();
  }
  if (mergeError) std::rethrow_exception(mergeError);
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
