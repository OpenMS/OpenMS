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

#include <boost/assign.hpp>
#include <cassert>
#include <cstring>
#include <exception>
#include <functional>
#include <cfloat>
#include <cmath>
#include <unordered_set>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DataSet.h"
#include "Globals.h"
#include "MassHandler.h"
#include "Normalizer.h"
#include "PosteriorEstimator.h"
#include "Scores.h"
#include "SetHandler.h"
#include "ssl.h"
#include "IsotonicPEP.h"
#include "SortLikeStd.h"
#include "OpenMPTasks.h"  // OpenMS: PERCOLATOR_OMP_TASKS, PERCOLATOR_TASK_PRAGMA

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// OpenMS: the smallest block of PSMs a task scores (Scores::calcScores)
[[maybe_unused]] constexpr std::ptrdiff_t kMinScoreTaskSize = 2048;
// OpenMS: the smallest block of PSMs a task of getScoreLabelPairs, calcQvals or normalizeScores
// handles
[[maybe_unused]] constexpr std::ptrdiff_t kMinPairTaskSize = 4 * kMinScoreTaskSize;

// OpenMS: the grain size of a taskloop over n elements: at least minGrain, and at most kMaxLoopTasks
// tasks. (libgomp runs all tasks of a taskloop in the encountering thread, one after the other, if
// its team's queue would then hold more than 64 tasks per thread; the grain size does not change
// what an element computes.)
[[maybe_unused]] constexpr std::ptrdiff_t kMaxLoopTasks = 64;
[[maybe_unused]] static std::ptrdiff_t loopGrain(std::ptrdiff_t n, std::ptrdiff_t minGrain) {
  return std::max(minGrain, (n + kMaxLoopTasks - 1) / kMaxLoopTasks);
}

// OpenMS: from this many PSMs on, the validations of the (Cpos, Cneg) pairs (countPositives) score
// the PSMs and sort as tasks of the enclosing OpenMP parallel region (in which the cross-validation
// runs them), and so does getFeatureRows
[[maybe_unused]] constexpr std::ptrdiff_t kLargeSet = std::ptrdiff_t(1) << 17;

// OpenMS: whether ScoreHolder's operator< and operator> (score, scan, expMass, label) are a
// strict weak ordering on [first, last): no score or mass is NaN.
static bool comparableScores(std::vector<ScoreHolder>::const_iterator first,
                             std::vector<ScoreHolder>::const_iterator last) {
  return std::none_of(first, last, [](const ScoreHolder& sh) {
    return std::isnan(sh.score) || std::isnan(sh.pPSM->expMass);
  });
}

// OpenMS: some sorts of ScoreHolders sort keys instead: what the comparison compares, and the
// position of the ScoreHolder. If the comparisons of the keys have the outcomes of the
// comparisons of their ScoreHolders, a sort of the keys (the same algorithm) makes the same
// moves as a sort of the ScoreHolders: it leaves the same permutation, also of ScoreHolders that
// compare equal. permuteScores then puts scores[keys[k].index] at position k.
// OpenMS: inside an OpenMP parallel region, large sets are permuted by its tasks (blocks of
// positions).
template <typename Key>
static void permuteScores(std::vector<ScoreHolder>& scores, const std::vector<Key>& keys) {
#ifdef _OPENMP
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(keys.size());
  if (n >= 4 * kMinScoreTaskSize && keys.size() == scores.size() && omp_in_parallel()) {
    std::vector<ScoreHolder> permuted(keys.size());
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(n, 4 * kMinScoreTaskSize)) shared(scores, keys, permuted))
    for (std::ptrdiff_t k = 0; k < n; ++k) {
      permuted[static_cast<std::size_t>(k)] = scores[keys[static_cast<std::size_t>(k)].index];
    }
    scores.swap(permuted);
    return;
  }
#endif
  std::vector<ScoreHolder> permuted;
  permuted.reserve(scores.size());
  for (const Key& key : keys) permuted.push_back(scores[key.index]);
  scores.swap(permuted);
}

// OpenMS: operator> of ScoreHolder (ScoreHolder.cpp), as a comparison that the sorts below can
// inline (same outcomes, so the sorts make the same moves)
struct ScoreHolderGreater {
  bool operator()(const ScoreHolder& one, const ScoreHolder& other) const {
    return (one.score > other.score) ||
           (one.score == other.score && one.pPSM->scan > other.pPSM->scan) ||
           (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
            one.pPSM->expMass > other.pPSM->expMass) ||
           (one.score == other.score && one.pPSM->scan == other.pPSM->scan &&
            one.pPSM->expMass == other.pPSM->expMass && one.label > other.label);
  }
};

// OpenMS: sorts the scores like std::sort, with the tasks of an enclosing OpenMP parallel
// region (see SortLikeStd.h). Scores that are no strict weak ordering get std::sort itself.
template <typename Less>
static void sortScoresLikeStdSort(std::vector<ScoreHolder>& scores, Less less) {
  if (comparableScores(scores.begin(), scores.end())) {
    sortLikeStdSort(scores, less);
  } else {
    std::sort(scores.begin(), scores.end(), less);
  }
}

void Scores::merge(std::vector<Scores>& sv,
                   double fdr,
                   bool skipNormalizeScores,
                   std::vector<std::vector<double> >& all_w) {
  scores_.clear();
  std::vector<Scores>::iterator cvBinScores = sv.begin();
  std::vector<std::vector<double> >::iterator weights = all_w.begin();
  for (; cvBinScores != sv.end(); cvBinScores++, weights++) {
    sortScoresLikeStdSort(cvBinScores->scores_, ScoreHolderGreater());
    cvBinScores->checkSeparationAndSetPi0();
    cvBinScores->calcQvals(fdr);
    if (!skipNormalizeScores) {
      cvBinScores->normalizeScores(fdr, *weights);
    }
    copy(cvBinScores->begin(), cvBinScores->end(), back_inserter(scores_));
  }
  postMergeStep();
}

void Scores::postMergeStep() {
  sortScoresLikeStdSort(scores_, ScoreHolderGreater());
  totalNumberOfDecoys_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isDecoy)));
  totalNumberOfTargets_ = static_cast<unsigned int>(
      count_if(scores_.begin(), scores_.end(), mem_fn(&ScoreHolder::isTarget)));
  targetDecoySizeRatio_ =
      totalNumberOfTargets_ / max(1.0, (double)totalNumberOfDecoys_);
  checkSeparationAndSetPi0();
}

double Scores::calcScore(const double* feat,
                         const std::vector<double>& w) const {
  std::size_t ix = FeatureNames::getNumFeatures();
  double score = w[ix];
  for (; ix--;) {
    score += feat[ix] * w[ix];
  }
  return score;
}

void Scores::scoreAndAddPSM(ScoreHolder& sh,
                            const std::vector<double>& rawWeights,
                            FeatureMemoryPool& featurePool) {
  const unsigned int numFeatures =
      static_cast<unsigned int>(FeatureNames::getNumFeatures());

  for (unsigned int j = 0; j < numFeatures; j++) {
    sh.score += sh.pPSM->features[j] * rawWeights[j];
  }
  sh.score += rawWeights[numFeatures];

  featurePool.deallocate(sh.pPSM->features);
  sh.pPSM->deleteRetentionFeatures();

  if (sh.isTarget()) {
    ++totalNumberOfTargets_;
  } else if (sh.isDecoy()) {
    ++totalNumberOfDecoys_;
  }

  if (!sh.isTarget() && !sh.isDecoy()) {
    std::cerr << "Warning: the PSM " << sh.pPSM->getId()
              << " has a label not in {1,-1} and will be ignored." << std::endl;
    PSMDescription::deletePtr(sh.pPSM);
  } else {
    scores_.push_back(sh);
  }
}

void Scores::print(LabelType label, std::ostream& os) {
 std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  os << "PSMId\t";
  if (PSMDescription::hasSpectrumFileName()) {
    os << "filename\t";
  }
  os << "score\tq-value\tposterior_error_prob\tpeptide\tproteinIds";
  if (is_output_rt_) {
    os << "\trt\n";
  } else {
    os << "\n";
  }

  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->label == label) {
      std::ostringstream out;
      scoreIt->pPSM->printProteins(out);
      ResultHolder rh(scoreIt->score, scoreIt->q, scoreIt->pep,
                      scoreIt->pPSM->getId(), scoreIt->pPSM->getFullPeptide(),
                      out.str(), scoreIt->pPSM->getSpectrumFileName());
      if (is_output_rt_) {
        rh.retentionTime = scoreIt->pPSM->getRetentionTime();
        rh.outputRT =
            true;  // Ideally this should be set once outside this loop, but
                   // ResultHolder needs to be initialized before.
      }

      os << rh << std::endl;
    }
  }
}

void Scores::populateWithPSMs(SetHandler& setHandler) {
  scores_.clear();
  setHandler.populateScoresWithPSMs(scores_, LabelType::TARGET);
  setHandler.populateScoresWithPSMs(scores_, LabelType::DECOY);
  totalNumberOfTargets_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::TARGET));
  totalNumberOfDecoys_ =
      static_cast<unsigned int>(setHandler.getSizeFromLabel(LabelType::DECOY));
  targetDecoySizeRatio_ = (double)totalNumberOfTargets_ / totalNumberOfDecoys_;

  if (VERB > 1) {
    cerr << "Train/test set contains " << totalNumberOfTargets_
         << " positives and " << totalNumberOfDecoys_
         << " negatives, size ratio=" << targetDecoySizeRatio_
         << " and pi0=" << pi0_ << endl;
  }

  if (totalNumberOfTargets_ == 0) {
    ostringstream oss;
    oss << "Error: no target PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  if (totalNumberOfDecoys_ == 0) {
    ostringstream oss;
    oss << "Error: no decoy PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
    } else {
      throw MyException(oss.str());
    }
  }

  // check for the minimum recommended number of positive and negative hits
  if (totalNumberOfTargets_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of positive samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
  if (totalNumberOfDecoys_ <= (unsigned)(FeatureNames::getNumFeatures() * 5)) {
    std::cerr << "Warning : the number of negative samples read is too small "
                 "to perform a correct classification.\n"
              << std::endl;
  }
}

/**
 * Divides the PSMs from pin file into xval_fold cross-validation sets based on
 * their spectrum scan number
 * @param train vector containing the training sets of PSMs
 * @param test vector containing the test sets of PSMs
 * @param xval_fold number of folds in train and test
 * @param featurePool pool for storing features to optimize memory access in SVM
 * @param decoyFractionTraining fraction of decoys to use for training, will
 * be 1 unless RESET algorithm is used
 */
void Scores::createXvalSetsBySpectrum(std::vector<Scores>& train,
                                      std::vector<Scores>& test,
                                      unsigned int xval_fold,
                                      FeatureMemoryPool& featurePool,
                                      double decoyFractionTraining,
                                      unsigned int decoysPerTarget,
                                      const std::function<unsigned int(const PSMDescription&)>* fixedFold) {
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold, Scores(usePi0_));
  test.resize(xval_fold, Scores(usePi0_));
  size_t approxTestSize = scores_.size() / xval_fold + 1;
  size_t approxTrainSize = (approxTestSize < scores_.size()) ? scores_.size() - approxTestSize : scores_.size();
#ifdef _OPENMP
  const bool fillInTasks = xval_fold > 1u && omp_in_parallel();
#else
  const bool fillInTasks = false;
#endif
  // OpenMS: filled by tasks (below), each fold reserves exactly its sizes there (here, a fold of whole spectra could
  // outgrow the estimate and reallocate once)
  for (unsigned int i = 0; i < xval_fold && !fillInTasks; ++i) {
    train[i].reserve(approxTrainSize);
    test[i].reserve(approxTestSize);
  }
  // remain keeps track of residual space available in each fold
  std::vector<int> remain(xval_fold);
  // set values for remain: initially each fold is assigned (tot number of
  // scores_ / tot number of folds)
  int fold = static_cast<int>(xval_fold), ix = static_cast<int>(scores_.size());
  while (fold--) {
    remain[static_cast<std::size_t>(fold)] = ix / (fold + 1);
    ix -= remain[static_cast<std::size_t>(fold)];
  }

  // OpenMS: sorted like std::sort(scores_.begin(), scores_.end(), OrderScanHash()) (see
  // SortLikeStd.h), on (hash, position) keys: the hashes are those OrderScanHash compares, so the
  // comparisons have the same outcomes and the sort makes the same moves (see permuteScores),
  // without hashing two PSMs for every comparison
  {
    struct HashKey {
      size_t hash;
      std::size_t index;
    };
    std::vector<HashKey> keys(scores_.size());
    const std::ptrdiff_t numScores = static_cast<std::ptrdiff_t>(scores_.size());
    // (blocks of PSMs as the tasks of an enclosing parallel region, if any)
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numScores, 4 * kMinScoreTaskSize)) shared(keys) if(numScores >= 8 * kMinScoreTaskSize))
    for (std::ptrdiff_t i = 0; i < numScores; ++i) {
      const PSMDescription* psm = scores_[static_cast<std::size_t>(i)].pPSM;
      const size_t hash = fast_uint_hash(psm->specFileNr) ^ (fast_uint_hash(psm->scan) << 1);
      keys[static_cast<std::size_t>(i)] = {hash, static_cast<std::size_t>(i)};
    }
    sortLikeStdSort(keys, [](const HashKey& x, const HashKey& y) { return x.hash < y.hash; });
    permuteScores(scores_, keys);
  }

  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  // put scores into the folds; choose a fold (at random) and change it only
  // when scores from a new spectra are encountered
  // OpenMS: the folds and labels are drawn here, in the order of the PSMs as before; inside an
  // OpenMP parallel region, the folds are filled afterwards, each by a task of its own, in the
  // same order (see below)
  std::vector<unsigned int> foldOfPSM;
  std::vector<LabelType> labelOfPSM;
  if (fillInTasks) {
    foldOfPSM.reserve(scores_.size());
    labelOfPSM.reserve(scores_.size());
  }
  unsigned int previousSpectrum = scores_.begin()->pPSM->scan;
  // OpenMS: with fixedFold (out-of-fold capped training), each PSM's fold is given (the caller
  // keeps the PSMs of a spectrum in one fold); no fold is drawn
  const bool foldsFixed = fixedFold != nullptr && static_cast<bool>(*fixedFold);
  size_t randIndex = foldsFixed ? 0u : PseudoRandom::lcg_rand() % xval_fold;
  for (std::vector<ScoreHolder>::iterator it = scores_.begin();
       it != scores_.end(); ++it) {
    const unsigned int curScan = (*it).pPSM->scan;
    ScoreHolder sh = (*it);

    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    if (foldsFixed) {
      randIndex = (*fixedFold)(*sh.pPSM);
      if (randIndex >= xval_fold) {
        throw MyException("Error: a fixed cross-validation fold is out of range.\n");
      }
    } else if (previousSpectrum != curScan) {
      randIndex = PseudoRandom::lcg_rand() % xval_fold;
      // allow only indexes of folds that are non-full
      while (remain[randIndex] <= 0) {
        randIndex = PseudoRandom::lcg_rand() % xval_fold;
      }
    }

    // if we use multiple folds with RESET, assign 1-decoyFractionTraining as
    // pseudo targets.
    if (xval_fold > 1u && decoyFractionTraining < 1.0 &&
        sh.label == LabelType::DECOY &&
        PseudoRandom::lcg_uniform_rand() > decoyFractionTraining) {
      // From Algorithm S3 of the percolator-RESET supplementary material
      // decoyFractionTraining - the probability of assigning a decoy to the
      // training set
      sh.label = LabelType::PSEUDO_TARGET;
    }

    // insert
    if (fillInTasks) {
      foldOfPSM.push_back(static_cast<unsigned int>(randIndex));
      labelOfPSM.push_back(sh.label);
    } else {
      for (unsigned int i = 0; i < xval_fold; ++i) {
        if (i == randIndex) {
          test[i].addScoreHolder(sh);
        }
        if (i != randIndex || xval_fold == 1) {
          train[i].addScoreHolder(sh);
        }
      }
    }
    // update number of free position for used fold
    --remain[randIndex];
    // set previous spectrum to current one for next iteration
    previousSpectrum = curScan;
  }
  if (fillInTasks) {
    // OpenMS: fold i gets the PSMs drawn for it (test) and the others (train), in the order of
    // the PSMs, as the loop above would insert them
    std::vector<std::exception_ptr> fillErrors(xval_fold);
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(train, test, foldOfPSM, labelOfPSM, fillErrors))
    for (int i = 0; i < static_cast<int>(xval_fold); ++i) {
      try {
        const std::size_t testSize = static_cast<std::size_t>(
            std::count(foldOfPSM.begin(), foldOfPSM.end(), static_cast<unsigned int>(i)));
        test[static_cast<std::size_t>(i)].reserve(testSize);
        train[static_cast<std::size_t>(i)].reserve(scores_.size() - testSize);
        for (std::size_t k = 0; k < scores_.size(); ++k) {
          ScoreHolder sh = scores_[k];
          sh.label = labelOfPSM[k];
          if (foldOfPSM[k] == static_cast<unsigned int>(i)) {
            test[static_cast<std::size_t>(i)].addScoreHolder(sh);
          } else {
            train[static_cast<std::size_t>(i)].addScoreHolder(sh);
          }
        }
      } catch (...) {
        fillErrors[static_cast<std::size_t>(i)] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : fillErrors) {
      if (error) std::rethrow_exception(error);
    }
  }

  // without RESET: decoysPerTarget = 1 and decoyFractionTraining = 1, then
  // factor = 1 and trainNullTargetWinProb = 1/2 (and decoyFactor = 1)
  //
  // with RESET 1 decoy: decoysPerTarget = 1 and decoyFractionTraining = 0.5,
  // then factor = 1.5 and trainNullTargetWinProb = 3/4 (and decoyFactor = 3)
  //
  // with RESET 2 decoys: decoysPerTarget = 2 and decoyFractionTraining = 0.5,
  // then factor = 2 and trainNullTargetWinProb = 2/3 (and decoyFactor = 2)
  double factor =
      (decoysPerTarget + 1) - decoyFractionTraining * decoysPerTarget;
  double trainNullTargetWinProb = factor / (1.0 + decoysPerTarget);

  // calculate ratios of target over decoy for train and test set
  for (unsigned int i = 0; i < xval_fold; ++i) {
    train[i].recalculateSizes();
    test[i].recalculateSizes();

    train[i].setNullTargetWinProb(trainNullTargetWinProb);
    test[i].setNullTargetWinProb(trainNullTargetWinProb);
  }

  if (featurePool.isInitialized()) {
    // OpenMS: the rows of the test sets' PSMs in the order the loop below moves them to (the
    // first row, the second, ...); FeatureMemoryPool::moveRowsInOrder moves them there directly
    // if they are the pool's rows, each once (then the loop leaves the same rows and pointers)
    std::vector<double**> rowsInOrder;
    rowsInOrder.reserve(scores_.size());
    for (unsigned int i = 0; i < xval_fold; ++i) {
      for (const bool isTarget : {true, false}) {
        for (ScoreHolder& sh : test[i]) {
          if (sh.isTarget() == isTarget) rowsInOrder.push_back(&sh.pPSM->features);
        }
      }
    }
    if (featurePool.moveRowsInOrder(rowsInOrder, FeatureNames::getNumFeatures())) return;
    boost::unordered_map<double*, double*> movedAddresses;
    size_t idx = 0;
    for (unsigned int i = 0; i < xval_fold; ++i) {
      bool isTarget = true;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
      isTarget = false;
      test[i].reorderFeatureRows(featurePool, isTarget, movedAddresses, idx);
    }
  }
}

void Scores::recalculateSizes() {
  totalNumberOfTargets_ = 0;
  totalNumberOfDecoys_ = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget()) {
      ++totalNumberOfTargets_;
    } else {
      ++totalNumberOfDecoys_;
    }
  }
  targetDecoySizeRatio_ = totalNumberOfTargets_ / (double)totalNumberOfDecoys_;
}

void Scores::reorderFeatureRows(
    FeatureMemoryPool& featurePool,
    bool isTarget,
    boost::unordered_map<double*, double*>& movedAddresses,
    size_t& idx) {
  size_t numFeatures = FeatureNames::getNumFeatures();
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() == isTarget) {
      double* newAddress =
          featurePool.addressFromIdx(static_cast<unsigned int>(idx++));
      double* oldAddress = scoreIt->pPSM->features;
      while (movedAddresses.find(oldAddress) != movedAddresses.end()) {
        oldAddress = movedAddresses[oldAddress];
      }
      std::swap_ranges(oldAddress, oldAddress + numFeatures, newAddress);
      scoreIt->pPSM->features = newAddress;
      if (oldAddress != newAddress) {
        movedAddresses[newAddress] = oldAddress;
      }
    }
  }
}

/**
 * Linear rescaling of SVM scores such that q=fdr = 0 and the median decoy = -1.
 * @param fdr fdr-cutoff to use, typically --trainFDR
 * @param weights SVM weights used to compute SVM scores
 */
void Scores::normalizeScores(double fdr, std::vector<double>& weights) {
  unsigned int medianIndex = std::max(0u, totalNumberOfDecoys_ / 2u),
               decoys = 0u;

  if (scores_.size() == 0) {
    ostringstream oss;
    oss << "Error: no scored PSMs were provided.\n";
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error."
           << std::endl;
      return;
    } else {
      throw MyException(oss.str());
    }
  }

  std::vector<ScoreHolder>::iterator it = scores_.begin();
  double fdrScore = it->score;
  double medianDecoyScore = fdrScore + 1.0;

  for (; it != scores_.end(); ++it) {
    if (it->q < fdr)
      fdrScore = it->score;
    if (it->isDecoy()) {
      if (++decoys == medianIndex) {
        medianDecoyScore = it->score;
        break;
      }
    }
  }

  double diff = fdrScore - medianDecoyScore;
  if (diff <= 0) {
    ostringstream oss;
    oss << "Error: median decoy score <= score at " << fdr * 100
        << "% FDR. Cannot rescale scores to merge cross validation bins, try "
        << "lowering --trainFDR.\n";
    if (NO_TERMINATE) {
      cerr << oss.str()
           << "No-terminate flag set: apply offset such that median "
           << "decoy has score -1, but skipping rescaling." << std::endl;
      diff = 1.0;
      fdrScore += 1.0;
    } else {
      throw MyException(oss.str());
    }
  }

#ifdef _OPENMP
  // OpenMS: inside an OpenMP parallel region, large sets by its tasks (blocks of PSMs)
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(scores_.size());
  if (n >= 2 * kMinPairTaskSize && omp_in_parallel()) {
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(n, kMinPairTaskSize)) shared(fdrScore, diff))
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
      sh.score -= fdrScore;
      sh.score /= diff;
    }
    Normalizer::endScoreNormalizeWeights(weights, weights, fdrScore, diff);
    return;
  }
#endif
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    scoreIt->score -= fdrScore;
    scoreIt->score /= diff;
  }
  Normalizer::endScoreNormalizeWeights(weights, weights, fdrScore, diff);
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * and calculate q values for each PSM
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcScoresAndQvals(std::vector<double>& w,
                               double fdr,
                               bool skipDecoysPlusOne) {
  calcScores(w);
  return calcQvals(fdr, skipDecoysPlusOne);
}

/**
 * OpenMS: the number of positives calcScoresAndQvals(w, fdr, skipDecoysPlusOne) returns,
 * without changing the scores (CrossValidation validates the (Cpos, Cneg) pairs of a fold
 * with it, concurrently, on the same test set). The count depends on the order of the scores
 * and on the labels of each group of equal scores, not on the order within a group
 * (calcQvals treats ties as one group), so the (score, label) pairs are sorted by score
 * only. With a NaN score (no strict weak ordering), calcScoresAndQvals on a copy.
 */
int Scores::countPositives(const std::vector<double>& w, double fdr,
                           bool skipDecoysPlusOne) const {
  std::vector<pair<double, bool> > combined;
  combined.reserve(scores_.size());
  for (const ScoreHolder& sh : scores_) {
    combined.emplace_back(calcScore(sh.pPSM->features, w), sh.isTarget());
  }
  return countPositivesOfScores(combined, w, fdr, skipDecoysPlusOne);
}

/**
 * OpenMS: the feature rows of the PSMs (pPSM->features) with their isTarget(), sorted by address:
 * countPositives(w, fdr, skipDecoysPlusOne, rows) reads the rows in the order they are in memory
 * instead of following each PSM's pointers in the order of the scores. The rows stay valid as
 * long as the set holds the same PSMs with the same labels, in any order.
 */
void Scores::getFeatureRows(std::vector<pair<const double*, bool> >& rows) const {
  rows.clear();
  rows.reserve(scores_.size());
  for (const ScoreHolder& sh : scores_) {
    rows.emplace_back(sh.pPSM->features, sh.isTarget());
  }
  const auto byAddress = [](const pair<const double*, bool>& a, const pair<const double*, bool>& b) {
    return std::less<const double*>()(a.first, b.first);
  };
  if (static_cast<std::ptrdiff_t>(rows.size()) >= kLargeSet) {
    sortLikeStdSort(rows, byAddress);
  } else {
    std::sort(rows.begin(), rows.end(), byAddress);
  }
}

/**
 * OpenMS: countPositives(w, fdr, skipDecoysPlusOne) with the rows of getFeatureRows of this set
 * (the count does not depend on the order of the scores, see countPositives).
 */
int Scores::countPositives(const std::vector<double>& w, double fdr, bool skipDecoysPlusOne,
                           const std::vector<pair<const double*, bool> >& rows,
                           std::vector<double>* rowScores) const {
  std::vector<pair<double, bool> > combined;
  const std::ptrdiff_t numRows = static_cast<std::ptrdiff_t>(rows.size());
  if (rowScores != nullptr) rowScores->resize(rows.size());
  double* const rs = rowScores != nullptr ? rowScores->data() : nullptr;
#ifdef _OPENMP
  if (numRows >= kLargeSet && omp_in_parallel()) {
    combined.resize(rows.size());
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(numRows / 64) shared(combined, rows, w))
    for (std::ptrdiff_t i = 0; i < numRows; ++i) {
      const pair<const double*, bool>& row = rows[static_cast<std::size_t>(i)];
      const double score = calcScore(row.first, w);
      combined[static_cast<std::size_t>(i)] = {score, row.second};
      if (rs != nullptr) rs[i] = score;
    }
  } else
#endif
  {
    combined.reserve(rows.size());
    for (std::ptrdiff_t i = 0; i < numRows; ++i) {
      const pair<const double*, bool>& row = rows[static_cast<std::size_t>(i)];
      const double score = calcScore(row.first, w);
      combined.emplace_back(score, row.second);
      if (rs != nullptr) rs[i] = score;
    }
  }
  return countPositivesOfScores(combined, w, fdr, skipDecoysPlusOne);
}

/**
 * OpenMS (CV-A2): calcScoresAndQvals(w, fdr, skipDecoysPlusOne) where the scores of w were computed
 * already, by countPositives(w, ..., rows, &rowScores) on the rows of getFeatureRows of this set:
 * each PSM gets the score of its row (calcScore of the same row and weights, the same bits), found by
 * its address; then the sort and the q-values as calcScoresAndQvals. A PSM whose row is not in rows
 * (the set changed) makes it score all PSMs with w instead.
 */
int Scores::calcScoresOfRowsAndQvals(std::vector<double>& w,
                                     const std::vector<pair<const double*, bool> >& rows,
                                     const std::vector<double>& rowScores, double fdr,
                                     bool skipDecoysPlusOne) {
  if (rowScores.size() != rows.size() || rows.size() != scores_.size()) {
    return calcScoresAndQvals(w, fdr, skipDecoysPlusOne);
  }
  const auto byAddress = [](const pair<const double*, bool>& a, const double* b) {
    return std::less<const double*>()(a.first, b);
  };
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(scores_.size());
  std::vector<char> found;
  const auto assign = [&](std::ptrdiff_t b, std::ptrdiff_t e) {
    bool ok = true;
    for (std::ptrdiff_t i = b; i < e; ++i) {
      ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
      const auto it = std::lower_bound(rows.begin(), rows.end(), sh.pPSM->features, byAddress);
      if (it == rows.end() || it->first != sh.pPSM->features) {
        ok = false;
        continue;
      }
      sh.score = rowScores[static_cast<std::size_t>(it - rows.begin())];
    }
    return ok;
  };
#ifdef _OPENMP
  if (n >= 2 * kMinPairTaskSize && omp_in_parallel()) {
    const std::ptrdiff_t blockSize = loopGrain(n, kMinPairTaskSize);
    const std::ptrdiff_t numBlocks = (n + blockSize - 1) / blockSize;
    found.assign(static_cast<std::size_t>(numBlocks), 1);
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(found, assign))
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) {
      found[static_cast<std::size_t>(block)] = assign(block * blockSize, std::min(n, (block + 1) * blockSize));
    }
  } else
#endif
  {
    found.assign(1, assign(0, n));
  }
  if (!std::all_of(found.begin(), found.end(), [](char x) { return x != 0; })) {
    return calcScoresAndQvals(w, fdr, skipDecoysPlusOne);
  }
  sortScores();
  return calcQvals(fdr, skipDecoysPlusOne);
}

// OpenMS: countPositives from the (score, isTarget) pairs of the PSMs, in any order
int Scores::countPositivesOfScores(std::vector<pair<double, bool> >& combined,
                                   const std::vector<double>& w, double fdr,
                                   bool skipDecoysPlusOne) const {
  if (std::any_of(combined.begin(), combined.end(),
                  [](const pair<double, bool>& x) { return std::isnan(x.first); })) {
    Scores copy = *this;
    std::vector<double> weights = w;
    return copy.calcScoresAndQvals(weights, fdr, skipDecoysPlusOne);
  }
  int numPos = 0;
  if (countPositivesOfTop(combined, fdr, skipDecoysPlusOne, numPos)) {
    return numPos;
  }
  std::sort(combined.begin(), combined.end(),
            [](const pair<double, bool>& a, const pair<double, bool>& b) {
              return a.first > b.first;
            });
  return countPositivesOfPairs(combined, fdr, skipDecoysPlusOne);
}

/**
 * OpenMS: the count of countPositivesOfPairs (combined sorted by score, no NaN score) from the
 * PSMs with the highest scores only, where that gives the same count; then true (combined is
 * permuted), otherwise false (only permuted).
 * With pi0 >= 1 (finite) and a finite decoy factor >= 0, getQValues makes no mix-max correction:
 * at the end of each group of equal scores, FDR = (decoys * pi0 + 0) / max(1, targets) *
 * decoyFactor (decoys and targets counted from the top, decoys from 1 without skipDecoysPlusOne),
 * and the q-value of a PSM is the minimum of min(FDR, 1) at the ends of its group and the groups
 * after it. So the positives are the targets up to the end of the last group with min(FDR, 1) <
 * fdr. The groups of a top part (all PSMs scoring at least as high as the k-th highest, sorted)
 * are walked as in getQValues; a group after the top part ends with at least the top part's
 * decoys and at most all targets, and each operation of the expression is monotone, so its FDR
 * is at least that of (decoys of the top part, all targets). If min(that, 1) >= fdr, no later
 * group counts, and the top part is sorted and walked. Otherwise k grows; once the top part is
 * all PSMs, false (the caller sorts them).
 */
bool Scores::countPositivesOfTop(std::vector<pair<double, bool> >& combined, double fdr,
                                 bool skipDecoysPlusOne, int& numPos) const {
  const double pi0 = pi0_;
  const double decoyFactor = nullTargetWinProb_ / (1.0 - nullTargetWinProb_);
  // (with fdr > 1 or NaN, no bound can exclude the later groups)
  if (!(pi0 >= 1.0) || !std::isfinite(pi0) || !(decoyFactor >= 0.0) ||
      !std::isfinite(decoyFactor) || !(fdr <= 1.0)) {
    return false;
  }
  const auto greater = [](const pair<double, bool>& a, const pair<double, bool>& b) {
    return a.first > b.first;
  };
  const std::size_t n = combined.size();
  int numTargets = 0;
  for (const pair<double, bool>& x : combined) {
    if (x.second) ++numTargets;
  }
  for (std::size_t k = std::max<std::size_t>(4096, n / 8); k < n; k *= 4) {
    std::nth_element(combined.begin(), combined.begin() + static_cast<std::ptrdiff_t>(k - 1),
                     combined.end(), greater);
    // the top part: the k highest scores and all PSMs scoring equal to the k-th (whole groups)
    const double kth = combined[k - 1].first;
    const auto topEnd = std::partition(
        combined.begin() + static_cast<std::ptrdiff_t>(k), combined.end(),
        [kth](const pair<double, bool>& x) { return x.first == kth; });
    if (topEnd == combined.end()) break;
    // the bound for the groups after the top part first (from its decoys, in any order)
    int topDecoys = skipDecoysPlusOne ? 0 : 1;
    for (auto it = combined.begin(); it != topEnd; ++it) {
      if (!it->second) ++topDecoys;
    }
    const double tailFdr =
        (topDecoys * pi0 + 0.0) / (double)((std::max)(1, numTargets)) * decoyFactor;
    if (!((std::min)(tailFdr, 1.0) >= fdr)) continue;
    if (topEnd - combined.begin() >= kLargeSet) {
      sortLikeStdSort(combined.data(), combined.data() + (topEnd - combined.begin()), greater);
    } else {
      std::sort(combined.begin(), topEnd, greater);
    }
    // the loop of getQValues on the top part
    PosteriorEstimator::setNegative(true);  // as countPositivesOfPairs
    int n_z_ge_w = skipDecoysPlusOne ? 0 : 1, n_w_ge_w = 0;
    int count = 0;
    for (auto it = combined.begin(); it != topEnd; ++it) {
      if (it->second) {
        ++n_w_ge_w;
      } else {
        ++n_z_ge_w;
      }
      if (it + 1 == topEnd || it->first != (it + 1)->first) {
        const double groupFdr =
            (n_z_ge_w * pi0 + 0.0) / (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
        if ((std::min)(groupFdr, 1.0) < fdr) count = n_w_ge_w;
      }
    }
    numPos = count;
    return true;
  }
  return false;
}

/**
 * Calculates the SVM cost/score of each PSM and sorts them
 * @param w normal vector used for SVM cost
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return 0 on successfull execution
 */
int Scores::calcScores(std::vector<double>& w) {
  // OpenMS: inside an OpenMP parallel region (CrossValidation), the PSMs are scored by the tasks
  // of the region, in blocks; each score is computed as before, by calcScore.
  const std::ptrdiff_t numScores = static_cast<std::ptrdiff_t>(scores_.size());
#ifdef _OPENMP
  if (numScores >= 2 * kMinScoreTaskSize && omp_in_parallel()) {
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numScores, kMinScoreTaskSize)) shared(w))
    for (std::ptrdiff_t i = 0; i < numScores; ++i) {
      ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
      sh.score = calcScore(sh.pPSM->features, w);
    }
  } else
#endif
  {
    std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
    for (; scoreIt != scores_.end(); ++scoreIt) {
      scoreIt->score = calcScore(scoreIt->pPSM->features, w);
    }
  }
  return sortScores();
}

/**
 * OpenMS: calcScoresAndQvals(w, fdr, skipDecoysPlusOne) for scores that calcScores(w) set
 * already and that have not changed since (in any order): the PSMs are not scored again.
 */
int Scores::sortScoresAndCalcQvals(double fdr, bool skipDecoysPlusOne) {
  // OpenMS: the sort is skipped where it cannot change what any later step sees (CrossValidation sorts
  // the training scores that the step before sorted): see sortedUpToInterchangeable
  // (the check is skipped after it has failed once for this set: on tie-heavy sets it costs as much as
  // part of the sort and does not succeed)
  if (VERB > 3 || sortCertificateFailed_ || !comparableScores(scores_.begin(), scores_.end()) ||
      !sortedUpToInterchangeable()) {
    if (VERB <= 3) sortCertificateFailed_ = true;
    sortScores();
  }
  return calcQvals(fdr, skipDecoysPlusOne);
}

// OpenMS: whether every two neighbours are either strictly ordered (ScoreHolderGreater) or
// interchangeable: bit-identical score, equal label, scan, expMass, spectrum file and bit-identical feature rows
// (only the PSM objects differ). A sort of such an order can only permute holders within runs of
// interchangeable ones (the equalities are transitive), so every later consumer of the training set
// (q-values by position, the selection of the best positive per spectrum with the same moves and
// value-identical representatives, the SVM rows in the same order, the counts) sees the same values
// at every position as after the sort; the identity of a training holder never reaches the output.
// (With NaN scores or masses the caller sorts: comparableScores.)
bool Scores::sortedUpToInterchangeable() const {
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(scores_.size());
  const ScoreHolderGreater greater;
  const std::size_t featureBytes = FeatureNames::getNumFeatures() * sizeof(double);
  const auto interchangeable = [featureBytes](const ScoreHolder& a, const ScoreHolder& b) {
    return std::memcmp(&a.score, &b.score, sizeof(double)) == 0 && a.label == b.label && a.pPSM->scan == b.pPSM->scan &&
           a.pPSM->expMass == b.pPSM->expMass && a.pPSM->specFileNr == b.pPSM->specFileNr &&
           (a.pPSM->features == b.pPSM->features ||
            std::memcmp(a.pPSM->features, b.pPSM->features, featureBytes) == 0);
  };
  const auto ordered = [&](std::ptrdiff_t b, std::ptrdiff_t e) {
    for (std::ptrdiff_t i = b; i + 1 < e; ++i) {
      const ScoreHolder& x = scores_[static_cast<std::size_t>(i)];
      const ScoreHolder& y = scores_[static_cast<std::size_t>(i + 1)];
      if (!greater(x, y) && !interchangeable(x, y)) return false;
    }
    return true;
  };
#ifdef _OPENMP
  if (n >= 2 * kMinPairTaskSize && omp_in_parallel()) {
    const std::ptrdiff_t blockSize = loopGrain(n, kMinPairTaskSize);
    const std::ptrdiff_t numBlocks = (n + blockSize - 1) / blockSize;
    std::vector<char> ok(static_cast<std::size_t>(numBlocks), 1);
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(ok, ordered))
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) {
      // a block and its link to the next one
      ok[static_cast<std::size_t>(block)] = ordered(block * blockSize, std::min(n, (block + 1) * blockSize + 1));
    }
    return std::all_of(ok.begin(), ok.end(), [](char x) { return x != 0; });
  }
#endif
  return ordered(0, n);
}

// the sort of calcScores
int Scores::sortScores() {
  std::size_t ix;
  sortScoresLikeStdSort(scores_, ScoreHolderGreater());
  if (VERB > 3) {
    if (scores_.size() >= 10) {
      cerr << "10 best scores and labels" << endl;
      for (ix = 0; ix < 10; ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
      cerr << "10 worst scores and labels" << endl;
      for (ix = scores_.size() - 10; ix < scores_.size(); ix++) {
        cerr << scores_[ix].score << " " << scores_[ix].label << endl;
      }
    } else {
      cerr << "Too few scores to display top and bottom PSMs ("
           << scores_.size() << " scores found)." << endl;
    }
  }
  return 0;
}

void Scores::getScoreLabelPairs(std::vector<pair<double, bool> >& combined) {
  combined.resize(scores_.size());
#ifdef _OPENMP
  // OpenMS: inside an OpenMP parallel region, large sets by its tasks (blocks of PSMs)
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(scores_.size());
  if (n >= 2 * kMinPairTaskSize && omp_in_parallel()) {
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(n, kMinPairTaskSize)) shared(combined))
    for (std::ptrdiff_t i = 0; i < n; ++i) {
      combined[static_cast<std::size_t>(i)] = scores_[static_cast<std::size_t>(i)].toPair();
    }
    return;
  }
#endif
  transform(scores_.begin(), scores_.end(), combined.begin(),
            mem_fn(&ScoreHolder::toPair));
}

/**
 * Calculates the q-value for each psm in scores_: the q-value is the minimal
 * FDR of any set that includes the particular psm
 * @param fdr FDR threshold specified by user (default 0.01)
 * @return number of true positives
 */
int Scores::calcQvals(double fdr, bool skipDecoysPlusOne) {
  std::vector<pair<double, bool> > combined;
  return calcQvals(fdr, skipDecoysPlusOne, combined);
}

int Scores::calcQvals(double fdr, bool skipDecoysPlusOne,
                      std::vector<pair<double, bool> >& combined) {
  assert(totalNumberOfDecoys_ + totalNumberOfTargets_ == size());

  getScoreLabelPairs(combined);

  std::vector<double> qvals;
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  // OpenMS: with the tasks of an enclosing OpenMP parallel region where it can (see there)
  if (!PosteriorEstimator::getQValuesInTasks(pi0_, combined, qvals, skipDecoysPlusOne,
                                             nullTargetWinProb_)) {
    PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne,
                                   nullTargetWinProb_);
  }

  // set q-values and count number of positives
#ifdef _OPENMP
  // OpenMS: inside an OpenMP parallel region, large sets by its tasks (blocks of PSMs; the counts
  // of the blocks are added up afterwards)
  const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(qvals.size());
  if (n >= 2 * kMinPairTaskSize && qvals.size() == scores_.size() && omp_in_parallel()) {
    const std::ptrdiff_t blockSize = loopGrain(n, kMinPairTaskSize);
    const std::ptrdiff_t numBlocks = (n + blockSize - 1) / blockSize;
    std::vector<int> positivesOfBlock(static_cast<std::size_t>(numBlocks), 0);
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(qvals, positivesOfBlock, fdr))
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) {
      const std::ptrdiff_t last = std::min(n, (block + 1) * blockSize);
      int positives = 0;
      for (std::ptrdiff_t i = block * blockSize; i < last; ++i) {
        ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
        sh.q = qvals[static_cast<std::size_t>(i)];
        if (sh.q < fdr && sh.isTarget()) ++positives;
      }
      positivesOfBlock[static_cast<std::size_t>(block)] = positives;
    }
    int numPos = 0;
    for (const int positives : positivesOfBlock) numPos += positives;
    return numPos;
  }
#endif
  std::vector<double>::const_iterator qIt = qvals.begin();
  std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();

  int numPos = 0;
  for (; qIt != qvals.end(); ++qIt, ++scoreIt) {
    scoreIt->q = *qIt;
    if (scoreIt->q < fdr && scoreIt->isTarget())
      ++numPos;
  }

  return numPos;
}

void Scores::generateNegativeTrainingSet(AlgIn& data, const double cneg) {
  std::size_t ix2 = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isDecoy()) {
      data.vals[ix2] = scoreIt->pPSM->features;
      data.Y[ix2] = -1;
      // data.C[ix2] = cneg;
      ix2++;
    }
  }
  data.negatives = static_cast<int>(ix2);
}

void Scores::generatePositiveTrainingSet(AlgIn& data,
                                         const double fdr,
                                         const double cpos,
                                         const bool trainBestPositive) {
  std::size_t ix2 = static_cast<std::size_t>(data.negatives);
  int p = 0;

  // OpenMS: the targets that pass the FDR, in the order of scores_.
  std::vector<ScoreHolder> positives;
  for (const ScoreHolder& sh : scores_) {
    if (sh.isTarget() && sh.q <= fdr) {
      positives.push_back(sh);
    }
  }
  // OpenMS: with trainBestPositive, only the best-scoring of these per spectrum
  // (file, scan). A higher score never has a higher q-value, so if any target
  // of a spectrum passes the FDR, a best-scoring one does. Upstream sorted
  // scores_ itself by file, scan and label only and kept an arbitrary PSM per
  // spectrum; and std::unique overwrote the PSMs it dropped from scores_, which
  // CrossValidation::mergeCpCnPairs trains on again (nested-xval-bins > 1).
  // stable_sort: equal scores keep the order of OrderScanLabel, on every
  // standard library.
  if (trainBestPositive) {
    // OpenMS: sorted like std::sort (see SortLikeStd.h); std::sort itself with a NaN score
    sortScoresLikeStdSort(positives, OrderScanLabel());
    positives.erase(
        std::unique(positives.begin(), positives.end(), UniqueScanLabel()),
        positives.end());
    std::stable_sort(positives.begin(), positives.end(),
                     greater<ScoreHolder>());
  }

  for (const ScoreHolder& sh : positives) {
    data.vals[ix2] = sh.pPSM->features;
    data.Y[ix2] = 1;
    // data.C[ix2] = cpos;
    ix2++;
    ++p;
  }
  data.positives = p;
  data.m = static_cast<int>(ix2);
}

/**
 * OpenMS: what generateNegativeTrainingSet(data, 1.0) and then
 * generatePositiveTrainingSet(data, fdr, 1.0, trainBestPositive) give on this set, without
 * copying or sorting the PSMs; returns false (and changes nothing) where that is not certain.
 * The negatives are the decoys in their order. The positives are the targets with q <= fdr;
 * with trainBestPositive, only the first of each spectrum (file, scan) by OrderScanLabel (the
 * highest score, then the smallest PSM id), which generatePositiveTrainingSet's sort by
 * OrderScanLabel and std::unique keep, and then sorted by score (stable_sort):
 *  - the PSMs of a spectrum are found by sorting (key, position) pairs, a total order, and the
 *    first of each by OrderScanLabel itself; if two of them are equivalent for it (equal score
 *    and PSM id), which one the sort leaves first is not certain: false;
 *  - if each positive scores strictly higher (operator>) than the next positive here, no two
 *    positives are equivalent, so every sort by score leaves them in this order. Otherwise, and
 *    with a NaN score or mass, false.
 */
bool Scores::generateTrainingSetsOfSortedSet(AlgIn& data, const double fdr,
                                             const bool trainBestPositive) const {
  // kept[i]: whether the PSM scores_[i] is kept (with trainBestPositive: the first of its spectrum)
  std::vector<char> kept(scores_.size(), 1);
  if (trainBestPositive && !scores_.empty()) {
    if (!comparableScores(scores_.begin(), scores_.end()) || scores_.size() >= (std::size_t(1) << 32)) return false;
    // OpenMS (CV-A1, BP-PAR): one 64-bit key (file << 32 | scan) per candidate (a target with q <= fdr), kNone for
    // the other PSMs, with the position; in blocks as tasks inside a parallel region
    struct PackedKey {
      std::uint64_t key;
      std::uint32_t index;
    };
    constexpr std::uint64_t kNone = ~std::uint64_t(0);
    const std::ptrdiff_t sn = static_cast<std::ptrdiff_t>(scores_.size());
    std::ptrdiff_t blockSize = sn;
#ifdef _OPENMP
    const bool inTasks = sn >= 2 * kMinPairTaskSize && omp_in_parallel();
    if (inTasks) blockSize = loopGrain(sn, kMinPairTaskSize);
#endif
    const std::ptrdiff_t numBlocks = (sn + blockSize - 1) / blockSize;
    std::vector<PackedKey> packed(scores_.size());
    // per block: whether every key fits (file < 2^31: no key is kNone), whether a spectrum's first is certain
    std::vector<char> blockOk(static_cast<std::size_t>(numBlocks), 1);
    const auto fill = [&](std::ptrdiff_t block) {
      const std::ptrdiff_t b = block * blockSize, e = std::min(sn, (block + 1) * blockSize);
      for (std::ptrdiff_t i = b; i < e; ++i) {
        const ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
        std::uint64_t key = kNone;
        if (sh.isTarget() && sh.q <= fdr) {
          if (!(sh.pPSM->specFileNr < (1u << 31))) blockOk[static_cast<std::size_t>(block)] = 0;
          key = (std::uint64_t(sh.pPSM->specFileNr) << 32) | std::uint64_t(sh.pPSM->scan);
        }
        packed[static_cast<std::size_t>(i)] = {key, static_cast<std::uint32_t>(i)};
      }
    };
    // the first by OrderScanLabel of each spectrum whose run of keys starts in the block (the run may end in a later
    // block); the others are not kept
    const auto pick = [&](std::ptrdiff_t block) {
      const std::ptrdiff_t b = block * blockSize, e = std::min(sn, (block + 1) * blockSize);
      for (std::ptrdiff_t k = b; k < e; ++k) {
        const std::uint64_t key = packed[static_cast<std::size_t>(k)].key;
        if (key == kNone) break;  // (sorted last)
        if (k > 0 && packed[static_cast<std::size_t>(k - 1)].key == key) continue;
        std::size_t first = packed[static_cast<std::size_t>(k)].index;
        bool certain = true;
        for (std::ptrdiff_t m = k + 1; m < sn && packed[static_cast<std::size_t>(m)].key == key; ++m) {
          const std::size_t other = packed[static_cast<std::size_t>(m)].index;
          if (OrderScanLabel()(scores_[other], scores_[first])) {
            kept[first] = 0;
            first = other;
            certain = true;
          } else {
            kept[other] = 0;
            if (!OrderScanLabel()(scores_[first], scores_[other])) certain = false;
          }
        }
        if (!certain) blockOk[static_cast<std::size_t>(block)] = 0;
      }
    };
#ifdef _OPENMP
    if (inTasks) {
      PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(fill))
      for (std::ptrdiff_t block = 0; block < numBlocks; ++block) fill(block);
    } else
#endif
    {
      for (std::ptrdiff_t block = 0; block < numBlocks; ++block) fill(block);
    }
    if (std::any_of(blockOk.begin(), blockOk.end(), [](char x) { return x == 0; })) return false;
    sortLikeStdSort(packed, [](const PackedKey& x, const PackedKey& y) {
      return x.key < y.key || (x.key == y.key && x.index < y.index);
    });
#ifdef _OPENMP
    if (inTasks) {
      PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(pick))
      for (std::ptrdiff_t block = 0; block < numBlocks; ++block) pick(block);
    } else
#endif
    {
      for (std::ptrdiff_t block = 0; block < numBlocks; ++block) pick(block);
    }
    if (std::any_of(blockOk.begin(), blockOk.end(), [](char x) { return x == 0; })) return false;
  }
  // OpenMS (SETUP-PAR): the positives and the decoys of blocks of PSMs (as tasks inside a parallel region),
  // then the blocks in order: the same positives, the same check, the same rows in the same order
  const std::ptrdiff_t total = static_cast<std::ptrdiff_t>(scores_.size());
  std::ptrdiff_t blockSize = total > 0 ? total : 1;
#ifdef _OPENMP
  const bool inTasks = total >= 2 * kMinPairTaskSize && omp_in_parallel();
  if (inTasks) blockSize = loopGrain(total, kMinPairTaskSize);
#endif
  const std::ptrdiff_t numBlocks = (total + blockSize - 1) / blockSize;
  std::vector<std::vector<std::size_t> > positivesOfBlock(static_cast<std::size_t>(numBlocks));
  std::vector<std::size_t> decoysOfBlock(static_cast<std::size_t>(numBlocks), 0);
  std::vector<char> orderedBlock(static_cast<std::size_t>(numBlocks), 1);
  const auto collect = [&](std::ptrdiff_t block) {
    std::vector<std::size_t>& pos = positivesOfBlock[static_cast<std::size_t>(block)];
    std::size_t decoys = 0;
    const std::ptrdiff_t last = std::min(total, (block + 1) * blockSize);
    for (std::ptrdiff_t j = block * blockSize; j < last; ++j) {
      const std::size_t i = static_cast<std::size_t>(j);
      const ScoreHolder& sh = scores_[i];
      if (sh.isDecoy()) ++decoys;
      if (kept[i] && sh.isTarget() && sh.q <= fdr) {
        if (trainBestPositive && !pos.empty() && !(scores_[pos.back()] > sh)) {
          orderedBlock[static_cast<std::size_t>(block)] = 0;
        }
        pos.push_back(i);
      }
    }
    decoysOfBlock[static_cast<std::size_t>(block)] = decoys;
  };
#ifdef _OPENMP
  if (inTasks) {
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(collect))
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) collect(block);
  } else
#endif
  {
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) collect(block);
  }
  // the check of each positive against the one before it: inside the blocks, and across them
  std::size_t numNeg = 0, numPos = 0;
  const std::size_t* prev = nullptr;
  for (std::ptrdiff_t block = 0; block < numBlocks; ++block) {
    const std::vector<std::size_t>& pos = positivesOfBlock[static_cast<std::size_t>(block)];
    if (!orderedBlock[static_cast<std::size_t>(block)]) return false;
    if (trainBestPositive && prev != nullptr && !pos.empty() && !(scores_[*prev] > scores_[pos.front()])) {
      return false;
    }
    if (!pos.empty()) prev = &pos.back();
    numNeg += decoysOfBlock[static_cast<std::size_t>(block)];
    numPos += pos.size();
  }
  // the rows: the decoys in order, then the positives in order, at the offsets of their blocks
  std::vector<std::size_t> negOffset(static_cast<std::size_t>(numBlocks)), posOffset(static_cast<std::size_t>(numBlocks));
  for (std::size_t block = 0, neg = 0, pos = numNeg; block < static_cast<std::size_t>(numBlocks); ++block) {
    negOffset[block] = neg;
    posOffset[block] = pos;
    neg += decoysOfBlock[block];
    pos += positivesOfBlock[block].size();
  }
  const auto emit = [&](std::ptrdiff_t block) {
    std::size_t ix2 = negOffset[static_cast<std::size_t>(block)];
    const std::ptrdiff_t last = std::min(total, (block + 1) * blockSize);
    for (std::ptrdiff_t j = block * blockSize; j < last; ++j) {
      const ScoreHolder& sh = scores_[static_cast<std::size_t>(j)];
      if (sh.isDecoy()) {
        data.vals[ix2] = sh.pPSM->features;
        data.Y[ix2] = -1;
        ix2++;
      }
    }
    ix2 = posOffset[static_cast<std::size_t>(block)];
    for (std::size_t i : positivesOfBlock[static_cast<std::size_t>(block)]) {
      data.vals[ix2] = scores_[i].pPSM->features;
      data.Y[ix2] = 1;
      ix2++;
    }
  };
#ifdef _OPENMP
  if (inTasks) {
    PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(emit))
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) emit(block);
  } else
#endif
  {
    for (std::ptrdiff_t block = 0; block < numBlocks; ++block) emit(block);
  }
  data.negatives = static_cast<int>(numNeg);
  data.positives = static_cast<int>(numPos);
  data.m = static_cast<int>(numNeg + numPos);
  return true;
}

void Scores::weedOutRedundant() {
  std::map<std::string, unsigned int> peptideSpecCounts;
  double specCountQvalThreshold = -1.0;
  weedOutRedundant(peptideSpecCounts, specCountQvalThreshold);
}

/**
 * Routine that sees to that only unique peptides are kept (used for analysis
 * on peptide-fdr rather than psm-fdr)
 */
void Scores::weedOutRedundant(
    std::map<std::string, unsigned int>& peptideSpecCounts,
    double specCountQvalThreshold) {
  // lexicographically order the scores_ (based on peptides names,labels and
  // scores)
  std::sort(scores_.begin(), scores_.end(), lexicOrderProb());

  /*
   * much simpler version but it does not fill up the peptide-PSM map:
   * scores_.erase(std::unique(scores_.begin(), scores_.end(), mycmp),
   * scores_.end());
   */

  std::string previousPeptide = "";
  LabelType previousLabel = LabelType::UNDEFINED;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx) {
    std::string currentPeptide = scores_[idx].pPSM->getPeptideSequence();
    LabelType currentLabel = scores_[idx].label;
    if (currentPeptide != previousPeptide || currentLabel != previousLabel) {
      // insert as a new score
      scores_[lastWrittenIdx++] = scores_[idx];
      previousPeptide = currentPeptide;
      previousLabel = currentLabel;
    }
    // append the psm
    peptidePsmMap_[scores_[lastWrittenIdx - 1].pPSM].push_back(
        scores_[idx].pPSM);
    if (specCountQvalThreshold > 0.0 &&
        scores_[idx].q < specCountQvalThreshold) {
      ++peptideSpecCounts[currentPeptide];
    }
  }
  scores_.resize(lastWrittenIdx);
  postMergeStep();
}

/**
 * Routine that sees to that only unique spectra are kept for TDC
 */
void Scores::weedOutRedundantTDC() {
  // order the scores (based on spectra id and score)
  // OpenMS: like std::sort, with the tasks of an enclosing OpenMP parallel region (see
  // SortLikeStd.h), if the comparison is a strict weak ordering (no NaN mass or score)
  if (comparableScores(scores_.begin(), scores_.end())) {
    sortLikeStdSort(scores_, OrderScanMassCharge());
  } else {
    std::sort(scores_.begin(), scores_.end(), OrderScanMassCharge());
  }
  scores_.erase(
      std::unique(scores_.begin(), scores_.end(), UniqueScanMassCharge()),
      scores_.end());

  /* does not actually release memory because of memory fragmentation
  double previousExpMass = 0.0;
  unsigned int previousScan = 0u;
  size_t lastWrittenIdx = 0u;
  for (size_t idx = 0u; idx < scores_.size(); ++idx){
    double currentExpMass = scores_.at(idx).pPSM->expMass;
    int currentScan = scores_.at(idx).pPSM->scan;
    if (currentExpMass != previousExpMass || currentScan != previousScan) {
      // insert as a new score
      scores_.at(lastWrittenIdx++).swap(scores_.at(idx));
      previousScan = currentScan;
      previousExpMass = currentExpMass;
    } else {
      PSMDescription::deletePtr(scores_.at(idx).pPSM);
    }
  }
  scores_.resize(lastWrittenIdx);
  */
  postMergeStep();
}

/**
 * Routine that sees to that only 1 target and 1 decoy spectra are kept for
 * mix-max when using multiple hits per spectrum and separate searches
 */
void Scores::weedOutRedundantMixMax() {
  // order the scores (based on spectra id and score)
  std::sort(scores_.begin(), scores_.end(), OrderScanMassLabelCharge());
  scores_.erase(
      std::unique(scores_.begin(), scores_.end(), UniqueScanMassLabelCharge()),
      scores_.end());

  postMergeStep();
}

// OpenMS: the number of positives calcQvals(fdr, skipDecoysPlusOne) returns, without setting
// the q-values; combined holds (score, isTarget) in the order calcQvals would see.
int Scores::countPositivesOfPairs(const std::vector<pair<double, bool> >& combined,
                               double fdr, bool skipDecoysPlusOne) const {
  std::vector<double> qvals;
  return countPositivesOfPairs(combined, fdr, skipDecoysPlusOne, qvals);
}

int Scores::countPositivesOfPairs(const std::vector<pair<double, bool> >& combined,
                                  double fdr, bool skipDecoysPlusOne,
                                  std::vector<double>& qvals) const {
  // Without the mix-max correction (pi0 not below 1), getQValues gives every PSM of a group of equal scores the
  // minimum of min(FDR, 1) at the ends of its group and of the groups after it (FDR = decoys * pi0 / max(1, targets)
  // * decoyFactor, counted from the top as below). So q < fdr exactly for the PSMs up to the end of the last group
  // with min(FDR, 1) < fdr, and the count is the number of targets up to there: walked without the q-values
  // (as countPositivesOfTop). A NaN FDR (which the minimum does not order) uses getQValues. Only where doubles are
  // evaluated in double precision (FLT_EVAL_METHOD 0, e.g. SSE2 / AArch64): getQValues rounds every q-value to double
  // when it stores it, which an x87 build with excess precision would not do in registers.
#if defined(FLT_EVAL_METHOD) && FLT_EVAL_METHOD == 0
  if (!(pi0_ < 1.0)) {
    const double pi0 = pi0_;
    const double decoyFactor = nullTargetWinProb_ / (1.0 - nullTargetWinProb_);
    PosteriorEstimator::setNegative(true);  // as below
    int n_z_ge_w = skipDecoysPlusOne ? 0 : 1, n_w_ge_w = 0;
    int numPos = 0;
    bool nanFdr = false;
    const std::size_t n = combined.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (combined[i].second) {
        ++n_w_ge_w;
      } else {
        ++n_z_ge_w;
      }
      if (i + 1 == n || combined[i].first != combined[i + 1].first) {
        const double groupFdr = (n_z_ge_w * pi0 + 0.0) / (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
        const double q = (std::min)(groupFdr, 1.0);
        if (std::isnan(q)) {
          nanFdr = true;
          break;
        }
        if (q < fdr) numPos = n_w_ge_w;
      }
    }
    if (!nanFdr) return numPos;
  }
#endif
  // (one q-value per pair: room for all of them at once, instead of growing while they are added)
  qvals.clear();
  qvals.reserve(combined.size());
  PosteriorEstimator::setNegative(true);  // also get q-values for decoys
  PosteriorEstimator::getQValues(pi0_, combined, qvals, skipDecoysPlusOne,
                                 nullTargetWinProb_);
  int numPos = 0;
  for (std::size_t i = 0; i < qvals.size(); ++i) {
    if (qvals[i] < fdr && combined[i].second) ++numPos;
  }
  return numPos;
}

int Scores::getInitDirection(const double initialSelectionFdr,
                             std::vector<double>& direction,
                             std::ostream& report) {
  int bestPositives = -1;
  int bestFeature = -1;
  bool lowBest = false;

  // for determining the initial direction, the decoys+1 in the FDR estimates
  // is too restrictive for small datasets
  bool skipDecoysPlusOne = true;

  const unsigned int numFeatures =
      static_cast<unsigned int>(FeatureNames::getNumFeatures());
  // positives[2 * featNo + i]: in forward (i = 0) and backward (i = 1) direction
  std::vector<int> positives(2 * numFeatures);

  // OpenMS: inside an OpenMP parallel region (CrossValidation), the scores are sorted by every
  // feature in turn as below (sorting starts from the order the previous feature left, which
  // decides the order of PSMs that compare equal, and the order after the last feature is where
  // the training starts from), with the tasks of the region (see SortLikeStd.h), on keys (see
  // permuteScores): the feature value, the rank of the PSM's (scan, expMass, label) among the
  // PSMs (what operator< compares after the score), and the position. The PSMs are put into the
  // final order (and get the scores of the last feature) once, at the end. After each sort, a task counts the
  // positives of the feature on the (score, label) pairs in that order: the q-values, and so
  // the count, depend on the order of the scores and on the labels of each group of equal
  // scores, not on the order within a group (calcQvals treats ties as one group). The q-values
  // are set for the last feature, as by the loop below. A NaN feature value or mass (no strict
  // weak ordering) uses the loop below.
  bool counted = false;
#ifdef _OPENMP
  const bool anyNaN = std::any_of(scores_.begin(), scores_.end(), [numFeatures](const ScoreHolder& sh) {
    return std::isnan(sh.pPSM->expMass) ||
           std::any_of(sh.pPSM->features, sh.pPSM->features + numFeatures,
                       [](double x) { return std::isnan(x); });
  });
  const std::size_t n = scores_.size();
  if (numFeatures > 0 && !anyNaN && n > 0 && n <= std::numeric_limits<unsigned int>::max() &&
      omp_in_parallel()) {
    std::exception_ptr chainError;
    std::vector<std::exception_ptr> errors(numFeatures);
    struct ChainKey {
      double score;
      unsigned int tieRank;
      unsigned int index;
    };
    // OpenMS: the buffers of the counts (the n (score, label) pairs in the order of a feature's sorted keys, and
    // the q-value scratch), taken from a list and put back when a count is done, so the counts of later features
    // reuse memory that is already in use; released at the end of this call. pairsOfFeature hands a feature's pairs
    // from the chain to its count task.
    struct CountBuffers {
      std::vector<pair<double, bool> > pairs;
      std::vector<double> qvals;
    };
    std::mutex buffersMutex;
    std::vector<std::unique_ptr<CountBuffers> > freeCountBuffers;
    std::vector<std::unique_ptr<CountBuffers> > pairsOfFeature(numFeatures);
    std::vector<char> isTarget(n);
    // OpenMS: the copies and moves of keys of a large set below are tasks of blocks of this size
    // (each key is written once, so the order of the blocks does not matter)
    constexpr std::ptrdiff_t kMinChainTaskSize = 16384;
    const std::ptrdiff_t numKeys = static_cast<std::ptrdiff_t>(n);
    PERCOLATOR_TASK_PRAGMA(omp taskgroup)
    {
      try {
        // the rank of (scan, expMass, label) in the lexicographic order of operator<
        struct TieKey {
          unsigned int scan;
          double expMass;
          LabelType label;
          unsigned int index;
        };
        const auto tieLess = [](const TieKey& x, const TieKey& y) {
          if (x.scan != y.scan) return x.scan < y.scan;
          if (x.expMass != y.expMass) return x.expMass < y.expMass;
          return x.label < y.label;
        };
        std::vector<TieKey> ties;
        ties.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
          ties.push_back({scores_[i].pPSM->scan, scores_[i].pPSM->expMass, scores_[i].label,
                          static_cast<unsigned int>(i)});
        }
        sortLikeStdSort(ties, tieLess);  // (any sort: only the ranks are used)
        std::vector<unsigned int> tieRank(n);
        unsigned int rank = 0;
        for (std::size_t k = 0; k < n; ++k) {
          if (k > 0 && tieLess(ties[k - 1], ties[k])) ++rank;
          tieRank[ties[k].index] = rank;
        }
        // the feature values by feature, and the labels (blocks of PSMs as tasks; every value is
        // written there, so the array is not zero-filled first)
        // (a block of PSMs per task, column by column: consecutive writes)
        const std::unique_ptr<double[]> columns(new double[static_cast<std::size_t>(numFeatures) * n]);
        adviseHugePages(columns.get(), static_cast<std::size_t>(numFeatures) * n * sizeof(double));  // (ssl.h)
        const std::ptrdiff_t blockSize = loopGrain(numKeys, kMinScoreTaskSize);
        const std::ptrdiff_t numBlocks = (numKeys + blockSize - 1) / blockSize;
        PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(columns, isTarget))
        for (std::ptrdiff_t block = 0; block < numBlocks; ++block) {
          const std::size_t first = static_cast<std::size_t>(block * blockSize);
          const std::size_t last = std::min(n, first + static_cast<std::size_t>(blockSize));
          for (std::size_t i = first; i < last; ++i) {
            isTarget[i] = scores_[i].isTarget();
          }
          for (unsigned int featNo = 0; featNo < numFeatures; ++featNo) {
            double* column = &columns[featNo * n];
            for (std::size_t i = first; i < last; ++i) {
              column[i] = scores_[i].pPSM->features[featNo];
            }
          }
        }
        std::vector<ChainKey> keys;
        keys.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
          keys.push_back({0.0, tieRank[i], static_cast<unsigned int>(i)});
        }
        // OpenMS: packed keys. A feature with at most kMaxPackedCodes distinct values is sorted on 64-bit words
        // (code of the value << (tie bits + index bits) | tieRank << index bits | index), compared on their bits above
        // the index: the code is the rank of the value among the feature's distinct values (equal values, -0.0 and
        // +0.0 included, get equal codes, a smaller value a smaller code), so every comparison of the sort has the
        // outcome of (score, tieRank) and the sort makes the same moves on half the bytes. Afterwards the keys get their
        // own values back. (No NaN here; the other features are sorted as ChainKeys.)
        constexpr std::size_t kMaxPackedCodes = std::size_t(1) << 16;
        const auto bitsFor = [](std::uint64_t maxValue) {
          unsigned int b = 0;
          while (b < 64 && (maxValue >> b) != 0) ++b;
          return b;
        };
        const unsigned int indexBits = bitsFor(n - 1);
        const unsigned int tieBits = bitsFor(rank);
        const bool packable = indexBits + tieBits + 16 <= 64;
        std::vector<std::vector<double> > valuesOfFeature(numFeatures);  // sorted distinct values; empty = not packed
        if (packable) {
          PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(1) shared(valuesOfFeature, columns))
          for (unsigned int f = 0; f < numFeatures; ++f) {
            // (no exception may leave a task: if an allocation fails, the feature is sorted as ChainKeys, which makes
            // the same moves; int14 review r3)
            try {
              const double* column = &columns[static_cast<std::size_t>(f) * n];
              std::unordered_set<double> distinct;
              bool few = true;
              for (std::size_t i = 0; i < n && few; ++i) {
                distinct.insert(column[i]);
                few = distinct.size() <= kMaxPackedCodes;
              }
              if (few) {
                std::vector<double> values(distinct.begin(), distinct.end());
                std::sort(values.begin(), values.end());
                valuesOfFeature[f].swap(values);
              }
            } catch (...) {
              valuesOfFeature[f].clear();
            }
          }
        }
        std::vector<std::uint64_t> packed;
        for (unsigned int featNo = 0; featNo < numFeatures; featNo++) {
          const double* column = &columns[featNo * n];
          if (!valuesOfFeature[featNo].empty()) {
            const std::vector<double>& values = valuesOfFeature[featNo];
            const unsigned int shift = tieBits + indexBits;
            packed.resize(n);
            PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys, packed, values))
            for (std::ptrdiff_t i = 0; i < numKeys; ++i) {
              const ChainKey& key = keys[static_cast<std::size_t>(i)];
              const std::uint64_t code = static_cast<std::uint64_t>(
                  std::lower_bound(values.begin(), values.end(), column[key.index]) - values.begin());
              packed[static_cast<std::size_t>(i)] =
                  (code << shift) | (static_cast<std::uint64_t>(key.tieRank) << indexBits) | key.index;
            }
            sortLikeStdSort(packed, [indexBits](std::uint64_t x, std::uint64_t y) {
              return (x >> indexBits) < (y >> indexBits);
            });
            const std::uint64_t indexMask = (std::uint64_t(1) << indexBits) - 1;
            PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys, packed, tieRank))
            for (std::ptrdiff_t i = 0; i < numKeys; ++i) {
              const unsigned int index = static_cast<unsigned int>(packed[static_cast<std::size_t>(i)] & indexMask);
              keys[static_cast<std::size_t>(i)] = {column[index], tieRank[index], index};
            }
          } else {
          PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys))
          for (std::ptrdiff_t i = 0; i < numKeys; ++i) {
            ChainKey& key = keys[static_cast<std::size_t>(i)];
            key.score = column[key.index];
          }
          // operator< of ScoreHolder
          sortLikeStdSort(keys, [](const ChainKey& x, const ChainKey& y) {
            return x.score < y.score || (x.score == y.score && x.tieRank < y.tieRank);
          });
          }
          {
            // OpenMS: the (score, label) pairs of the count, in the order of the sorted keys, written by tasks from
            // the keys directly (a pooled buffer keeps its n elements and is overwritten)
            std::unique_ptr<CountBuffers> buffers;
            {
              std::lock_guard<std::mutex> lock(buffersMutex);
              if (!freeCountBuffers.empty()) {
                buffers = std::move(freeCountBuffers.back());
                freeCountBuffers.pop_back();
              }
            }
            if (!buffers) {
              buffers.reset(new CountBuffers());
              buffers->pairs.reserve(n);
              adviseHugePages(buffers->pairs.data(), n * sizeof(pair<double, bool>));
              buffers->pairs.resize(n);
            }
            pair<double, bool>* pairs = buffers->pairs.data();
            PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys, isTarget))
            for (std::ptrdiff_t i = 0; i < numKeys; ++i) {
              const ChainKey& key = keys[static_cast<std::size_t>(i)];
              pairs[i] = pair<double, bool>(key.score, isTarget[key.index] != 0);
            }
            pairsOfFeature[featNo] = std::move(buffers);
          }
          PERCOLATOR_TASK_PRAGMA(omp task firstprivate(featNo) shared(pairsOfFeature, positives, errors, initialSelectionFdr, skipDecoysPlusOne, buffersMutex, freeCountBuffers))
          try {
            std::unique_ptr<CountBuffers> buffers = std::move(pairsOfFeature[featNo]);
            std::vector<pair<double, bool> >& pairs = buffers->pairs;
            positives[2 * featNo] = countPositivesOfPairs(pairs, initialSelectionFdr, skipDecoysPlusOne, buffers->qvals);
            reverse(pairs.begin(), pairs.end());
            positives[2 * featNo + 1] = countPositivesOfPairs(pairs, initialSelectionFdr, skipDecoysPlusOne, buffers->qvals);
            std::lock_guard<std::mutex> lock(buffersMutex);
            freeCountBuffers.push_back(std::move(buffers));
          } catch (...) {
            errors[featNo] = std::current_exception();
          }
          // reversed (the pairs of positions as tasks)
          PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys))
          for (std::ptrdiff_t i = 0; i < numKeys / 2; ++i) {
            std::swap(keys[static_cast<std::size_t>(i)], keys[static_cast<std::size_t>(numKeys - 1 - i)]);
          }
        }
        // the PSMs in the order of the keys (as permuteScores), with the scores of the last feature
        // (blocks as tasks)
        std::vector<ScoreHolder> permuted(n);
        PERCOLATOR_TASK_PRAGMA(omp taskloop grainsize(loopGrain(numKeys, kMinChainTaskSize)) shared(keys, permuted))
        for (std::ptrdiff_t i = 0; i < numKeys; ++i) {
          const ChainKey& key = keys[static_cast<std::size_t>(i)];
          permuted[static_cast<std::size_t>(i)] = scores_[key.index];
          permuted[static_cast<std::size_t>(i)].score = key.score;
        }
        scores_.swap(permuted);
        std::vector<pair<double, bool> > combined;
        calcQvals(initialSelectionFdr, skipDecoysPlusOne, combined);
      } catch (...) {
        chainError = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors) {
      if (error) std::rethrow_exception(error);
    }
    if (chainError) std::rethrow_exception(chainError);
    counted = true;
  }
#endif
  if (!counted) {
    // pre-allocate once and reuse across all 2*numFeatures calcQvals calls
    std::vector<pair<double, bool> > combined;
    combined.reserve(scores_.size());

    for (unsigned int featNo = 0; featNo < numFeatures; featNo++) {
      for (std::vector<ScoreHolder>::iterator scoreIt = scores_.begin();
           scoreIt != scores_.end(); ++scoreIt) {
        scoreIt->score = scoreIt->pPSM->features[featNo];
      }
      sort(scores_.begin(), scores_.end());
      // check once in forward direction (i = 0, higher scores are better) and
      // once in backward direction (i = 1, lower scores are better)
      for (int i = 0; i < 2; i++) {
        if (i == 1) {
          reverse(scores_.begin(), scores_.end());
        }
        positives[2 * featNo + i] = calcQvals(initialSelectionFdr, skipDecoysPlusOne, combined);
      }
    }
  }
  for (unsigned int featNo = 0; featNo < numFeatures; featNo++) {
    for (int i = 0; i < 2; i++) {
      if (positives[2 * featNo + i] > bestPositives) {
        bestPositives = positives[2 * featNo + i];
        bestFeature = static_cast<int>(featNo);
        lowBest = (i == 0);
      }
    }
  }
  for (std::size_t ix = FeatureNames::getNumFeatures(); ix--;) {
    direction[ix] = 0;
  }

  if (bestPositives <= 0) {
    ostringstream oss;
    oss << "Error in the input data: cannot find an initial direction with "
        << "positive training examples. "
        << "Consider setting/raising the initial training FDR threshold "
           "(--train-fdr-initial)."
        << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str();
      std::cerr << "No-terminate flag set: setting initial direction to the "
                << "first feature and ignoring the error." << std::endl;
      bestFeature = 0;
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }

  if (bestFeature >= 0) {
    direction[static_cast<std::size_t>(bestFeature)] = (lowBest ? -1 : 1);
  }

  if (VERB > 1) {
    report << "Selected feature " << bestFeature + 1
         << " as initial direction. Could separate " << bestPositives
         << " training set positives with q<" << initialSelectionFdr
         << " in that direction." << endl;
  }
  return bestPositives;
}

void Scores::checkSeparationAndSetPi0() {
  std::vector<pair<double, bool> > combined;
  getScoreLabelPairs(combined);

  std::vector<double> pvals;
  PosteriorEstimator::getPValues(combined, pvals);

  pi0_ = 1.0;
  bool tooGoodSeparation = PosteriorEstimator::checkSeparation(pvals);
  if (tooGoodSeparation) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (usePi0_) {
        std::cerr
            << "No-terminate flag set: setting pi0 = 1 and ignoring error."
            << std::endl;
      } else {
        std::cerr << "No-terminate flag set: ignoring error." << std::endl;
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  } else if (usePi0_) {
    pi0_ = PosteriorEstimator::estimatePi0(pvals);
  }
}

void Scores::calcPep(const bool spline, const bool interp, const bool pava) {
    if (!spline) {
        // Derive PEPs from q-values so that average PEP ≈ q at every
        // threshold.  Use PAVA when --pava-pep is set, otherwise I-spline
        // for smooth output.
        const bool use_ispline = !pava;
        std::vector<double> target_q, sc;
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                target_q.push_back(sh.q);
                sc.push_back(sh.score);
            }
        }
        InferPEP reg(use_ispline);
        auto target_pep = interp
                            ? reg.qns_to_pep(target_q, sc)
                            : reg.q_to_pep(target_q);
        // Move PEPs to scoreholders. The PEPs are only defined for target,
        // We use interpolation for decoys.
        // Add elements avoiding overflow problems if last sh is a decoy
        target_pep.push_back(1.0); target_q.push_back(1.0);
        auto it_pep = target_pep.begin();
        auto it_q = target_q.begin();
        double l_q(0.0), l_pep(0.0);
        for (auto& sh : scores_) {
            if (sh.isTarget()) {
                sh.pep = *it_pep;
                // remember last (l_) pep and q for interpolation
                l_pep = *it_pep;
                l_q = *it_q;
                it_pep++; it_q++;
            } else {
                double pep = reg.interpolate(sh.q,l_q,*it_q,l_pep,*it_pep);
                sh.pep = std::max(0.0, std::min(1.0, pep));
            }
        }
    } else {
        std::vector<pair<double, bool> > combined;
        getScoreLabelPairs(combined);

        std::vector<double> peps;
        // Logistic regression on the data
        PosteriorEstimator::estimatePEP(combined, usePi0_, pi0_, peps, true);
        for (size_t ix = 0; ix < scores_.size(); ix++) {
            scores_[ix].pep = peps[ix];
        }
    }
}

unsigned Scores::getQvaluesBelowLevel(double level) {
  unsigned hits = 0;
  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != scores_.end(); ++scoreIt) {
    if (scoreIt->isTarget() && scoreIt->q < level) {
      hits++;
    }
  }
  return hits;
}

void Scores::setUsePi0(bool usePi0) {
  usePi0_ = usePi0;
}

}}}  // namespace OpenMS::Internal::Percolator
