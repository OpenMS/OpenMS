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
#include <exception>
#include <functional>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

// OpenMS: the smallest block of PSMs a task scores (Scores::calcScores)
[[maybe_unused]] constexpr std::ptrdiff_t kMinScoreTaskSize = 2048;

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
template <typename Key>
static void permuteScores(std::vector<ScoreHolder>& scores, const std::vector<Key>& keys) {
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
                                      unsigned int decoysPerTarget) {
  // set the number of cross validation folds for train and test to xval_fold
  train.resize(xval_fold, Scores(usePi0_));
  test.resize(xval_fold, Scores(usePi0_));
  size_t approxTestSize = scores_.size() / xval_fold + 1;
  size_t approxTrainSize = (approxTestSize < scores_.size()) ? scores_.size() - approxTestSize : scores_.size();
  for (unsigned int i = 0; i < xval_fold; ++i) {
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
    std::vector<HashKey> keys;
    keys.reserve(scores_.size());
    for (std::size_t i = 0; i < scores_.size(); ++i) {
      const PSMDescription* psm = scores_[i].pPSM;
      const size_t hash = fast_uint_hash(psm->specFileNr) ^ (fast_uint_hash(psm->scan) << 1);
      keys.push_back({hash, i});
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
#ifdef _OPENMP
  const bool fillInTasks = xval_fold > 1u && omp_in_parallel();
#else
  const bool fillInTasks = false;
#endif
  if (fillInTasks) {
    foldOfPSM.reserve(scores_.size());
    labelOfPSM.reserve(scores_.size());
  }
  unsigned int previousSpectrum = scores_.begin()->pPSM->scan;
  size_t randIndex = PseudoRandom::lcg_rand() % xval_fold;
  for (std::vector<ScoreHolder>::iterator it = scores_.begin();
       it != scores_.end(); ++it) {
    const unsigned int curScan = (*it).pPSM->scan;
    ScoreHolder sh = (*it);

    // if current score is from a different spectra than the one encountered in
    // the previous iteration, choose new fold
    if (previousSpectrum != curScan) {
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
#pragma omp taskloop grainsize(1) shared(train, test, foldOfPSM, labelOfPSM, fillErrors)
    for (int i = 0; i < static_cast<int>(xval_fold); ++i) {
      try {
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
  std::sort(rows.begin(), rows.end(),
            [](const pair<const double*, bool>& a, const pair<const double*, bool>& b) {
              return std::less<const double*>()(a.first, b.first);
            });
}

/**
 * OpenMS: countPositives(w, fdr, skipDecoysPlusOne) with the rows of getFeatureRows of this set
 * (the count does not depend on the order of the scores, see countPositives).
 */
int Scores::countPositives(const std::vector<double>& w, double fdr, bool skipDecoysPlusOne,
                           const std::vector<pair<const double*, bool> >& rows) const {
  std::vector<pair<double, bool> > combined;
  combined.reserve(rows.size());
  for (const pair<const double*, bool>& row : rows) {
    combined.emplace_back(calcScore(row.first, w), row.second);
  }
  return countPositivesOfScores(combined, w, fdr, skipDecoysPlusOne);
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
    std::sort(combined.begin(), topEnd, greater);
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
#pragma omp taskloop grainsize(kMinScoreTaskSize) shared(w)
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
  sortScores();
  return calcQvals(fdr, skipDecoysPlusOne);
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

  std::vector<ScoreHolder>::iterator lastUniqueIt = scores_.end();
  if (trainBestPositive) {
    // OpenMS: sorted like std::sort (see SortLikeStd.h); file and scan numbers and labels are
    // no NaN
    sortLikeStdSort(scores_, OrderScanLabel());
    lastUniqueIt =
        std::unique(scores_.begin(), scores_.end(), UniqueScanLabel());
    if (!comparableScores(scores_.begin(), lastUniqueIt)) {
      std::sort(scores_.begin(), lastUniqueIt, greater<ScoreHolder>());
    } else if (lastUniqueIt != scores_.begin()) {
      sortLikeStdSort(&*scores_.begin(), &*scores_.begin() + (lastUniqueIt - scores_.begin()),
                      greater<ScoreHolder>());
    }
  }

  std::vector<ScoreHolder>::const_iterator scoreIt = scores_.begin();
  for (; scoreIt != lastUniqueIt; ++scoreIt) {
    if (scoreIt->isTarget()) {
      if (scoreIt->q <= fdr) {
        data.vals[ix2] = scoreIt->pPSM->features;
        data.Y[ix2] = 1;
        // data.C[ix2] = cpos;
        ix2++;
        ++p;
      }
    }
  }
  data.positives = p;
  data.m = static_cast<int>(ix2);
}

/**
 * OpenMS: what generateNegativeTrainingSet(data, 1.0) and then
 * generatePositiveTrainingSet(data, fdr, 1.0, trainBestPositive) give on a copy of this set,
 * without the copy and without sorting the PSMs; returns false (and changes nothing) where that
 * is not certain. The negatives are the decoys in their order. With trainBestPositive, the copy
 * is sorted by spectrum and label (OrderScanLabel), keeps the first PSM of each group
 * (std::unique) and sorts the kept PSMs by score again; the positives are the kept targets with
 * q <= fdr, in that order:
 *  - which PSM of a group comes first depends on how the sort orders the PSMs of the group, so
 *    the (spectrum, label) keys are sorted from the order of the PSMs here, like std::sort
 *    (sortLikeStdSort, the sort generatePositiveTrainingSet uses): the comparisons of the keys
 *    have the outcomes of OrderScanLabel's, so the sort makes the same moves;
 *  - if each positive scores strictly higher (operator>) than the next positive here, no two
 *    positives are equivalent, so every sort by score leaves them in this order. Otherwise, and
 *    with a NaN score or mass, false.
 */
bool Scores::generateTrainingSetsOfSortedSet(AlgIn& data, const double fdr,
                                             const bool trainBestPositive) const {
  // kept[i]: whether the PSM scores_[i] is kept (the first of its group)
  std::vector<char> kept(scores_.size(), 1);
  if (trainBestPositive && !scores_.empty()) {
    if (!comparableScores(scores_.begin(), scores_.end())) return false;
    struct ScanLabelKey {
      unsigned int specFileNr, scan;
      LabelType label;
      std::size_t index;
    };
    std::vector<ScanLabelKey> keys;
    keys.reserve(scores_.size());
    for (std::size_t i = 0; i < scores_.size(); ++i) {
      keys.push_back({scores_[i].pPSM->specFileNr, scores_[i].pPSM->scan, scores_[i].label, i});
    }
    // the comparison of OrderScanLabel, on the keys
    sortLikeStdSort(keys, [](const ScanLabelKey& x, const ScanLabelKey& y) {
      if (x.specFileNr != y.specFileNr) return x.specFileNr < y.specFileNr;
      if (x.scan != y.scan) return x.scan < y.scan;
      return x.label > y.label;
    });
    // std::unique with UniqueScanLabel: keeps the first of each run of equal keys
    for (std::size_t k = 1; k < keys.size(); ++k) {
      if (keys[k].specFileNr == keys[k - 1].specFileNr && keys[k].scan == keys[k - 1].scan &&
          keys[k].label == keys[k - 1].label) {
        kept[keys[k].index] = 0;
      }
    }
  }
  std::vector<std::size_t> positives;
  for (std::size_t i = 0; i < scores_.size(); ++i) {
    const ScoreHolder& sh = scores_[i];
    if (kept[i] && sh.isTarget() && sh.q <= fdr) {
      if (trainBestPositive && !positives.empty() && !(scores_[positives.back()] > sh)) {
        return false;
      }
      positives.push_back(i);
    }
  }

  std::size_t ix2 = 0;
  for (const ScoreHolder& sh : scores_) {
    if (sh.isDecoy()) {
      data.vals[ix2] = sh.pPSM->features;
      data.Y[ix2] = -1;
      ix2++;
    }
  }
  data.negatives = static_cast<int>(ix2);
  for (std::size_t i : positives) {
    data.vals[ix2] = scores_[i].pPSM->features;
    data.Y[ix2] = 1;
    ix2++;
  }
  data.positives = static_cast<int>(positives.size());
  data.m = static_cast<int>(ix2);
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
    // the keys after the sort of each feature (until its count), and the labels
    std::vector<std::vector<ChainKey> > keysOfFeature(numFeatures);
    std::vector<char> isTarget(n);
#pragma omp taskgroup
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
        // the feature values by feature, and the labels (blocks of PSMs as tasks)
        std::vector<double> columns(static_cast<std::size_t>(numFeatures) * n);
        const std::ptrdiff_t numScores = static_cast<std::ptrdiff_t>(n);
#pragma omp taskloop grainsize(kMinScoreTaskSize) shared(columns, isTarget)
        for (std::ptrdiff_t i = 0; i < numScores; ++i) {
          const ScoreHolder& sh = scores_[static_cast<std::size_t>(i)];
          for (unsigned int featNo = 0; featNo < numFeatures; ++featNo) {
            columns[featNo * n + static_cast<std::size_t>(i)] = sh.pPSM->features[featNo];
          }
          isTarget[static_cast<std::size_t>(i)] = sh.isTarget();
        }
        std::vector<ChainKey> keys;
        keys.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
          keys.push_back({0.0, tieRank[i], static_cast<unsigned int>(i)});
        }
        for (unsigned int featNo = 0; featNo < numFeatures; featNo++) {
          const double* column = &columns[featNo * n];
          for (ChainKey& key : keys) key.score = column[key.index];
          // operator< of ScoreHolder
          sortLikeStdSort(keys, [](const ChainKey& x, const ChainKey& y) {
            return x.score < y.score || (x.score == y.score && x.tieRank < y.tieRank);
          });
          keysOfFeature[featNo] = keys;
#pragma omp task firstprivate(featNo) shared(keysOfFeature, isTarget, positives, errors, initialSelectionFdr, skipDecoysPlusOne)
          try {
            std::vector<pair<double, bool> > pairs;
            pairs.reserve(keysOfFeature[featNo].size());
            for (const ChainKey& key : keysOfFeature[featNo]) {
              pairs.emplace_back(key.score, isTarget[key.index] != 0);
            }
            std::vector<ChainKey>().swap(keysOfFeature[featNo]);
            positives[2 * featNo] = countPositivesOfPairs(pairs, initialSelectionFdr, skipDecoysPlusOne);
            reverse(pairs.begin(), pairs.end());
            positives[2 * featNo + 1] = countPositivesOfPairs(pairs, initialSelectionFdr, skipDecoysPlusOne);
          } catch (...) {
            errors[featNo] = std::current_exception();
          }
          reverse(keys.begin(), keys.end());
        }
        permuteScores(scores_, keys);
        for (std::size_t i = 0; i < n; ++i) scores_[i].score = keys[i].score;
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
