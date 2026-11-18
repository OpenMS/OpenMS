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
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "DataSet.h"
#include "Scores.h"
#include "Normalizer.h"
#include "SanityCheck.h"
#include "Globals.h"
#include "PosteriorEstimator.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;


namespace {
// why a merged count finds no usable list, see Scores::countMergedTargets
const char* const kUnusableListReasons =
    "the scores of a test bin cannot be rescaled, or targets and decoys "
    "separate too well in a test bin (test, --static), in the merged list or "
    "in its best PSM per spectrum (-Y) or per peptide";
}  // namespace

void SanityCheck::checkAndSetDefaultDir() {
  if (!initDefaultDir_ && initDefaultDirName_.size() > 0) {
    int sign = 1;
    if (initDefaultDirName_[0] == '-') {
      initDefaultDirName_.erase(0,1);
      sign = -1;
    }
    initDefaultDir_ = sign * DataSet::getFeatureNames().getFeatureNumber(initDefaultDirName_);
    if (initDefaultDir_ == 0) {
      ostringstream temp;
      temp << "ERROR: Initial direction feature name \"" << initDefaultDirName_ << "\" not found" << std::endl;
      throw MyException(temp.str());
    }
  }
}

int SanityCheck::getInitDirection(vector<Scores>& testset,
                                  vector<Scores>& trainset,
                                  const Normalizer* pNorm,
                                  vector<vector<double> > & w,
                                  double test_fdr,
                                  double initial_train_fdr) {
  test_fdr_= test_fdr;
  initial_train_fdr_ = initial_train_fdr;
  if (initWeightFN_.size() > 0) {
    vector<double> ww(FeatureNames::getNumFeatures() + 1);
    ifstream weightStream(initWeightFN_.data(), ios::in);
    if (weightStream.is_open()) {
      readWeights(weightStream, ww);
      weightStream.close();
      assert(pNorm);
      pNorm->normalizeweight(ww, w[0]);
      for (size_t set = 1; set < w.size(); ++set) {
        copy(w[0].begin(), w[0].end(), w[set].begin());
      }
      fileWeights_ = w[0];
    } else {
      std::cerr << "WARNING: Could not find weights input file " << initWeightFN_
                << ". Using default weights instead." << std::endl;
      getDefaultDirection(trainset, w);
    }
  } else {
    getDefaultDirection(trainset, w);
  }
  initPositives_ = 0;
  for (size_t set = 0; set < w.size(); ++set) {
    initPositives_ += testset[set].calcScoresAndQvals(w[set], test_fdr);
  }
  // The pooled check compares the trained direction with the initial one
  // counted in the same way: per fold, merged and pooled. The bins are
  // still scored above, so that both checks leave them in the same state.
  if (usePooledCount(w.size())) {
    initPositives_ = std::max(
        0, Scores::countMergedTargets(testset, trainset, w, plan_.scheme,
                                      plan_.fdr, plan_.skipNormalizeScores,
                                      test_fdr, plan_.lists));
  }
  // The second chance counts the vectors of a reset only when it runs.
  // Without --init-weights, they are the initial ones, the default vectors,
  // which a reset chooses again on the same training sets; with --static,
  // the vector of --init-weights, which a reset restores. Training changes
  // w, so they are kept for it. Otherwise, under --init-weights, the reset
  // chooses per fold when the check fails (chooseResetVectors).
  initWeights_.clear();
  if (useSecondChance(w.size())) {
    initWeights_ = w;
  }
  return initPositives_;
}

void SanityCheck::getDefaultDirection(const vector<Scores>& trainset,
                                      vector<vector<double> >& w) const {
  for (size_t set = 0; set < w.size(); ++set) {
    getDefaultDirection(trainset[set], set, w[set], false);
  }
}

bool SanityCheck::getDefaultDirection(const Scores& trainset, size_t set,
                                      vector<double>& w, bool quiet) const {
  //If I have not been given a initial direction
  if (!initDefaultDir_) {
    //if the default weights from the input file are not present
    if (defaultWeights_.size() == 0) {
      // Set init direction to be the most discriminative direction
      if (quiet) {
        return trainset.findInitDirection(initial_train_fdr_, w) > 0;
      }
      if (VERB > 1) std::cerr << "Split " << set + 1 << ":\t";
      return trainset.getInitDirection(initial_train_fdr_, w) > 0;
    }
    // I want to assign the default vector that is present in the input file
    for (size_t ix = 0; ix < w.size(); ix++) {
      w[ix] = 0;
      if (ix < defaultWeights_.size()) {
        w[ix] = defaultWeights_[ix];
      }
    }
    return true;
  }
  for (size_t ix = 0; ix < w.size(); ix++) {
    w[ix] = 0;
  }
  w[static_cast<std::size_t>(abs(initDefaultDir_) - 1)] =
      (initDefaultDir_ < 0 ? -1 : 1);
  return true;
}

void SanityCheck::chooseResetVectors(const vector<Scores>& trainset) {
  resetWeights_.assign(trainset.size(), fileWeights_);
  resetCounts_.assign(trainset.size(), std::make_pair(0, 0));
  for (size_t set = 0; set < trainset.size(); ++set) {
    Scores::FeatureRows rows;
    trainset[set].getFeatureRows(rows);
    Scores::CountBuffers buffers;
    // Both vectors are counted as getInitDirection counts the features:
    // the targets with q < the initial training FDR in the training set,
    // without the decoy + 1 of the FDR estimate. For the single feature of
    // the default direction this is the count that chose it.
    const bool skipDecoysPlusOne = true;
    const int fileCount = trainset[set].countTargetsBelowFdr(
        rows, fileWeights_, initial_train_fdr_, skipDecoysPlusOne, buffers);
    vector<double> defaultWeights(fileWeights_.size(), 0.0);
    int defaultCount = 0;
    if (getDefaultDirection(trainset[set], set, defaultWeights, true)) {
      defaultCount = trainset[set].countTargetsBelowFdr(
          rows, defaultWeights, initial_train_fdr_, skipDecoysPlusOne,
          buffers);
    }
    resetCounts_[set] = std::make_pair(fileCount, defaultCount);
    if (defaultCount > fileCount) {
      resetWeights_[set] = defaultWeights;
    }
  }
}

bool SanityCheck::wrongSignOfDefaultDirection(
    const vector<vector<double> >& w) const {
  if (!initDefaultDir_) {
    return false;
  }
  for (size_t set = 0; set < w.size(); ++set) {
    if (w[set][static_cast<std::size_t>(abs(initDefaultDir_) - 1)] *
            initDefaultDir_ <= 0) {
      return true;
    }
  }
  return false;
}

bool SanityCheck::secondChance(const vector<Scores>& testset,
                               const vector<Scores>& trainset,
                               const vector<vector<double> >& w,
                               int& trainedCount) const {
  // the anchors of the training sets exist for test bins of any size
  trainedCount = Scores::countMergedTargets(
      testset, trainset, w, MergeNormalization::TRAIN_DECOYS, plan_.fdr,
      plan_.skipNormalizeScores, test_fdr_, plan_.lists);
  // The fallback that the second chance replaces, counted as the run would
  // report it: the vectors of a reset (each fold's default vector, or under
  // --init-weights per fold its vector or the default vector, whichever
  // separates more of the training set), merged with the run's scheme
  // (test up to three folds) and counted at the test FDR. A fallback whose
  // list cannot be used would end the run and counts 0.
  int resetCount = 0;
  const vector<vector<double> >& resetWeights = resetVectors();
  if (resetWeights.size() == w.size()) {
    resetCount = std::max(
        0, Scores::countMergedTargets(testset, trainset, resetWeights,
                                      plan_.scheme, plan_.fdr,
                                      plan_.skipNormalizeScores, test_fdr_,
                                      plan_.lists));
  }
  const bool countsKeep = trainedCount > 0 && trainedCount >= resetCount;
  // A run that keeps the trained SVMs is merged with the plan's
  // secondChanceScheme. Another scheme than train-decoys, given by
  // --merge-normalization, can end that run where train-decoys does not:
  // test checks the separation of each bin, and the scores of test, train
  // and none are merged on another scale. Its list therefore has to be
  // usable as well; the count of train-decoys still decides. --static
  // merges the raw scores, which the count above has already checked.
  const bool otherKeptScheme =
      !plan_.skipNormalizeScores &&
      plan_.secondChanceScheme != MergeNormalization::TRAIN_DECOYS;
  const bool unusableAsKept =
      countsKeep && otherKeptScheme &&
      Scores::countMergedTargets(testset, trainset, w,
                                 plan_.secondChanceScheme, plan_.fdr,
                                 plan_.skipNormalizeScores, test_fdr_,
                                 plan_.lists) < 0;
  // the sign check of --default-direction applies as after a passed count
  const bool wrongSign =
      countsKeep && !unusableAsKept && wrongSignOfDefaultDirection(w);
  const bool keep = countsKeep && !unusableAsKept && !wrongSign;
  cerr << "Second chance on the merged out-of-fold list"
       << (plan_.skipNormalizeScores ? "" : ", rescaled by the training decoys")
       << ": ";
  if (trainedCount < 0) {
    cerr << "the trained SVMs give no usable list (" << kUnusableListReasons
         << ")";
  } else {
    cerr << trainedCount << " targets with q<" << test_fdr_
         << " for the trained SVMs, " << resetCount << " for the "
         << (fileWeights_.empty()
                 ? "default vectors"
                 : (choosesResetVectors() ? "--init-weights or default vectors"
                                          : "--init-weights vector"))
         << " of a reset";
    if (!plan_.skipNormalizeScores) {
      cerr << " merged with " << mergeNormalizationName(plan_.scheme);
    }
    if (unusableAsKept) {
      cerr << "; merged with "
           << mergeNormalizationName(plan_.secondChanceScheme)
           << ", as a run that keeps them would be, the trained SVMs give no "
           << "usable list (" << kUnusableListReasons << ")";
    }
  }
  if (overRule_) {
    cerr << "; --override keeps the trained SVMs in any case." << endl;
  } else if (keep) {
    cerr << "; the trained SVMs are kept";
    if (!plan_.skipNormalizeScores) {
      cerr << ", merged with "
           << mergeNormalizationName(plan_.secondChanceScheme);
    }
    cerr << "." << endl;
  } else if (wrongSign) {
    cerr << "; the trained SVMs are not kept, as the weight of the main "
         << "scoring direction (--default-direction) has the wrong sign."
         << endl;
  } else {
    cerr << "; the trained SVMs are not kept." << endl;
  }
  return keep;
}

SanityCheck::DirectionOutcome SanityCheck::validateDirection(
    vector<Scores>& testset,
    vector<Scores>& trainset,
    vector<vector<double> >& w) {
  resetWeights_.clear();
  resetCounts_.clear();
  const bool pooled = usePooledCount(w.size());
  int overFDR = 0;
  if (pooled) {
    // on copies: the test bins are only scored once the direction is decided
    overFDR = Scores::countMergedTargets(testset, trainset, w, plan_.scheme,
                                         plan_.fdr, plan_.skipNormalizeScores,
                                         test_fdr_, plan_.lists);
    if (overFDR < 0) {
      cerr << "The trained weights give no usable merged out-of-fold list: "
           << kUnusableListReasons << endl;
      resetDirection(trainset, w);
      return DirectionOutcome::FAILED;
    }
    if (VERB > 1) {
      cerr << "Found " << overFDR << " targets with q<" << test_fdr_
           << " in the merged out-of-fold list, " << initPositives_
           << " in the initial direction" << endl;
    }
  } else {
    for (size_t set = 0; set < w.size(); ++set) {
      overFDR += testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  bool keptBySecondChance = false;
  if (overFDR <= 0 || initPositives_ > overFDR) {
    if (overFDR <= 0) {
      cerr << "No targets found with q<" << test_fdr_ << endl;
    } else {
      cerr << "Less identifications (" << overFDR << " vs " << initPositives_
          << ") after percolator processing than before processing" << endl;
    }
    // Where the per-bin count would reset, the second chance judges the
    // trained SVMs on the merged out-of-fold list, which small test bins
    // do not limit. It decides nothing under --override.
    const bool secondChanceTried = useSecondChance(w.size());
    // under --init-weights the second chance compares with the vectors that
    // a reset would choose
    if (choosesResetVectors() && (secondChanceTried || !overRule_)) {
      chooseResetVectors(trainset);
    }
    int pooledCount = 0;
    keptBySecondChance = secondChanceTried &&
                         secondChance(testset, trainset, w, pooledCount) &&
                         !overRule_;
    if (!keptBySecondChance) {
      resetDirection(trainset, w);
      if (overFDR <= 0 && !overRule_) {
        // not for an unusable merged list (-1), whose reason is printed
        if (pooled || (secondChanceTried && pooledCount == 0)) {
          cerr << "Hint: q<" << test_fdr_ << " needs about "
               << 1.0 / test_fdr_ << " targets above the best-scoring decoy "
               << "of the merged list; a larger --testFDR is an alternative."
               << endl;
        } else if (w.size() > 1u && !secondChanceTried) {
          cerr << "Hint: each of the " << w.size() << " test folds may be too "
               << "small to reach q<" << test_fdr_ << ". --num-folds 1 trains "
               << "and evaluates on the same PSMs," << endl
               << "so its q values and PEPs are optimistic, and wrong for "
               << "inputs without signal; a larger --testFDR is an "
               << "alternative." << endl;
        }
      }
      return DirectionOutcome::FAILED;
    }
  }
  if (wrongSignOfDefaultDirection(w)) {
    if (VERB > 1) {
      cerr << "Warning, wrong sign of the weight for main scoring direction"
           << endl;
    }
    resetDirection(trainset, w);
    return DirectionOutcome::FAILED;
  }
  if (pooled) {
    // the direction is kept: score the test bins with their folds' weights
    // for the merge, as the count of the bins check does
    for (size_t set = 0; set < w.size(); ++set) {
      testset[set].calcScoresAndQvals(w[set], test_fdr_);
    }
  }
  // the per-bin count has already scored the test bins with w
  return keptBySecondChance ? DirectionOutcome::KEPT_BY_SECOND_CHANCE
                            : DirectionOutcome::KEPT;
}

/**
 * Reads the raw weights of the first split from a file written with -w:
 * comment lines, then per split a line of feature names followed by m0, a
 * line of normalized weights and a line of raw weights. The weights are
 * used by position, so the names, if present, must be the feature names of
 * the input, and the raw weights line must hold one finite number per
 * feature plus m0, with at least one feature weight that is not 0.
 */
void SanityCheck::readWeights(istream& weightStream, vector<double>& w) {
  const size_t numWeights = FeatureNames::getNumFeatures() + 1;
  std::string line, namesLine;
  unsigned int lineNr = 0;
  bool foundWeights = false;
  // Find first line with weights, but skip this and get second line containing raw features
  while (std::getline(weightStream, line)) {
    ++lineNr;
    if (line[0] == '-' || (line[0] >= '0' && line[0] <= '9')) {
      foundWeights = true;
      break;
    }
    if (!line.empty() && line[0] != '#') namesLine = line;
  }
  std::string rawLine;
  if (!foundWeights || !std::getline(weightStream, rawLine)) {
    throw MyException("ERROR: The weights file " + initWeightFN_ +
                      " has no line of normalized weights followed by a "
                      "line of raw weights.\n");
  }
  if (!namesLine.empty()) {
    std::vector<std::string> names;
    std::istringstream nameStream(namesLine);
    for (std::string name; std::getline(nameStream, name, '\t');) {
      names.push_back(name);
    }
    while (!names.empty() && names.back().find_first_not_of(" \r") ==
                                 std::string::npos) {
      names.pop_back();  // a trailing tab or carriage return
    }
    if (!names.empty()) {
      size_t ix = names.back().find_last_not_of(" \r");
      names.back().erase(ix + 1);
    }
    // the names as -w writes them for this input
    std::vector<std::string> expected;
    std::istringstream expectedStream(
        DataSet::getFeatureNames().getFeatureNames() + "\tm0");
    for (std::string name; std::getline(expectedStream, name, '\t');) {
      expected.push_back(name);
    }
    for (size_t ix = 0; ix < std::max(names.size(), expected.size()); ++ix) {
      if (ix < names.size() && ix < expected.size() &&
          names[ix] == expected[ix]) {
        continue;
      }
      ostringstream oss;
      oss << "ERROR: The weights file " << initWeightFN_ << " has "
          << names.size() << " columns, but the input has "
          << expected.size() << " (its features and m0). Column " << ix + 1
          << " is " << (ix < names.size() ? "\"" + names[ix] + "\"" :
                                            std::string("missing"))
          << " in the weights file and "
          << (ix < expected.size() ? "\"" + expected[ix] + "\"" :
                                     std::string("missing"))
          << " in the input. Weights are used by the position of the "
          << "feature, so the file must come from an input with the same "
          << "features in the same order." << std::endl;
      throw MyException(oss.str());
    }
  }
  std::istringstream rawStream(rawLine);
  std::vector<double> raw;
  for (double value; rawStream >> value;) raw.push_back(value);
  bool allFinite = std::all_of(raw.begin(), raw.end(),
                               [](double x) { return std::isfinite(x); });
  if (!rawStream.eof() || raw.size() != numWeights || !allFinite) {
    ostringstream oss;
    oss << "ERROR: Line " << lineNr + 1 << " of the weights file "
        << initWeightFN_ << " should hold " << numWeights
        << " finite raw weights, one per feature and m0, but ";
    if (!rawStream.eof() || !allFinite) {
      oss << "holds a value that is not a finite number";
    } else {
      oss << "holds " << raw.size();
    }
    oss << "." << std::endl;
    throw MyException(oss.str());
  }
  // Feature weights that are all 0 score every PSM alike, whatever m0 is:
  // they cannot order the PSMs, and the run would end later, where the
  // scores of the cross-validation bins cannot be rescaled.
  if (numWeights > 1u && std::all_of(raw.begin(), raw.end() - 1,
                                     [](double x) { return x == 0.0; })) {
    ostringstream oss;
    oss << "ERROR: Line " << lineNr + 1 << " of the weights file "
        << initWeightFN_ << " gives every feature the weight 0, so these "
        << "weights cannot order the PSMs." << std::endl;
    throw MyException(oss.str());
  }
  std::copy(raw.begin(), raw.end(), w.begin());
  cerr << "Read weights from file" << endl;
  for (unsigned int ix = 0; ix < FeatureNames::getNumFeatures() + 1; ix++) {
    cerr << w[ix] << "\t";
  }
  cerr << endl;
}

void SanityCheck::resetDirection(vector<Scores>& trainset,
                                 vector<vector<double> >& w) {
  if (!overRule_ && choosesResetVectors()) {
    // per fold the vector given with -W or the default that it replaced,
    // chosen on the fold's training set
    if (resetWeights_.size() != w.size()) {
      chooseResetVectors(trainset);
    }
    cerr << "Resetting score vector, in each fold to the weights file or to "
         << "the default vector, whichever gives more targets with q<"
         << initial_train_fdr_ << " in the fold's training set (the weights "
         << "file on ties; targets of the weights file vs the default "
         << "vector):";
    for (size_t set = 0; set < w.size(); ++set) {
      const std::pair<int, int>& counts = resetCounts_[set];
      cerr << (set > 0u ? "," : "") << " fold " << set + 1 << " "
           << (counts.second > counts.first ? "default vector"
                                            : "weights file")
           << " (" << counts.first << " vs " << counts.second << ")";
      w[set] = resetWeights_[set];
    }
    cerr << ". Use --override flag to prevent this." << endl;
    resetWeights_.clear();
    resetCounts_.clear();
  } else if (!overRule_ && !fileWeights_.empty()) {
    // --static: the vector given with -W is the model
    cerr << "Resetting score vector, using the initial weights from the "
         << "weights file. Use --override flag to prevent this." << endl;
    for (size_t set = 0; set < w.size(); ++set) {
      w[set] = fileWeights_;
    }
  } else if (!overRule_) {
    cerr << "Resetting score vector, using default vector. Use --override flag to prevent this." << endl;
    getDefaultDirection(trainset, w);
  }
}

}}}  // namespace OpenMS::Internal::Percolator
