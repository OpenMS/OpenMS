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
#ifndef SANITYCHECK_H_
#define SANITYCHECK_H_

#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "Scores.h"

namespace OpenMS { namespace Internal { namespace Percolator {


class Normalizer;

/*
 * SanityCheck sets the initial direction of the SVM, either from an input
 * weight file, a default direction given by the user or in the input file,
 * or the single most discriminative feature, and checks that the trained
 * direction is not worse than the initial one: by the targets with q < test
 * FDR summed over the test bins, each with its own q values, or by those of
 * the merged out-of-fold list, see DirectionCheck and Caller's
 * --direction-check. Where the trained direction is worse, the folds fall
 * back (resetDirection) to the default direction chosen on each fold's
 * training set or, with a weights file, per fold to its weights or that
 * default direction, whichever separates more targets of the training set
 * (see chooseResetVectors), unless --override keeps the trained vectors.
 */
class SanityCheck {
 public:
  // OpenMS: validateDirection is virtual, so that OpenMS' wrapper can check each fold on its own training set
  // (out-of-fold capped training)
  virtual ~SanityCheck() = default;
  // How validateDirection counts the targets with q < test FDR. With one
  // bin, as with RESET, all three are the same.
  enum class DirectionCheck {
    BINS,           // within each test bin, with the bin's own q values
    SECOND_CHANCE,  // as BINS; where that count would reset, the trained
                    // SVMs are kept if the merged out-of-fold list,
                    // rescaled by the training decoys, can be used (see
                    // MergePlan::lists) and holds targets with
                    // q < test FDR and no fewer than the reset would
                    // report: the folds' default vectors, or under
                    // --init-weights the vectors of chooseResetVectors,
                    // merged with the scheme of the plan; their list merged
                    // with the plan's secondChanceScheme has to be usable
                    // too
    POOLED          // in the merged out-of-fold list only, merged with the
                    // scheme of the plan
  };
  // What validateDirection decided about the trained SVMs
  enum class DirectionOutcome {
    KEPT,                    // the count of the check passed
    KEPT_BY_SECOND_CHANCE,   // the per-bin count failed, the second chance
                             // kept them: the bins are merged with the
                             // plan's secondChanceScheme
    FAILED                   // not accepted: the folds were reset, or under
                             // --override keep their trained vectors
  };

  void readWeights(std::istream& weightStream, std::vector<double>& w);
  void checkAndSetDefaultDir();
  int getInitDirection(std::vector<Scores>& testset,
                       std::vector<Scores>& trainset, const Normalizer* pNorm,
                       std::vector<std::vector<double> >& w, double test_fdr,
                       double initial_train_fdr);
  virtual DirectionOutcome validateDirection(std::vector<Scores>& testset,
                                     std::vector<Scores>& trainset,
                                     std::vector<std::vector<double> >& w);
  void resetDirection(std::vector<Scores>& trainset,
                      std::vector<std::vector<double> >& w);

  void setInitWeightFN(const std::string& fn) { initWeightFN_ = fn; }
  void setInitDefaultDirName(const std::string& dirName) {
    initDefaultDirName_ = dirName;
  }
  void setOverrule(bool orl) { overRule_ = orl; }
  void setDirectionCheck(DirectionCheck check) { directionCheck_ = check; }
  DirectionCheck directionCheck() const { return directionCheck_; }
  // How the test bins are merged, which the merged counts repeat, see
  // MergePlan, Scores::merge and Scores::countMergedTargets: the pooled
  // check counts with the plan's scheme; the second chance counts the
  // trained SVMs with train-decoys and the fallback with the plan's scheme.
  void setMergePlan(const MergePlan& plan) { plan_ = plan; }
  const MergePlan& mergePlan() const { return plan_; }
  void addDefaultWeights(const std::vector<double>& defaultWeights) {
    defaultWeights_ = defaultWeights;
  }
  const std::vector<double>& getDefaultWeights() const {
    return defaultWeights_;
  }
  bool concatenatedSearch() const { return concatenatedSearch_; }
  void setConcatenatedSearch(const bool s) { concatenatedSearch_ = s; }

 protected:
  void getDefaultDirection(const std::vector<Scores>& trainset,
                           std::vector<std::vector<double> >& w) const;
  // The default direction of one fold, as getDefaultDirection chooses it:
  // the feature of --default-direction, the default weights of the input or
  // the single feature that separates most targets of the fold's training
  // set. quiet leaves out the messages of the latter and its error where no
  // feature gives a target; w is then all 0 and the result false.
  bool getDefaultDirection(const Scores& trainset, std::size_t set,
                           std::vector<double>& w, bool quiet) const;
  // Whether a reset chooses per fold between the vector of --init-weights
  // and the fold's default direction. Not with --static
  // (skipNormalizeScores), whose given weights are the model: a reset
  // restores them.
  bool choosesResetVectors() const {
    return !fileWeights_.empty() && !plan_.skipNormalizeScores;
  }
  // The vectors of a reset under --init-weights (choosesResetVectors): per
  // fold, the vector of the weights file or the fold's default direction,
  // whichever gives more targets with q < the initial training FDR in the
  // fold's training set, counted as getInitDirection counts a feature; the
  // weights file on ties, and where no default direction gives a target.
  // Sets resetWeights_ and resetCounts_; leaves the bins and the random
  // number generator as they are.
  void chooseResetVectors(const std::vector<Scores>& trainset);
  // The vectors that resetDirection gives the folds: under --init-weights
  // those of chooseResetVectors, which validateDirection calls before the
  // second chance; otherwise the initial ones, which getInitDirection keeps
  // for the second chance: each fold's default vector, which a reset
  // chooses again on the same training set, or with --static the given
  // weights in every fold.
  const std::vector<std::vector<double> >& resetVectors() const {
    return choosesResetVectors() ? resetWeights_ : initWeights_;
  }
  // whether the targets are counted in the merged list of these bins
  bool usePooledCount(std::size_t numBins) const {
    return directionCheck_ == DirectionCheck::POOLED && numBins > 1u;
  }
  // whether a failed per-bin count of these bins gets the second chance
  bool useSecondChance(std::size_t numBins) const {
    return directionCheck_ == DirectionCheck::SECOND_CHANCE && numBins > 1u;
  }
  // whether the weight of the --default-direction feature has the wrong
  // sign in one of the folds' vectors w
  bool wrongSignOfDefaultDirection(
      const std::vector<std::vector<double> >& w) const;
  // The second chance: the targets with q < test FDR in the merged
  // out-of-fold list with train-decoys for the trained weights w, against
  // those of the fallback as the run would report it (resetVectors, merged
  // with the scheme of the plan), and whether they keep w, which also
  // needs the right sign of the --default-direction feature. A list that
  // cannot be used, also by the lists of the plan, counts -1 for w
  // and 0 for the fallback. w is not kept either where its list merged with
  // the plan's secondChanceScheme, which the run then uses,
  // cannot be used. Leaves the bins, the weights and the random number
  // generator as they are. Prints its decision.
  bool secondChance(const std::vector<Scores>& testset,
                    const std::vector<Scores>& trainset,
                    const std::vector<std::vector<double> >& w,
                    int& trainedCount) const;

  int initPositives_ = 0;
  double test_fdr_ = 0.01, initial_train_fdr_ = 0.01;
  bool overRule_ = false;
  std::string initWeightFN_;
  std::string initDefaultDirName_;
  // Default Direction, 0=do not use, positive integer = feature number,
  // negative integer = lower score better
  int initDefaultDir_ = 0;
  // initial direction given in the input file
  std::vector<double> defaultWeights_;
  // normalized initial weights read from the weights file (-W), if any
  std::vector<double> fileWeights_;
  // input from concatenated search, i.e. not separate target and decoy searches
  bool concatenatedSearch_ = true;
  // the direction check and the merge that its merged counts repeat
  DirectionCheck directionCheck_ = DirectionCheck::BINS;
  MergePlan plan_;
  // the second chance: the initial direction of each fold, which is also
  // the vector a reset gives it without --init-weights
  std::vector<std::vector<double> > initWeights_;
  // the vectors of a reset under --init-weights and, per fold, the targets
  // of the weights file and of the default direction (chooseResetVectors)
  std::vector<std::vector<double> > resetWeights_;
  std::vector<std::pair<int, int> > resetCounts_;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*SANITYCHECK_H_*/
