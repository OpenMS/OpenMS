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
#ifndef SCORES_H_
#define SCORES_H_

#ifndef WIN32
#include <stdint.h>
#endif

#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "LabelType.h"
#include "Normalizer.h"
#include "PSMDescription.h"
#include "PseudoRandom.h"
#include "ScoreHolder.h"

namespace OpenMS { namespace Internal { namespace Percolator {


class SetHandler;
class AlgIn;

/*
 * How Scores::merge puts the SVM scores of the cross-validation bins on a
 * common scale before merging them, see Caller's --merge-normalization.
 */
enum class MergeNormalization {
  TEST,          // the score of the last PSM with q < fdr and the median
                 // decoy of each test bin become 0 and -1
  TRAIN,         // as TEST, with both scores taken from the fold's training
                 // set scored with the fold's weights
  TRAIN_DECOYS,  // (score - mean) / sd of the decoys of the fold's training
                 // set scored with the fold's weights
  NONE           // the raw SVM scores
};

// The name of the scheme as --merge-normalization takes it: "test",
// "train", "train-decoys" or "none".
const char* mergeNormalizationName(MergeNormalization scheme);

/*
 * The linear map score -> (score - offset) / scale that puts the SVM scores
 * of a cross-validation bin on the common scale. A scale <= 0 means that the
 * anchors of the map do not exist, e.g. for too few decoys.
 */
struct ScoreScale {
  double offset = 0.0;
  double scale = 0.0;
};

/*
 * The lists that the run makes from the merged list, in this order, before
 * it reports them: the best-scoring PSM of each spectrum (-Y,
 * Scores::weedOutRedundantTDC) and, from the PSMs left, the best-scoring PSM
 * of each peptide (the peptide level, Scores::weedOutRedundant). Each of
 * them ends the run, as the merged list does, when its targets and decoys
 * separate too well (postMergeStep).
 */
struct ReportedLists {
  bool bestPsmPerSpectrum = false;
  bool bestPsmPerPeptide = false;
};

/*
 * How the cross-validation bins of a run are merged, for the merge itself
 * (CrossValidation::postIterationProcessing, Scores::merge) and for the
 * direction check, whose merged counts repeat it (SanityCheck,
 * Scores::countMergedTargets). Caller::run builds the one value.
 */
struct MergePlan {
  // --merge-normalization, or the default of the number of folds
  // (CrossValidation::defaultMergeNormalization): the scheme of the merge,
  // of the pooled count and of the fallback that the second chance counts
  MergeNormalization scheme = MergeNormalization::TEST;
  // The scheme of the merge where the second chance of the direction check
  // kept the trained SVMs: train-decoys, as it counts them, unless a scheme
  // was given, whose list the second chance then checks as well
  MergeNormalization secondChanceScheme = MergeNormalization::TRAIN_DECOYS;
  // --trainFDR as CrossValidation resolved it, which sets it: the anchors of
  // test and train are taken at q = fdr
  double fdr = 0.01;
  // --static: the raw scores are merged, whatever the scheme
  bool skipNormalizeScores = false;
  // the lists the run makes from the merged one (-Y, the peptide level),
  // which the merged counts also find unusable when they separate too well
  ReportedLists lists;

  // The scheme the bins are merged with. --static rescales nothing and
  // ignores the scheme, after a second chance too: its bins are merged as
  // with test, which also makes -N rescale the scores of all PSMs.
  MergeNormalization usedScheme(bool keptBySecondChance) const {
    if (skipNormalizeScores) {
      return MergeNormalization::TEST;
    }
    return keptBySecondChance ? secondChanceScheme : scheme;
  }
};

/*
 * Scores is a container of ScoreHolders that allows you to do a sorted merge
 * of vectors of ScoreHolder.
 *
 * Here are some useful abbreviations:
 * FDR - False Discovery Rate
 * Pi0 - prior probability of null hypothesis
 * TDC - Target Decoy Competition
 *
 */
class Scores {
 public:
  using iterator = std::vector<ScoreHolder>::iterator;
  using const_iterator = std::vector<ScoreHolder>::const_iterator;

  Scores(bool usePi0)
      : usePi0_(usePi0),
        pi0_(1.0),
        targetDecoySizeRatio_(1.0),
        nullTargetWinProb_(0.5),
        totalNumberOfDecoys_(0),
        totalNumberOfTargets_(0) {}

  // Replaces the PSMs of this set by those of the cross-validation bins sv,
  // scored with the weights all_w (after a failed direction check under -N,
  // a reset or --override, CrossValidation scores them all with all_w[0]
  // for TEST), and rescales each bin's scores
  // and weights as the scheme says: TEST takes the anchors from the bin's
  // own scores, TRAIN and TRAIN_DECOYS take binScales[bin] (see
  // calcScoreScale), NONE keeps the raw scores. The result is sorted by
  // score, and the bins are emptied.
  void merge(std::vector<Scores>& sv,
             double fdr,
             bool skipNormalizeScores,
             std::vector<std::vector<double> >& all_w,
             MergeNormalization scheme = MergeNormalization::TEST,
             const std::vector<ScoreScale>& binScales =
                 std::vector<ScoreScale>());
  // tooGoodSeparationIsError = false sets pi0 = 1 without a message when the
  // separation is too good; returns whether it is
  bool postMergeStep(bool tooGoodSeparationIsError = true);
  // The number of targets with q < testFdr in the list that merge would give
  // for the test bins scored with the weights w, one per bin: each bin is
  // rescaled as the scheme says, with the anchors of TRAIN and TRAIN_DECOYS
  // taken from trainBins[bin] scored with w[bin] (see calcScoreScale), and
  // the merged list gets its pi0 and q values as the final list does
  // (postMergeStep, calcQvals). The count is that of the merged PSMs; the
  // lists made from them are not counted. Works on copies, so neither the
  // bins nor the weights change, and leaves the global random number
  // generator in the state it found.
  // skipNormalizeScores (--static) merges the raw scores, which merge then
  // treats as with TEST.
  // Returns -1, a list that cannot be used, when the scores of one of
  // several bins cannot be rescaled, which would end the merge, or when
  // targets and decoys separate too well in a bin whose pi0 the merge
  // estimates (TEST, skipNormalizeScores), in the merged list or in one of
  // the lists the run makes from it (lists), which would end the run at
  // their postMergeStep.
  static int countMergedTargets(const std::vector<Scores>& testBins,
                                const std::vector<Scores>& trainBins,
                                const std::vector<std::vector<double> >& w,
                                MergeNormalization scheme,
                                double fdr,
                                bool skipNormalizeScores,
                                double testFdr,
                                ReportedLists lists = ReportedLists());

  // The rescaling of train and train-decoys for each bin, taken from the
  // training set of its fold scored with the fold's weights w[bin]
  // (calcScoreScale): what merge takes as binScales. Empty for the other
  // schemes and with skipNormalizeScores. onCopies leaves the training sets
  // as they are and prints nothing.
  static std::vector<ScoreScale> calcBinScales(
      std::vector<Scores>& trainBins,
      const std::vector<std::vector<double> >& w,
      MergeNormalization scheme,
      double fdr,
      bool skipNormalizeScores,
      bool onCopies = false);

  iterator begin() { return scores_.begin(); }
  iterator end() { return scores_.end(); }
  const_iterator begin() const { return scores_.begin(); }
  const_iterator end() const { return scores_.end(); }

  double calcScore(const double* features, const std::vector<double>& w) const;
  void scoreAndAddPSM(ScoreHolder& sh,
                      const std::vector<double>& rawWeights,
                      FeatureMemoryPool& featurePool);
  int calcScoresAndQvals(const vector<double>& w,
                         double fdr,
                         bool skipDecoysPlusOne = false);
  int calcScores(const vector<double>& w);
  int calcQvals(double fdr, bool skipDecoysPlusOne = false);
  // whether a sort by decreasing score would leave the PSMs as they are
  bool isSortedByScore() const;
  // The feature row and label of a PSM in one word: the address of its
  // features with the lowest bit set for a target (feature rows are
  // aligned to at least 2 bytes). See getFeatureRows.
  typedef std::uintptr_t FeatureRow;
  typedef std::vector<FeatureRow> FeatureRows;
  static inline FeatureRow makeFeatureRow(const double* features,
                                          bool isTarget) {
    return reinterpret_cast<std::uintptr_t>(features) |
           static_cast<std::uintptr_t>(isTarget);
  }
  static inline const double* rowFeatures(FeatureRow row) {
    return reinterpret_cast<const double*>(row & ~std::uintptr_t(1));
  }
  static inline bool rowIsTarget(FeatureRow row) { return (row & 1u) != 0; }
  void getFeatureRows(FeatureRows& rows) const;
  // Buffers of the counts below, kept by a thread for all its counts so
  // that they are allocated once
  struct CountBuffers {
    std::vector<double> values;       // one value per row
    std::vector<double> decoyValues;  // the values of the decoys
    // the PSMs that are sorted: those ahead of the cut-off decoy
    std::vector<std::pair<double, bool> > high, low;
  };
  int countTargetsBelowFdr(const FeatureRows& rows,
                           const std::vector<double>& w,
                           double fdr,
                           bool skipDecoysPlusOne,
                           CountBuffers& buffers) const;
  // The numbers of targets with q < fdr that calcQvals would count, with
  // pi0 = 1 and the given null target win probability, after sorting the
  // PSMs by decreasing value (highCount) and by increasing value
  // (lowCount); either may be null. buffers.values holds one value per
  // row of rows, which gives the labels. Only the PSMs ahead of the decoy
  // from which on no cut-off can reach q < fdr are sorted. Returns false
  // without counting if a value is NaN, or if the null target win
  // probability is not in [0, 1).
  static bool countTargetsBySelection(const FeatureRows& rows,
                                      CountBuffers& buffers,
                                      double fdr,
                                      bool skipDecoysPlusOne,
                                      double nullTargetWinProb,
                                      int* highCount,
                                      int* lowCount);
  // PEPs by P-IRLS (spline), or by I-spline/PAVA (pava) regression of either
  // target q values (default) or target-decoy labels (tdc) on the score
  // (interpol) or the rank
  void calcPep(const bool spline = false, const bool interpol = true,
               const bool pava = false, const bool tdc = false);
 
  void populateWithPSMs(SetHandler& setHandler);

  int getInitDirection(const double initialSelectionFdr,
                       std::vector<double>& direction) const;
  // The direction of getInitDirection, the single feature that gives the
  // most targets with q < initialSelectionFdr, without its message and
  // without its error where no feature gives a target: direction is then
  // all 0. Returns the number of those targets.
  int findInitDirection(const double initialSelectionFdr,
                        std::vector<double>& direction) const;
  // OpenMS: with fixedFold (a function that is set), the fold of each PSM is given instead of drawn; the caller keeps
  // the PSMs of a spectrum in one fold (out-of-fold capped training)
  void createXvalSetsBySpectrum(std::vector<Scores>& train,
                                std::vector<Scores>& test,
                                unsigned int xval_fold,
                                double decoyFractionTraining = 1.0,
                                unsigned int decoysPerTarget = 1u,
                                const std::function<unsigned int(const PSMDescription&)>* fixedFold = nullptr);
  // number of distinct spectra (spectrum file and scan number)
  std::size_t getNumSpectra() const;

  void generatePositiveTrainingSet(AlgIn& data,
                                   const double fdr,
                                   const double cpos,
                                   const bool trainBestPositive) const;
  void generateNegativeTrainingSet(AlgIn& data, const double cneg);
  void reserveTrainingSet(AlgIn& data, const double fdr) const;

  void recalculateSizes();
  // singleBin: the scores are the only cross-validation bin, so they need
  // no common scale, and they are kept as they are when they cannot be
  // rescaled
  void normalizeScores(double fdr, std::vector<double>& weights,
                       bool singleBin = false);
  // The rescaling of the scheme taken from these scores, which have to be
  // sorted, with q values for TEST and TRAIN: the score of the last PSM with
  // q < fdr above the median decoy becomes 0 and the median decoy -1 (TEST,
  // TRAIN), or the decoys get mean 0 and standard deviation 1
  // (TRAIN_DECOYS). NONE gives the identity.
  ScoreScale getScoreScale(MergeNormalization scheme, double fdr) const;
  // Scores these PSMs with w, sorts them and returns getScoreScale(scheme,
  // fdr), with the q values computed as merge does for a test bin, except
  // that a too good separation gives pi0 = 1 instead of an error. Used to
  // take the anchors of TRAIN and TRAIN_DECOYS from a training set.
  ScoreScale calcScoreScale(const std::vector<double>& w,
                            MergeNormalization scheme,
                            double fdr);
  // Rescales the scores and the weights by s. When s.scale <= 0 this is an
  // error, except with --no-terminate and for the only bin (singleBin),
  // which keeps its raw scores and weights as normalizeScores does; scheme
  // and fdr only word the message.
  void applyScoreScale(const ScoreScale& s,
                       std::vector<double>& weights,
                       bool singleBin,
                       MergeNormalization scheme,
                       double fdr);

  void weedOutRedundant();
  void weedOutRedundant(std::map<std::string, unsigned int>& peptideSpecCounts,
                        double specCountQvalThreshold);
  void weedOutRedundantTDC();


  void print(LabelType label, std::ostream& os = std::cout);

  inline double getPi0() const { return pi0_; }
  inline double getTargetDecoySizeRatio() const {
    return targetDecoySizeRatio_;
  }
  inline unsigned int size() const {
    return totalNumberOfTargets_ + totalNumberOfDecoys_;
  }
  inline unsigned int posSize() const { return totalNumberOfTargets_; }
  inline unsigned int negSize() const { return totalNumberOfDecoys_; }

  inline void addScoreHolder(const ScoreHolder& sh) { scores_.push_back(sh); }
  // OpenMS: replaces the score holders (built elsewhere, e.g. in parallel); call recalculateSizes after
  inline void setScoreHolders(std::vector<ScoreHolder>&& scores) { scores_ = std::move(scores); }
  inline void reserve(size_t n) { scores_.reserve(n); }

  inline void setNullTargetWinProb(const double nullTargetWinProb) {
    nullTargetWinProb_ = nullTargetWinProb;
  }

  // removes all PSMs and releases their memory
  void reset() {
    std::vector<ScoreHolder>().swap(scores_);
    totalNumberOfTargets_ = 0;
    totalNumberOfDecoys_ = 0;
  }
  void setUsePi0(bool usePi0);
  inline void setOutputRT(const bool is_output_rt) {
    is_output_rt_ = is_output_rt;
  }

 protected:
  bool usePi0_;

  double pi0_;
  double targetDecoySizeRatio_, nullTargetWinProb_;
  unsigned int totalNumberOfDecoys_, totalNumberOfTargets_;

  std::vector<ScoreHolder> scores_;

  void getScoreLabelPairs(std::vector<std::pair<double, bool> >& combined);
  // What ends the merge of the cross-validation bins, see mergeBins
  enum class MergeStatus {
    OK,
    MISSING_SCALES,           // train, train-decoys: no scale for each bin
    BIN_SEPARATES_TOO_WELL,   // a bin whose pi0 the merge estimates
    CANNOT_RESCALE,           // the anchors of a bin do not exist
    LIST_SEPARATES_TOO_WELL   // the merged list
  };
  // the one merge behind merge and countMergedTargets (Scores.cpp)
  template <class GetBin>
  MergeStatus mergeBins(std::size_t numBins,
                        std::size_t numScores,
                        GetBin getBin,
                        std::vector<std::vector<double> >& weights,
                        double fdr,
                        bool skipNormalizeScores,
                        MergeNormalization scheme,
                        const std::vector<ScoreScale>& binScales,
                        bool report);
  // Whether there are no PSMs, which normalizeScores reports as an error
  // (ignored with a message under --no-terminate)
  bool noScoredPsms() const;
  // the p value of each target, in ascending order, from the decoys that
  // score at least as high; the scores have to be sorted
  void getTargetPValues(std::vector<double>& pvals) const;
  // the selections of weedOutRedundant and weedOutRedundantTDC, which then
  // call postMergeStep
  void keepBestPsmPerPeptide(
      std::map<std::string, unsigned int>& peptideSpecCounts,
      double specCountQvalThreshold);
  // The same with the given hash of the peptides, which groups the PSMs:
  // for the unit tests, which force peptides to share a hash value.
  typedef std::uint64_t (*PeptideHash)(std::string_view peptide);
  void keepBestPsmPerPeptide(
      std::map<std::string, unsigned int>& peptideSpecCounts,
      double specCountQvalThreshold, PeptideHash hash);
  void keepBestPsmPerSpectrum();
  // Sorts the scores and returns whether their targets and decoys separate
  // too well, which postMergeStep reports as an error; estimates no pi0, so
  // it draws no random numbers.
  bool separatesTooWell();
  int countTargetsBelowFdr(
      const std::vector<std::pair<double, bool> >& combined,
      double fdr,
      bool skipDecoysPlusOne) const;
  // tooGoodSeparationIsError = false sets pi0 = 1 without a message when the
  // separation is too good; returns whether it is
  bool checkSeparationAndSetPi0(bool tooGoodSeparationIsError = true);
  // The feature and direction of getInitDirection and its number of targets
  // (-1 without features, with bestFeature -1)
  int bestInitDirection(double initialSelectionFdr, int& bestFeature,
                        bool& lowBest) const;
  bool is_output_rt_ = false;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*SCORES_H_*/
