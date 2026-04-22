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
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>

#include "FeatureMemoryPool.h"
#include "FeatureNames.h"
#include "LabelType.h"
#include "Normalizer.h"
#include "PSMDescription.h"
#include "PseudoRandom.h"
#include "ScoreHolder.h"

#include <boost/unordered/unordered_map.hpp>

class SetHandler;
class AlgIn;

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

  void merge(std::vector<Scores>& sv,
             double fdr,
             bool skipNormalizeScores,
             std::vector<std::vector<double> >& all_w);
  void postMergeStep();

  iterator begin() { return scores_.begin(); }
  iterator end() { return scores_.end(); }
  const_iterator begin() const { return scores_.begin(); }
  const_iterator end() const { return scores_.end(); }

  double calcScore(const double* features, const std::vector<double>& w) const;
  void scoreAndAddPSM(ScoreHolder& sh,
                      const std::vector<double>& rawWeights,
                      FeatureMemoryPool& featurePool);
  int calcScoresAndQvals(vector<double>& w,
                         double fdr,
                         bool skipDecoysPlusOne = false);
  int calcScores(vector<double>& w);
  int calcQvals(double fdr, bool skipDecoysPlusOne = false);
  void calcPep(const bool spline = false, const bool interpol = false, const bool from_q = false);
 
  void populateWithPSMs(SetHandler& setHandler);

  int getInitDirection(const double initialSelectionFdr,
                       std::vector<double>& direction);
  void createXvalSetsBySpectrum(std::vector<Scores>& train,
                                std::vector<Scores>& test,
                                unsigned int xval_fold,
                                FeatureMemoryPool& featurePool,
                                double decoyFractionTraining = 1.0,
                                unsigned int decoysPerTarget = 1u);

  void generatePositiveTrainingSet(AlgIn& data,
                                   const double fdr,
                                   const double cpos,
                                   const bool trainBestPositive);
  void generateNegativeTrainingSet(AlgIn& data, const double cneg);

  void recalculateSizes();
  void normalizeScores(double fdr, std::vector<double>& weights);

  void weedOutRedundant();
  void weedOutRedundant(std::map<std::string, unsigned int>& peptideSpecCounts,
                        double specCountQvalThreshold);
  void weedOutRedundantTDC();
  void weedOutRedundantMixMax();

  void printRetentionTime(ostream& outs, double fdr);
  unsigned getQvaluesBelowLevel(double level);

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

  inline void setNullTargetWinProb(const double nullTargetWinProb) {
    nullTargetWinProb_ = nullTargetWinProb;
  }

  std::vector<PSMDescription*>& getPsms(PSMDescription* pPSM) {
    return peptidePsmMap_[pPSM];
  }

  void reset() {
    scores_.clear();
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
  std::map<PSMDescription*, std::vector<PSMDescription*> > peptidePsmMap_;

  void reorderFeatureRows(
      FeatureMemoryPool& featurePool,
      bool isTarget,
      boost::unordered_map<double*, double*>& movedAddresses,
      size_t& idx);
  void getScoreLabelPairs(std::vector<std::pair<double, bool> >& combined);
  void checkSeparationAndSetPi0();
  bool is_output_rt_ = false;
};

#endif /*SCORES_H_*/
