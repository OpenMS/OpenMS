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

#ifndef POSTERIORESTIMATOR_H_
#define POSTERIORESTIMATOR_H_
#include <vector>
#include <string>
#include <utility>
#include <cfloat>

#include "LogisticRegression.h"
#include "PseudoRandom.h"

namespace OpenMS { namespace Internal { namespace Percolator {



class PosteriorEstimator {
 public:
  PosteriorEstimator(){};
  virtual ~PosteriorEstimator(){};
  string greeter();
  int run();
  // PEP = oddsFactor * (fitted decoy odds), capped at 1; oddsFactor is the
  // expected number of incorrect targets per decoy (2 with RESET).
  // reversed: combined is sorted by increasing score, lower scores are better
  static void estimatePEP(std::vector<std::pair<double, bool> >& combined,
          bool usePi0, double pi0, std::vector<double>& peps,
		      bool include_negative = false, double oddsFactor = 1.0,
		      bool reversed = false);
  static void getPValues(const std::vector<std::pair<double, bool> >& combined,
                         std::vector<double>& p);
  // q holds one value per target, or per entry with includeNegatives
  static void getQValues(double pi0,
                         const std::vector<std::pair<double, bool> >& combined,
                         std::vector<double>& q, bool skipDecoysPlusOne = false,
                         double nullTargetWinProb = 0.5,
                         bool includeNegatives = false);
  static void getQValuesFromP(double pi0, const std::vector<double>& p,
                              std::vector<double>& q);
  static void getQValuesFromPEP(const std::vector<double>& pep,
                              std::vector<double>& q);
  // For p value input with the decoys in the result (-d): the q values of
  // all entries of combined, the targets with those of getQValuesFromP.
  static std::vector<double> getAllQValuesFromP(double pi0,
      const std::vector<double>& p, const std::vector<double>& targetQ,
      const std::vector<std::pair<double, bool> >& combined);
  static bool checkSeparation(std::vector<double>& p);
  static double estimatePi0(std::vector<double>& p,
                            const unsigned int numBoot = 100);
  // Sets the process-wide settings of the estimation, which the options of
  // qvality and earlier runs change, back to their initial values.
  static void resetToDefaults();
  static void setUsePi0(bool usePi0) {
    usePi0_ = usePi0;
  }
 protected:
  void finishStandalone(std::vector<std::pair<double, bool> >& combined,
                        const std::vector<double>& peps,
                        const std::vector<double>& p, double pi0);
  static void getMixMaxCounts(const std::vector<std::pair<double, bool> >& combined,
                       std::vector<double>& h_w_le_z,
                       std::vector<double>& h_z_le_z);

  static void estimate(std::vector<std::pair<double, bool> >& combined,
                       LogisticRegression& lr, bool usePi0, double pi0,
                       bool reversed);
  // PEPs by I-spline or PAVA regression, from the target q values or, with
  // tdcPep_, from the decoy rate; same entries as finishStandalone reports
  void estimateMonotonePEP(const std::vector<std::pair<double, bool> >& combined,
                           const std::vector<double>& p, double pi0,
                           std::vector<double>& peps) const;
  // Reads the whitespace-separated scores of a qvality input file, and
  // throws if the file cannot be opened, holds no score, or holds something
  // other than a finite number.
  static void readScores(const std::string& fileName, const char* role,
                         std::vector<double>& scores);
  static void binData(const std::vector<std::pair<double, bool> >& combined,
                      double pi0, std::vector<double>& medians,
                      std::vector<double>& negatives,
                      std::vector<double>& sizes);

  // used for standalone execution
  std::string targetFile, decoyFile;
  static bool pvalInput, competition, usePi0_;
  // lower scores are better (-r, p value input); report the decoys too (-d)
  bool reversed_ = false, includeNegativesInResult_ = false;
  std::string resultFileName;
  // PEP method: I-spline by default, PAVA, or P-IRLS; regression on the
  // target q values by default or on the target-decoy labels (tdcPep_)
  bool irlsPep_ = false, pavaPep_ = false, tdcPep_ = false, rankPep_ = false;
  bool nothingToRun_ = false;  // --version, --help: run does nothing
};

}}}  // namespace OpenMS::Internal::Percolator

#endif /*POSTERIORESTIMATOR_H_*/
