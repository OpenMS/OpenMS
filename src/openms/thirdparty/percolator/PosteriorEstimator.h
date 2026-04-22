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


class PosteriorEstimator {
 public:
  PosteriorEstimator(){};
  virtual ~PosteriorEstimator(){};
  bool parseOptions(int argc, char** argv);
  string greeter();
  int run();
  static void estimatePEP(std::vector<std::pair<double, bool> >& combined,
          bool usePi0, double pi0, std::vector<double>& peps,
		      bool include_negative = false);
  static void getPValues(const std::vector<std::pair<double, bool> >& combined,
                         std::vector<double>& p);
  static void getQValues(double pi0,
                         const std::vector<std::pair<double, bool> >& combined,
                         std::vector<double>& q, bool skipDecoysPlusOne = false,
                         double nullTargetWinProb = 0.5);
  static void getQValuesFromP(double pi0, const std::vector<double>& p,
                              std::vector<double>& q);
  static void getQValuesFromPEP(const std::vector<double>& pep,
                              std::vector<double>& q);
  static bool checkSeparation(std::vector<double>& p);
  static double estimatePi0(std::vector<double>& p,
                            const unsigned int numBoot = 100);
  static void setReversed(bool status) {
	  reversed = status;
  }
  static void setNegative(bool negative) {
    includeNegativesInResult = negative;
  }
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
                       LogisticRegression& lr, bool usePi0, double pi0);
  static void binData(const std::vector<std::pair<double, bool> >& combined,
                      double pi0, std::vector<double>& medians,
                      std::vector<double>& negatives,
                      std::vector<double>& sizes);

  // used for standalone execution
  std::string targetFile, decoyFile;
  static bool reversed, pvalInput, includeNegativesInResult, competition, usePi0_;
  std::string resultFileName;
};

#endif /*POSTERIORESTIMATOR_H_*/
