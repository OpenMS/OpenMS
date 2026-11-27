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

#include<cmath>
#include<vector>
#include<utility>
#include<cstdlib>
#include<limits.h>
#include<fstream>
#include<sstream>
#include<iterator>
#include<algorithm>
#include<numeric>
#include<functional>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Version.h"
#include "PosteriorEstimator.h"
#include "Transform.h"
#include "Globals.h"
#include "IsotonicPEP.h"
#include "FdrWalk.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;


static const int kDefaultNoIntervals = 500;
static int noIntervals = kDefaultNoIntervals;
static unsigned int numLambda = 100;
static double maxLambda = 0.5;

bool PosteriorEstimator::pvalInput = false;
bool PosteriorEstimator::competition = false;
bool PosteriorEstimator::usePi0_ = true;

void PosteriorEstimator::resetToDefaults() {
  pvalInput = false;
  competition = false;
  usePi0_ = true;
  noIntervals = kDefaultNoIntervals;
}

pair<double, bool> make_my_pair(double d, bool b) {
  return make_pair(d, b);
}

class IsDecoy {
  public:
    bool operator()(const pair<double, bool>& aPair) {
      return !aPair.second;
    }
};

bool isMixed(const pair<double, bool>& aPair) {
  return aPair.second;
}

template<class T> void bootstrap(const vector<T>& in, vector<T>& out,
                                 size_t max_size = 1000) {
  out.clear();
  double n = static_cast<double>(in.size());
  size_t num_draw = min(in.size(), max_size);
  for (size_t ix = 0; ix < num_draw; ++ix) {
    size_t draw = (size_t)(PseudoRandom::lcg_uniform_rand() * n);
    out.push_back(in[draw]);
  }
  // sort in desending order
  sort(out.begin(), out.end());
}

double mymin(double a, double b) {
  return a > b ? b : a;
}

void PosteriorEstimator::estimatePEP(vector<pair<double, bool> >& combined,
    bool usePi0, double pi0, vector<double>& peps, bool include_negative,
    double oddsFactor, bool reversed) {
  // Logistic regression on the data
  LogisticRegression lr;
  estimate(combined, lr, usePi0, pi0, reversed);
  vector<double> xvals(0);
  vector<pair<double, bool> >::const_iterator elem = combined.begin();
  for (; elem != combined.end(); ++elem) {
    if (elem->second || include_negative) { // target PSM
      xvals.push_back(elem->first);
    }
  }
  lr.predict(xvals, peps);
#define OUTPUT_DEBUG_FILES
#undef OUTPUT_DEBUG_FILES
#ifdef OUTPUT_DEBUG_FILES
  ofstream drFile("decoyRate.all", ios::out), xvalFile("xvals.all", ios::out);
  ostream_iterator<double> drIt(drFile, "\n"), xvalIt(xvalFile, "\n");
  copy(peps.begin(), peps.end(), drIt);
  copy(xvals.begin(), xvals.end(), xvalIt);
#endif
  if (peps.size() == 0) return;
  
  double top =
      min(1.0, oddsFactor * exp(*max_element(peps.begin(), peps.end())));
  vector<double>::iterator pep = peps.begin();
  bool crap = false;
  for (; pep != peps.end(); ++pep) {
    if (crap) {
      *pep = top;
      continue;
    }
    *pep = oddsFactor * exp(*pep);
    if (*pep >= top) {
      *pep = top;
      crap = true;
    }
  }
  partial_sum(peps.rbegin(), peps.rend(), peps.rbegin(), mymin);
}

void PosteriorEstimator::estimate(vector<pair<double, bool> >& combined,
    LogisticRegression& lr, bool usePi0, double pi0, bool reversed) {
  // switch sorting order if we do not use mix-max
  if (!reversed && !usePi0) {
    reverse(combined.begin(), combined.end());
  }
  vector<double> medians, negatives, sizes;
  binData(combined, pi0, medians, negatives, sizes);
  if (medians.size() < 2) {
    ostringstream oss;
    oss << "ERROR: Only 1 bin available for PEP estimation, "
        << "no distinguishing feature present." << std::endl;
    if (NO_TERMINATE) {
      cerr << oss.str() << "No-terminate flag set: ignoring error and giving "
          << "all scores the same PEP." << std::endl;
      // The bin is padded to 4 below, with the same counts, and the fit is
      // flat. The noise of 1e-20 added to the scores before broke no ties
      // unless the scores were close to 0, and ordered tied scores by their
      // label. The draw keeps the later random numbers as they were.
      PseudoRandom::lcg_uniform_rand();
    } else {
      throw MyException(oss.str());
    }
  }

  if (medians.size() < 4) {
    ostringstream oss;
    oss << "ERROR: Fewer than 4 bins available for PEP estimation, "
        << "cannot perform quadratic spline." << std::endl;
    if (NO_TERMINATE) {
      if (medians.size() > 0) {
        cerr << oss.str() << "No-terminate flag set: ignoring error and adding "
             << "additional bins to get 4 bins." << std::endl;
      } else {
        cerr << oss.str() << "No-terminate flag set: ignoring error and skipping "
             << "PEP calculation." << std::endl;
        return;
      }

      // Copies of the last bin at distinct scores beyond it, in the order of
      // the bins: equal scores give intervals of length 0 in the spline.
      PseudoRandom::lcg_uniform_rand();  // as above
      double step = 1e-6 * max(1.0, max(fabs(medians.back()),
          fabs(medians.back() - medians.front())));
      if (medians.front() > medians.back()) {
        step = -step;
      }
      while (medians.size() < 4) {
        medians.push_back(medians.back() + step);
        negatives.push_back(negatives.back());
        sizes.push_back(sizes.back());
      }
    } else {
      throw MyException(oss.str());
    }
  }
  
  // the IRLS implementation requires the medians to be in ascending order
  if (medians.size() > 0 && medians.front() > medians.back()) {
    reverse(medians.begin(), medians.end());
    reverse(negatives.begin(), negatives.end());
    reverse(sizes.begin(), sizes.end());
  }
  
  lr.setData(medians, negatives, sizes);
  lr.roughnessPenaltyIRLS();
  // restore sorting order
  if (!reversed && !usePi0) {
    reverse(combined.begin(), combined.end());
  }
}

// Estimates q-values and prints
void PosteriorEstimator::finishStandalone(
    vector<pair<double, bool> >& combined, const vector<double>& peps,
    const vector<double>& p, double pi0) {
  vector<double> q(0), xvals(0);
  if (pvalInput) {
    getQValuesFromP(pi0, p, q);
    if (includeNegativesInResult_) {
      // one q value per line, as for score input: the decoys are printed too
      q = getAllQValuesFromP(pi0, p, q, combined);
    }
  } else {
    getQValues(pi0, combined, q, false, 0.5, includeNegativesInResult_);
  }
  vector<pair<double, bool> >::const_iterator elem = combined.begin();
  for (; elem != combined.end(); ++elem)
  {
    if (!includeNegativesInResult_ && elem->second)
    {
      xvals.push_back(elem->first);
    }
    else if(includeNegativesInResult_)
    {
      xvals.push_back(elem->first);
    }
  }
  vector<double>::iterator xval = xvals.begin();
  vector<double>::const_iterator qv = q.begin(), pep = peps.begin();
  // '\n' instead of endl: a flush per line costs one write() call per line
  if (resultFileName.empty()) {
    cout << "Score\tPEP\tq-value\n";
    for (; xval != xvals.end(); ++xval, ++pep, ++qv)
    {
      cout << *xval << "\t" << *pep << "\t" << *qv << '\n';
    }
    flushStandardOutput();
  } else {
    ofstream resultstream(resultFileName.c_str());
    resultstream << "Score\tPEP\tq-value\n";
    for (; xval != xvals.end(); ++xval, ++pep, ++qv)
    {
      resultstream << *xval << "\t" << *pep << "\t" << *qv << '\n';
    }
    closeOutputFile(resultstream, resultFileName);
  }
}

/*
 * If pi0 == 1.0 this is equal to the "traditional" binning
 */
void PosteriorEstimator::binData(const vector<pair<double, bool> >& combined,
    double pi0, vector<double>& medians, vector<double> & negatives,
    vector<double> & sizes) {
  std::vector<double> h_w_le_z, h_z_le_z; // N_{w<=z} and N_{z<=z}
  getMixMaxCounts(combined, h_w_le_z, h_z_le_z);
  
  double estPx_lt_zj = 0.0;
  double E_f1_mod_run_tot = 0.0;
  
  int binsLeft = noIntervals - 1;
  double targetedBinSize = max(static_cast<double>(combined.size()) / (double)(noIntervals), 1.0);
  
  std::vector<pair<double, bool> >::const_iterator myPair = combined.begin();
  int n_z_ge_w = 0, sum_n_z_ge_w = 0; // N_{z>=w} in bin and total
  int decoyQueue = 0, psmsInBin = 0, binStartIdx = 0;
  for (; myPair != combined.end(); ++myPair) {
    if (!(myPair->second)) { // decoy PSM
      ++n_z_ge_w;
      ++decoyQueue;
    }
    ++psmsInBin;
    
    // handles ties
    if (myPair+1 == combined.end() || myPair->first != (myPair+1)->first) {
      if (pi0 < 1.0 && decoyQueue > 0) {
        int j = static_cast<int>(h_w_le_z.size()) - sum_n_z_ge_w - n_z_ge_w;
        int cnt_w = static_cast<int>(h_w_le_z.at(static_cast<std::size_t>(j)));
        int cnt_z = static_cast<int>(h_z_le_z.at(static_cast<std::size_t>(j)));
        estPx_lt_zj = (double)(cnt_w - pi0*cnt_z) / ((1.0 - pi0)*cnt_z);
        estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
        estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
        E_f1_mod_run_tot += decoyQueue * estPx_lt_zj * (1.0 - pi0);
      }
      decoyQueue = 0;
      
      if (static_cast<int>(combined.size()) - binStartIdx - psmsInBin <= binsLeft * targetedBinSize) {
        double median = combined.at(static_cast<std::size_t>(binStartIdx + psmsInBin / 2)).first;
        double numNegatives = n_z_ge_w * pi0 + E_f1_mod_run_tot;
        double numPsmsCorrected = psmsInBin - n_z_ge_w + numNegatives;
        
        if (medians.size() > 0 && *(medians.rbegin()) == median) {
          *(negatives.rbegin()) += numNegatives;
          *(sizes.rbegin()) += numPsmsCorrected;
        } else {
          medians.push_back(median);
          negatives.push_back(numNegatives);
          sizes.push_back(numPsmsCorrected);
        }
        
        if (VERB > 4) {
          std::cerr << "Median = " << median << ", Num psms = " << psmsInBin 
                    << ", Num psms corrected = " << numPsmsCorrected
                    << ", Num decoys = " << n_z_ge_w 
                    << ", Num negatives = " << numNegatives << std::endl;
        }
        binStartIdx += psmsInBin;
        --binsLeft;
        
        psmsInBin = 0;
        E_f1_mod_run_tot = 0.0;
        sum_n_z_ge_w += n_z_ge_w;
        n_z_ge_w = 0;
      }
    }
  }
}

/**
 * The mix-max counts of binData. Assumes that scores are sorted in
 * descending order
 */
void PosteriorEstimator::getMixMaxCounts(const vector<pair<double, bool> >& combined,
    std::vector<double>& h_w_le_z, std::vector<double>& h_z_le_z) {
  int cnt_z = 0, cnt_w = 0, queue = 0;
  std::vector<pair<double, bool> >::const_reverse_iterator myPairRev = combined.rbegin();
  for ( ; myPairRev != combined.rend(); ++myPairRev) {
    if (myPairRev->second) {
      ++cnt_w; // target PSM
    } else {
      ++cnt_z; // decoy PSM
      ++queue;
    }
    
    // handles ties
    if (myPairRev+1 == combined.rend() || myPairRev->first != (myPairRev+1)->first) {
      for (int i = 0; i < queue; ++i) {
        h_w_le_z.push_back(static_cast<double>(cnt_w));
        h_z_le_z.push_back(static_cast<double>(cnt_z));
      }
      queue = 0;
    }
  }
}

/**
 * The q values of the entries of combined, which have to be sorted with the
 * best score first: one per target, or one per entry with includeNegatives.
 * The FDR estimates are those of walkFdrEstimates; with pi0 == 1.0 this is
 * the "traditional" q-value calculation.
 */
void PosteriorEstimator::getQValues(double pi0,
    const vector<pair<double, bool> >& combined, vector<double>& q,
    bool skipDecoysPlusOne, double nullTargetWinProb, bool includeNegatives) {
  double decoyFactor = nullTargetWinProb / (1.0 - nullTargetWinProb);
  walkFdrEstimates(
      combined.size(),
      [&combined](std::size_t ix) { return combined[ix].second; },
      [&combined](std::size_t ix) {
        return combined[ix].first == combined[ix + 1].first;
      },
      pi0, decoyFactor, skipDecoysPlusOne,
      [&q, includeNegatives](std::size_t begin, std::size_t end,
                             int numGroupTargets, int, double fdr) {
        const std::size_t numValues =
            includeNegatives ? end - begin
                             : static_cast<std::size_t>(numGroupTargets);
        for (std::size_t i = 0; i < numValues; ++i) {
          q.push_back((std::min)(fdr, 1.0));
        }
      });
  // Convert the FDRs into q-values.
  partial_sum(q.rbegin(), q.rend(), q.rbegin(), mymin);
}

void PosteriorEstimator::getQValuesFromP(double pi0,
                                         const vector<double>& p, vector<double> & q) {
	double m = (double)p.size();
	int nP = 1;
	// assuming combined sorted in decending order
	for (vector<double>::const_iterator myP = p.begin(); myP != p.end(); ++myP, ++nP) {
		q.push_back((*myP * m * pi0) / (double)nP);
	}
	partial_sum(q.rbegin(), q.rend(), q.rbegin(), mymin);
}

/**
 * The q values of all entries of combined for p value input, where the
 * targets are the p values and the decoys are the quantiles of the uniform
 * distribution that run() adds. A target keeps its value of
 * getQValuesFromP (targetQ). A decoy gets the same estimate at its own p
 * value d, pi0 * m * d / #{p <= d}, limited by the q value of the next
 * target, so that the values do not decrease along the list. Both p and
 * combined have to be sorted in ascending order.
 */
vector<double> PosteriorEstimator::getAllQValuesFromP(double pi0,
    const vector<double>& p, const vector<double>& targetQ,
    const vector<pair<double, bool> >& combined) {
  const double m = static_cast<double>(p.size());
  vector<double> q;
  q.reserve(combined.size());
  size_t numTargets = 0;  // targets written so far
  size_t numLe = 0;       // p values <= the current decoy
  for (const auto& elem : combined) {
    if (elem.second) {
      // a target of combined is a p value, so targetQ has a value for it
      q.push_back(targetQ.at(numTargets++));
      continue;
    }
    while (numLe < p.size() && p[numLe] <= elem.first) {
      ++numLe;
    }
    double fdr = pi0 * m * elem.first /
                 static_cast<double>((std::max)(numLe, static_cast<size_t>(1)));
    if (numLe < targetQ.size()) {
      fdr = (std::min)(fdr, targetQ[numLe]);
    }
    q.push_back(fdr);
  }
  return q;
}

// assumes pep are sorted in ascending order
void PosteriorEstimator::getQValuesFromPEP(const vector<double>& pep, vector<double>& q) {
	int nP = 1;
	double sum = 0.0;
	for (vector<double>::const_iterator myP = pep.begin(); myP != pep.end(); ++myP, ++nP) {
	  sum += *myP;
	  q.push_back(sum / (double)nP);
	}
	partial_sum(q.rbegin(), q.rend(), q.rbegin(), mymin);
}

void PosteriorEstimator::getPValues(const vector<pair<double, bool> >& combined,
                                    vector<double>& p) {
  // assuming combined sorted in best hit first order
  vector<pair<double, bool> >::const_iterator myPair = combined.begin();
  size_t nDecoys = 1, posSame = 0, negSame = 0;
  for ( ; myPair != combined.end(); ++myPair) {
    if (myPair->second) { // isTarget
      ++posSame;
    } else { // isDecoy
      ++negSame;
    }
    if (myPair+1 == combined.end() || myPair->first != (myPair+1)->first) {
      for (size_t ix = 0; ix < posSame; ++ix) {
        p.push_back(static_cast<double>(nDecoys) + 
          static_cast<double>(negSame * (ix + 1)) / (double)(posSame + 1) );
      }
      nDecoys += negSame;
      negSame = 0;
      posSame = 0;
    }
  }
  // p values sorted in ascending order
transform(p.begin(), p.end(), p.begin(), 
          [nDecoys](double value) { return value / static_cast<double>(nDecoys); });
}

bool PosteriorEstimator::checkSeparation(std::vector<double>& p) {
  double minLambda = (1.0 / numLambda) * maxLambda;
  // Find the index of the first element in p that is < lambda.
  // N.B. Assumes p is sorted in ascending order.
  std::vector<double>::iterator start = lower_bound(p.begin(), p.end(), minLambda);
  // Calculates the difference in index between start and end
  size_t Wl = static_cast<std::size_t>(distance(start, p.end()));
  return (Wl == 0u);
}

/*
 * Described in Storey, "A direct approach to false discovery rates."
 * JRSS 2002.
 */
double PosteriorEstimator::estimatePi0(vector<double>& p,
                                       const unsigned int numBoot) {
  vector<double> lambdas, pi0s;
  vector<double>::iterator start;
  size_t n = p.size();
  // Calculate pi0 for different values for lambda
  // N.B. numLambda and maxLambda are global variables.
  for (unsigned int ix = 0; ix <= numLambda; ++ix) {
    double lambda = ((ix + 1) / (double)numLambda) * maxLambda;
    // Find the index of the first element in p that is < lambda.
    // N.B. Assumes p is sorted in ascending order.
    start = lower_bound(p.begin(), p.end(), lambda);
    // Calculates the difference in index between start and end
    double Wl = (double)distance(start, p.end());
    double pi0 = Wl / static_cast<double>(n) / (1. - lambda);
    if (pi0 > 0.0) {
      lambdas.push_back(lambda);
      pi0s.push_back(pi0);
    }
  }
  
  if (pi0s.size() == 0) {
    ostringstream oss;
    oss << "Error in the input data: too good separation between target "
        << "and decoy PSMs.\n";
    if (NO_TERMINATE) {
      cerr << oss.str();
      if (usePi0_) {
        std::cerr << "No-terminate flag set: setting pi0 = 1 and ignoring error." << std::endl;
        return 1.0;
      } else {
        throw MyException(oss.str() + "Cannot continue without pi0 estimate. Terminating.\n");
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }
  double minPi0 = *min_element(pi0s.begin(), pi0s.end());
  
  vector<double> pBoot, mse;
  // Initialize the vector mse with zeroes.
  fill_n(back_inserter(mse), pi0s.size(), 0.0);
  // Examine which lambda level that is most stable under bootstrap
  for (unsigned int boot = 0; boot < numBoot; ++boot) {
    // Create an array of bootstrapped p-values, and sort in ascending order.
    bootstrap<double> (p, pBoot);
    n = pBoot.size();
    for (unsigned int ix = 0; ix < lambdas.size(); ++ix) {
      start = lower_bound(pBoot.begin(), pBoot.end(), lambdas[ix]);
      double Wl = (double)distance(start, pBoot.end());
      double pi0Boot = Wl / static_cast<double>(n) / (1. - lambdas[ix]);
      // Estimated mean-squared error.
      mse[ix] += (pi0Boot - minPi0) * (pi0Boot - minPi0);
    }
  }
  // Which index did the iterator get?
  unsigned int minIx = static_cast<unsigned int>(distance(mse.begin(), 
                                min_element(mse.begin(), mse.end())));
  double pi0 = max(min(pi0s[minIx], 1.0), 0.0);
  return pi0;
}

void PosteriorEstimator::readScores(const std::string& fileName,
                                    const char* role,
                                    std::vector<double>& scores) {
  ifstream in(fileName.c_str(), ios::in);
  if (!in.is_open()) {
    throw MyException("ERROR: Could not open the " + std::string(role) +
                      " file " + fileName + ".\n");
  }
  // strtod on each whitespace-separated token gives the same values as
  // operator>>, which also converts with strtod, but reports the token at
  // which reading stops instead of ending the list there
  std::string token;
  while (in >> token) {
    char* end = nullptr;
    double value = std::strtod(token.c_str(), &end);
    if (end == token.c_str() || *end != '\0' || !std::isfinite(value)) {
      ostringstream oss;
      oss << "ERROR: Value " << scores.size() + 1 << " of the " << role
          << " file " << fileName << " is not a finite number: \"" << token
          << "\"." << std::endl;
      throw MyException(oss.str());
    }
    scores.push_back(value);
  }
  if (in.bad()) {
    throw MyException("ERROR: Could not read the " + std::string(role) +
                      " file " + fileName + ".\n");
  }
  if (scores.empty()) {
    throw MyException("ERROR: The " + std::string(role) + " file " +
                      fileName + " contains no values.\n");
  }
}

int PosteriorEstimator::run() {
  if (nothingToRun_) {
    return true;  // parseOptions printed the version or the help
  }
  // Merge a labeled version of the two lists into a combined list
  vector<pair<double, bool> > combined;
  vector<double> pvals;
  if (!pvalInput) {
      vector<double> targetScores, decoyScores;
      readScores(targetFile, "target score", targetScores);
      readScores(decoyFile, "null score", decoyScores);
      std::transform(targetScores.begin(), targetScores.end(),
                     std::back_inserter(combined),
                     [](double value) { return make_my_pair(value, true); });
      size_t targetSize = combined.size();
      std::transform(decoyScores.begin(), decoyScores.end(),
                     std::back_inserter(combined),
                     [](double value) { return make_my_pair(value, false); });

      if (VERB > 0) {
            std::cerr << "Read " << targetSize << " target scores and "
                      << (combined.size() - targetSize) << " decoy scores" << std::endl;
      }
  } else {
      readScores(targetFile, "p value", pvals);
      std::sort(pvals.begin(), pvals.end());
      std::transform(pvals.begin(), pvals.end(), std::back_inserter(combined),
                       [](double value) { return make_my_pair(value, true); });

      size_t nDec = pvals.size();
      double step = 1.0 / 2.0 / static_cast<double>(nDec);
      for (size_t ix = 0; ix < nDec; ++ix) {
          combined.push_back(make_my_pair(step * static_cast<double>(1 + 2 * ix), false));
      }
      reversed_ = true;

      if (VERB > 0) {
            std::cerr << "Read " << pvals.size() << " statistics" << std::endl;
      }
  }
  if (reversed_) {
    if (VERB > 0) {
      cerr << "Reversing all scores" << endl;
    }
  }
  // std::sort here and above, not parallelSort: -0 and 0 in the input compare
  // equal but print differently, and only a sort on one thread puts them in
  // an order that does not depend on the number of threads.
  if (reversed_) // sorting in ascending order
  {
    sort(combined.begin(), combined.end());
  }
  else	// sorting in decending order
  {
    sort(combined.begin(), combined.end(), greater<pair<double, bool> > ());
  }
  if (!pvalInput) {
    getPValues(combined, pvals);
  }
  vector<double> peps;
  
  double pi0 = 1.0;
  if (usePi0_) {
    pi0 = estimatePi0(pvals);
    if (pi0 < 0) { //NOTE there was an error
      return false;
    }
    if (VERB > 1) {
      std::cerr << "Selecting pi_0=" << pi0 << std::endl;
    }
  }
  
  if (irlsPep_) {
    // Logistic regression on the data
    estimatePEP(combined, usePi0_, pi0, peps, includeNegativesInResult_, 1.0,
                reversed_);
  } else {
    estimateMonotonePEP(combined, pvals, pi0, peps);
  }
  finishStandalone(combined, peps, pvals, pi0);

  return true;
}

void PosteriorEstimator::estimateMonotonePEP(
    const vector<pair<double, bool> >& combined, const vector<double>& p,
    double pi0, vector<double>& peps) const {
  // combined is sorted best first; the regressions expect larger x = better
  const double sign = reversed_ ? -1.0 : 1.0;
  InferPEP reg(!pavaPep_);
  peps.clear();
  if (tdcPep_ && !pvalInput) {
    // decoy rate p(x) and PEP = p/(1-p); with separate target and decoy
    // searches the odds are scaled by pi0 and the ratio of targets to decoys
    vector<double> isDecoy, x;
    double numTargets = 0.0, numDecoys = 0.0;
    for (const auto& elem : combined) {
      isDecoy.push_back(elem.second ? 0.0 : 1.0);
      x.push_back(sign * elem.first);
      (elem.second ? numTargets : numDecoys) += 1.0;
    }
    const double oddsFactor =
        usePi0_ ? pi0 * numTargets / std::max(1.0, numDecoys) : 1.0;
    vector<double> all = rankPep_ ? reg.tdc_to_pep(isDecoy, {}, oddsFactor)
                                  : reg.tdc_to_pep(isDecoy, x, oddsFactor);
    for (size_t ix = 0; ix < combined.size(); ++ix) {
      if (combined[ix].second || includeNegativesInResult_) {
        peps.push_back(all[ix]);
      }
    }
    return;
  }
  if (tdcPep_ && VERB > 0) {
    cerr << "Warning: --tdc-pep needs target and null scores, using the "
         << "q values of the targets instead." << endl;
  }
  // raw PEPs from the target q values, n q_n - (n-1) q_{n-1}, regressed on
  // the score or the rank
  vector<double> q, targetQ, targetX;
  if (pvalInput) {
    getQValuesFromP(pi0, p, q);
    if (includeNegativesInResult_) {
      q = getAllQValuesFromP(pi0, p, q, combined);
    }
  } else {
    getQValues(pi0, combined, q, false, 0.5, includeNegativesInResult_);
  }
  // q holds one value per target, or per entry with includeNegativesInResult_
  const bool allEntries = includeNegativesInResult_;
  for (size_t ix = 0, iq = 0; ix < combined.size(); ++ix) {
    if (!combined[ix].second && !allEntries) continue;
    if (combined[ix].second && iq < q.size()) {
      targetQ.push_back(q[iq]);
      targetX.push_back(sign * combined[ix].first);
    }
    ++iq;
  }
  vector<double> targetPep = rankPep_ ? reg.q_to_pep(targetQ)
                                      : reg.qns_to_pep(targetQ, targetX);
  if (!includeNegativesInResult_) {
    peps = targetPep;
    return;
  }
  // decoys get PEPs interpolated by their q values between the targets
  // around them, as in Percolator
  targetPep.push_back(1.0);
  targetQ.push_back(1.0);
  size_t it = 0;
  double lastQ = 0.0, lastPep = 0.0;
  for (size_t ix = 0; ix < combined.size(); ++ix) {
    const double qv = (allEntries && ix < q.size()) ? q[ix] : lastQ;
    if (combined[ix].second) {
      peps.push_back(targetPep[it]);
      lastQ = targetQ[it];
      lastPep = targetPep[it];
      ++it;
    } else {
      peps.push_back(reg.interpolate(qv, lastQ, targetQ[it], lastPep,
                                     targetPep[it]));
    }
  }
}

string PosteriorEstimator::greeter() {
  ostringstream oss;
  oss << "qvality version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Distributed under Apache 2.0 License" << endl;
  oss << "Written by Lukas Käll (lukas.kall@cbr.su.se) in the" << endl;
  oss << "Department of Genome Sciences at the University of Washington."
      << endl;
  return oss.str();
}

}}}  // namespace OpenMS::Internal::Percolator
