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

#include "PosteriorEstimator.h"
#include "Transform.h"
#include "Globals.h"
#include "Version.h"
#include <exception>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;

static int noIntervals = 500;
static unsigned int numLambda = 100;
static double maxLambda = 0.5;

bool PosteriorEstimator::reversed = false;
bool PosteriorEstimator::pvalInput = false;
bool PosteriorEstimator::competition = false;
std::atomic<bool> PosteriorEstimator::includeNegativesInResult{false};
bool PosteriorEstimator::usePi0_ = true;

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

// OpenMS: bootstrap with the numbers drawn from the caller's state of the generator
template<class T> void bootstrapFrom(const vector<T>& in, vector<T>& out, uint64_t& state,
                                     size_t max_size = 1000) {
  out.clear();
  double n = static_cast<double>(in.size());
  size_t num_draw = min(in.size(), max_size);
  for (size_t ix = 0; ix < num_draw; ++ix) {
    size_t draw = (size_t)(PseudoRandom::lcg_uniform_rand(state) * n);
    out.push_back(in[draw]);
  }
  // sort in desending order
  sort(out.begin(), out.end());
}

double mymin(double a, double b) {
  return a > b ? b : a;
}

void PosteriorEstimator::estimatePEP(vector<pair<double, bool> >& combined,
    bool usePi0, double pi0, vector<double>& peps, bool include_negative) {
  // Logistic regression on the data
  LogisticRegression lr;
  estimate(combined, lr, usePi0, pi0);
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
  
  double top = min(1.0, exp(*max_element(peps.begin(), peps.end())));
  vector<double>::iterator pep = peps.begin();
  bool crap = false;
  for (; pep != peps.end(); ++pep) {
    if (crap) {
      *pep = top;
      continue;
    }
    *pep = exp(*pep);
    if (*pep >= top) {
      *pep = top;
      crap = true;
    }
  }
  partial_sum(peps.rbegin(), peps.rend(), peps.rbegin(), mymin);
}

void PosteriorEstimator::estimate(vector<pair<double, bool> >& combined,
    LogisticRegression& lr, bool usePi0, double pi0) {
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
      cerr << oss.str() << "No-terminate flag set: ignoring error and adding "
          << "random noise to scores to break ties for PEP calculation." << std::endl;

      medians.clear();
      negatives.clear();
      sizes.clear();

      double random_offset = PseudoRandom::lcg_uniform_rand() * 1e-20;
      int i = 0;
      vector<pair<double, bool> >::iterator elem = combined.begin();
      for (; elem != combined.end(); ++elem, ++i) {
        elem->first += i*random_offset;
      }
      binData(combined, pi0, medians, negatives, sizes);
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

      // OpenMS: LogisticRegression::setData pads to 4 bins, in the spline's
      // transformed coordinates. Upstream added i * lcg_uniform_rand() * 1e-20
      // here, which is lost to rounding unless the median is ~0, and the padding
      // could change the transform chosen from the medians. The random draw is
      // kept, so the random stream stays the same.
      PseudoRandom::lcg_uniform_rand();
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
  } else {
    getQValues(pi0, combined, q);
  }
  vector<pair<double, bool> >::const_iterator elem = combined.begin();
  for (; elem != combined.end(); ++elem)
  {
    if (!includeNegativesInResult && elem->second)
    {
      xvals.push_back(elem->first);
    }
    else if(includeNegativesInResult)
    {
      xvals.push_back(elem->first);
    }
  }
  vector<double>::iterator xval = xvals.begin();
  vector<double>::const_iterator qv = q.begin(), pep = peps.begin();
  if (resultFileName.empty()) {
    cout << "Score\tPEP\tq-value" << endl;
    for (; xval != xvals.end(); ++xval, ++pep, ++qv)
    {
      cout << *xval << "\t" << *pep << "\t" << *qv << endl;
    }
  } else {
    ofstream resultstream(resultFileName.c_str());
    resultstream << "Score\tPEP\tq-value" << endl;
    for (; xval != xvals.end(); ++xval, ++pep, ++qv)
    {
      resultstream << *xval << "\t" << *pep << "\t" << *qv << endl;
    }
    resultstream.close();
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
 * Assumes that scores are sorted in descending order
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
 * This is a reimplementation of 
 *   Crux/src/app/AssignConfidenceApplication.cpp::compute_decoy_qvalues_mixmax 
 * Which itself was a reimplementation of Uri Keich's code written in R.
 *
 * Assumes that scores are sorted in descending order
 * 
 * If pi0 == 1.0 this is equal to the "traditional" q-value calculation
 */
void PosteriorEstimator::getQValues(double pi0, 
    const vector<pair<double, bool> >& combined, vector<double>& q,
    bool skipDecoysPlusOne, double nullTargetWinProb) {
  double decoyFactor = nullTargetWinProb / (1.0 - nullTargetWinProb);

  std::vector<double> h_w_le_z, h_z_le_z; // N_{w<=z} and N_{z<=z}
  if (pi0 < 1.0) {
    getMixMaxCounts(combined, h_w_le_z, h_z_le_z);
  }

  double estPx_lt_zj = 0.0;
  double E_f1_mod_run_tot = 0.0;
  double fdr = 0.0;

  int n_z_ge_w = 1, n_w_ge_w = 0; // N_{z>=w} and N_{w>=w}
  if (skipDecoysPlusOne) n_z_ge_w = 0;
  
  std::vector<pair<double, bool> >::const_iterator myPair = combined.begin();
  int decoyQueue = 0, targetQueue = 0; // handles ties
  for ( ; myPair != combined.end(); ++myPair) {
    if (myPair->second) { 
      ++n_w_ge_w; // target PSM
      ++targetQueue;
    } else {
      ++n_z_ge_w; // decoy PSM
      ++decoyQueue;
    }
    
    // handles ties
    if (myPair+1 == combined.end() || myPair->first != (myPair+1)->first) {
      if (pi0 < 1.0 && decoyQueue > 0) {
        int j = static_cast<int>(h_w_le_z.size()) - (n_z_ge_w - 1);
        int cnt_w = static_cast<int>(h_w_le_z.at(static_cast<std::size_t>(j)));
        int cnt_z = static_cast<int>(h_z_le_z.at(static_cast<std::size_t>(j)));
        estPx_lt_zj = (double)(cnt_w - pi0*cnt_z) / ((1.0 - pi0)*cnt_z);
        estPx_lt_zj = estPx_lt_zj > 1 ? 1 : estPx_lt_zj;
        estPx_lt_zj = estPx_lt_zj < 0 ? 0 : estPx_lt_zj;
        E_f1_mod_run_tot += decoyQueue * estPx_lt_zj * (1.0 - pi0);
        if (VERB > 4) {
          std::cerr << "Mix-max num negatives correction: "
            << (1.0-pi0) * n_z_ge_w << " vs. " << E_f1_mod_run_tot << std::endl;
        }
      }
      
      if (includeNegativesInResult) {
        targetQueue += decoyQueue;
      }
      fdr = (n_z_ge_w * pi0 + E_f1_mod_run_tot) / (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
      for (int i = 0; i < targetQueue; ++i) {
        q.push_back((std::min)(fdr, 1.0));
      }
      decoyQueue = 0;
      targetQueue = 0;
    }
  }
  // Convert the FDRs into q-values.
  partial_sum(q.rbegin(), q.rend(), q.rbegin(), mymin);
}

// OpenMS: the smallest block of PSMs a task of getQValuesInTasks handles
[[maybe_unused]] static const size_t kMinQValueBlock = 4096;

/**
 * OpenMS: getQValues(pi0, combined, q, skipDecoysPlusOne, nullTargetWinProb) with the tasks of
 * an enclosing OpenMP parallel region, for a finite pi0 >= 1 (then getQValues makes no mix-max
 * correction), a finite decoy factor and an empty q; otherwise returns false and does nothing.
 * The PSMs are split into blocks of whole groups of equal scores (at fixed positions, not
 * depending on the number of threads). The targets and decoys of the blocks are counted first,
 * so each block computes the FDRs of its groups from the same counts as the loop of getQValues,
 * by the same expression. The q-values (the minimum FDR from a position to the end: the FDRs
 * are numbers in [0, 1], so the order of the minimisation does not matter) are the minima
 * within each block, and then the minimum of these and the q-value at the start of the next
 * block.
 */
bool PosteriorEstimator::getQValuesInTasks(double pi0,
    const vector<pair<double, bool> >& combined, vector<double>& q,
    bool skipDecoysPlusOne, double nullTargetWinProb) {
#ifdef _OPENMP
  const size_t n = combined.size();
  const double decoyFactor = nullTargetWinProb / (1.0 - nullTargetWinProb);
  if (!(pi0 >= 1.0) || !std::isfinite(pi0) || !(decoyFactor >= 0.0) ||
      !std::isfinite(decoyFactor) || !q.empty() || n < 2 * kMinQValueBlock ||
      n > static_cast<size_t>(INT_MAX) || !omp_in_parallel()) {
    return false;
  }
  const bool withNegatives = includeNegativesInResult;
  // the blocks [starts[b], starts[b + 1]): from multiples of kMinQValueBlock to the next group
  vector<size_t> starts(1, 0);
  for (size_t at = kMinQValueBlock; at < n; at += kMinQValueBlock) {
    size_t start = (std::max)(at, starts.back() + 1);
    while (start < n && combined[start - 1].first == combined[start].first) ++start;
    if (start < n) starts.push_back(start);
  }
  starts.push_back(n);
  const int numBlocks = static_cast<int>(starts.size() - 1);
  vector<int> targets(static_cast<size_t>(numBlocks)), decoys(static_cast<size_t>(numBlocks));
#pragma omp taskloop grainsize(1) shared(combined, starts, targets, decoys)
  for (int b = 0; b < numBlocks; ++b) {
    int t = 0, d = 0;
    for (size_t i = starts[static_cast<size_t>(b)]; i < starts[static_cast<size_t>(b) + 1]; ++i) {
      if (combined[i].second) ++t; else ++d;
    }
    targets[static_cast<size_t>(b)] = t;
    decoys[static_cast<size_t>(b)] = d;
  }
  // the counts before each block, and the position of its first q-value
  vector<int> zBefore(static_cast<size_t>(numBlocks)), wBefore(static_cast<size_t>(numBlocks));
  vector<size_t> qBefore(static_cast<size_t>(numBlocks) + 1);
  int n_z = skipDecoysPlusOne ? 0 : 1, n_w = 0;
  qBefore[0] = 0;
  for (size_t b = 0; b < static_cast<size_t>(numBlocks); ++b) {
    zBefore[b] = n_z;
    wBefore[b] = n_w;
    n_z += decoys[b];
    n_w += targets[b];
    qBefore[b + 1] = qBefore[b] + static_cast<size_t>(targets[b]) +
                     (withNegatives ? static_cast<size_t>(decoys[b]) : 0u);
  }
  q.resize(qBefore[static_cast<size_t>(numBlocks)]);
#pragma omp taskloop grainsize(1) shared(combined, starts, zBefore, wBefore, qBefore, q)
  for (int b = 0; b < numBlocks; ++b) {
    // the loop of getQValues on the block, without the mix-max correction (pi0 >= 1)
    const double E_f1_mod_run_tot = 0.0;
    int n_z_ge_w = zBefore[static_cast<size_t>(b)], n_w_ge_w = wBefore[static_cast<size_t>(b)];
    int decoyQueue = 0, targetQueue = 0;
    size_t out = qBefore[static_cast<size_t>(b)];
    const size_t end = starts[static_cast<size_t>(b) + 1];
    for (size_t i = starts[static_cast<size_t>(b)]; i < end; ++i) {
      if (combined[i].second) {
        ++n_w_ge_w;
        ++targetQueue;
      } else {
        ++n_z_ge_w;
        ++decoyQueue;
      }
      if (i + 1 == end || combined[i].first != combined[i + 1].first) {
        if (withNegatives) {
          targetQueue += decoyQueue;
        }
        double fdr = (n_z_ge_w * pi0 + E_f1_mod_run_tot) / (double)((std::max)(1, n_w_ge_w)) * decoyFactor;
        for (int k = 0; k < targetQueue; ++k) {
          q[out++] = (std::min)(fdr, 1.0);
        }
        decoyQueue = 0;
        targetQueue = 0;
      }
    }
    // the minima within the block
    for (size_t k = qBefore[static_cast<size_t>(b) + 1]; k-- > qBefore[static_cast<size_t>(b)] + 1;) {
      q[k - 1] = mymin(q[k], q[k - 1]);
    }
  }
  // the minimum of the q-values after each block
  vector<double> after(static_cast<size_t>(numBlocks), 2.0);  // above every FDR (none after the last block)
  for (size_t b = static_cast<size_t>(numBlocks) - 1; b-- > 0;) {
    const bool empty = qBefore[b + 1] == qBefore[b + 2];
    after[b] = empty ? after[b + 1] : mymin(after[b + 1], q[qBefore[b + 1]]);
  }
#pragma omp taskloop grainsize(1) shared(after, qBefore, q)
  for (int b = 0; b < numBlocks - 1; ++b) {
    for (size_t k = qBefore[static_cast<size_t>(b)]; k < qBefore[static_cast<size_t>(b) + 1]; ++k) {
      q[k] = mymin(after[static_cast<size_t>(b)], q[k]);
    }
  }
  return true;
#else
  (void)pi0; (void)combined; (void)q; (void)skipDecoysPlusOne; (void)nullTargetWinProb;
  return false;
#endif
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
        std::cerr << "No-terminate flag set: ignoring error." << std::endl;
      }
    } else {
      throw MyException(oss.str() + "Terminating.\n");
    }
  }
  double minPi0 = *min_element(pi0s.begin(), pi0s.end());
  
  vector<double> pBoot, mse;
  // Initialize the vector mse with zeroes.
  fill_n(back_inserter(mse), pi0s.size(), 0.0);
  // OpenMS: inside an OpenMP parallel region (CrossValidation, Scores::merge), the bootstraps
  // are tasks: each draws its numbers from the generator's state at its first draw (the
  // bootstraps draw min(n, 1000) numbers each, see bootstrap; PseudoRandom::advance), the
  // errors are added up in the order of the bootstraps as below, and the generator ends in
  // the state the loop below leaves.
  bool booted = false;
#ifdef _OPENMP
  if (numBoot > 1 && !p.empty() && omp_in_parallel()) {
    const uint64_t numDraw = min(p.size(), static_cast<size_t>(1000));
    const uint64_t seed = PseudoRandom::getSeed();
    const size_t numLambdas = lambdas.size();
    vector<double> pi0Boots(static_cast<size_t>(numBoot) * numLambdas);
    vector<std::exception_ptr> errors(numBoot);
#pragma omp taskloop grainsize(4) shared(p, lambdas, pi0Boots, errors)
    for (int boot = 0; boot < static_cast<int>(numBoot); ++boot) {
      try {
        uint64_t state = PseudoRandom::advance(seed, static_cast<uint64_t>(boot) * numDraw);
        vector<double> pBootOfTask;
        bootstrapFrom<double>(p, pBootOfTask, state);
        const size_t nBoot = pBootOfTask.size();
        for (size_t ix = 0; ix < numLambdas; ++ix) {
          vector<double>::const_iterator startOfTask =
              lower_bound(pBootOfTask.begin(), pBootOfTask.end(), lambdas[ix]);
          double Wl = (double)distance(startOfTask, pBootOfTask.cend());
          pi0Boots[static_cast<size_t>(boot) * numLambdas + ix] =
              Wl / static_cast<double>(nBoot) / (1. - lambdas[ix]);
        }
      } catch (...) {
        errors[static_cast<size_t>(boot)] = std::current_exception();
      }
    }
    for (const std::exception_ptr& error : errors) {
      if (error) std::rethrow_exception(error);
    }
    for (unsigned int boot = 0; boot < numBoot; ++boot) {
      for (unsigned int ix = 0; ix < lambdas.size(); ++ix) {
        double pi0Boot = pi0Boots[static_cast<size_t>(boot) * numLambdas + ix];
        // Estimated mean-squared error.
        mse[ix] += (pi0Boot - minPi0) * (pi0Boot - minPi0);
      }
    }
    PseudoRandom::setSeed(static_cast<unsigned long>(
        PseudoRandom::advance(seed, static_cast<uint64_t>(numBoot) * numDraw)));
    booted = true;
  }
#endif
  // Examine which lambda level that is most stable under bootstrap
  for (unsigned int boot = 0; !booted && boot < numBoot; ++boot) {
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

int PosteriorEstimator::run() {
  ifstream target(targetFile.c_str(), ios::in), decoy(decoyFile.c_str(),
                                                      ios::in);
  istream_iterator<double> tarIt(target), decIt(decoy), tarEnd, decEnd;
  // Merge a labeled version of the two lists into a combined list
  vector<pair<double, bool> > combined;
  vector<double> pvals;
  if (!pvalInput) {
      std::transform(tarIt, tarEnd, std::back_inserter(combined),
                       [](double value) { return make_my_pair(value, true); });
      size_t targetSize = combined.size();
      std::transform(decIt, decEnd, std::back_inserter(combined),
                       [](double value) { return make_my_pair(value, false); });

      if (VERB > 0) {
            std::cerr << "Read " << targetSize << " target scores and "
                      << (combined.size() - targetSize) << " decoy scores" << std::endl;
      }
  } else {
      std::copy(tarIt, tarEnd, std::back_inserter(pvals));
      std::sort(pvals.begin(), pvals.end());
      std::transform(pvals.begin(), pvals.end(), std::back_inserter(combined),
                       [](double value) { return make_my_pair(value, true); });

      size_t nDec = pvals.size();
      double step = 1.0 / 2.0 / static_cast<double>(nDec);
      for (size_t ix = 0; ix < nDec; ++ix) {
          combined.push_back(make_my_pair(step * static_cast<double>(1 + 2 * ix), false));
      }
      reversed = true;

      if (VERB > 0) {
            std::cerr << "Read " << pvals.size() << " statistics" << std::endl;
      }
  }
  if (reversed) {
    if (VERB > 0) {
      cerr << "Reversing all scores" << endl;
    }
  }
  if (reversed) // sorting in ascending order
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
  
  // Logistic regression on the data
  estimatePEP(combined, usePi0_, pi0, peps,includeNegativesInResult);
  finishStandalone(combined, peps, pvals, pi0);

  return true;
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

