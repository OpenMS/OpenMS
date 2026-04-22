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
using namespace std;

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "Option.h"
#include "Version.h"
#include "PosteriorEstimator.h"
#include "Transform.h"
#include "Globals.h"

static int noIntervals = 500;
static unsigned int numLambda = 100;
static double maxLambda = 0.5;

bool PosteriorEstimator::reversed = false;
bool PosteriorEstimator::pvalInput = false;
bool PosteriorEstimator::competition = false;
bool PosteriorEstimator::includeNegativesInResult = false;
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

      double random_offset = PseudoRandom::lcg_uniform_rand() * 1e-20;
      for (int i = medians.size(); i < 4; ++i) {
        medians.push_back(medians.back() + i*random_offset);
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

bool PosteriorEstimator::parseOptions(int argc, char** argv) {
  // init
  ostringstream intro;
  intro << greeter() << endl;
  intro << "Usage:" << endl;
  intro << "   qvality [options] target_file null_file" << endl << "or"
      << endl;
  intro << "   qvality [options] pvalue_file" << endl << endl;
  intro
      << "target_file and null_file are files containing scores from a mixed model"
      << endl;
  intro
      << "and a null model, each score separated with whitespace or line feed."
      << endl;
  intro
      << "Alternatively, accuate p-value could be provided in a single file pvalue_file."
      << endl << endl;
  intro << "The method reports score, PEP, and q-value calculated directly from scores" << endl; 
  intro << "Convergence of the method can be checked by comparing the average PEP of identifications" << endl;
  intro << "above threshold with the q-value" << endl; 
  CommandLineParser cmd(intro.str());
  // finally parse and handle return codes (display help etc...)
  cmd.defineOption("v",
                   "verbose",
                   "Set verbosity of output: 0=no processing info, 5=all, default is 2",
                   "level");
  cmd.defineOption("s",
                   "epsilon-step",
                   "The relative step size used as treshhold before cross validation error is calculated",
                   "value");
  cmd.defineOption("n",
                   "number-of-bins",
                   "The number of spline knots used when interpolating spline function. Default is 500.",
                   "bins");
  cmd.defineOption("c",
                   "epsilon-cross-validation",
                   "The relative crossvalidation step size used as treshhold before ending the iterations",
                   "value");
  cmd.defineOption("r",
                   "reverse",
                   "Indicating that the scoring mechanism is reversed, i.e., that low scores are better than higher scores",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("o",
                   "output-file",
                   "Output results to file instead of stdout",
                   "file");
  cmd.defineOption("Y",
                   "tdc-input",
                   "Turns off the pi0 correction for search results from a concatenated database.",
                   "",
                   TRUE_IF_SET);
  cmd.defineOption("d",
                   "include-negative",
                   "Include negative hits (decoy) probabilities in the results",
		    "",
		    TRUE_IF_SET);

  cmd.parseArgs(argc, argv);
  if (cmd.isOptionSet("verbose")) {
    Globals::getInstance()->setVerbose(cmd.getInt("verbose", 0, 10));
  }
  if (cmd.isOptionSet("number-of-bins")) {
    noIntervals = cmd.getInt("number-of-bins", 1, INT_MAX);
  }
  if (cmd.isOptionSet("epsilon-cross-validation")) {
    BaseSpline::convergeEpsilon = cmd.getDouble("epsilon-cross-validation", 0.0, 1.0);
  }
  if (cmd.isOptionSet("epsilon-step")) {
    BaseSpline::stepEpsilon = cmd.getDouble("epsilon-step", 0.0, 1.0);
  }
  if (cmd.isOptionSet("output-file")) {
    resultFileName = cmd.options["output-file"];
  }
  if (cmd.isOptionSet("reverse")) {
    PosteriorEstimator::setReversed(true);
  }
  if (cmd.isOptionSet("tdc-input")) {
    PosteriorEstimator::setUsePi0(false);
  }
  if (cmd.isOptionSet("include-negative")) {
    PosteriorEstimator::setNegative(true);
  }
  if (cmd.arguments.size() > 2) {
    cerr << "Too many arguments given" << endl;
    cmd.help();
  }
  if (cmd.arguments.size() == 0) {
    cerr << "No arguments given" << endl;
    cmd.help();
  }
  targetFile = cmd.arguments[0];
  if (cmd.arguments.size() == 2) {
    decoyFile = cmd.arguments[1];
  } else {
    PosteriorEstimator::setReversed(true);
    pvalInput = true;
  }
  return true;
}

