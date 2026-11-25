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
#include<iterator>
#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<fstream>
#include<sstream>
#include<cmath>
#include<iostream>

#include "Globals.h"
#include "LogisticRegression.h"

namespace OpenMS { namespace Internal { namespace Percolator {

using namespace std;


const double LogisticRegression::gRange = 35.0;

void invlogit(double& out, double in) {
  double e = exp(in);
  out = e / (1 + e);
}

double invlogit(double in) {
  double e = exp(in);
  return e / (1 + e);
}

double logit(double p) {
  return log(p / (1 - p));
}

// Distinct scores can give equal knots after the log or logit transform of
// BaseSpline::setData (two scores near 1e9 one ulp apart). The interval
// between equal knots has length 0, and reinschStep cannot take a step.
// Such knots are merged. The spline needs 4 knots: if fewer would remain,
// nothing is merged, the step fails, and roughnessPenaltyIRLS prints its
// warning. Knots that are only close together are not merged.
void LogisticRegression::mergeEqualKnots(std::vector<double>& xx,
                                         std::vector<double>& yy,
                                         std::vector<double>& mm) {
  const std::size_t n = xx.size();
  std::size_t distinct = n > 0 ? 1 : 0;
  for (std::size_t i = 1; i < n; ++i) {
    distinct += xx[i] != xx[i - 1];
  }
  if (distinct == n || distinct < 4) {
    return;
  }
  std::size_t k = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (xx[i] == xx[k]) {
      yy[k] += yy[i];
      mm[k] += mm[i];
    } else {
      ++k;
      xx[k] = xx[i];
      yy[k] = yy[i];
      mm[k] = mm[i];
    }
  }
  xx.resize(distinct);
  yy.resize(distinct);
  mm.resize(distinct);
  if (VERB > 2) {
    cerr << "Merged " << n - distinct << " of " << n << " bins for the PEP "
         << "fit: their knots are equal." << endl;
  }
}

// Knots that are not equal but much closer together than the knots around
// them leave the system of each step ill-conditioned also in double-double.
// The fit can then end at its limit with the warning, or converge, and the
// PEPs can depend on the last bits of exp and log or be wrong (24 qvality
// scores 1 ulp apart at 1.5 and 400 others: best PEP 6.3e-16, against 3.8e-4
// with the 24 merged; protein scores 0 and 1e-20 among 300: PEPs up to 0.85
// off those with the two scores 3e-6 apart). A rule that merged such
// knots would have to tell the knots that the fit separates from those that
// it cannot, and merging changes the PEPs. They are not merged; a warning is
// printed instead. The smallest ratio of a gap to the median of the gaps
// around it was 4.0e-6 on 150 synthetic sets and 3.5e-7 in qvality inputs
// with two scores 1e-8 apart; at 7.3e-11 (two of 300 protein scores 1e-12
// apart) the PEPs differed by up to 5.4e-6 with and without the FMA variants
// of exp and log. The larger middle gap is taken so that a run of nearly
// equal knots counts at its ends; with the smaller one, runs of 7 or more
// knots 1 ulp apart would not count. Where all knots form runs of 3 or more
// knots 1 ulp apart, the gaps around each pair are small too, and few or no
// pairs count (none on 18 qvality inputs with runs of 4 to 8 or 12, where
// the fit did not converge and printed its warning).
static const double kNearTie = 1e-9;
static const std::size_t kNearTieWindow = 5;

std::size_t LogisticRegression::countNearTies(const std::vector<double>& xx) {
  const std::size_t n = xx.size();
  std::size_t count = 0;
  std::vector<double> gaps;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    gaps.clear();
    for (std::size_t j = i - min(i, kNearTieWindow);
         j < min(n - 1, i + 1 + kNearTieWindow); ++j) {
      const double gap = xx[j + 1] - xx[j];
      if (j != i && !std::isnan(gap)) {
        gaps.push_back(gap);
      }
    }
    if (gaps.empty()) {
      continue;
    }
    const std::size_t mid = gaps.size() / 2;
    std::nth_element(gaps.begin(), gaps.begin() + mid, gaps.end());
    count += xx[i + 1] - xx[i] < kNearTie * gaps[mid];
  }
  return count;
}

void LogisticRegression::setData(const std::vector<double>& xx,
                                 const std::vector<double>& yy,
                                 const std::vector<double>& mm) {
  BaseSpline::setData(xx);
  y = yy;
  m = mm;
  mergeEqualKnots(x, y, m);
  const std::size_t nearTies = countNearTies(x);
  if (nearTies > 0 && VERB > 0) {
    cerr << "Warning: " << nearTies << " pairs of scores in the PEP fit are "
         << "nearly equal; the PEPs may be inaccurate." << endl;
  }
}

void LogisticRegression::limitg() {
  for (int ix = gnew.numberEntries(); ix--;) {
    gnew.packedReplace(ix, min(gRange, max(-gRange, gnew[ix])));
    assert(std::isfinite(gnew[ix]));
  }
}

void LogisticRegression::calcPZW() {

  for (std::size_t ix = static_cast<std::size_t>(z.numberEntries()); ix--;) {
    assert(std::isfinite(g[ix]));
    double e = exp(g[ix]);
    assert(std::isfinite(e));
    double epsilon = 1e-15;
    p[ix] = min(max(e / (1 + e), epsilon), 1
        - epsilon);
    assert(std::isfinite(p[ix]));
    w.packedReplace(static_cast<int>(ix), max(m[ix] * p[ix] * (1 - p[ix]), epsilon));
    assert(std::isfinite(w[ix]));
    z.packedReplace(static_cast<int>(ix), min(gRange, max(-gRange, g[ix] +
        (y[ix] - p[ix] * m[ix]) / w[ix])));
    assert(std::isfinite(z[ix]));
  }
}

void LogisticRegression::initg() {
  BaseSpline::initg();
  int n = static_cast<int>(x.size());
  p.resize(n);
  gnew = PackedVector(n);
  for (std::size_t ix = static_cast<std::size_t>(g.size()); ix--;) {
    double p = (y[ix] + 0.05) / (m[ix] + 0.1);
    gnew.packedReplace(static_cast<int>(ix), log(p / (1 - p)));
    assert(std::isfinite(p));
    assert(std::isfinite(g[ix]));
  }
#define OUTPUT_DEBUG_FILES
#undef OUTPUT_DEBUG_FILES
#ifdef OUTPUT_DEBUG_FILES
  ofstream drFile("decoyRate.bins", ios::out), xvalFile("xvals.bins", ios::out);
  ostream_iterator<double> xvalIt(xvalFile, "\n");
  copy(x.begin(), x.end(), xvalIt);
  for (size_t yix = 0; yix < y.size(); ++yix) {
    drFile << y[yix] / m[yix] << endl;
  }
  drFile.close();
#endif
}

}}}  // namespace OpenMS::Internal::Percolator
