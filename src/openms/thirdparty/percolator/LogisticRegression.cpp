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
#include<iterator>
#include<vector>
#include<algorithm>
#include<numeric>
#include<functional>
#include<fstream>
#include<sstream>

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

void LogisticRegression::limitg() {
  for (int ix = gnew.numberEntries(); ix--;) {
    gnew.packedReplace(ix, min(gRange, max(-gRange, gnew[ix])));
    assert(isfinite(gnew[ix]));
  }
}

void LogisticRegression::limitgamma() {
  for (int ix = gamma.numberEntries(); ix--;) {
    gamma.packedReplace(ix, min(gRange, max(-gRange, gamma[ix])));
    assert(isfinite(gamma[ix]));
  }
}

void LogisticRegression::calcPZW() {

  for (std::size_t ix = static_cast<std::size_t>(z.numberEntries()); ix--;) {
    assert(isfinite(g[ix]));
    double e = exp(g[ix]);
    assert(isfinite(e));
    double epsilon = 1e-15;
    p[ix] = min(max(e / (1 + e), epsilon), 1
        - epsilon);
    assert(isfinite(p[ix]));
    w.packedReplace(static_cast<int>(ix), max(m[ix] * p[ix] * (1 - p[ix]), epsilon));
    assert(isfinite(w[ix]));
    z.packedReplace(static_cast<int>(ix), min(gRange, max(-gRange, g[ix] +
        (y[ix] - p[ix] * m[ix]) / w[ix])));
    assert(isfinite(z[ix]));
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
    assert(isfinite(p));
    assert(isfinite(g[ix]));
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


// OpenMS: the spline needs at least 4 strictly increasing x (initiateQR), and
// keeps their distances in sparse vectors that treat |dx| <= 1e-15 as 0. Merge
// the points it cannot tell apart and pad to 4 points spaced like the data.
// Both happen in the transformed coordinates the spline fits in, which
// setData chose from the real data, so the padding cannot change the
// transform. Expects x in ascending order.
void LogisticRegression::mergeTiesAndPad() {
  std::size_t n = 0;
  for (std::size_t ix = 0; ix < x.size(); ++ix) {
    if (n > 0 && x[ix] - x[n - 1] <= 1e-12 * max(1.0, fabs(x[ix]))) {
      y[n - 1] += y[ix];
      m[n - 1] += m[ix];
    } else {
      x[n] = x[ix];
      y[n] = y[ix];
      m[n] = m[ix];
      ++n;
    }
  }
  x.resize(n);
  y.resize(n);
  m.resize(n);
  if (n == 0) {
    return;
  }
  double offset = n > 1 ? x.back() - x.front() : max(1.0, fabs(x.back()));
  for (std::size_t ix = n; ix < 4; ++ix) {
    x.push_back(x.back() + static_cast<double>(ix) * offset);
    y.push_back(y.back());
    m.push_back(m.back());
  }
}

}}}  // namespace OpenMS::Internal::Percolator
