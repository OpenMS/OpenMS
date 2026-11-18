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
#include<cmath>
#include "BaseSpline.h"
#include "Globals.h"

namespace OpenMS { namespace Internal { namespace Percolator {


using namespace std;

class SplinePredictor {
    BaseSpline* bs;
  public:
    SplinePredictor(BaseSpline* b) {
      bs = b;
    }
    double operator()(double x) {
      return bs->predict(x);
    }
};

double BaseSpline::convergeEpsilon = 1e-4;
double BaseSpline::stepEpsilon = 1e-8;
double BaseSpline::weightSlope = 1e1;
double BaseSpline::scaleAlpha = 1;

double BaseSpline::splineEval(double xx) {
  xx = transf(xx);
  std::size_t n = x.size();
  vector<double>::iterator left, right = lower_bound(x.begin(),
                                                     x.end(),
                                                     xx);
  if (right == x.end()) {
    double derl = (g[n - 1] - g[n - 2]) / (x[n - 1] - x[n - 2])
        + (x[n - 1] - x[n - 2]) / 6 * gamma[n - 3];
    double gx = g[n - 1] + (xx - x[n - 1]) * derl;
    return gx;
  }
  int rix = static_cast<int>(right - x.begin());
  if (*right == xx) {
    return g[rix];
  }
  if (rix > 0) {
    left = right;
    left--;
    double dr = *right - xx;
    double dl = xx - *left;
    double gamr = (rix < static_cast<int>(n - 1) ? gamma[rix - 1] : 0.0);
    double gaml = (rix > 1 ? gamma[rix - 1 - 1] : 0.0);
    double h = *right - *left;
    double gx = (dl * g[rix] + dr * g[rix - 1]) / h - dl * dr / 6 * ((1.0
        + dl / h) * gamr + (1.0 + dr / h) * gaml);
    return gx;
  }
  // if (rix==0)
  double derr = (g[1] - g[0]) / (x[1] - x[0]) - (x[1] - x[0]) / 6
      * gamma[0];
  double gx = g[0] - (x[0] - xx) * derr;
  return gx;
}

static double tao = 2 / (1 + sqrt(5.0)); // inverse of golden section

namespace {

// Double-double arithmetic: a value is the unevaluated sum hi + lo with
// |lo| <= ulp(hi) / 2, about 32 significant digits (Dekker 1971; Hida, Li
// and Bailey 2001). The error-free transforms assume IEEE 754 double
// arithmetic rounded to nearest, without extended precision and without
// contraction of a * b + c into a fused multiply-add: this file is compiled
// with -ffp-contract=off, and must not be compiled with -ffast-math.
// std::fma gives the exact error of a product on every platform.
struct DD {
  double hi, lo;
  DD(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};

// s + e == a + b exactly, if |a| >= |b|
inline DD fastTwoSum(double a, double b) {
  double s = a + b;
  return DD(s, b - (s - a));
}

// s + e == a + b exactly
inline DD twoSum(double a, double b) {
  double s = a + b;
  double bb = s - a;
  return DD(s, (a - (s - bb)) + (b - bb));
}

inline DD operator+(const DD& a, const DD& b) {
  DD s = twoSum(a.hi, b.hi);
  DD t = twoSum(a.lo, b.lo);
  s.lo += t.hi;
  s = fastTwoSum(s.hi, s.lo);
  s.lo += t.lo;
  return fastTwoSum(s.hi, s.lo);
}

inline DD operator-(const DD& a) { return DD(-a.hi, -a.lo); }

inline DD operator-(const DD& a, const DD& b) { return a + (-b); }

inline DD operator*(const DD& a, const DD& b) {
  double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  return fastTwoSum(p, e);
}

inline DD operator/(const DD& a, const DD& b) {
  double q1 = a.hi / b.hi;
  DD r = a - DD(q1) * b;
  double q2 = r.hi / b.hi;
  r = r - DD(q2) * b;
  double q3 = r.hi / b.hi;
  return fastTwoSum(q1, q2) + DD(q3);
}

}  // namespace

// The step is computed in double-double. Two knots close together (1e-7 to
// 1e-5 apart on a score range of a few units, as when each entry is its own
// bin) give M = R + Q' alpha W^-1 Q a condition number of 1e10 or more,
// also after diagonal scaling. In double precision the error of the solve
// (median 2e-4 in g on a synthetic set, 7e-2 in the worst one) was then
// larger than the tolerance on the step: the iteration ran to the cap of 20
// steps, and g depended on the last bits of exp and log. In double-double
// the error stays below the rounding of g, and Newton's method converges in
// a few steps.
bool BaseSpline::reinschStep(double alpha) {
  const int n = static_cast<int>(x.size());
  const int k = n - 2;
  // the intervals h_j = x[j + 1] - x[j], exact in double-double; equal
  // knots give no step
  std::vector<DD> h(n - 1), inv(n - 1);
  for (int j = 0; j < n - 1; ++j) {
    h[j] = twoSum(x[j + 1], -x[j]);
    if (!(h[j].hi > 0.0)) {
      return false;
    }
    inv[j] = DD(1.0) / h[j];
    if (!std::isfinite(inv[j].hi)) {
      return false;
    }
  }
  std::vector<DD> aw(n);  // alpha / w_i
  for (int i = 0; i < n; ++i) {
    aw[i] = DD(alpha) / DD(w[i]);
  }
  // Column j of Q holds 1/h_j, -1/h_j - 1/h_(j+1) and 1/h_(j+1) in the
  // rows j, j + 1 and j + 2. M is symmetric and pentadiagonal: m0 holds its
  // diagonal, m1 and m2 the first and second superdiagonals. gam starts as
  // Q' z.
  std::vector<DD> m0(k), m1(k), m2(k), gam(k);
  for (int j = 0; j < k; ++j) {
    const DD& a = inv[j];
    const DD& c = inv[j + 1];
    DD b = -(a + c);
    m0[j] = (h[j] + h[j + 1]) / DD(3.0) + aw[j] * a * a + aw[j + 1] * b * b
        + aw[j + 2] * c * c;
    if (j + 1 < k) {
      DD bNext = -(inv[j + 1] + inv[j + 2]);
      m1[j] = h[j + 1] / DD(6.0) + aw[j + 1] * b * inv[j + 1]
          + aw[j + 2] * c * bNext;
    }
    if (j + 2 < k) {
      m2[j] = aw[j + 2] * c * inv[j + 2];
    }
    gam[j] = a * DD(z[j]) + b * DD(z[j + 1]) + c * DD(z[j + 2]);
  }
  // M = L D L' with unit lower triangular L of bandwidth 2, without
  // pivoting, as M is positive definite; l1[i] = L(i, i-1), l2[i] = L(i, i-2)
  std::vector<DD> d(k), l1(k), l2(k);
  for (int i = 0; i < k; ++i) {
    DD di = m0[i];
    if (i >= 2) {
      l2[i] = m2[i - 2] / d[i - 2];
      di = di - l2[i] * l2[i] * d[i - 2];
    }
    if (i >= 1) {
      DD s = m1[i - 1];
      if (i >= 2) {
        s = s - l2[i] * d[i - 2] * l1[i - 1];
      }
      l1[i] = s / d[i - 1];
      di = di - l1[i] * l1[i] * d[i - 1];
    }
    d[i] = di;
  }
  for (int i = 1; i < k; ++i) {
    gam[i] = gam[i] - l1[i] * gam[i - 1];
    if (i >= 2) {
      gam[i] = gam[i] - l2[i] * gam[i - 2];
    }
  }
  for (int i = 0; i < k; ++i) {
    gam[i] = gam[i] / d[i];
  }
  for (int i = k - 2; i >= 0; --i) {
    gam[i] = gam[i] - l1[i + 1] * gam[i + 1];
    if (i + 2 < k) {
      gam[i] = gam[i] - l2[i + 2] * gam[i + 2];
    }
  }
  // gnew = z - alpha W^-1 Q gamma; row i of Q has entries in the columns
  // i - 2, i - 1 and i
  std::vector<double> gammaNew(k), gNew(n);
  for (int j = 0; j < k; ++j) {
    gammaNew[j] = gam[j].hi;
    if (!std::isfinite(gammaNew[j])) {
      return false;
    }
  }
  for (int i = 0; i < n; ++i) {
    DD qg;
    if (i < k) {
      qg = qg + inv[i] * gam[i];
    }
    if (i >= 1 && i - 1 < k) {
      qg = qg - (inv[i - 1] + inv[i]) * gam[i - 1];
    }
    if (i >= 2) {
      qg = qg + inv[i - 1] * gam[i - 2];
    }
    gNew[i] = (DD(z[i]) - aw[i] * qg).hi;
    if (!std::isfinite(gNew[i])) {
      return false;
    }
  }
  gamma = PackedVector(k);
  for (int j = 0; j < k; ++j) {
    gamma.packedReplace(j, gammaNew[j]);
  }
  gnew = PackedVector(n);
  for (int i = 0; i < n; ++i) {
    gnew.packedReplace(i, gNew[i]);
  }
  return true;
}

void BaseSpline::roughnessPenaltyIRLS() {
  // at least 4 points are needed for the cubic smoothing spline
  assert(x.size() >= 4);
  initg();
  double p1 = 1 - tao;
  double p2 = tao;
  double alpha = alphaLinearSearchBA(0.0,
                          1.0,
                          p1,
                          p2,
                          evaluateSlope(-scaleAlpha*log(p1)),
                          evaluateSlope(-scaleAlpha*log(p2)));
  if (VERB > 2) {
    cerr << "Alpha selected to be " << alpha << endl;
  }
  if (!iterativeReweightedLeastSquares(alpha) && VERB > 0) {
    cerr << "Warning: the P-IRLS fit of the PEPs did not converge ("
         << x.size() << " bins); the PEPs may be inaccurate." << endl;
  }
}

bool BaseSpline::iterativeReweightedLeastSquares(double alpha) {
  double step = 0.0;
  int iter = 0;
  unsigned int n = static_cast<unsigned int>(x.size());
  do {
    g = gnew;
    calcPZW();
    if (!reinschStep(alpha)) {
      // no step: keep g (limitg would turn a NaN into -35)
      gnew = g;
      return false;
    }
    limitg();
    PackedVector difference = g.packedSubtract(gnew);
    step = packedNorm(difference) / n;
    if (VERB > 3) {
      cerr << "step size:" << step << endl;
    }
  } while ((step > stepEpsilon || step < 0.0) && (++iter < 20));
  g = gnew;
  // false also for a step that is not a number
  return step >= 0.0 && step <= stepEpsilon;
}

double BaseSpline::alphaLinearSearchBA(double min_p,
                                       double max_p,
                                       double p1, double p2,
                                       double cv1, double cv2) {
  // Minimize Slope score
  // Use neg log of 0<p<1 so that we allow for searches 0<alpha<inf
  double oldCV = 0.0;
  if (cv2 < cv1) {
    // keep point 2
    min_p = p1;
    p1 = p2;
    p2 = min_p + tao * (max_p - min_p);
    oldCV = cv1;
    cv1 = cv2;
    cv2 = evaluateSlope(-scaleAlpha*log(p2));
    if (VERB > 3) {
      cerr << "New point with alpha=" << -scaleAlpha*log(p2) << ", giving slopeScore=" << cv2 << endl;
    }
  } else {
    // keep point 1
    max_p = p2;
    p2 = p1;
    p1 = min_p + (1 - tao) * (max_p - min_p);
    oldCV = cv2;
    cv2 = cv1;
    cv1 = evaluateSlope(-scaleAlpha*log(p1));
    if (VERB > 3) {
      cerr << "New point with alpha=" << -scaleAlpha*log(p1) << ", giving slopeScore=" << cv1 << endl;
    }
  }
  if ((oldCV - min(cv1, cv2)) / oldCV < 1e-5 || (abs(p2 - p1) < 1e-10)) {
    return (cv1 < cv2 ? -scaleAlpha*log(p1) : -scaleAlpha*log(p2));
  }
  return alphaLinearSearchBA(min_p, max_p, p1, p2, cv1, cv2);
}

double BaseSpline::evaluateSlope(double alpha) {
  // Calculate a spline for current alpha
  iterativeReweightedLeastSquares(alpha);
  // Find highest point (we only want to evaluate things to the right of that point)
  int n = g.numberEntries();
  int mixg=1; // Ignore 0 and n-1
  double maxg = g[mixg];
  for (int ix=mixg;ix<n-1;++ix) {
    assert(ix==g.index(ix)); //This should be a filled vector
    if (g[ix]>=maxg) {
      maxg = g[ix];
      mixg = ix;
    }
  }
  double maxSlope = -10e6;
  int slopeix=-1;
  for (int ix=mixg+1;ix<n-2; ++ix) {
    double slope=g[ix-1]-g[ix];
    if (slope>maxSlope) {
      maxSlope = slope;
      slopeix = ix;
    }
  }
  // Now score the fit based on a linear combination between
  // The bump area and alpha

  if (VERB > 3) {
    cerr << "mixg=" << mixg << ", maxg=" << maxg << ", maxBA=" << maxSlope << " at ix=" << slopeix << ", alpha=" << alpha << endl;
  }
  return maxSlope*weightSlope + alpha;
}


void BaseSpline::predict(const vector<double>& xx, vector<double>& predict) {
  predict.clear();
  transform(xx.begin(),
            xx.end(),
            back_inserter(predict),
            SplinePredictor(this));
}

void BaseSpline::setData(const vector<double>& xx) {
  x.clear();
  double minV = *min_element(xx.begin(), xx.end());
  double maxV = *max_element(xx.begin(), xx.end());
  if (minV >= 0.0 && maxV <= 1.0) {
    if (VERB > 1) {
      cerr << "Logit transforming all scores prior to PEP calculation"
        << endl;
    }
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20,
                       maxV < 1.0 ? 0.0 : 1e-10,
                       true);
  } else if (minV >= 0.0) {
    if (VERB > 1) {
      cerr << "Log transforming all scores prior to PEP calculation"
          << endl;
    }
    transf = Transform(minV > 0.0 ? 0.0 : 1e-20, 0.0, false, true);
  }
  transform(xx.begin(), xx.end(), back_inserter(x), transf);
}

}}}  // namespace OpenMS::Internal::Percolator
