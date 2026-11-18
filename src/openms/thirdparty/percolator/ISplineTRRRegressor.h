#ifndef ISPLINE_TRR_REGRESSOR_H_
#define ISPLINE_TRR_REGRESSOR_H_

#include "MonotoneRegressor.h"
#include "MyException.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

namespace OpenMS { namespace Internal { namespace Percolator {


// Small dense linear algebra for the p x p problems of the I-spline fit
// (p is about 200), on std::vector.
typedef std::vector<double> Vec;

struct DenseMatrix {
  DenseMatrix() : nRows(0), nCols(0) {}
  DenseMatrix(int rows, int cols)
      : nRows(rows), nCols(cols),
        values(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0.0) {}
  int rows() const { return nRows; }
  int cols() const { return nCols; }
  double& operator()(int i, int j) {
    return values[static_cast<size_t>(i) * nCols + j];
  }
  double operator()(int i, int j) const {
    return values[static_cast<size_t>(i) * nCols + j];
  }
  int nRows, nCols;
  std::vector<double> values;  // row-major
};

inline double dot(const Vec& a, const Vec& b) {
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

inline Vec multiply(const DenseMatrix& A, const Vec& x) {
  Vec y(A.rows(), 0.0);
  for (int i = 0; i < A.rows(); ++i) {
    const double* row = &A.values[static_cast<size_t>(i) * A.cols()];
    double s = 0.0;
    for (int j = 0; j < A.cols(); ++j) s += row[j] * x[j];
    y[i] = s;
  }
  return y;
}

// Utility: project x to [lb, ub]
inline void projectToBox(Vec& x, const Vec& lb, const Vec& ub) {
  for (size_t i = 0; i < x.size(); ++i) {
    if (x[i] < lb[i]) x[i] = lb[i];
    else if (x[i] > ub[i]) x[i] = ub[i];
  }
}


// Quadratic objective model:
//   f(x) = x^T G x - 2 h^T x + c
//   g(x) = 2 (G x - h)
struct QuadraticModel {
  DenseMatrix G;
  Vec h;
  double c;
};

// Solves x'Gx - 2h'x = min subject to lb <= x <= ub for a symmetric
// positive definite G with a primal active-set method. Each iteration
// minimises exactly over the free variables (Cholesky), moves towards that
// minimum until a bound blocks, or, when nothing blocks, releases the bound
// variable whose gradient points furthest into the box. For a strictly
// convex problem this ends at the exact optimum after a finite number of
// iterations. converged is false if the iteration limit (maxIters, or
// 10 p + 100 if that is 0) was reached first or G is not positive definite.
struct BoxQPResult {
  Vec x;
  double fval;
  int iters;
  bool converged;
};

// Solves G_FF y = rhs by Cholesky; returns false if G_FF is not positive
// definite
inline bool solveOnSubset(const DenseMatrix& G, const std::vector<int>& F,
                          const Vec& rhs, Vec& y) {
  const int m = static_cast<int>(F.size());
  DenseMatrix L(m, m);
  for (int j = 0; j < m; ++j) {
    double s = G(F[j], F[j]);
    for (int k = 0; k < j; ++k) s -= L(j, k) * L(j, k);
    if (!(s > 0.0)) return false;
    const double ljj = std::sqrt(s);
    L(j, j) = ljj;
    for (int i = j + 1; i < m; ++i) {
      double t = G(F[i], F[j]);
      for (int k = 0; k < j; ++k) t -= L(i, k) * L(j, k);
      L(i, j) = t / ljj;
    }
  }
  y.assign(m, 0.0);
  for (int i = 0; i < m; ++i) {
    double t = rhs[i];
    for (int k = 0; k < i; ++k) t -= L(i, k) * y[k];
    y[i] = t / L(i, i);
  }
  for (int i = m - 1; i >= 0; --i) {
    double t = y[i];
    for (int k = i + 1; k < m; ++k) t -= L(k, i) * y[k];
    y[i] = t / L(i, i);
  }
  return true;
}

inline BoxQPResult solveBoxQP(const QuadraticModel& q, const Vec& lb,
                              const Vec& ub, Vec x, int maxIterations = 0) {
  const int p = static_cast<int>(x.size());
  projectToBox(x, lb, ub);
  std::vector<char> atBound(p);
  for (int i = 0; i < p; ++i) atBound[i] = (x[i] == lb[i] || x[i] == ub[i]);
  const int maxIters = maxIterations > 0 ? maxIterations : 10 * p + 100;
  int it = 0;
  bool converged = false;
  for (; it < maxIters; ++it) {
    std::vector<int> F;
    for (int i = 0; i < p; ++i) {
      if (!atBound[i]) F.push_back(i);
    }
    // minimum over the free variables, with the bound variables fixed
    Vec rhs(F.size());
    for (size_t a = 0; a < F.size(); ++a) {
      double t = q.h[F[a]];
      for (int j = 0; j < p; ++j) {
        if (atBound[j]) t -= q.G(F[a], j) * x[j];
      }
      rhs[a] = t;
    }
    Vec xF;
    if (!solveOnSubset(q.G, F, rhs, xF)) break;
    // move towards it until the first bound blocks
    double alpha = 1.0;
    int blocking = -1;
    for (size_t a = 0; a < F.size(); ++a) {
      const int i = F[a];
      const double step = xF[a] - x[i];
      if (step < 0.0 && xF[a] < lb[i]) {
        const double r = (lb[i] - x[i]) / step;
        if (r < alpha) { alpha = r; blocking = i; }
      } else if (step > 0.0 && xF[a] > ub[i]) {
        const double r = (ub[i] - x[i]) / step;
        if (r < alpha) { alpha = r; blocking = i; }
      }
    }
    for (size_t a = 0; a < F.size(); ++a) {
      const int i = F[a];
      x[i] = (alpha == 1.0) ? xF[a] : x[i] + alpha * (xF[a] - x[i]);
    }
    if (blocking >= 0) {
      x[blocking] = (x[blocking] - lb[blocking] < ub[blocking] - x[blocking])
                        ? lb[blocking] : ub[blocking];
      atBound[blocking] = 1;
      continue;
    }
    // optimal on the free set; check the bound variables
    const Vec Gx = multiply(q.G, x);
    int release = -1;
    double mostInward = 1e-12;
    for (int i = 0; i < p; ++i) {
      if (!atBound[i]) continue;
      const double g = 2.0 * (Gx[i] - q.h[i]);
      const double inward = (x[i] == lb[i]) ? -g : g;
      if (inward > mostInward) {
        mostInward = inward;
        release = i;
      }
    }
    if (release < 0) {
      converged = true;
      break;
    }
    atBound[release] = 0;
  }
  const Vec Gx = multiply(q.G, x);
  const double fval = dot(x, Gx) - 2.0 * dot(q.h, x) + q.c;
  return {x, fval, it, converged};
}

namespace ispline_detail {

// negated: x holds the negated scores, for the error message
inline std::vector<double> normalize_to_unit_interval(const std::vector<double>& x,
                                                      bool negated = false) {
  if (x.empty()) return {};

  auto mm = std::minmax_element(x.begin(), x.end());
  const double xmin = *mm.first;
  const double xmax = *mm.second;
  const double span = xmax - xmin;
  // Scores near 1e308 and -1e308 have an infinite span, and (x - xmin) / span
  // is then inf / inf for the largest of them: every PEP came out as NaN.
  // The comparison is also false for a NaN.
  if (!(span <= std::numeric_limits<double>::max())) {
    std::ostringstream oss;
    oss << "ERROR: Cannot estimate PEPs by I-spline regression: the scores "
        << "range from " << (negated ? -xmax : xmin) << " to "
        << (negated ? -xmin : xmax) << ", and their span is not "
        << "a finite number. Check the input for scores that are too large "
        << "in magnitude." << std::endl;
    throw MyException(oss.str());
  }
  if (!(span > 0.0)) {
    return std::vector<double>(x.size(), 0.5);
  }

  std::vector<double> out(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    out[i] = (x[i] - xmin) / span;
  }
  return out;
}

static constexpr int kMaxDegree = 31;

// Evaluates the nonzero I-spline values at x.
//
// I-spline columns are cumulative tail sums of B-splines:
//   I_j(x) = sum_{m >= j} B_m(x),  j = 1, ..., n_bspline-1
// For x in knot span k:
//   I_j(x) = 1             if j <= k-d
//          = cum[j-(k-d)]  if k-d < j <= k
//          = 0             if j > k
// so a row of the design matrix is a run of ones, d partial values and
// zeros. Returns k and fills cum[0..d].
//
// B-splines are evaluated via the de Boor triangular recursion (O(d^2)
// per data point, independent of the number of basis functions).
//
// The knot vector must be a clamped B-spline knot vector: degree+1
// copies of the left boundary, sorted internal knots, degree+1 copies
// of the right boundary.
inline int ispline_span_values(double x, int d, const std::vector<double>& knots,
                               double* cum) {
  const int n_knots = static_cast<int>(knots.size());
  const double* t = knots.data();

  // Find knot span k: t[k] <= x < t[k+1], clamped to valid range [d, n_knots-d-2]
  int k;
  if (x >= t[n_knots - d - 1]) {
    k = n_knots - d - 2;
  } else if (x <= t[d]) {
    k = d;
  } else {
    int lo = d, hi = n_knots - d - 1;
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (t[mid + 1] <= x) lo = mid + 1;
      else hi = mid;
    }
    k = lo;
  }

  // De Boor triangular recursion: compute B_{k-d}(x), ..., B_k(x)
  assert(d <= kMaxDegree && "ispline_degree must be <= 31");
  double B[kMaxDegree + 1];
  B[0] = 1.0;
  for (int j = 1; j <= d; ++j) {
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      double left_val  = x - t[k + 1 - j + r];
      double right_val = t[k + 1 + r] - x;
      double denom = left_val + right_val;
      double temp = (denom > 0.0) ? B[r] / denom : 0.0;
      B[r] = saved + right_val * temp;
      saved = left_val * temp;
    }
    B[j] = saved;
  }
  // B[r] == B_{k-d+r, d}(x)  for r = 0, ..., d

  // Cumulative sum from right: cum[r] = sum_{s=r}^{d} B[s]
  cum[d] = B[d];
  for (int r = d - 1; r >= 0; --r) {
    cum[r] = cum[r + 1] + B[r];
  }
  return k;
}

} // namespace ispline_detail

// Build the dense I-spline design matrix (one row per score), with an
// optional intercept column first. The fit itself does not form this
// matrix; it is kept for inspection and tests.
inline DenseMatrix build_ispline_design(
    const std::vector<double>& scores,
    int degree,
    const std::vector<double>& knots,
    bool include_intercept)
{
  const int n = static_cast<int>(scores.size());
  const int n_bspline = static_cast<int>(knots.size()) - degree - 1;
  const int d = degree;

  if (n_bspline <= 1) {
    DenseMatrix X(n, 1);
    std::fill(X.values.begin(), X.values.end(), 1.0);
    return X;
  }

  const int col0 = include_intercept ? 1 : 0;
  DenseMatrix X(n, col0 + n_bspline - 1);
  double cum[ispline_detail::kMaxDegree + 1];
  for (int i = 0; i < n; ++i) {
    if (include_intercept) X(i, 0) = 1.0;
    const int k = ispline_detail::ispline_span_values(scores[i], d, knots, cum);
    for (int j = 1; j < n_bspline; ++j) {
      X(i, col0 + j - 1) = (j <= k - d) ? 1.0 : (j > k) ? 0.0 : cum[j - (k - d)];
    }
  }
  return X;
}


// Build a clamped B-spline knot vector with degree+1 repeated boundary
// knots and internal knots placed at data quantiles.
static std::vector<double> make_default_knots(const std::vector<double>& scores, int degree) {
  const int n = (int)scores.size();
  if (n == 0) return {};

  // the scores usually arrive sorted, then the copy and sort are not needed
  std::vector<double> sortedCopy;
  const std::vector<double>* sortedPtr = &scores;
  if (!std::is_sorted(scores.begin(), scores.end())) {
    sortedCopy = scores;
    std::sort(sortedCopy.begin(), sortedCopy.end());
    sortedPtr = &sortedCopy;
  }
  const std::vector<double>& sorted = *sortedPtr;

  double lo = sorted.front();
  double hi = sorted.back();
  if (hi <= lo) hi = lo + 1e-6;

  int num_internal = std::min(200, (int)std::sqrt(n));

  // Place internal knots at quantiles, skipping duplicates and boundaries
  std::vector<double> internal;
  internal.reserve(num_internal);
  double prev = lo;
  for (int k = 1; k <= num_internal; ++k) {
    double q = (double)k / (num_internal + 1);
    size_t idx = std::min((size_t)(q * (n - 1)), (size_t)(n - 1));
    double val = sorted[idx];
    if (val > lo + 1e-12 && val < hi - 1e-12 && val > prev + 1e-12) {
      internal.push_back(val);
      prev = val;
    }
  }

  // Clamped knot vector: (degree+1) copies of lo, internal knots, (degree+1) copies of hi
  const int order = degree + 1;
  std::vector<double> knots;
  knots.reserve(2 * order + (int)internal.size());
  for (int i = 0; i < order; ++i) knots.push_back(lo);
  for (double v : internal) knots.push_back(v);
  for (int i = 0; i < order; ++i) knots.push_back(hi);

  return knots;
}

class ISplineTRRRegressor final : public MonotoneRegressor {
public:
  using MonotoneRegressor::MonotoneRegressor;

  std::vector<double> fit_xy(const std::vector<double>& x,
                             const std::vector<double>& y,
                             double clip_lo = 0.0, double clip_hi = 1.0) override {
    assert(x.size() == y.size());
    const int n = (int)x.size();

    std::vector<double> x_work;
    const std::vector<double>* x_ptr = &x;
    if (params_.y_decreasing_in_x) {
      x_work = x;
      for (auto& v : x_work) v = -v;
      x_ptr = &x_work;
    }
    const std::vector<double> x_norm = ispline_detail::normalize_to_unit_interval(
        *x_ptr, params_.y_decreasing_in_x);
    x_work.clear();
    x_work.shrink_to_fit();

    auto knots = params_.knots;
    if (knots.size() <= 0) {
      knots = make_default_knots(x_norm, params_.ispline_degree);
    }
    const int d = params_.ispline_degree;
    const int n_bspline = (int)knots.size() - d - 1;
    const bool splines = n_bspline > 1;
    // columns: optional intercept, then the I-splines
    const int col0 = (splines && params_.include_intercept) ? 1 : 0;
    const int p = splines ? col0 + n_bspline - 1 : 1;

    // Normal equation terms of the objective
    // (1/n)||X beta - y||^2 + lambda_ridge ||beta||^2 + lambda_smooth ||D beta_spline||^2
    // where D is the second-difference matrix on the I-spline coefficients.
    // The design matrix X is not formed: a row in knot span k is a run of
    // ones followed by d partial values, so X'X and X'y are accumulated per
    // span and assembled afterwards.
    QuadraticModel q;
    q.G = DenseMatrix(p, p);
    q.h.assign(p, 0.0);
    double yy = 0.0;
    double cum[ispline_detail::kMaxDegree + 1];
    if (!splines) {
      // X is a single column of ones
      double ySum = 0.0;
      for (int i = 0; i < n; ++i) {
        ySum += y[i];
        yy += y[i] * y[i];
      }
      q.G(0, 0) = n;
      q.h[0] = ySum;
    } else {
      // per span k: count, sum of y, sums of cum[t], cum[t]*cum[u] and
      // y*cum[t] for t, u = 1..d
      const int numSpans = (int)knots.size();
      const int dd = d + 1;
      std::vector<double> count(numSpans, 0.0), ySum(numSpans, 0.0);
      std::vector<double> cumSum(static_cast<size_t>(numSpans) * dd, 0.0);
      std::vector<double> yCumSum(static_cast<size_t>(numSpans) * dd, 0.0);
      std::vector<double> cumCumSum(static_cast<size_t>(numSpans) * dd * dd, 0.0);
      for (int i = 0; i < n; ++i) {
        const int k = ispline_detail::ispline_span_values(x_norm[i], d, knots, cum);
        count[k] += 1.0;
        ySum[k] += y[i];
        yy += y[i] * y[i];
        double* cs = &cumSum[static_cast<size_t>(k) * dd];
        double* ycs = &yCumSum[static_cast<size_t>(k) * dd];
        double* ccs = &cumCumSum[static_cast<size_t>(k) * dd * dd];
        for (int t = 1; t <= d; ++t) {
          cs[t] += cum[t];
          ycs[t] += y[i] * cum[t];
          for (int u = t; u <= d; ++u) ccs[t * dd + u] += cum[t] * cum[u];
        }
      }
      for (int k = d; k < numSpans; ++k) {
        if (count[k] == 0.0) continue;
        const int a = k - d;  // number of leading I-spline columns equal to 1
        // columns that are 1: the intercept and I-splines 1..a
        std::vector<int> ones;
        if (col0) ones.push_back(0);
        for (int j = 1; j <= a; ++j) ones.push_back(col0 + j - 1);
        const double* cs = &cumSum[static_cast<size_t>(k) * dd];
        const double* ycs = &yCumSum[static_cast<size_t>(k) * dd];
        const double* ccs = &cumCumSum[static_cast<size_t>(k) * dd * dd];
        for (int r : ones) {
          q.h[r] += ySum[k];
          for (int c : ones) q.G(r, c) += count[k];
          for (int t = 1; t <= d; ++t) {
            const int c = col0 + a + t - 1;
            q.G(r, c) += cs[t];
            q.G(c, r) += cs[t];
          }
        }
        for (int t = 1; t <= d; ++t) {
          const int r = col0 + a + t - 1;
          q.h[r] += ycs[t];
          for (int u = t; u <= d; ++u) {
            const int c = col0 + a + u - 1;
            q.G(r, c) += ccs[t * dd + u];
            if (c != r) q.G(c, r) += ccs[t * dd + u];
          }
        }
      }
    }
    const double inv_n = 1.0 / std::max(1, n);
    for (double& v : q.G.values) v *= inv_n;
    for (double& v : q.h) v *= inv_n;
    q.c = inv_n * yy;
    const double lambda = std::max(0.0, params_.ridge_lambda);
    if (lambda > 0.0) {
      for (int j = 0; j < p; ++j) q.G(j, j) += lambda;
    }

    // Second-difference smoothing penalty on I-spline coefficients.
    // Penalizes curvature (gamma_{j+2} - 2*gamma_{j+1} + gamma_j)^2,
    // encouraging smooth transitions instead of step functions.
    const int col0_s = params_.include_intercept ? 1 : 0;
    const int n_ispline = p - col0_s;
    const double smooth_lambda = std::max(0.0, params_.smooth_lambda);
    if (smooth_lambda > 0.0 && n_ispline > 2) {
      // add smooth_lambda * D'D, with rows (1, -2, 1) of D
      const double coef[3] = {1.0, -2.0, 1.0};
      for (int di = 0; di < n_ispline - 2; ++di) {
        for (int r = 0; r < 3; ++r) {
          for (int c = 0; c < 3; ++c) {
            q.G(col0_s + di + r, col0_s + di + c) += smooth_lambda * coef[r] * coef[c];
          }
        }
      }
    }

    // Bounds: non-negative except possibly intercept
    Vec lb(p, 0.0);
    Vec ub(p, std::numeric_limits<double>::infinity());
    if (params_.intercept_col >= 0 && params_.intercept_col < p)
      lb[params_.intercept_col] = -std::numeric_limits<double>::infinity();

    BoxQPResult sol =
        solveBoxQP(q, lb, ub, Vec(p, 0.0), params_.max_solver_iters);
    // The solver stops without the optimum at its iteration limit, or when
    // the normal equations are not positive definite. The coefficients are
    // within their bounds, so the fit is monotone, but it is not the
    // least-squares fit.
    if (!sol.converged) {
      std::cerr << "Warning: the I-spline regression of the PEPs did not "
                << "converge (" << sol.iters << " iterations, " << p
                << " coefficients); the PEPs may be inaccurate." << std::endl;
    }

    // Predict + clamp; prefix[a] is the sum of the first a I-spline
    // coefficients, so a row costs O(d)
    std::vector<double> yhat(n);
    if (!splines) {
      for (int i = 0; i < n; ++i) yhat[i] = std::min(clip_hi, std::max(clip_lo, sol.x[0]));
      return yhat;
    }
    Vec prefix(n_bspline, 0.0);
    for (int j = 1; j < n_bspline; ++j) prefix[j] = prefix[j - 1] + sol.x[col0 + j - 1];
    const double intercept = col0 ? sol.x[0] : 0.0;
    for (int i = 0; i < n; ++i) {
      const int k = ispline_detail::ispline_span_values(x_norm[i], d, knots, cum);
      const int a = k - d;
      double v = intercept + prefix[a];
      for (int t = 1; t <= d; ++t) v += cum[t] * sol.x[col0 + a + t - 1];
      if (v < clip_lo) v = clip_lo;
      if (v > clip_hi) v = clip_hi;
      yhat[i] = v;
    }
    return yhat;
  }

  std::vector<double> fit_y(const std::vector<double>& y,
                            double clip_lo = 0.0, double clip_hi = 1.0) override {
    // No x provided: use rank as a x-value.
    std::vector<double> x(y.size());
    std::iota(x.rbegin(), x.rend(), 1.0);  // x[0]=N, ..., x[N-1]=1
    return fit_xy(x, y, clip_lo, clip_hi);
  }

};

}}}  // namespace OpenMS::Internal::Percolator

#endif /* ISPLINE_TRR_REGRESSOR_H_ */
