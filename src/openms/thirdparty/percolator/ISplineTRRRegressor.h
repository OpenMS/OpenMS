#pragma once

#include "MonotoneRegressor.h"
#include <iostream> // for std::cerr, std::endl
#include <Eigen/Dense>
#include <limits>
#include <numeric>
#include <cassert>

namespace OpenMS { namespace Internal { namespace Percolator {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Utility: project x to [lb, ub]
inline void projectToBox(VectorXd& x, const VectorXd& lb, const VectorXd& ub) {
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < lb[i]) x[i] = lb[i];
    else if (x[i] > ub[i]) x[i] = ub[i];
  }
}


// Quadratic objective model:
//   f(x) = x^T G x - 2 h^T x + c
//   g(x) = 2 (G x - h)
struct QuadraticModel {
  MatrixXd G;
  VectorXd h;
  double c;
};

struct CostGrad {
  double f;
  VectorXd g;
};

inline CostGrad costGrad(const QuadraticModel& q, const VectorXd& x) {
  const VectorXd Gx = q.G * x;
  CostGrad cg;
  cg.f = x.dot(Gx) - 2.0 * q.h.dot(x) + q.c;
  cg.g = 2.0 * (Gx - q.h);
  return cg;
}

// Identify the "free" set: variables not at a bound OR where the gradient points interior
inline std::vector<int> freeSet(const VectorXd& x, const VectorXd& g,
                                const VectorXd& lb, const VectorXd& ub, double tol = 1e-12) {
  std::vector<int> F;
  F.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    bool atL = std::abs(x[i] - lb[i]) <= tol, atU = std::abs(x[i] - ub[i]) <= tol;
    bool moveIn = (!atL && !atU) || (atL && g[i] < 0.0) || (atU && g[i] > 0.0);
    if (moveIn) F.push_back(i);
  }
  return F;
}

// Steihaug–Toint truncated CG: solve (H)p = -g on free vars, with ||p|| <= Delta
struct TRStep { VectorXd p; double pred_red; bool hit_boundary; };

inline TRStep steihaugCG(const MatrixXd& H,
                  const VectorXd& x,
                  const VectorXd& g,
                  const std::vector<int>& F,
                  double Delta) {
  // Work only on free coordinates via views.
  auto Hmul = [&](const VectorXd& vF) {
    // Build a full vector v with zeros outside F, apply H, then extract F.
    VectorXd v = VectorXd::Zero(x.size());
    for (int k = 0; k < (int)F.size(); ++k) v[F[k]] = vF[k];
    VectorXd Hv = H * v;
    VectorXd out(vF.size());
    for (int k = 0; k < (int)F.size(); ++k) out[k] = Hv[F[k]];
    return out;
  };
  // Gradient on free set: gF
  VectorXd gF(F.size());
  for (int k = 0; k < (int)F.size(); ++k) gF[k] = g[F[k]];

  // CG init
  VectorXd pF = VectorXd::Zero(F.size());
  VectorXd rCG = -gF;                  // residual
  VectorXd d = rCG;                    // search direction
  double rr = rCG.squaredNorm();
  bool hit_boundary = false;

  const double eps = 1e-12;
  const int max_cg_iters = static_cast<int>(F.size());
  for (int cg_iter = 0; cg_iter < max_cg_iters && rr > eps; ++cg_iter) {
    VectorXd Hd = Hmul(d);
    double dHd = d.dot(Hd);
    if (dHd <= 0) {
      // Non-positive curvature, step to boundary
      // Solve ||pF + tau d|| = Delta  for tau > 0
      double a = d.squaredNorm();
      double bq = 2.0 * pF.dot(d);
      double c = pF.squaredNorm() - Delta * Delta;
      double tau = (-bq + std::sqrt(std::max(0.0, bq*bq - 4*a*c))) / (2*a);
      pF = pF + tau * d;
      hit_boundary = true;
      break;
    }
    double alpha = rr / dHd;
    VectorXd pTrial = pF + alpha * d;
    if (pTrial.norm() >= Delta) {
      // Intersect with trust-region boundary
      double a = d.squaredNorm();
      double bq = 2.0 * pF.dot(d);
      double c = pF.squaredNorm() - Delta * Delta;
      double tau = (-bq + std::sqrt(std::max(0.0, bq*bq - 4*a*c))) / (2*a);
      pF = pF + tau * d;
      hit_boundary = true;
      break;
    }
    // Regular CG step
    pF = pTrial;
    VectorXd rCG_new = rCG - alpha * Hd;
    double rr_new = rCG_new.squaredNorm();
    double beta = rr_new / rr;
    d = rCG_new + beta * d;
    rCG = rCG_new;
    rr = rr_new;
  }

  // Predicted reduction: -g_F^T pF - 0.5 pF^T H pF
  VectorXd HpF = Hmul(pF);
  double pred = -gF.dot(pF) - 0.5 * pF.dot(HpF);

  // Lift pF back to full space
  VectorXd p = VectorXd::Zero(x.size());
  for (int k = 0; k < (int)F.size(); ++k) p[F[k]] = pF[k];
  return {p, pred, hit_boundary};
}

// Reflect/clamp to bounds along step p; return alpha in (0,1] that hits first bound
inline double stepToBox(const VectorXd& x, const VectorXd& p,
                        const VectorXd& lb, const VectorXd& ub) {
  double alpha = 1.0;
  for (int i = 0; i < x.size(); ++i) {
    if (p[i] > 0.0) alpha = std::min(alpha, (ub[i] - x[i]) / p[i]);
    else if (p[i] < 0.0) alpha = std::min(alpha, (lb[i] - x[i]) / p[i]);
  }
  return std::max(0.0, alpha);
}

// Top-level TRR solver
struct TRROpts {
  double Delta0 = 1.0;
  double DeltaMax = 1e6;
  double eta = 1e-3;     // accept if rho >= eta
  double gtol = 1e-8;    // projected gradient tolerance (inf-norm)
  int max_iters = 200;
};

struct TRRResult {
  VectorXd x;
  double fval;
  int iters;
  bool converged;
};

inline TRRResult trr_boxed_qp(const QuadraticModel& q,
                       const VectorXd& lb,
                       const VectorXd& ub,
                       VectorXd x0,
                       const TRROpts& opts = {}) {
  VectorXd x = x0; projectToBox(x, lb, ub);
  double Delta = opts.Delta0;
  const MatrixXd H = 2.0 * q.G;

  auto proj_grad_inf = [&](const VectorXd& x, const VectorXd& g) {
    double normInf = 0.0;
    for (int i = 0; i < x.size(); ++i) {
      double gi = g[i];
      if      (x[i] <= lb[i] + 1e-12) gi = std::min(0.0, gi); // cannot go below lb
      else if (x[i] >= ub[i] - 1e-12) gi = std::max(0.0, gi); // cannot go above ub
      normInf = std::max(normInf, std::abs(gi));
    }
    return normInf;
  };

  for (int it = 0; it < opts.max_iters; ++it) {
    CostGrad cg = costGrad(q, x);
    if (proj_grad_inf(x, cg.g) < opts.gtol)
      return {x, cg.f, it, true};

    // Free set and TR subproblem
    auto F = freeSet(x, cg.g, lb, ub);
    if (F.empty()) return {x, cg.f, it, true}; // stuck at KKT on bounds

    TRStep st = steihaugCG(H, x, cg.g, F, Delta);
    if (st.p.norm() == 0.0) return {x, cg.f, it, true};

    // Respect bounds by clamping along p
    double alphaBox = stepToBox(x, st.p, lb, ub);
    VectorXd x_trial = x + alphaBox * st.p;

    // Actual & predicted reduction
    double f_old = cg.f;
    double f_new = costGrad(q, x_trial).f;
    double ared = f_old - f_new;
    VectorXd s = alphaBox * st.p;
    double pred = -(cg.g.dot(s) + 0.5 * s.dot(H * s));
    double rho = (pred > 0.0) ? (ared / pred) : -std::numeric_limits<double>::infinity();

    // Update trust region
    if (rho < 0.25) Delta *= 0.25;
    else if (rho > 0.75 && std::abs(alphaBox - 1.0) < 1e-12 && st.hit_boundary)
      Delta = std::min(2.0 * Delta, opts.DeltaMax);

    // Accept/reject
    if (rho >= opts.eta) x = x_trial; // accept
    // else reject and keep x, with smaller Delta
  }

  CostGrad cg_end = costGrad(q, x);
  return {x, cg_end.f, opts.max_iters, false};
}

namespace ispline_detail {

inline std::vector<double> normalize_to_unit_interval(const std::vector<double>& x) {
  if (x.empty()) return {};

  auto mm = std::minmax_element(x.begin(), x.end());
  const double xmin = *mm.first;
  const double xmax = *mm.second;
  const double span = xmax - xmin;
  if (!(span > 0.0)) {
    return std::vector<double>(x.size(), 0.5);
  }

  std::vector<double> out(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    out[i] = (x[i] - xmin) / span;
  }
  return out;
}

} // namespace ispline_detail

// Build I-spline design matrix from proper B-spline basis.
//
// I-spline columns are cumulative tail sums of B-splines:
//   I_j(x) = sum_{m >= j} B_m(x),  j = 1, ..., n_bspline-1
// With clamped (repeated) boundary knots, each I_j is a smooth monotone
// function from 0 to 1 spanning multiple knot intervals.  Non-negative
// coefficients on these columns guarantee a monotone non-decreasing fit.
//
// B-splines are evaluated via the de Boor triangular recursion (O(d^2)
// per data point, independent of the number of basis functions).
//
// The knot vector must be a clamped B-spline knot vector: degree+1
// copies of the left boundary, sorted internal knots, degree+1 copies
// of the right boundary.
inline Eigen::MatrixXd build_ispline_design(
    const std::vector<double>& scores,
    int degree,
    const std::vector<double>& knots,
    bool include_intercept)
{
  const int n = static_cast<int>(scores.size());
  const int n_knots = static_cast<int>(knots.size());
  const int n_bspline = n_knots - degree - 1;
  const int d = degree;

  if (n_bspline <= 1) {
    Eigen::MatrixXd X(n, 1);
    X.setOnes();
    return X;
  }

  const int n_ispline = n_bspline - 1;
  const int p = (include_intercept ? 1 : 0) + n_ispline;
  Eigen::MatrixXd X = Eigen::MatrixXd::Zero(n, p);

  int col0 = 0;
  if (include_intercept) {
    X.col(0).setOnes();
    col0 = 1;
  }

  const double* t = knots.data();

  for (int i = 0; i < n; ++i) {
    double x = scores[i];

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
    static constexpr int kMaxDegree = 31;
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
    double cum[kMaxDegree + 1];
    cum[d] = B[d];
    for (int r = d - 1; r >= 0; --r) {
      cum[r] = cum[r + 1] + B[r];
    }

    // Fill I-spline columns.
    // I_j(x) = sum_{m >= j} B_m(x)  for j = 1 .. n_bspline-1
    //   = 1           if j <= k-d    (all non-zero B-splines included)
    //   = cum[j-(k-d)] if k-d < j <= k (partial sum of non-zero B-splines)
    //   = 0           if j > k        (no non-zero B-splines included)
    for (int j = 1; j < n_bspline; ++j) {
      double val;
      if (j <= k - d) {
        val = 1.0;
      } else if (j > k) {
        val = 0.0;
      } else {
        val = cum[j - (k - d)];
      }
      X(i, col0 + j - 1) = val;
    }
  }

  return X;
}


// Build a clamped B-spline knot vector with degree+1 repeated boundary
// knots and internal knots placed at data quantiles.
static std::vector<double> make_default_knots(const std::vector<double>& scores, int degree) {
  const int n = (int)scores.size();
  if (n == 0) return {};

  std::vector<double> sorted = scores;
  std::sort(sorted.begin(), sorted.end());

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
    const std::vector<double> x_norm = ispline_detail::normalize_to_unit_interval(*x_ptr);

    auto knots = params_.knots;
    if (knots.size() <= 0) {
      knots = make_default_knots(x_norm, params_.ispline_degree);
    }
    // Build X
    Eigen::MatrixXd X = build_ispline_design(
        x_norm, params_.ispline_degree, knots, params_.include_intercept);
    const int p = (int)X.cols();

    // Precompute normal-equation terms for objective:
    // (1/n)||X beta - y||^2 + lambda_ridge ||beta||^2 + lambda_smooth ||D beta_spline||^2
    // where D is the second-difference matrix on the I-spline coefficients.
    const Eigen::Map<const Eigen::VectorXd> y_vec(y.data(), n);
    const double inv_n = 1.0 / std::max(1, n);
    QuadraticModel q;
    q.G = inv_n * (X.transpose() * X);
    q.h = inv_n * (X.transpose() * y_vec);
    q.c = inv_n * y_vec.squaredNorm();
    const double lambda = std::max(0.0, params_.ridge_lambda);
    if (lambda > 0.0) {
      q.G.diagonal().array() += lambda;
    }

    // Second-difference smoothing penalty on I-spline coefficients.
    // Penalizes curvature (gamma_{j+2} - 2*gamma_{j+1} + gamma_j)^2,
    // encouraging smooth transitions instead of step functions.
    const int col0_s = params_.include_intercept ? 1 : 0;
    const int n_ispline = p - col0_s;
    const double smooth_lambda = std::max(0.0, params_.smooth_lambda);
    if (smooth_lambda > 0.0 && n_ispline > 2) {
      const int nd = n_ispline - 2;
      Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nd, p);
      for (int di = 0; di < nd; ++di) {
        D(di, col0_s + di)     =  1.0;
        D(di, col0_s + di + 1) = -2.0;
        D(di, col0_s + di + 2) =  1.0;
      }
      q.G += smooth_lambda * (D.transpose() * D);
    }

    // Bounds: non-negative except possibly intercept
    Eigen::VectorXd lb = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd ub = Eigen::VectorXd::Constant(p, std::numeric_limits<double>::infinity());
    if (params_.intercept_col >= 0 && params_.intercept_col < p)
      lb[params_.intercept_col] = -std::numeric_limits<double>::infinity();

    Eigen::VectorXd x0 = Eigen::VectorXd::Zero(p);
    TRROpts opts; opts.Delta0 = 1.0; opts.gtol = 1e-8; opts.max_iters = 1000;
    TRRResult sol = trr_boxed_qp(q, lb, ub, x0, opts);
    // Predict + clamp
    std::vector<double> yhat(n);
    for (int i = 0; i < n; ++i) {
      double v = X.row(i).dot(sol.x);
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
