// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
#include <OpenMS/MATH/MISC/SmoothingSpline.h>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace OpenMS
{
namespace Math
{
// small helper: solve linear system A x = b using Gaussian elimination with partial pivoting
static bool solve_linear_system(std::vector<double> A, std::vector<double> b, std::vector<double>& x, int n)
{
  const double EPS = 1e-18;
  x.assign(n, 0.0);
  // augment and perform elimination in-place on A (row-major)
  for (int i = 0; i < n; ++i)
  {
    // pivot
    int piv = i;
    double maxv = std::fabs(A[i * n + i]);
    for (int r = i + 1; r < n; ++r)
    {
      double v = std::fabs(A[r * n + i]);
      if (v > maxv) { maxv = v; piv = r; }
    }
    if (maxv < EPS) return false; // singular
    if (piv != i)
    {
      for (int c = i; c < n; ++c) std::swap(A[i * n + c], A[piv * n + c]);
      std::swap(b[i], b[piv]);
    }
    double diag = A[i * n + i];
    for (int r = i + 1; r < n; ++r)
    {
      double fac = A[r * n + i] / diag;
      for (int c = i; c < n; ++c) A[r * n + c] -= fac * A[i * n + c];
      b[r] -= fac * b[i];
    }
  }
  // back substitution
  for (int i = n - 1; i >= 0; --i)
  {
    double s = b[i];
    for (int c = i + 1; c < n; ++c) s -= A[i * n + c] * x[c];
    double denom = A[i * n + i];
    if (std::fabs(denom) < EPS) return false;
    x[i] = s / denom;
  }
  return true;
}

// compute DtD for second-difference operator D (rows m = n-2)
static void build_DtD(int n, std::vector<double>& DtD)
{
  DtD.assign(n * n, 0.0);
  if (n < 3) return;
  int m = n - 2;
  for (int r = 0; r < m; ++r)
  {
    int a = r;
    int b = r + 1;
    int c = r + 2;
    // add outer product of [1, -2, 1]
    DtD[a * n + a] += 1.0;
    DtD[a * n + b] += -2.0;
    DtD[a * n + c] += 1.0;

    DtD[b * n + a] += -2.0;
    DtD[b * n + b] += 4.0;
    DtD[b * n + c] += -2.0;

    DtD[c * n + a] += 1.0;
    DtD[c * n + b] += -2.0;
    DtD[c * n + c] += 1.0;
  }
}

// compute smoothed c = (I + lambda DtD)^{-1} y and optionally return trace(M^{-1})
static bool solve_smoothed(const std::vector<double>& DtD,
                           const std::vector<double>& y,
                           double lambda,
                           std::vector<double>& c,
                           double* trace_out = nullptr)
{
  int n = static_cast<int>(y.size());
  std::vector<double> M(n * n, 0.0);
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j) M[i * n + j] = lambda * DtD[i * n + j];
    M[i * n + i] += 1.0;
  }
  // solve M c = y
  bool ok = solve_linear_system(M, y, c, n);
  if (!ok) return false;

  if (trace_out)
  {
    // compute diagonal of inverse by solving M X = I and summing diagonal
    std::vector<double> col(n);
    double trace = 0.0;
    std::vector<double> e(n);
    for (int k = 0; k < n; ++k)
    {
      std::fill(e.begin(), e.end(), 0.0);
      e[k] = 1.0;
      std::vector<double> sol;
      bool ok2 = solve_linear_system(M, e, sol, n);
      if (!ok2) return false;
      trace += sol[k];
    }
    *trace_out = trace;
  }
  return true;
}

SmoothingSpline::SmoothingSpline(const std::vector<double>& xs,
                                 const std::vector<double>& ys,
                                 double lambda)
  : xs_(xs), smoothed_y_(), spline_(std::vector<double>(1,0.0), std::vector<double>(1,0.0))
{
  ok_ = false;
  const int n = static_cast<int>(xs.size());
  if (n < 2 || static_cast<int>(ys.size()) != n) return;

  // For very small n fall back to interpolation (no smoothing)
  if (n < 4)
  {
    smoothed_y_ = ys;
    try { spline_ = CubicSpline2d(xs_, smoothed_y_); ok_ = true; } catch (...) { ok_ = false; }
    return;
  }

  // build DtD (discrete second-difference penalty)
  std::vector<double> DtD;
  build_DtD(n, DtD);

  // scale DtD to account for non-unit knot spacing. The discrete second
  // difference approximates h^2 f''; the continuous roughness penalty
  // integral scales like sum((f'')^2) * h. Combining these gives an
  // effective factor ~ 1 / h^3 for the discrete DtD matrix. Use the
  // average knot spacing as an approximation.
  double h = (xs_.back() - xs_.front()) / static_cast<double>(n - 1);
  if (h > 0.0)
  {
    double scale = 1.0 / (h * h * h);
    for (int i = 0; i < n * n; ++i) DtD[i] *= scale;
  }

  // choose lambda by GCV if requested (lambda <= 0)
  double lambda_opt = lambda;
  if (lambda_opt <= 0.0)
  {
    // grid search on log-space
    const int GRID = 120;
    double best_l = 0.0;
    double best_gcv = std::numeric_limits<double>::infinity();
    for (int i = 0; i < GRID; ++i)
    {
      double logl = -12.0 + 24.0 * static_cast<double>(i) / static_cast<double>(GRID - 1);
      double l = std::pow(10.0, logl);
      std::vector<double> ctmp;
      double trace = 0.0;
      bool ok = solve_smoothed(DtD, ys, l, ctmp, &trace);
      if (!ok) continue;
      // compute RSS
      double rss = 0.0;
      for (int k = 0; k < n; ++k) { double r = ys[k] - ctmp[k]; rss += r * r; }
      double denom = static_cast<double>(n) - trace;
      if (denom <= 0.0) continue;
      double gcv = rss / (denom * denom);
      if (gcv < best_gcv) { best_gcv = gcv; best_l = l; }
    }
    if (best_gcv < std::numeric_limits<double>::infinity()) lambda_opt = best_l; else lambda_opt = 1.0;
  }

  // final solve with lambda_opt
  std::vector<double> c;
  bool ok = solve_smoothed(DtD, ys, lambda_opt, c, nullptr);
  if (!ok)
  {
    // fallback to original values
    smoothed_y_ = ys;
    try { spline_ = CubicSpline2d(xs_, smoothed_y_); ok_ = true; } catch (...) { ok_ = false; }
    return;
  }

  smoothed_y_ = c;
  try { spline_ = CubicSpline2d(xs_, smoothed_y_); ok_ = true; } catch (...) { ok_ = false; }
}


double SmoothingSpline::eval(double x) const
{
  if (!ok_) return std::numeric_limits<double>::quiet_NaN();
  try { return spline_.eval(x); } catch (...) { return std::numeric_limits<double>::quiet_NaN(); }
}

} // namespace Math
} // namespace OpenMS
