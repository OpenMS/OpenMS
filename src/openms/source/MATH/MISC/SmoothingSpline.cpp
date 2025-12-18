// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
#include <OpenMS/MATH/MISC/SmoothingSpline.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{
namespace Math
{
// small helper: solve linear system A x = b using Gaussian elimination with partial pivoting
static bool solve_linear_system(std::vector<double> A, std::vector<double> b, std::vector<double>& x, int n)
{
  // slightly relaxed pivot tolerance to avoid false singular detections
  // make EPS a bit laxer to tolerate near-singular matrices from DtD scaling
  const double EPS = 1e-08;
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
  // small ridge added to diagonal to stabilise ill-conditioned systems
  // increased to improve numerical robustness
  const double RIDGE = 1e-4;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j) M[i * n + j] = lambda * DtD[i * n + j];
    M[i * n + i] += 1.0 + RIDGE;
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
      if (!ok2)
      {
        // If any column solve fails, consider this lambda invalid for GCV
        return false;
      }
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

  // Quick early-entry debug so we can detect whether the constructor is
  // invoked at all in the failing unit-test runs. This uses the same env
  // gates as the later diagnostics.
  {
    const char* dbg_top = std::getenv("OPENMS_TEST_VERBOSE");
    const char* dbg_force_top = std::getenv("OPENMS_FORCE_SMOOTH_DEBUG");
    const char* dbg_log_top = std::getenv("OPENMS_DEBUG_LOG");
    if ((!dbg_top || dbg_top[0] == '\0') && dbg_force_top && dbg_force_top[0] != '\0') dbg_top = dbg_force_top;
    if ((!dbg_top || dbg_top[0] == '\0') && dbg_log_top && dbg_log_top[0] != '\0') dbg_top = dbg_log_top;
    if (dbg_top && dbg_top[0] != '\0')
    {
      // unbuffered syscall for strace and immediate stderr print
      char tmpbuf[128];
      int tmplen = std::snprintf(tmpbuf, sizeof(tmpbuf), "SmoothingSpline ctor ENTRY n=%d\n", n);
      if (tmplen > 0) (void)write(STDERR_FILENO, tmpbuf, static_cast<size_t>(tmplen));
      std::cerr << "SmoothingSpline ctor ENTRY n=" << n << std::endl;
    }
  }

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
  // candidate list for GCV diagnostics (declared outer so debug printing can access it)
  struct Candidate { double logl; double l; double gcv; double rss; double trace; };
  std::vector<Candidate> candidates;
  if (lambda_opt <= 0.0)
  {
    // grid search on log-space. widen the search and increase resolution
    const int GRID = 200;
    double best_l = 0.0;
    double best_gcv = std::numeric_limits<double>::infinity();
    for (int i = 0; i < GRID; ++i)
    {
      // search log10(lambda) in [-10, +4]
      double logl = -10.0 + 14.0 * static_cast<double>(i) / static_cast<double>(GRID - 1);
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
      candidates.push_back(Candidate{logl, l, gcv, rss, trace});
      if (gcv < best_gcv) { best_gcv = gcv; best_l = l; }
    }
    if (!candidates.empty())
    {
      // choose best candidate
      lambda_opt = best_l;
    }
    else
    {
      lambda_opt = 1.0;
    }
  }
  // optional debug output to help diagnose GCV/solve issues
  const char* dbg_env = std::getenv("OPENMS_TEST_VERBOSE");
  // allow forced smoothing debug via OPENMS_FORCE_SMOOTH_DEBUG in case
  // the standard test variable is not set (helps during CI/debug runs)
  const char* dbg_force = std::getenv("OPENMS_FORCE_SMOOTH_DEBUG");
  // also accept the project's OpenMS debug env var so instrumentation
  // integrates with the normal logging switches used in this repo
  const char* dbg_log = std::getenv("OPENMS_DEBUG_LOG");
  if ((!dbg_env || dbg_env[0] == '\0') && dbg_force && dbg_force[0] != '\0') dbg_env = dbg_force;
  if ((!dbg_env || dbg_env[0] == '\0') && dbg_log && dbg_log[0] != '\0') dbg_env = dbg_log;
  // If forced debug is present, also write a tiny marker to a temp file so
  // we can detect runtime execution even if stderr is redirected or
  // suppressed by the test harness.
  if ((dbg_force && dbg_force[0] != '\0') || (dbg_log && dbg_log[0] != '\0'))
  {
    // Use the project File utility so the marker file is created in the
    // configured temp directory (OPENMS_TMPDIR or system temp). This keeps
    // instrumentation consistent with the project's file handling.
    OpenMS::String marker = File::getTempDirectory();
    marker.ensureLastChar('/');
    marker += "SmoothingSpline_instrumentation.log";

    FILE* f = std::fopen(marker.c_str(), "a");
    if (f)
    {
      std::fprintf(f, "SmoothingSpline ctor entered (n=%d)\n", n);
      std::fclose(f);
    }

    // Also emit an unbuffered write to stderr to ensure visibility when the
    // test harness redirects or buffers stdio. Use a single syscall so it
    // appears in strace output reliably.
    {
      char buf[256];
      int len = std::snprintf(buf, sizeof(buf), "SmoothingSpline ctor (unbuf) n=%d marker=%s\n", n, marker.c_str());
      if (len > 0) (void)write(STDERR_FILENO, buf, static_cast<size_t>(len));
    }

    // Also print via std::cerr so the debug information appears on stderr
    // (the test harness captures stderr and shows it in logs).
    std::cerr << "SmoothingSpline ctor entered (n=" << n << ", marker=" << marker << ")" << std::endl;
  }
  if (dbg_env && dbg_env[0] != '\0')
  {
    std::cerr << "SmoothingSpline debug: n=" << n << " lambda_opt=" << lambda_opt << "\n";
    if (!candidates.empty())
    {
      // sort candidates by gcv ascending and print the top ones
      std::sort(candidates.begin(), candidates.end(), [](const Candidate& a, const Candidate& b){ return a.gcv < b.gcv; });
      std::cerr << "SmoothingSpline debug: top GCV candidates (log10(lambda), lambda, gcv, rss, trace):\n";
      int limit = std::min(static_cast<int>(candidates.size()), 10);
      for (int i = 0; i < limit; ++i)
      {
        const Candidate& cnd = candidates[i];
        std::cerr << "  " << cnd.logl << ", " << cnd.l << ", " << cnd.gcv << ", " << cnd.rss << ", " << cnd.trace << "\n";
      }
    }
  }

  // final solve with lambda_opt
  std::vector<double> c;
  bool ok = solve_smoothed(DtD, ys, lambda_opt, c, nullptr);
  if (!ok)
  {
    // fallback to original values
    smoothed_y_ = ys;
    try { spline_ = CubicSpline2d(xs_, smoothed_y_); ok_ = true; }
    catch (...) { ok_ = false; }
    if (dbg_env && dbg_env[0] != '\0')
    {
      std::cerr << "SmoothingSpline debug: final solve failed, falling back to raw ys (n=" << n << ")\n";
    }
    return;
  }

  smoothed_y_ = c;
  try { spline_ = CubicSpline2d(xs_, smoothed_y_); ok_ = true; } catch (...) { ok_ = false; }
  if (dbg_env && dbg_env[0] != '\0')
  {
    std::cerr << "SmoothingSpline debug: solved OK, smoothed_y = [";
    for (int i = 0; i < n; ++i)
    {
      if (i) std::cerr << ", ";
      std::cerr << smoothed_y_[i];
    }
    std::cerr << "]\n";
    // predicted value at the rightmost knot
    std::cerr << "SmoothingSpline debug: pred_at_max_x = " << smoothed_y_.back() << " (x=" << xs_.back() << ")\n";
  }
}


double SmoothingSpline::eval(double x) const
{
  if (!ok_) return std::numeric_limits<double>::quiet_NaN();
  try { return spline_.eval(x); } catch (...) { return std::numeric_limits<double>::quiet_NaN(); }
}

} // namespace Math
} // namespace OpenMS
