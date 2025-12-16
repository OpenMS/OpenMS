// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <complex>
#include <cmath>

// Use the evergreen umbrella header so the vendor headers are included
// consistently inside the `evergreen` namespace and their helper
// macros are defined. Then include the FFT submodule *inside* the
// `evergreen` namespace so its template and type definitions live in
// `evergreen::` as the rest of the code expects.
#include <Evergreen/evergreen.hpp>
namespace evergreen {
#include <FFT/FFT.hpp>
}

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <iterator>
#include <utility>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <boost/math/distributions/normal.hpp>

namespace OpenMS
{
namespace Math
{
// Helper: stable argsort ascending
static std::vector<std::size_t> argsort_asc(const std::vector<double>& a)
{
  std::vector<std::size_t> idx(a.size());
  for (std::size_t i = 0; i < a.size(); ++i) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(), [&](std::size_t i, std::size_t j){ return a[i] < a[j]; });
  return idx;
}

std::vector<double> qvalue(const std::vector<double>& p_values, double pi0, bool pfdr)
{
  const std::size_t m_total = p_values.size();
  if (m_total == 0) return std::vector<double>();

  // mark finite entries
  std::vector<char> rm_na(m_total, 0);
  std::vector<double> p;
  p.reserve(m_total);
  for (std::size_t i = 0; i < m_total; ++i)
  {
    if (std::isfinite(p_values[i])) { rm_na[i] = 1; p.push_back(p_values[i]); }
  }

  std::vector<double> qvals_out(m_total, std::numeric_limits<double>::quiet_NaN());

  const std::size_t m = p.size();
  if (m == 0) return qvals_out;

  // validate
  for (auto vv : p) if (vv < 0.0 || vv > 1.0) throw std::invalid_argument("qvalue: p-values not in [0,1]");
  if (pi0 < 0.0 || pi0 > 1.0) throw std::invalid_argument("qvalue: pi0 not in [0,1]");

  // argsort ascending
  std::vector<std::size_t> u = argsort_asc(p);

  // ranks (max)
  std::vector<double> v = RankData::rankdata<double>(p, RankData::Method::Max, RankData::NaNPolicy::Propagate);

  std::vector<double> qvals(m);
  if (pfdr)
  {
    for (std::size_t i = 0; i < m; ++i)
    {
      double denom = v[i] * (1.0 - std::pow(1.0 - p[i], static_cast<double>(m)));
      if (denom == 0.0) qvals[i] = std::numeric_limits<double>::infinity();
      else qvals[i] = (pi0 * static_cast<double>(m) * p[i]) / denom;
    }
  }
  else
  {
    for (std::size_t i = 0; i < m; ++i)
    {
      double denom = v[i];
      if (denom == 0.0) qvals[i] = std::numeric_limits<double>::infinity();
      else qvals[i] = (pi0 * static_cast<double>(m) * p[i]) / denom;
    }
  }

  // enforce monotonicity over sorted indices u
  // qvals[u[m-1]] = min(qvals[u[m-1]], 1)
  std::size_t last = u[m-1];
  qvals[last] = std::min(qvals[last], 1.0);
  if (qvals[last] != qvals[last]) qvals[last] = 1.0; // NaN guard

  for (std::size_t ii = m - 1; ii-- > 0; )
  {
    std::size_t idx = u[ii];
    std::size_t idx_next = u[ii+1];
    qvals[idx] = std::min(qvals[idx], qvals[idx_next]);
  }

  // place back into output
  std::size_t k = 0;
  for (std::size_t i = 0; i < m_total; ++i)
  {
    if (rm_na[i]) { qvals_out[i] = qvals[k++]; }
  }

  return qvals_out;
}

// helper percentile: simple empirical percentile (nearest-rank)
static double percentile(const std::vector<double>& a, double p)
{
  if (a.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double> cp = a;
  std::sort(cp.begin(), cp.end());
  double pos = p * (static_cast<double>(cp.size()) - 1.0);
  std::size_t idx = static_cast<std::size_t>(std::floor(pos + 0.5));
  if (idx >= cp.size()) idx = cp.size() - 1;
  return cp[idx];
}

Pi0Result pi0est(const std::vector<double>& p_values, const std::vector<double>& lambda_, const std::string& pi0_method, int smooth_df, bool smooth_log_pi0)
{
  // filter finite
  std::vector<double> p;
  for (double v : p_values) if (std::isfinite(v)) p.push_back(v);
  const std::size_t m = p.size();
  if (m == 0) throw std::invalid_argument("pi0est: no finite p-values provided");

  // default lambda range if empty
  std::vector<double> lambda_v = lambda_;
  if (lambda_v.empty())
  {
    for (double l = 0.05; l < 1.0 - 1e-12; l += 0.05) lambda_v.push_back(l);
  }

  const std::size_t ll = lambda_v.size();
  if (ll == 0) throw std::invalid_argument("pi0est: empty lambda");

  // validations
  for (double v : p) if (v < 0.0 || v > 1.0) throw std::invalid_argument("pi0est: p-values not in [0,1]");
  for (double l : lambda_v) if (l < 0.0 || l >= 1.0) throw std::invalid_argument("pi0est: lambda must be in [0,1)");

  Pi0Result res;
  res.lambda_ = lambda_v;

  if (ll == 1)
  {
    double l = lambda_v[0];
    double frac = 0.0;
    for (double v : p) if (v >= l) frac += 1.0;
    frac /= static_cast<double>(m);
    double pi0 = frac / (1.0 - l);
    res.pi0_lambda = {std::min(pi0, 1.0)};
    res.pi0 = std::min(pi0, 1.0);
    res.pi0_smooth = false;
    return res;
  }

  // compute pi0 for each lambda
  std::vector<double> pi0s;
  pi0s.reserve(ll);
  for (double l : lambda_v)
  {
    double cnt = 0.0;
    for (double v : p) if (v >= l) cnt += 1.0;
    double pi0 = (cnt / static_cast<double>(m)) / (1.0 - l);
    pi0s.push_back(pi0);
  }
  res.pi0_lambda = pi0s;

  if (ll == 1)
  {
    // already handled above, but keep defensive
    return res;
  }

  // Choose method: 'smoother' (default) or 'bootstrap'
  if (pi0_method == "smoother")
  {
    // If too few lambda values are provided for a smoothing fit, fall back
    // to the naive (non-smoothed) minimum pi0 across the grid. This keeps
    // behaviour compatible with existing callers that may pass small
    // lambda vectors (see unit tests).
    if (ll < 4)
    {
      res.pi0 = std::min(*std::min_element(pi0s.begin(), pi0s.end()), 1.0);
      res.pi0_smooth = false;
      return res;
    }

    // Optionally work on log(pi0)
    std::vector<double> y = pi0s;
    if (smooth_log_pi0)
    {
      for (double &v : y) v = (v > 0.0) ? std::log(v) : -std::numeric_limits<double>::infinity();
    }

    // Fit a small-degree polynomial (up to cubic) by least squares to approximate
    // the smoother used by the R/qvalue reference implementation. Use a
    // centered-and-scaled Vandermonde basis for numerical stability which tends
    // to reduce small numeric discrepancies versus an unscaled fit.
    const int deg_i = std::min(3, std::max(1, smooth_df));
    const std::size_t deg = static_cast<std::size_t>(deg_i);
    const std::size_t D = std::min<std::size_t>(deg + 1, ll);

    // compute centering/scale for lambda values
    double lambda_min = *std::min_element(lambda_v.begin(), lambda_v.end());
    double lambda_max = *std::max_element(lambda_v.begin(), lambda_v.end());
    double lambda_mean = 0.5 * (lambda_min + lambda_max);
    double lambda_scale = lambda_max - lambda_min;
    if (!(lambda_scale > 0.0)) lambda_scale = 1.0;

    // Build normal equations A coeffs = B using scaled x = (x - mean)/scale
    std::vector<std::vector<double>> A(D, std::vector<double>(D, 0.0));
    std::vector<double> B(D, 0.0);
    for (std::size_t i = 0; i < ll; ++i)
    {
      double x = (lambda_v[i] - lambda_mean) / lambda_scale;
      double xp = 1.0;
      std::vector<double> xv(D);
      for (std::size_t j = 0; j < D; ++j) { xv[j] = xp; xp *= x; }
      for (std::size_t r = 0; r < D; ++r)
      {
        for (std::size_t c = 0; c < D; ++c) A[r][c] += xv[r] * xv[c];
        B[r] += xv[r] * y[i];
      }
    }

    // Solve via Gaussian elimination on augmented matrix (D small <=4)
    std::vector<std::vector<double>> M(D, std::vector<double>(D + 1, 0.0));
    for (std::size_t i = 0; i < D; ++i)
    {
      for (std::size_t j = 0; j < D; ++j) M[i][j] = A[i][j];
      M[i][D] = B[i];
    }
    bool solve_ok = true;
    for (std::size_t i = 0; i < D; ++i)
    {
      // pivot selection
      std::size_t piv = i;
      for (std::size_t r = i + 1; r < D; ++r) if (std::fabs(M[r][i]) > std::fabs(M[piv][i])) piv = r;
      if (std::fabs(M[piv][i]) < 1e-16) { solve_ok = false; break; }
      if (piv != i) std::swap(M[piv], M[i]);
      double diag = M[i][i];
      for (std::size_t c = i; c < D + 1; ++c) M[i][c] /= diag;
      for (std::size_t r = 0; r < D; ++r)
      {
        if (r == i) continue;
        double fac = M[r][i];
        for (std::size_t c = i; c < D + 1; ++c) M[r][c] -= fac * M[i][c];
      }
    }

    if (solve_ok)
    {
      std::vector<double> coeffs(D, 0.0);
      for (std::size_t i = 0; i < D; ++i) coeffs[i] = M[i][D];
  // evaluate polynomial at scaled max lambda
  double max_lambda = *std::max_element(lambda_v.begin(), lambda_v.end());
  double x_eval = (max_lambda - lambda_mean) / lambda_scale;
  double xpow = 1.0;
  double pred = 0.0;
  for (std::size_t j = 0; j < D; ++j) { pred += coeffs[j] * xpow; xpow *= x_eval; }
      if (smooth_log_pi0)
      {
        if (pred <= -1e300 || !std::isfinite(pred)) pred = *std::min_element(pi0s.begin(), pi0s.end());
        else pred = std::exp(pred);
      }
      if (!std::isfinite(pred)) pred = *std::min_element(pi0s.begin(), pi0s.end());
      res.pi0 = std::min(std::max(pred, 0.0), 1.0);
      res.pi0_smooth = true;
      return res;
    }

    // fallback to naive min pi0
    res.pi0 = std::min(*std::min_element(pi0s.begin(), pi0s.end()), 1.0);
    res.pi0_smooth = false;
    return res;
  }
  else if (pi0_method == "bootstrap")
  {
    // bootstrap variant as used by pyprophet/qvalue
    double minpi0 = percentile(pi0s, 0.1);
    std::vector<double> W;
    W.reserve(ll);
    for (double l : lambda_v)
    {
      double cnt = 0.0;
      for (double v : p) if (v >= l) cnt += 1.0;
      W.push_back(cnt);
    }
    std::vector<double> mse(ll, 0.0);
    for (std::size_t i = 0; i < ll; ++i)
    {
      double l = lambda_v[i];
      double w = W[i];
      double term1 = 0.0;
      double denom = static_cast<double>(m) * static_cast<double>(m) * (1.0 - l) * (1.0 - l);
      if (denom > 0.0) term1 = (w / denom) * (1.0 - w / static_cast<double>(m));
      double term2 = (pi0s[i] - minpi0) * (pi0s[i] - minpi0);
      mse[i] = term1 + term2;
    }
    std::size_t argmin = 0;
    double best = mse[0];
    for (std::size_t i = 1; i < mse.size(); ++i) if (mse[i] < best) { best = mse[i]; argmin = i; }
    res.pi0 = std::min(pi0s[argmin], 1.0);
    res.pi0_smooth = false;
    return res;
  }
  else
  {
    throw std::invalid_argument("pi0est: pi0_method must be 'smoother' or 'bootstrap'");
  }
}

std::vector<double> pnorm(const std::vector<double>& stat, const std::vector<double>& stat0)
{
  const std::size_t n = stat.size();
  std::vector<double> out(n, std::numeric_limits<double>::quiet_NaN());

  // require stat0 non-empty
  if (stat0.empty()) throw std::invalid_argument("pnorm: stat0 must be non-empty");

  // filter finite entries in stat0
  std::vector<double> s0;
  s0.reserve(stat0.size());
  for (double v : stat0) if (std::isfinite(v)) s0.push_back(v);
  const std::size_t m = s0.size();
  if (m == 0) throw std::invalid_argument("pnorm: stat0 contains no finite values");

  // compute mean
  double sum = 0.0;
  for (double v : s0) sum += v;
  double mu = sum / static_cast<double>(m);

  // compute sample standard deviation (ddof=1) when possible
  double var = 0.0;
  if (m > 1)
  {
    for (double v : s0)
    {
      double d = v - mu;
      var += d * d;
    }
    var /= static_cast<double>(m - 1);
  }
  else
  {
    var = 0.0;
  }
  double sigma = std::sqrt(var);

  const double sqrt2 = std::sqrt(2.0);

  // If sigma == 0, treat distribution as degenerate at mu: P(X>stat_i) = 1 if stat_i < mu, 0 otherwise
  const bool degenerate = !(sigma > 0.0);

  for (std::size_t i = 0; i < n; ++i)
  {
    double v = stat[i];
    if (!std::isfinite(v)) { out[i] = std::numeric_limits<double>::quiet_NaN(); continue; }
    if (degenerate)
    {
      out[i] = (v < mu) ? 1.0 : 0.0;
      continue;
    }
    double z = (v - mu) / sigma;
    double cdf = 0.5 * (1.0 + std::erf(z / sqrt2));
    double tail = 1.0 - cdf;
    out[i] = tail;
  }

  return out;
}

double bw_nrd0(const std::vector<double>& x)
{
  // filter finite
  std::vector<double> xf;
  xf.reserve(x.size());
  for (double v : x) if (std::isfinite(v)) xf.push_back(v);
  const std::size_t n = xf.size();
  if (n < 2) return 0.0;

  // sample standard deviation (ddof=1)
  double mean = 0.0;
  for (double v : xf) mean += v;
  mean /= static_cast<double>(n);
  double var = 0.0;
  for (double v : xf) { double d = v - mean; var += d * d; }
  var /= static_cast<double>(n - 1);
  double sd = std::sqrt(var);

  // compute percentiles using existing Math::quantile (Type 7, numpy-like linear interpolation)
  std::vector<double> cp = xf;
  std::sort(cp.begin(), cp.end());
  double q25 = Math::quantile(cp.begin(), cp.end(), 0.25);
  double q75 = Math::quantile(cp.begin(), cp.end(), 0.75);
  double iqr = q75 - q25;

  double lo = std::min(sd, iqr / 1.34);
  // fallbacks to match pyprophet/numpy implementation
  if (!(lo > 0.0))
  {
    if (sd > 0.0) lo = sd;
    else if (!cp.empty()) lo = std::abs(cp[0]);
    else lo = 1.0;
  }

  double bw = 0.9 * lo * std::pow(static_cast<double>(n), -0.2); // Silverman variant used by statsmodels/pyprophet
  return bw;
}

std::vector<double> linbin(const std::vector<double>& x, double xmin, double xmax, std::size_t nbins, const std::vector<double>* weights)
{
  if (nbins == 0) throw std::invalid_argument("linbin: nbins must be > 0");
  std::vector<double> bins(nbins, 0.0);
  if (x.empty()) return bins;
  if (!(xmax > xmin)) throw std::invalid_argument("linbin: xmax must be > xmin");

  const double width = (xmax - xmin) / static_cast<double>(nbins);
  if (!(width > 0.0)) return bins;

  const bool use_weights = (weights != nullptr && weights->size() == x.size());

  for (std::size_t i = 0; i < x.size(); ++i)
  {
    double v = x[i];
    if (!std::isfinite(v)) continue;
    if (v < xmin || v > xmax) continue; // ignore outside range; include xmax in last bin
    std::size_t idx = static_cast<std::size_t>(std::floor((v - xmin) / width));
    if (idx >= nbins) idx = nbins - 1; // edge case when v == xmax
    double w = 1.0;
    if (use_weights) w = (*weights)[i];
    bins[idx] += w;
  }

  return bins;
}

// forward using evergreen real FFT (Munro-packed / rfft-like ordering)
// Matches numpy.fft.rfft semantics (no forward 1/N scaling).
std::vector<double> forrt(const std::vector<double>& X, std::size_t M)
{
  if (M == 0) M = X.size();
  if (M == 0) return std::vector<double>();

  const std::size_t half = M / 2;
  const std::size_t out_len = M;

  // create a 1-D tensor and copy (zero-pad if needed)
  evergreen::Tensor<double> t({M});
  for (std::size_t i = 0; i < M; ++i)
  {
    double v = 0.0;
    if (i < X.size() && std::isfinite(X[i])) v = X[i];
    t[i] = v;
  }

  // use DIF-based real FFT provided by evergreen (packed real format)
  // Use FRESHLY_ZERO_PADDED=true for the single-axis case to match optimized path
  evergreen::Tensor<evergreen::cpx> packed = evergreen::real_fft<evergreen::DIF, false, false, true>(t);

  const std::size_t ny = packed.flat_size(); // should be half + 1
  std::vector<double> out(out_len, 0.0);

  // store real parts
  for (std::size_t k = 0; k < ny; ++k) out[k] = packed[k].r;
  // store imag parts for k=1..half-1
  for (std::size_t k = 1; k < ny - 1; ++k) out[half + k] = packed[k].i;

  return out;
}

// inverse using evergreen real_ifft on Munro-packed spectrum
std::vector<double> revrt(const std::vector<double>& Xp, std::size_t M)
{
  if (M == 0) M = Xp.size();
  if (M == 0) return std::vector<double>();

  const std::size_t half = M / 2;
  const std::size_t ny = half + 1;

  if (Xp.size() < M) throw std::invalid_argument("revrt: input length must equal M");

  // build packed complex tensor
  evergreen::Tensor<evergreen::cpx> packed({ny});
  for (std::size_t k = 0; k < ny; ++k)
  {
    double re = Xp[k];
    double im = 0.0;
    if (k > 0 && k < ny - 1)
    {
      im = Xp[half + k];
    }
  packed[k] = evergreen::cpx{re, im};
  }

  // use evergreen real_ifft which returns a real tensor of length M
  evergreen::Tensor<double> rec = evergreen::real_ifft<evergreen::DIF, false, false>(packed);

  std::vector<double> out(M, 0.0);
  const std::size_t ncopy = std::min(rec.flat_size(), M);
  for (std::size_t i = 0; i < ncopy; ++i) out[i] = rec[i];
  return out;
}

// Munro-packed Silverman Gaussian kernel FFT (matches statsmodels implementation)
std::vector<double> silverman_kernel_fft(double bw, std::size_t M, double RANGE)
{
  if (M == 0) return std::vector<double>();
  const std::size_t half = M / 2;
  // J = 0..M/2
  std::vector<double> J(half + 1);
  for (std::size_t k = 0; k <= half; ++k) J[k] = static_cast<double>(k);
  const double FAC1 = 2.0 * std::pow(M_PI * bw / RANGE, 2);
  std::vector<double> FAC(half + 1);
  for (std::size_t k = 0; k <= half; ++k)
  {
    double j = J[k];
    double jfac = j * j * FAC1;
    double BC = 1.0 - (j * (M_PI / static_cast<double>(M))) * (j * (M_PI / static_cast<double>(M))) / 3.0;
    if (!(BC > 0.0)) BC = std::numeric_limits<double>::min();
    FAC[k] = std::exp(-jfac) / BC;
  }

  // mirror to Munro-packed length M: [FAC[0..half], FAC[1..half-1]]
  std::vector<double> out(M, 0.0);
  // copy FAC[0..half]
  for (std::size_t k = 0; k <= half; ++k) out[k] = FAC[k];
  // mirror FAC[1..half-1]
  if (half >= 2)
  {
    for (std::size_t k = 1; k < half; ++k) out[half + k] = FAC[k];
  }
  return out;
}

// helper: next power of two >= v
static std::size_t next_pow2(std::size_t v)
{
  if (v == 0) return 1;
  std::size_t p = 1;
  while (p < v) p <<= 1;
  return p;
}

// fast linear binning matching statsmodels' endpoint rules (no weights support here)
static std::vector<double> fast_linbin(const std::vector<double>& x, double a, double b, std::size_t M)
{
  // linear binning with inclusive a, inclusive b into M bins (np.linspace(a,b,M))
  std::vector<double> bins(M, 0.0);
  if (x.empty() || !(b > a)) return bins;
  const double delta = (b - a) / static_cast<double>(M - 1); // linspace with M points has M-1 intervals
  for (double v : x)
  {
    if (!std::isfinite(v)) continue;
    if (v < a || v > b) continue;
    double pos = (v - a) / delta;
    std::size_t left = static_cast<std::size_t>(std::floor(pos));
    if (left >= M - 1)
    {
      bins[M - 1] += 1.0;
      continue;
    }
    double frac = pos - static_cast<double>(left);
    bins[left] += 1.0 - frac;
    bins[left + 1] += frac;
  }
  return bins;
}

// grid KDE via Silverman FFT + linear binning, mirrors statsmodels.kdensityfft
std::pair<std::vector<double>, std::vector<double>> grid_kde_fft(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
{
  std::vector<double> xx = x;
  // ensure column vector semantics
  // n
  const std::size_t n = xx.size();

  // M = 2^ceil(log2(max(gridsize, n, 512)))
  std::size_t target = std::max<std::size_t>(gridsize, std::max<std::size_t>(n, 512));
  std::size_t M = next_pow2(target);

  // compute a,b
  double a = 0.0, b = 0.0;
  if (n > 0)
  {
    auto mm = std::minmax_element(xx.begin(), xx.end());
    a = *mm.first - cut * bw;
    b = *mm.second + cut * bw;
  }
  else
  {
    a = -cut * bw;
    b = cut * bw;
  }
  // build grid (M points linspace a..b)
  std::vector<double> grid(M);
  double delta = (b - a) / static_cast<double>(M - 1);
  for (std::size_t i = 0; i < M; ++i) grid[i] = a + static_cast<double>(i) * delta;

  double RANGE = b - a;

  // linear binning normalized by (delta * n)
  std::vector<double> binned = fast_linbin(xx, a, b, M);
  if (n > 0)
  {
    for (double &v : binned) v /= (delta * static_cast<double>(n));
  }

  // forward Munro RFFT
  std::vector<double> Y = forrt(binned, M);

  // multiply by Silverman kernel FFT
  std::vector<double> K = silverman_kernel_fft(bw, M, RANGE);
  if (K.size() != Y.size()) throw std::runtime_error("grid_kde_fft: internal size mismatch");
  std::vector<double> Z(M);
  for (std::size_t i = 0; i < M; ++i) Z[i] = K[i] * Y[i];

  // inverse
  std::vector<double> dens = revrt(Z, M);

  // renormalize to enforce integral 1
  double sum = 0.0;
  for (double v : dens) sum += v;
  if (sum > 0.0 && delta > 0.0)
  {
    double scale = 1.0 / (sum * delta);
    for (double &v : dens) v *= scale;
  }

  return {dens, grid};
}

// Evaluate KDE at query points using FFT-grid method + cubic spline interpolation
std::vector<double> kde_fft_eval(const std::vector<double>& x, double bw, std::size_t gridsize, double cut)
{
  // build grid density
  auto pg = grid_kde_fft(x, bw, gridsize, cut);
  const std::vector<double>& dens = pg.first;
  const std::vector<double>& grid = pg.second;

  // build cubic spline (CubicSpline2d expects sorted x)
  OpenMS::CubicSpline2d spline(grid, dens);

  std::vector<double> out;
  out.reserve(x.size());
  for (double xi : x)
  {
    out.push_back(spline.eval(xi));
  }
  return out;
}

// Estimate local FDR / PEP from p-values using probit/logit transform and FFT-grid KDE
std::vector<double> lfdr(const std::vector<double>& p_values,
                         double pi0,
                         bool trunc,
                         bool monotone,
                         const std::string& transf,
                         double adj,
                         double eps,
                         std::size_t gridsize,
                         double cut)
{
  const std::size_t m_total = p_values.size();
  std::vector<double> lfdr_out(m_total, std::numeric_limits<double>::quiet_NaN());

  if (m_total == 0) return lfdr_out;

  // select finite entries
  std::vector<char> rm_na(m_total, 0);
  std::vector<double> p;
  p.reserve(m_total);
  for (std::size_t i = 0; i < m_total; ++i)
  {
    if (std::isfinite(p_values[i])) { rm_na[i] = 1; p.push_back(p_values[i]); }
  }

  const std::size_t m = p.size();
  if (m == 0)
  {
    return lfdr_out; // all NA
  }

  // validate ranges
  for (double pv : p) if (pv < 0.0 || pv > 1.0) throw std::invalid_argument("lfdr: p-values not in [0,1]");
  if (pi0 < 0.0 || pi0 > 1.0) throw std::invalid_argument("lfdr: pi0 not in [0,1]");

  // working vectors
  std::vector<double> x(m);
  std::vector<double> y;

  if (transf == "probit")
  {
    // clip p to avoid infinities
    for (std::size_t i = 0; i < m; ++i)
    {
      double pv = p[i];
      if (pv < eps) pv = eps;
      if (pv > 1.0 - eps) pv = 1.0 - eps;
      p[i] = pv; // replace clipped
    }

    // inverse normal (probit)
    boost::math::normal_distribution<double> nd(0.0, 1.0);
    for (std::size_t i = 0; i < m; ++i) x[i] = boost::math::quantile(nd, p[i]);

    double bw = bw_nrd0(x) * adj;
    y = kde_fft_eval(x, bw, gridsize, cut);

    // null density f0 = N(0,1)
    y.reserve(y.size());
    std::vector<double> f0(m);
    const double norm_const = 1.0 / std::sqrt(2.0 * M_PI);
    for (std::size_t i = 0; i < m; ++i)
    {
      f0[i] = norm_const * std::exp(-0.5 * x[i] * x[i]);
    }

    // lfdr = pi0 * f0 / y
    std::vector<double> lfdr_v(m, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < m; ++i)
    {
      if (!(y[i] > 0.0)) lfdr_v[i] = std::numeric_limits<double>::infinity();
      else lfdr_v[i] = pi0 * f0[i] / y[i];
    }
    y.swap(lfdr_v); // temporarily store lfdr in y variable
  }
  else if (transf == "logit")
  {
    for (std::size_t i = 0; i < m; ++i)
    {
      double pv = p[i];
      // use eps to stabilize
      x[i] = std::log((pv + eps) / (1.0 - pv + eps));
    }
    double bw = bw_nrd0(x) * adj;
    y = kde_fft_eval(x, bw, gridsize, cut);

    // dx = exp(x) / (1+exp(x))^2
    std::vector<double> dx(m);
    for (std::size_t i = 0; i < m; ++i)
    {
      double ex = std::exp(x[i]);
      double denom = (1.0 + ex);
      dx[i] = ex / (denom * denom);
    }

    std::vector<double> lfdr_v(m, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < m; ++i)
    {
      if (!(y[i] > 0.0)) lfdr_v[i] = std::numeric_limits<double>::infinity();
      else lfdr_v[i] = (pi0 * dx[i]) / y[i];
    }
    y.swap(lfdr_v);
  }
  else
  {
    throw std::invalid_argument("lfdr: invalid transf, expected 'probit' or 'logit'");
  }

  // y now holds lfdr for the finite entries in original order
  std::vector<double> lfdr_vec = y;

  // truncation
  if (trunc)
  {
    for (double &v : lfdr_vec) if (v > 1.0) v = 1.0;
  }

  // monotone enforcement in p: make lfdr non-decreasing in p
  if (monotone && m > 1)
  {
    // argsort p ascending
    std::vector<std::size_t> order = argsort_asc(p);
    // build sorted lfdr
    std::vector<double> lfdr_sorted(m);
    for (std::size_t i = 0; i < m; ++i) lfdr_sorted[i] = lfdr_vec[order[i]];
    // cumulative maximum (isotonic non-decreasing)
    for (std::size_t i = 1; i < m; ++i)
    {
      if (lfdr_sorted[i] < lfdr_sorted[i - 1]) lfdr_sorted[i] = lfdr_sorted[i - 1];
    }
    // compute min-rank mapping (rankdata 'min')
    std::vector<int> assigned(m, 0);
    std::vector<std::size_t> minrank(m, 0);
    for (std::size_t j = 0; j < m; ++j)
    {
      std::size_t idx = order[j];
      if (!assigned[idx])
      {
        minrank[idx] = j + 1; // 1-based
        assigned[idx] = 1;
      }
    }
    // map back: lfdr_mapped[i] = lfdr_sorted[minrank[i]-1]
    std::vector<double> lfdr_mapped(m, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < m; ++i)
    {
      std::size_t r = static_cast<std::size_t>(minrank[i]);
      if (r == 0) lfdr_mapped[i] = std::numeric_limits<double>::quiet_NaN();
      else lfdr_mapped[i] = lfdr_sorted[r - 1];
    }
    lfdr_vec.swap(lfdr_mapped);
  }

  // place back into output
  std::size_t k = 0;
  for (std::size_t i = 0; i < m_total; ++i)
  {
    if (rm_na[i]) { lfdr_out[i] = lfdr_vec[k++]; }
  }

  return lfdr_out;
}

} // namespace Math
} // namespace OpenMS
