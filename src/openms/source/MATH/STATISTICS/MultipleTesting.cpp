// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib> // for getenv
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

// Use the evergreen umbrella header so the vendor headers are included
// consistently inside the `evergreen` namespace and their helper
// macros are defined. Then include the FFT submodule *inside* the
// `evergreen` namespace so its template and type definitions live in
// `evergreen::` as the rest of the code expects.
#include <Evergreen/evergreen.hpp>
namespace evergreen {
#include <FFT/FFT.hpp>
}

#include <boost/math/distributions/normal.hpp>

#include <OpenMS/MATH/MISC/BSpline2d.h>
#include <OpenMS/MATH/MISC/BSplineSmoothingSpline.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/MATH/STATISTICS/KernelDensityEstimation.h>
#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

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

std::vector<double> qValue(const std::vector<double>& p_values, double pi0, bool pfdr)
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
  for (auto vv : p) if (vv < 0.0 || vv > 1.0) throw std::invalid_argument("qValue: p-values not in [0,1]");
  if (pi0 < 0.0 || pi0 > 1.0) throw std::invalid_argument("qValue: pi0 not in [0,1]");

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

Pi0Result pi0Est(const std::vector<double>& p_values, const std::vector<double>& lambda_, const std::string& pi0_method, int smooth_df, bool smooth_log_pi0)
{
  // filter finite
  std::vector<double> p;
  for (double v : p_values) if (std::isfinite(v)) p.push_back(v);
  const std::size_t m = p.size();
  if (m == 0) throw std::invalid_argument("pi0Est: no finite p-values provided");

  // default lambda range if empty
  std::vector<double> lambda_v = lambda_;
  if (lambda_v.empty())
  {
    for (double l = 0.05; l < 1.0 - 1e-12; l += 0.05) lambda_v.push_back(l);
  }

  const std::size_t ll = lambda_v.size();
  if (ll == 0) throw std::invalid_argument("pi0est: empty lambda");

  // validations
  for (double v : p) if (v < 0.0 || v > 1.0) throw std::invalid_argument("pi0Est: p-values not in [0,1]");
  for (double l : lambda_v) if (l < 0.0 || l >= 1.0) throw std::invalid_argument("pi0Est: lambda must be in [0,1)");

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
    // If too few lambda values for smoothing, fall back to minimum pi0
    if (ll < 4)
    {
      res.pi0 = std::min(*std::min_element(pi0s.begin(), pi0s.end()), 1.0);
      res.pi0_smooth = false;
      return res;
    }

    // Optionally work on log(pi0) to match PyProphet's smooth_log_pi0 behavior
    std::vector<double> y = pi0s;
    if (smooth_log_pi0)
    {
      for (double &v : y) v = (v > 0.0) ? std::log(v) : -std::numeric_limits<double>::infinity();
    }

    // Build a sorted unique mapping of lambda -> y (y may be log(pi0) if requested)
    std::map<double, double> xy;
    for (std::size_t i = 0; i < ll; ++i) xy[lambda_v[i]] = y[i];

    if (xy.size() < 2)
    {
      res.pi0 = std::min(*std::min_element(pi0s.begin(), pi0s.end()), 1.0);
      res.pi0_smooth = false;
      return res;
    }

    // Extract sorted vectors
    std::vector<double> xs; xs.reserve(xy.size());
    std::vector<double> ys; ys.reserve(xy.size());
    for (const auto &kv : xy) { xs.push_back(kv.first); ys.push_back(kv.second); }

    // Use BSplineSmoothingSpline which matches scipy's UnivariateSpline behavior
    // with automatic smoothing parameter selection (s = m - sqrt(2*m))
    double max_lambda = *std::max_element(lambda_v.begin(), lambda_v.end());
    OpenMS::Math::BSplineSmoothingSpline spl(xs, ys, -1.0, smooth_df);
    double pred = spl.eval(max_lambda);

    // If working in log-space, exponentiate the result
    if (smooth_log_pi0)
    {
      pred = std::exp(pred);
    }

    // Clamp to [0, 1] range
    res.pi0 = std::min(std::max(pred, 0.0), 1.0);
    res.pi0_smooth = true;
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

std::vector<double> pNorm(const std::vector<double>& stat, const std::vector<double>& stat0)
{
  const std::size_t n = stat.size();
  std::vector<double> out(n, std::numeric_limits<double>::quiet_NaN());

  // require stat0 non-empty
  if (stat0.empty()) throw std::invalid_argument("pNorm: stat0 must be non-empty");

  // filter finite entries in stat0
  std::vector<double> s0;
  s0.reserve(stat0.size());
  for (double v : stat0) if (std::isfinite(v)) s0.push_back(v);
  const std::size_t m = s0.size();
  if (m == 0) throw std::invalid_argument("pNorm: stat0 contains no finite values");

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

    double bw = bwNrd0(x) * adj;
    y = kdeFftEval(x, bw, gridsize, cut);

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
    double bw = bwNrd0(x) * adj;
    y = kdeFftEval(x, bw, gridsize, cut);

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
