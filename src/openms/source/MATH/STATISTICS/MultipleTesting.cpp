// SPDX-License-Identifier: BSD-3-Clause
// Implementation of qvalue and pi0est (bootstrap)

#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

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

Pi0Result pi0est(const std::vector<double>& p_values, const std::vector<double>& lambda_)
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

  // bootstrap selection: compute minpi0 at 10th percentile
  double minpi0 = percentile(pi0s, 0.1);

  std::vector<double> W;
  W.reserve(ll);
  for (double l : lambda_v)
  {
    double cnt = 0.0;
    for (double v : p) if (v >= l) cnt += 1.0;
    W.push_back(cnt);
  }

  std::vector<double> mse(ll);
  for (std::size_t i = 0; i < ll; ++i)
  {
    double w = W[i];
    double denom = static_cast<double>(m) * static_cast<double>(m) * std::pow(1.0 - lambda_v[i], 2);
    double term = 0.0;
    if (denom != 0.0) term = (w / denom) * (1.0 - w / static_cast<double>(m));
    double diff = (pi0s[i] - minpi0);
    mse[i] = term + diff * diff;
  }

  // pick index with minimal mse
  std::size_t argmin = 0;
  for (std::size_t i = 1; i < mse.size(); ++i) if (mse[i] < mse[argmin]) argmin = i;
  res.pi0 = std::min(pi0s[argmin], 1.0);
  res.pi0_smooth = false;
  return res;
}

} // namespace Math
} // namespace OpenMS
