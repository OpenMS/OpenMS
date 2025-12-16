// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <OpenMS/config.h>
#include <OpenMS/MATH/STATISTICS/RankData.h>

namespace OpenMS
{
namespace Math
{
struct OPENMS_DLLAPI Pi0Result
{
  double pi0 = 1.0;
  std::vector<double> pi0_lambda;
    std::vector<double> lambda_;
  bool pi0_smooth = false;
};

/// Compute model-based FDR estimates from posterior error probabilities
template <class T>
inline std::vector<double> compute_model_fdr(const std::vector<T>& data_in)
{
  using D = double;
  const std::size_t n = data_in.size();
  std::vector<double> fdr(n, std::numeric_limits<double>::quiet_NaN());

  if (n == 0) return fdr;

  auto is_nan_at = [&](std::size_t i) -> bool {
    if constexpr (std::is_floating_point<T>::value) return std::isnan(data_in[i]);
    else return false;
  };

  bool any_nan = false;
  if constexpr (std::is_floating_point<T>::value)
  {
    for (std::size_t i = 0; i < n; ++i) { if (is_nan_at(i)) { any_nan = true; break; } }
  }
  if (any_nan)
  {
    return fdr; // propagate
  }

  // argsort (stable)
  std::vector<std::size_t> order(n);
  for (std::size_t i = 0; i < n; ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t i, std::size_t j)
                   {
                     return static_cast<D>(data_in[i]) < static_cast<D>(data_in[j]);
                   });

  // sorted data
  std::vector<D> data_sorted(n);
  for (std::size_t i = 0; i < n; ++i) data_sorted[i] = static_cast<D>(data_in[order[i]]);

  // ranks for sorted data using 'max' tie method
  std::vector<double> ranks = RankData::rankdata<double>(data_sorted, RankData::Method::Max, RankData::NaNPolicy::Propagate);

  // cumulative sum of sorted data
  std::vector<D> cumsum(n);
  D acc = 0.0;
  for (std::size_t i = 0; i < n; ++i)
  {
    acc += data_sorted[i];
    cumsum[i] = acc;
  }

  // populate fdr in original order
  for (std::size_t i = 0; i < n; ++i)
  {
    const double r_d = ranks[i];
    if (std::isnan(r_d))
    {
      fdr[order[i]] = std::numeric_limits<double>::quiet_NaN();
      continue;
    }
    const std::size_t r_idx = static_cast<std::size_t>(static_cast<long long>(r_d) - 1);
    const double denom = r_d;
    const double numer = cumsum[std::min(r_idx, cumsum.size() - 1)];
    fdr[order[i]] = numer / denom;
  }

  return fdr;
}

/// Compute empirical p-values (bioconductor/qvalue empPvals translation)
template <class T>
inline std::vector<double> pemp(const std::vector<T>& stat, const std::vector<T>& stat0)
{
  using D = double;
  const std::size_t m = stat.size();
  const std::size_t m0 = stat0.size();
  if (m == 0 || m0 == 0) throw std::invalid_argument("pemp: input arrays must be non-empty");

  // concatenate
  std::vector<D> statc;
  statc.reserve(m + m0);
  for (auto&& v : stat) statc.push_back(static_cast<D>(v));
  for (auto&& v : stat0) statc.push_back(static_cast<D>(v));

  // v flags: True for stat, False for stat0
  std::vector<char> v;
  v.reserve(m + m0);
  for (std::size_t i = 0; i < m; ++i) v.push_back(1);
  for (std::size_t i = 0; i < m0; ++i) v.push_back(0);

  // argsort descending (stable)
  const std::size_t N = statc.size();
  std::vector<std::size_t> perm(N);
  for (std::size_t i = 0; i < N; ++i) perm[i] = i;
  std::stable_sort(perm.begin(), perm.end(), [&](std::size_t i, std::size_t j)
                   { return statc[i] > statc[j]; });

  // apply permutation to v
  std::vector<char> vperm(N);
  for (std::size_t i = 0; i < N; ++i) vperm[i] = v[perm[i]];

  // u: positions of True entries
  std::vector<std::size_t> u;
  for (std::size_t i = 0; i < N; ++i) if (vperm[i]) u.push_back(i);
  if (u.size() != m) throw std::runtime_error("pemp: internal error, unexpected u size");

  std::vector<double> p(m);
  for (std::size_t i = 0; i < m; ++i)
  {
    p[i] = static_cast<double>(static_cast<long long>(u[i]) - static_cast<long long>(i)) / static_cast<double>(m0);
  }

  // ranks: floor(rankdata(-stat)) - 1
  std::vector<D> neg_stat(m);
  for (std::size_t i = 0; i < m; ++i) neg_stat[i] = -static_cast<D>(stat[i]);
  std::vector<double> ranks = RankData::rankdata<double>(neg_stat, RankData::Method::Average, RankData::NaNPolicy::Propagate);

  std::vector<double> out(m);
  for (std::size_t i = 0; i < m; ++i)
  {
    double rf = std::floor(ranks[i]);
    std::size_t idx = static_cast<std::size_t>(static_cast<long long>(rf) - 1);
    if (idx >= p.size()) idx = p.size() - 1;
    out[i] = p[idx];
  }

  // enforce minimum 1/m0
  const double minp = 1.0 / static_cast<double>(m0);
  for (auto& vv : out) if (vv <= minp) vv = minp;

  return out;
}

// qvalue and pi0est are implemented in the .cpp file
OPENMS_DLLAPI std::vector<double> qvalue(const std::vector<double>& p_values, double pi0, bool pfdr = false);
OPENMS_DLLAPI Pi0Result pi0est(const std::vector<double>& p_values,
                 const std::vector<double>& lambda_ = std::vector<double>(),
                 const std::string& pi0_method = "smoother",
                 int smooth_df = 3,
                 bool smooth_log_pi0 = false);

/// Bandwidth selector using the "nrd0" (Silverman-ish) rule-of-thumb
OPENMS_DLLAPI double bw_nrd0(const std::vector<double>& x);

/// Bin data onto an equally spaced grid between xmin (inclusive) and xmax (exclusive)
/// Returns a vector of length nbins containing counts (or weighted counts if weights provided).
OPENMS_DLLAPI std::vector<double> linbin(const std::vector<double>& x, double xmin, double xmax, std::size_t nbins, const std::vector<double>* weights = nullptr);

/// Munro-packed real FFT (forward). Packs rfft(X,M)/M into length-M vector
/// format: [Re(Y_0..Y_{M/2}), Im(Y_1..Y_{M/2-1})]
OPENMS_DLLAPI std::vector<double> forrt(const std::vector<double>& X, std::size_t M = 0);

/// Inverse of `forrt`: reconstruct real signal (scaled by *M to invert forrt)
OPENMS_DLLAPI std::vector<double> revrt(const std::vector<double>& Xp, std::size_t M = 0);

/// Build Silverman FFT of Gaussian kernel on Munro frequency grid (Munro-packed)
/// See statsmodels.nonparametric.kdetools.silverman_transform for reference.
OPENMS_DLLAPI std::vector<double> silverman_kernel_fft(double bw, std::size_t M, double RANGE);

/// FFT-grid Gaussian KDE using Silverman/Munro FFT + linear binning.
/// Returns pair (density, grid) where grid is the M equally spaced points from a..b.
OPENMS_DLLAPI std::pair<std::vector<double>, std::vector<double>> grid_kde_fft(const std::vector<double>& x, double bw, std::size_t gridsize = 512, double cut = 3.0);

/// Evaluate KDE at query points using FFT-grid method + cubic spline interpolation
OPENMS_DLLAPI std::vector<double> kde_fft_eval(const std::vector<double>& x, double bw, std::size_t gridsize = 512, double cut = 3.0);

/// Estimate local FDR / posterior error probability (PEP) from p-values.
/// Supports "probit" and "logit" transforms and uses the FFT-grid
/// Silverman KDE implemented above. Mirrors the python/statsmodels
/// `lfdr` behavior used in qvalue-style local FDR estimation.
OPENMS_DLLAPI std::vector<double> lfdr(const std::vector<double>& p_values,
                                      double pi0,
                                      bool trunc = true,
                                      bool monotone = true,
                                      const std::string& transf = "probit",
                                      double adj = 1.5,
                                      double eps = 1e-8,
                                      std::size_t gridsize = 512,
                                      double cut = 3.0);

/// Compute tail probabilities under a normal distribution fitted to stat0
/// Returns P(X > stat_i) where X ~ N(mu, sigma^2) with mu/sigma estimated from stat0
OPENMS_DLLAPI std::vector<double> pnorm(const std::vector<double>& stat, const std::vector<double>& stat0);

} // namespace Math
} // namespace OpenMS
