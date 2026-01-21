// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/RankData.h>

namespace OpenMS
{
namespace Math
{
/**
  @brief Result of pi0 estimation for multiple testing correction.

  This structure contains the estimated proportion of true null hypotheses (pi0),
  which is a key parameter in FDR estimation. The pi0 value represents the proportion
  of hypotheses that are truly null (i.e., no real effect) among all tested hypotheses.

  The estimation uses the method described in:
  Storey JD and Tibshirani R. (2003) "Statistical significance for genome-wide experiments."
  PNAS 100: 9440-9445. doi: 10.1073/pnas.1530509100

  Fields:
    - pi0          : Estimated proportion of true null hypotheses (range: 0-1).
                     A value of 1.0 means all hypotheses are assumed null (most conservative).
    - pi0_lambda   : Vector of pi0 estimates at each lambda threshold value.
                     Used for diagnostic purposes and smoothing.
    - lambda_      : Vector of lambda threshold values used for estimation.
                     Typically ranges from 0 to 0.95 in steps.
    - pi0_smooth   : Boolean indicating whether smoothing (via spline fitting) was
                     successfully applied. If false, the minimum pi0 value is used
                     as a conservative fallback.

  @see pi0est() for the estimation procedure
*/
struct OPENMS_DLLAPI Pi0Result
{
  double pi0 = 1.0;
  std::vector<double> pi0_lambda;
  std::vector<double> lambda_;
  bool pi0_smooth = false;
};

/**
  @brief Statistical functions for multiple testing correction.

  Provides FDR estimation, q-value calculation, and local FDR computation
  based on the Storey-Tibshirani method. This struct contains static methods
  and enum classes for the various multiple testing correction procedures.

  References:
  - Storey JD and Tibshirani R. (2003) "Statistical significance for genome-wide experiments."
    PNAS 100: 9440-9445. doi: 10.1073/pnas.1530509100
  - Storey JD. (2002) "A direct approach to false discovery rates."
    J. R. Statist. Soc. B. 64(3): 479-498. doi: 10.1111/1467-9868.00346

  @ingroup MathFunctionsStatistics
*/
struct OPENMS_DLLAPI MultipleTesting
{
  /// Method for estimating proportion of true null hypotheses (pi0)
  enum class Pi0Method
  {
    Smoother,   ///< Spline smoothing method (default)
    Bootstrap   ///< Bootstrap resampling method
  };

  /// Transformation for local FDR estimation
  enum class LfdrTransform
  {
    Probit,     ///< Inverse normal (Gaussian) transformation (default)
    Logit       ///< Log-odds transformation
  };

  /// Convert Pi0Method enum to string representation
  static std::string pi0MethodToString(Pi0Method m);

  /// Convert string to Pi0Method enum (throws if invalid)
  static Pi0Method toPi0Method(const std::string& s);

  /// Convert LfdrTransform enum to string representation
  static std::string lfdrTransformToString(LfdrTransform t);

  /// Convert string to LfdrTransform enum (throws if invalid)
  static LfdrTransform toLfdrTransform(const std::string& s);

  /**
    @brief Calculate q-values (FDR-adjusted p-values) for multiple testing correction.

    Converts p-values to q-values using the Storey-Tibshirani method. A q-value is
    the minimum FDR at which a test is called significant.

    @param p_values Vector of p-values to be adjusted (must be in range [0,1])
    @param pi0 Proportion of true null hypotheses (typically estimated by pi0Est())
    @param pfdr If true, compute positive FDR; if false (default), compute regular FDR
    @return Vector of q-values corresponding to input p-values. NaN values are propagated.
  */
  static std::vector<double> qValue(const std::vector<double>& p_values, double pi0, bool pfdr = false);

  /**
    @brief Estimate the proportion of true null hypotheses (pi0).

    @param p_values Vector of p-values from hypothesis tests (range [0,1])
    @param lambda_ Vector of lambda threshold values for estimation. If empty, defaults to seq(0.05, 0.95, 0.05).
    @param method Method for pi0 estimation: Smoother (default) or Bootstrap
    @param smooth_df Degrees of freedom for smoothing spline (default: 3)
    @param smooth_log_pi0 If true, perform smoothing in log-space (default: false)
    @return Pi0Result structure with estimation results
  */
  static Pi0Result pi0Est(const std::vector<double>& p_values,
                          const std::vector<double>& lambda_ = std::vector<double>(),
                          Pi0Method method = Pi0Method::Smoother,
                          int smooth_df = 3,
                          bool smooth_log_pi0 = false);

  /**
    @brief Estimate local false discovery rate (local FDR) from p-values.

    @param p_values Vector of p-values from hypothesis tests
    @param pi0 Estimated proportion of true null hypotheses
    @param trunc If true, truncate lfdr values to [0,1] range (default true)
    @param monotone If true, enforce monotonicity constraint (default true)
    @param transf Transformation to apply: Probit (default) or Logit
    @param adj Bandwidth adjustment factor (default 1.5)
    @param eps Small constant to avoid division by zero (default 1e-8)
    @param gridsize Number of FFT grid points for KDE (default 512)
    @param cut Grid extension factor in units of bandwidth (default 3.0)
    @return Vector of local FDR values
  */
  static std::vector<double> lfdr(const std::vector<double>& p_values,
                                  double pi0,
                                  bool trunc = true,
                                  bool monotone = true,
                                  LfdrTransform transf = LfdrTransform::Probit,
                                  double adj = 1.5,
                                  double eps = 1e-8,
                                  std::size_t gridsize = 512,
                                  double cut = 3.0);

  /**
    @brief Compute tail probabilities under a fitted normal distribution.

    @param stat Vector of test statistics
    @param stat0 Vector of null statistics used to estimate N(mu, sigma^2)
    @return Vector of upper tail probabilities (p-values)
  */
  static std::vector<double> pNorm(const std::vector<double>& stat, const std::vector<double>& stat0);

  /**
    @brief Compute model-based FDR estimates from posterior error probabilities.

    @param data_in Vector of posterior error probabilities (PEPs)
    @return Vector of FDR estimates
  */
  template <class T>
  static std::vector<double> computeModelFDR(const std::vector<T>& data_in);

  /**
    @brief Compute empirical p-values from test statistics and null distribution.

    @param stat Vector of observed test statistics
    @param stat0 Vector of null test statistics
    @return Vector of empirical p-values
  */
  template <class T>
  static std::vector<double> pEmp(const std::vector<T>& stat, const std::vector<T>& stat0);
};

// Template implementation for MultipleTesting::computeModelFDR
template <class T>
inline std::vector<double> MultipleTesting::computeModelFDR(const std::vector<T>& data_in)
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
    const std::size_t r_idx = static_cast<std::size_t>(static_cast<Int64>(r_d) - 1);
    const double denom = r_d;
    const double numer = cumsum[std::min(r_idx, cumsum.size() - 1)];
    fdr[order[i]] = numer / denom;
  }

  return fdr;
}

// Template implementation for MultipleTesting::pEmp
template <class T>
inline std::vector<double> MultipleTesting::pEmp(const std::vector<T>& stat, const std::vector<T>& stat0)
{
  using D = double;
  const std::size_t m = stat.size();
  const std::size_t m0 = stat0.size();
  if (m == 0 || m0 == 0) throw std::invalid_argument("pEmp: input arrays must be non-empty");

  // concatenate
  std::vector<D> statc;
  statc.reserve(m + m0);
  for (auto v : stat) statc.push_back(static_cast<D>(v));
  for (auto v : stat0) statc.push_back(static_cast<D>(v));

  // v flags: True for stat, False for stat0
  std::vector<char> v;
  v.reserve(m + m0);
  v.insert(v.end(), m, 1);
  v.insert(v.end(), m0, 0);

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
    p[i] = static_cast<double>(static_cast<Int64>(u[i]) - static_cast<Int64>(i)) / static_cast<double>(m0);
  }

  // ranks: floor(rankdata(-stat)) - 1
  std::vector<D> neg_stat(m);
  for (std::size_t i = 0; i < m; ++i) neg_stat[i] = -static_cast<D>(stat[i]);
  std::vector<double> ranks = RankData::rankdata<double>(neg_stat, RankData::Method::Average, RankData::NaNPolicy::Propagate);

  std::vector<double> out(m);
  for (std::size_t i = 0; i < m; ++i)
  {
    double rf = std::floor(ranks[i]);
    std::size_t idx = static_cast<std::size_t>(static_cast<Int64>(rf) - 1);
    if (idx >= p.size()) idx = p.size() - 1;
    out[i] = p[idx];
  }

  // enforce minimum 1/m0
  const double minp = 1.0 / static_cast<double>(m0);
  for (auto& vv : out) if (vv <= minp) vv = minp;

  return out;
}

// -------------------------------------------------------------------------
// Backward-compatible free function wrappers
// -------------------------------------------------------------------------

/// @brief Backward-compatible wrapper for MultipleTesting::qValue
inline std::vector<double> qValue(const std::vector<double>& p_values, double pi0, bool pfdr = false)
{
  return MultipleTesting::qValue(p_values, pi0, pfdr);
}

/// @brief Backward-compatible wrapper for MultipleTesting::pi0Est (string-based API)
OPENMS_DLLAPI Pi0Result pi0Est(const std::vector<double>& p_values,
                               const std::vector<double>& lambda_ = std::vector<double>(),
                               const std::string& pi0_method = "smoother",
                               int smooth_df = 3,
                               bool smooth_log_pi0 = false);

/// @brief Backward-compatible wrapper for MultipleTesting::lfdr (string-based API)
OPENMS_DLLAPI std::vector<double> lfdr(const std::vector<double>& p_values,
                                       double pi0,
                                       bool trunc = true,
                                       bool monotone = true,
                                       const std::string& transf = "probit",
                                       double adj = 1.5,
                                       double eps = 1e-8,
                                       std::size_t gridsize = 512,
                                       double cut = 3.0);

/// @brief Backward-compatible wrapper for MultipleTesting::pNorm
inline std::vector<double> pNorm(const std::vector<double>& stat, const std::vector<double>& stat0)
{
  return MultipleTesting::pNorm(stat, stat0);
}

/// @brief Backward-compatible wrapper for MultipleTesting::computeModelFDR
template <class T>
inline std::vector<double> computeModelFDR(const std::vector<T>& data_in)
{
  return MultipleTesting::computeModelFDR(data_in);
}

/// @brief Backward-compatible wrapper for MultipleTesting::pEmp
template <class T>
inline std::vector<double> pEmp(const std::vector<T>& stat, const std::vector<T>& stat0)
{
  return MultipleTesting::pEmp(stat, stat0);
}

} // namespace Math
} // namespace OpenMS
