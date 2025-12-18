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
  @brief Compute model-based FDR estimates from posterior error probabilities.

  This function converts posterior error probabilities (PEPs) into false discovery
  rate (FDR) estimates using a model-based approach. The method follows the framework
  described in Kall et al. (2008).

  The FDR at rank i is calculated as:
    FDR(i) = (sum of PEPs up to rank i) / i

  This provides a direct relationship between PEPs and FDR without requiring
  permutation-based null hypothesis testing.

  Reference:
  Kall L, Storey JD, MacCoss MJ, Noble WS. (2008)
  "Posterior error probabilities and false discovery rates: two sides of the same coin."
  J Proteome Res. 7(1):40-4. doi: 10.1021/pr700739d

  @param data_in Vector of posterior error probabilities (PEPs), typically in range [0,1]
  @return Vector of FDR estimates corresponding to each input PEP value.
          NaN values in input are propagated to output.

  @note The function is templated to accept various numeric types but converts
        internally to double for computation.

  @ingroup MathFunctionsStatistics
*/
template <class T>
inline std::vector<double> computeModelFdr(const std::vector<T>& data_in)
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

/**
  @brief Compute empirical p-values from test statistics and null distribution.

  This function calculates empirical p-values by comparing observed test statistics
  against a null (permuted/randomized) distribution. The method is commonly used
  in permutation-based hypothesis testing.

  The empirical p-value for a statistic s is defined as:
    p(s) = (number of null statistics >= s) / (number of null statistics)

  This implementation follows the algorithm from the Bioconductor qvalue package
  (Bass et al. 2015) with rank-based tie handling for improved numerical stability.

  Reference:
  Bass AJ, Dabney A, Robinson D. (2015) "qvalue: Q-value estimation for false
  discovery rate control." R package. Bioconductor.
  http://bioconductor.org/packages/qvalue/

  @param stat Vector of observed test statistics (typically from real data). For example, this would be your target scores.
  @param stat0 Vector of null test statistics (from permutations or simulations). For example, this would be your decoy scores.
  @return Vector of empirical p-values, one for each element in @p stat.
          Values are bounded below by 1/m0 where m0 is the size of the null distribution.

  @exception Exception::InvalidArgument if either input vector is empty

  @note Assumes larger statistic values indicate stronger signals (one-sided test).
        For two-sided tests, use absolute values of statistics.

  @ingroup MathFunctionsStatistics
*/
template <class T>
inline std::vector<double> pEmp(const std::vector<T>& stat, const std::vector<T>& stat0)
{
  using D = double;
  const std::size_t m = stat.size();
  const std::size_t m0 = stat0.size();
  if (m == 0 || m0 == 0) throw std::invalid_argument("pEmp: input arrays must be non-empty");

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

/**
  @brief Calculate q-values (FDR-adjusted p-values) for multiple testing correction.

  Converts p-values to q-values using the Storey-Tibshirani method. A q-value is
  the minimum FDR at which a test is called significant. This provides a more
  interpretable measure than p-values when performing multiple hypothesis tests.

  The q-value for test i is defined as:
    q(i) = min_{t>=p(i)} FDR(t)

  where FDR(t) is the expected false discovery rate when calling all tests with
  p-value <= t significant.

  Reference:
  Storey JD. (2002) "A direct approach to false discovery rates."
  J. R. Statist. Soc. B. 64(3): 479-498. doi: 10.1111/1467-9868.00346

  @param p_values Vector of p-values to be adjusted (must be in range [0,1])
  @param pi0 Proportion of true null hypotheses (typically estimated by pi0est())
  @param pfdr If true, compute positive FDR; if false (default), compute regular FDR
  @return Vector of q-values corresponding to input p-values. NaN values are propagated.

  @exception Exception::InvalidArgument if p-values are outside [0,1] or pi0 is outside [0,1]

  @see pi0Est() for estimating pi0
  @ingroup MathFunctionsStatistics
*/
OPENMS_DLLAPI std::vector<double> qValue(const std::vector<double>& p_values, double pi0, bool pfdr = false);

/**
  @brief Estimate the proportion of true null hypotheses (pi0).

  This function estimates pi0, the proportion of truly null hypotheses among all
  tested hypotheses, which is a key parameter for FDR control. The estimation uses
  the method of Storey & Tibshirani (2003) with optional smoothing.

  The algorithm works by:
  1. Computing pi0 estimates at multiple lambda threshold values
  2. Optionally smoothing these estimates using a cubic spline
  3. Evaluating the smoothed curve at the maximum lambda value

  When @p smooth_log_pi0 is true, smoothing is performed in log-space to ensure
  positive pi0 estimates and better numerical stability, following the approach
  used in PyProphet (Reiter et al. 2011).

  References:
  - Storey JD and Tibshirani R. (2003) "Statistical significance for genome-wide experiments."
    PNAS 100: 9440-9445. doi: 10.1073/pnas.1530509100
  - Reiter L et al. (2011) "mProphet: automated data processing and statistical validation
    for large-scale SRM experiments." Nat Methods 8(5):430-5. doi: 10.1038/nmeth.1584

  @param p_values Vector of p-values from hypothesis tests (range [0,1])
  @param lambda_ Vector of lambda threshold values for estimation. If empty, defaults
                 to seq(0, 0.90, 0.05). Lambda values represent p-value thresholds.
  @param pi0_method Method for pi0 estimation: "smoother" (default) uses spline smoothing,
                    "bootstrap" uses bootstrap resampling 
  @param smooth_df Degrees of freedom for smoothing spline (default: 3 for cubic spline)
  @param smooth_log_pi0 If true, perform smoothing in log-space for improved stability
                        and to guarantee positive estimates (default: false)

  @return Pi0Result structure containing:
          - pi0: estimated proportion (range [0,1])
          - pi0_lambda: vector of pi0 estimates at each lambda
          - lambda_: lambda values used
          - pi0_smooth: whether smoothing succeeded (false means conservative minimum used)

  @note If smoothing fails or produces invalid results, the function returns the minimum
        pi0 estimate as a conservative fallback.

  @see Pi0Result for details on return structure
  @ingroup MathFunctionsStatistics
*/
OPENMS_DLLAPI Pi0Result pi0Est(const std::vector<double>& p_values,
                const std::vector<double>& lambda_ = std::vector<double>(),
                const std::string& pi0_method = "smoother",
                int smooth_df = 3,
                bool smooth_log_pi0 = false);

/**
  @brief Estimate local false discovery rate (local FDR) or posterior error probability (PEP).

  Computes local FDR values from p-values using kernel density estimation and the
  two-component mixture model framework. Local FDR represents the probability that
  a specific hypothesis is null given its test statistic, providing a more granular
  assessment than global FDR (q-values).

  The method transforms p-values to a more suitable scale (probit or logit), estimates
  the null and alternative density components using KDE, then computes:
    lfdr(z) = pi0 * f0(z) / f(z)

  where f0 is the null density, f is the mixture density, and z is the transformed value.

  The implementation follows the approach used in the R/Bioconductor qvalue package
  and Python statsmodels, using FFT-based KDE (Silverman's method) for efficiency.

  References:
  Efron B. (2004) "Large-scale simultaneous hypothesis testing: the choice of a null
  hypothesis." J Am Stat Assoc 99(465):96-104. DOI: 10.1198/016214504000000089

  Storey JD and Tibshirani R. (2003) "Statistical significance for genomewide studies."
  PNAS 100(16):9440-9445. DOI: 10.1073/pnas.1530509100

  @param p_values Vector of p-values from hypothesis tests
  @param pi0 Estimated proportion of true null hypotheses (from pi0est())
  @param trunc If true, truncate lfdr values to [0,1] range (default true)
  @param monotone If true, enforce monotonicity constraint (default true)
  @param transf Transformation to apply: "probit" (inverse normal) or "logit" (log-odds)
  @param adj Bandwidth adjustment factor (multiplied by automatic bandwidth selection)
  @param eps Small constant added to density estimates to avoid division by zero
  @param gridsize Number of FFT grid points for KDE (default 512)
  @param cut Grid extension factor in units of bandwidth (default 3.0)
  @return Vector of local FDR values, one for each input p-value

  @note The "probit" transform is generally preferred for p-values as it better stabilizes
        the variance across the support.

  @see pi0est() for estimating pi0
  @see grid_kde_fft() for the underlying KDE implementation

  @ingroup MathFunctionsStatistics
*/
OPENMS_DLLAPI std::vector<double> lfdr(const std::vector<double>& p_values,
                                      double pi0,
                                      bool trunc = true,
                                      bool monotone = true,
                                      const std::string& transf = "probit",
                                      double adj = 1.5,
                                      double eps = 1e-8,
                                      std::size_t gridsize = 512,
                                      double cut = 3.0);

/**
  @brief Compute tail probabilities under a fitted normal distribution.

  Estimates parameters (mean and standard deviation) from a null distribution (@p stat0),
  then computes upper tail probabilities P(X > stat_i) for each value in @p stat,
  where X ~ N(mu, sigma^2) with mu and sigma estimated from @p stat0.

  This is useful for computing empirical p-values when you have both target scores
  (@p stat) and a sample from the null distribution (@p stat0).

  @param stat Vector of test statistics for which to compute tail probabilities
  @param stat0 Vector of statistics from the null distribution used to estimate N(mu, sigma^2)
  @return Vector of upper tail probabilities (p-values), one for each value in @p stat

  @note Uses robust estimators (median, MAD) if the null distribution is suspected
        to be contaminated

  @ingroup MathFunctionsStatistics
*/
OPENMS_DLLAPI std::vector<double> pNorm(const std::vector<double>& stat, const std::vector<double>& stat0);

} // namespace Math
} // namespace OpenMS
