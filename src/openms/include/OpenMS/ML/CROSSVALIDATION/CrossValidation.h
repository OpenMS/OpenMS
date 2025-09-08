// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace OpenMS
{

/**
  @brief Lightweight K-fold / LOO cross-validation utilities and 1-D grid search.

  Provides:
    - `makeKFolds(n, K)`: deterministic round-robin fold assignment (LOO if K==n)
    - `gridSearch1D(...)`: evaluate a 1-D candidate grid via CV and pick the best

  Tie-breaking uses a tiny absolute tolerance and (optionally) prefers larger
  candidates (useful for smoother/regularized models).

  References:
    - Stone, M. (1974) Cross-Validatory Choice and Assessment of Statistical Predictions.
      J. Roy. Stat. Soc. B, 36(2):111–147.

  @see GridSearch

  @ingroup MachineLearning
*/
class CrossValidation
{
public:
  /**
    @brief Tie-breaking preference for equal (within tolerance) CV scores.

    - PreferLarger  : choose the larger candidate value on ties
    - PreferSmaller : choose the smaller candidate value on ties
    - PreferAny     : keep the first encountered (stable, no size preference)

    @ingroup MachineLearning
  */
  enum class CandidateTieBreak
  {
    PreferLarger,
    PreferSmaller,
    PreferAny
  };

  /**
    @brief Build @p K folds for indices [0, n).

    Deterministic round-robin assignment: fold(i) = i % K.
    For leave-one-out (LOO), use K = n.

    @param n Number of samples
    @param K Requested number of folds (clamped to [1, n])

    @exception Exception::InvalidValue if @p n == 0 or @p K == 0

    @ingroup MachineLearning
  */
  static std::vector<std::vector<Size>> makeKFolds(Size n, Size K)
  {
    if (n == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "n", String(n));
    }
    if (K == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "K", String(K));
    }
    if (K > n) K = n;

    std::vector<std::vector<Size>> folds(K);
    for (Size i = 0; i < n; ++i) folds[i % K].push_back(i);
    return folds;
  }

  /**
    @brief One-dimensional grid search with external cross-validation evaluation.

    Iterates candidates [@p cbegin, @p cend), calls @p train_eval(candidate, folds, abs_errs)
    to append absolute errors from all validation points, then scores them via
    @p score(abs_errs) (lower is better). Returns the best (candidate, score).

    Tie-breaking:
     - If |score - best_score| <= @p tie_tol, choose by @p prefer_larger (true → larger wins).

    @tparam CandIter  Random-access or forward iterator over candidate values
    @tparam TrainEval Callable of signature `void(const Cand&, const std::vector<std::vector<Size>>&, std::vector<double>&)`
    @tparam ScoreFn   Callable of signature `double(const std::vector<double>&)`

    @param cbegin         Begin iterator of candidate grid
    @param cend           End iterator of candidate grid
    @param folds          Fold index sets (e.g., from makeKFolds)
    @param train_eval     Callback: fit on train folds and append |error| for all held-out points
    @param score          Callback: convert accumulated errors to a scalar loss (lower is better)
    @param tie_tol        Absolute tolerance for tie detection (default: 1e-12)
    @param tie_break      Preference for ties (default: PreferLarger)

    @return (best_candidate, best_score)

    @exception Exception::InvalidRange if candidate range is empty

    @ingroup MachineLearning
  */
  template <typename CandIter, typename TrainEval, typename ScoreFn>
  static std::pair<typename std::iterator_traits<CandIter>::value_type, double>
  gridSearch1D(CandIter cbegin, CandIter cend,
               const std::vector<std::vector<Size>>& folds,
               TrainEval train_eval,
               ScoreFn score,
               double tie_tol = 1e-12,
               CandidateTieBreak tie_break = CandidateTieBreak::PreferLarger)
  {
    using CandT = typename std::iterator_traits<CandIter>::value_type;

    if (cbegin == cend)
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    CandT best_cand = *cbegin;
    double best_score = std::numeric_limits<double>::infinity();
    bool first = true;

    for (auto it = cbegin; it != cend; ++it)
    {
      const CandT cand = *it;

      std::vector<double> abs_errs;
      abs_errs.reserve(256); // grows as needed
      train_eval(cand, folds, abs_errs);

      const double s = score(abs_errs);

      // Prefer larger candidate on numerical ties (more stable smoothing, etc.)
      const bool better = (s < best_score - tie_tol);
      const bool tie    = (std::fabs(s - best_score) <= tie_tol);

      bool wins_on_tie = false;
      if (tie)
      {
        switch (tie_break)
        {
          case CandidateTieBreak::PreferLarger:  wins_on_tie = cand > best_cand; break;
          case CandidateTieBreak::PreferSmaller: wins_on_tie = cand < best_cand; break;
          case CandidateTieBreak::PreferAny:     wins_on_tie = false;            break;
        }
      }

      if (first || better || wins_on_tie)
      {
        best_cand  = cand;
        best_score = s;
        first = false;
      }
    }

    return {best_cand, best_score};
  }
};

} // namespace OpenMS
