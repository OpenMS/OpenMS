// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace OpenMS {
namespace Math {

  /**
    @brief Rank items (1-based) with SciPy-like tie handling.

    Replicates scipy.stats.rankdata(a, method=..., nan_policy=...), flattened/1D semantics.
    Tie methods:
      - average: average of ranks for ties
      - min: minimum rank for ties (aka "competition")
      - max: maximum rank for ties
      - dense: like "min" but next unique value increases rank by 1
      - ordinal: break ties by original order (stable), all ranks distinct

    NaN handling:
      - propagate: if any NaN in input, every output is NaN (matches SciPy 1.16.2 docs)
      - omit: ignore NaNs when ranking; NaN positions in output remain NaN
      - raise: throw std::invalid_argument if any NaN is present

    Ranks start at 1, matching SciPy.

    @tparam T Numeric type (float/double/integers supported; integers are internally promoted).
  */
  struct RankData
  {
    enum class Method
    {
      Average,
      Min,
      Max,
      Dense,
      Ordinal
    };

    enum class NaNPolicy
    {
      Propagate,
      Omit,
      Raise
    };

    template <class T>
    static std::vector<double> rankdata(const std::vector<T>& a,
                                        Method method = Method::Average,
                                        NaNPolicy nan_policy = NaNPolicy::Propagate)
    {
      using D = double;
      const std::size_t n = a.size();
      std::vector<double> ranks(n, std::numeric_limits<double>::quiet_NaN());

      if (n == 0) return ranks; // empty -> empty

      // Detect NaNs (only meaningful for floating types)
      auto is_nan_at = [&](std::size_t i) -> bool {
        if constexpr (std::is_floating_point<T>::value) {
          return std::isnan(a[i]);
        } else {
          (void)i; // suppress unused parameter warning
          return false;
        }
      };

      bool any_nan = false;
      if constexpr (std::is_floating_point<T>::value)
      {
        for (std::size_t i = 0; i < n; ++i) { if (is_nan_at(i)) { any_nan = true; break; } }
      }

      // Handle nan_policy
      if (any_nan)
      {
        switch (nan_policy)
        {
          case NaNPolicy::Propagate:
            // All-NaN output (already initialized as NaN)
            return ranks;
          case NaNPolicy::Raise:
            throw std::invalid_argument("RankData::rankdata: NaNs present but nan_policy=Raise.");
          case NaNPolicy::Omit:
            // Continue; we'll rank only non-NaN elements and keep NaN at NaN positions.
            break;
        }
      }

      // Collect indices to consider (omit NaNs if requested)
      std::vector<std::size_t> idx;
      idx.reserve(n);
      for (std::size_t i = 0; i < n; ++i)
      {
        if (nan_policy == NaNPolicy::Omit && is_nan_at(i)) continue;
        idx.push_back(i);
      }

      if (idx.empty())
      {
        // All NaN with 'omit' -> return NaN vector (already NaN)
        return ranks;
      }

      // Stable sort indices by value ascending to mimic SciPy's stable behavior
      // (important for 'ordinal' and tie handling)
      std::stable_sort(idx.begin(), idx.end(),
                       [&](std::size_t i, std::size_t j)
                       {
                         // Comparison uses promoted double to unify behavior across T
                         const D vi = static_cast<D>(a[i]);
                         const D vj = static_cast<D>(a[j]);
                         return vi < vj;
                       });

      // Helper to assign a constant rank to a tied block [lo, hi) in the sorted index array
      auto assign_block_rank = [&](std::size_t lo, std::size_t hi, double value)
      {
        for (std::size_t k = lo; k < hi; ++k)
        {
          ranks[idx[k]] = value;
        }
      };

      // Main logic per tie method
      switch (method)
      {
        case Method::Ordinal:
        {
          // Distinct ranks, breaking ties by original order due to stable sort.
          // Rank is position in sorted order + 1
          for (std::size_t pos = 0; pos < idx.size(); ++pos)
          {
            ranks[idx[pos]] = static_cast<double>(pos + 1);
          }
          break;
        }

        case Method::Dense:
        {
          // For each block of equal values, assign a compact rank 1..m
          std::size_t lo = 0;
          std::size_t dense_rank = 1;
          while (lo < idx.size())
          {
            std::size_t hi = lo + 1;
            const D v = static_cast<D>(a[idx[lo]]);
            while (hi < idx.size() && static_cast<D>(a[idx[hi]]) == v) ++hi;

            assign_block_rank(lo, hi, static_cast<double>(dense_rank));
            ++dense_rank;
            lo = hi;
          }
          break;
        }

        case Method::Min:
        case Method::Max:
        case Method::Average:
        {
          // For each block of equal values, compute min/max/average of standard ranks.
          // "Standard" rank of sorted position p is (p + 1) in 1-based indexing.
          std::size_t lo = 0;
          while (lo < idx.size())
          {
            std::size_t hi = lo + 1;
            const D v = static_cast<D>(a[idx[lo]]);
            while (hi < idx.size() && static_cast<D>(a[idx[hi]]) == v) ++hi;

            // 1-based ranks for this block span [lo+1, hi]
            const double r_min = static_cast<double>(lo + 1);
            const double r_max = static_cast<double>(hi);
            double r_value = 0.0;
            if (method == Method::Min)      r_value = r_min;
            else if (method == Method::Max) r_value = r_max;
            else                             r_value = 0.5 * (r_min + r_max); // average

            assign_block_rank(lo, hi, r_value);
            lo = hi;
          }
          break;
        }
      }

      return ranks;
    }

    // Concrete forwarders for Python bindings (avoid C++ templates on the Python side)
    static std::vector<double> rankdata_double(const std::vector<double>& a,
                                               Method m = Method::Average,
                                               NaNPolicy p = NaNPolicy::Propagate)
    { return rankdata<double>(a, m, p); }

    static std::vector<double> rankdata_float(const std::vector<float>& a,
                                              Method m = Method::Average,
                                              NaNPolicy p = NaNPolicy::Propagate)
    { return rankdata<float>(a, m, p); }

    static std::vector<double> rankdata_int(const std::vector<int>& a,
                                            Method m = Method::Average,
                                            NaNPolicy p = NaNPolicy::Propagate)
    { return rankdata<int>(a, m, p); }
  };

  /// Free function overloads (by-value) for double input
  inline std::vector<double>
  rankdata_double(std::vector<double> a,
                  RankData::Method m = RankData::Method::Average,
                  RankData::NaNPolicy p = RankData::NaNPolicy::Propagate)
  {
    return RankData::rankdata<double>(a, m, p);
  }

  /// Free function overloads (by-value) for float input
  inline std::vector<double>
  rankdata_float(std::vector<float> a,
                 RankData::Method m = RankData::Method::Average,
                 RankData::NaNPolicy p = RankData::NaNPolicy::Propagate)
  {
    return RankData::rankdata<float>(a, m, p);
  }

  /// Free function overloads (by-value) for int input
  inline std::vector<double>
  rankdata_int(std::vector<int> a,
               RankData::Method m = RankData::Method::Average,
               RankData::NaNPolicy p = RankData::NaNPolicy::Propagate)
  {
    return RankData::rankdata<int>(a, m, p);
  }

} // namespace Math
} // namespace OpenMS
