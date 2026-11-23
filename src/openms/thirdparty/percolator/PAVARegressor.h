#pragma once

#include "MonotoneRegressor.h"
#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>

namespace OpenMS { namespace Internal { namespace Percolator {


class PAVARegressor final : public MonotoneRegressor {
public:
  using MonotoneRegressor::MonotoneRegressor;

  std::vector<double> fit_xy(const std::vector<double>& x,
                             const std::vector<double>& y,
                             double clip_lo = 0.0, double clip_hi = 1.0) override {
    if (x.size() != y.size()) {
      throw std::invalid_argument("PAVARegressor::fit_xy: x and y size mismatch");
    }
    if (y.empty()) {
      return {};
    }

    // Sort observations by x, run 1D isotonic fit, then map predictions back.
    // If y is expected to decrease with x, enforce monotonicity on descending x.
    std::vector<size_t> order(y.size());
    std::iota(order.begin(), order.end(), 0);
    const bool descending_x = params_.y_decreasing_in_x;
    std::stable_sort(order.begin(), order.end(),
                     [&](size_t lhs, size_t rhs) {
                       if (x[lhs] == x[rhs]) {
                         return lhs < rhs;
                       }
                       return descending_x ? (x[lhs] > x[rhs]) : (x[lhs] < x[rhs]);
                     });

    // Observations with the same x form one block from the start: the fit
    // is a function of x, so they get one value, the mean of their y if no
    // neighbouring block is merged into it. Fed one by one in input order,
    // a run of equal x could be split, and e.g. PSMs with the same score
    // got different PEPs depending on where the decoys stood among them.
    std::vector<double> block_sum;
    std::vector<int> block_count;
    for (size_t i = 0; i < order.size(); ++i) {
      if (i > 0 && x[order[i]] == x[order[i - 1]]) {
        block_sum.back() += y[order[i]];
        ++block_count.back();
      } else {
        block_sum.push_back(y[order[i]]);
        block_count.push_back(1);
      }
    }
    const auto fit_sorted = pava_impl(block_sum, block_count, clip_lo, clip_hi);

    std::vector<double> out(y.size());
    for (size_t i = 0; i < order.size(); ++i) {
      out[order[i]] = fit_sorted[i];
    }
    return out;
  }

  std::vector<double> fit_y(const std::vector<double>& y,
                            double clip_lo = 0.0, double clip_hi = 1.0) override {
    // no x: every observation is a block of its own
    return pava_impl(y, std::vector<int>(y.size(), 1), clip_lo, clip_hi);
  }

private:
  // Pool-adjacent-violators on blocks given by the sum and the number of
  // their observations, in the order in which the fit must not decrease;
  // returns one value per observation.
  std::vector<double> pava_impl(const std::vector<double>& sums,
    const std::vector<int>& counts,
    const double min_value,
    const double max_value) const;
};


inline std::vector<double> PAVARegressor::pava_impl(
    const std::vector<double>& sums,
    const std::vector<int>& counts,
    const double min_value,
    const double max_value
) const
{
    assert(sums.size() == counts.size());
    const int n = static_cast<int>(sums.size());
    if (n == 0) {
        return {};
    }

    // We'll store blocks on a stack
    struct MergeBlock {
        double sum;
        int    count;
        double avg;
    };

    std::vector<MergeBlock> stack;
    stack.reserve(n);

    // 1. Left to right
    for (int i = 0; i < n; ++i) {
        MergeBlock newBlock { sums[i], counts[i], sums[i] / counts[i] };
        stack.push_back(newBlock);

        // 2. Merge while there's a violation of non-decreasing
        while (stack.size() > 1) {
            auto &top    = stack.back();
            auto &secTop = stack[stack.size() - 2];
            double mergedSum   = secTop.sum + top.sum;
            int    mergedCount = secTop.count + top.count;
            double mergedAvg   = mergedSum / mergedCount;

            // if (( secTop.avg > top.avg) || ( mergedAvg < min_value ) || ( mergedAvg > max_value )) {
            if ( secTop.avg > top.avg) {
                stack.pop_back();
                stack.pop_back();
                stack.push_back({ mergedSum, mergedCount, mergedAvg });
            } else {
                break;
            }
        }
    }

    // 3. Expand final solution
    std::vector<double> result;
    result.reserve(std::accumulate(counts.begin(), counts.end(), size_t(0)));

    for (auto &b : stack) {
        double val = std::max( std::min(b.avg, max_value), min_value);
        for (int c = 0; c < b.count; ++c) {
            result.push_back(val);
        }
    }
    return result;
}

}}}  // namespace OpenMS::Internal::Percolator
