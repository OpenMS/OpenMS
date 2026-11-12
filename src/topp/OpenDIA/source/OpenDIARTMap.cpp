// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIARTMap.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    /// Median of a nonempty range; even sizes average the two middle values.
    double median_(std::vector<double> values)
    {
      const std::size_t mid = values.size() / 2;
      std::nth_element(values.begin(), values.begin() + mid, values.end());
      const double upper = values[mid];
      if (values.size() % 2 == 1)
      {
        return upper;
      }
      return (*std::max_element(values.begin(), values.begin() + mid) + upper) / 2.0;
    }

    /// Knots from @p bins equal-count bins of anchors sorted by (library RT, run RT).
    RTMap binnedKnots_(const std::vector<std::pair<double, double>>& points, std::size_t bins)
    {
      RTMap map;
      const std::size_t n = points.size();
      for (std::size_t bin = 0; bin < bins; ++bin)
      {
        const std::size_t begin = bin * n / bins;
        const std::size_t end = (bin + 1) * n / bins;
        if (begin == end)
        {
          continue;
        }
        std::vector<double> x, y;
        for (std::size_t i = begin; i < end; ++i)
        {
          x.push_back(points[i].first);
          y.push_back(points[i].second);
        }
        map.library_rt.push_back(median_(std::move(x)));
        map.run_rt.push_back(median_(std::move(y)));
      }
      return map;
    }
  }

  double RTMap::operator()(double x) const
  {
    if (library_rt.empty())
    {
      throw std::logic_error("Empty RT map");
    }
    const std::size_t n = library_rt.size();
    // Outside the knots: linear extrapolation on the line through the end knot and the knot kEndSpan positions inward
    // (or the farthest knot when there are fewer). The slope is clamped at zero so the map stays nondecreasing. Holding
    // the end value constant instead puts every precursor beyond the knots at the end knot's run RT; on long-tailed
    // library RT ranges that is off by minutes at the gradient ends.
    const auto end_slope = [&](std::size_t end, std::size_t inner) {
      const double dx = library_rt[end] - library_rt[inner];
      return dx == 0.0 ? 0.0 : std::max((run_rt[end] - run_rt[inner]) / dx, 0.0);
    };
    const std::size_t span = std::min(kEndSpan, n - 1);
    if (x <= library_rt.front())
    {
      if (std::isfinite(low_slope))
      {
        return run_rt.front() + low_slope * (x - library_rt.front());
      }
      return run_rt.front() + (n > 1 ? end_slope(span, 0) : 0.0) * (x - library_rt.front());
    }
    if (x >= library_rt.back())
    {
      if (std::isfinite(high_slope))
      {
        return run_rt.back() + high_slope * (x - library_rt.back());
      }
      return run_rt.back() + (n > 1 ? end_slope(n - 1, n - 1 - span) : 0.0) * (x - library_rt.back());
    }
    // First knot strictly above x; the previous knot is at or below x.
    const std::size_t i = std::upper_bound(library_rt.begin(), library_rt.end(), x) - library_rt.begin();
    const double t = (x - library_rt[i - 1]) / (library_rt[i] - library_rt[i - 1]);
    return run_rt[i - 1] + t * (run_rt[i] - run_rt[i - 1]);
  }

  RTMap fitRTMap(std::vector<std::pair<double, double>> anchors, std::size_t bins)
  {
    if (anchors.size() < 20)
    {
      throw std::invalid_argument("RT calibration needs at least 20 anchors");
    }
    if (bins == 0)
    {
      throw std::invalid_argument("RT calibration needs at least one bin");
    }
    for (const auto& anchor : anchors)
    {
      if (!std::isfinite(anchor.first) || !std::isfinite(anchor.second))
      {
        throw std::invalid_argument("RT anchors must be finite");
      }
    }
    std::sort(anchors.begin(), anchors.end());
    RTMap map;
    // One fit plus up to three rejection refits.
    for (int round = 0; round < 4; ++round)
    {
      map = binnedKnots_(anchors, bins);
      std::vector<double> residuals;
      residuals.reserve(anchors.size());
      for (const auto& anchor : anchors)
      {
        residuals.push_back(std::fabs(anchor.second - map(anchor.first)));
      }
      map.sigma = 1.4826 * median_(residuals);
      map.anchors = anchors.size();
      const double cut = std::max(3.0 * map.sigma, 3.0);
      std::vector<std::pair<double, double>> kept;
      for (std::size_t i = 0; i < anchors.size(); ++i)
      {
        if (residuals[i] <= cut)
        {
          kept.push_back(anchors[i]);
        }
      }
      if (kept.size() == anchors.size())
      {
        break;
      }
      if (kept.size() < 20)
      {
        throw std::invalid_argument("RT calibration rejected too many anchors");
      }
      anchors = std::move(kept);
    }
    return map;
  }

  namespace
  {
    /// Theil-Sen slope (median pairwise slope) of points with distinct x; 0 if none.
    double theilSen_(const std::vector<std::pair<double, double>>& points)
    {
      std::vector<double> slopes;
      for (std::size_t i = 0; i < points.size(); ++i)
      {
        for (std::size_t j = i + 1; j < points.size(); ++j)
        {
          const double dx = points[j].first - points[i].first;
          if (dx > 0.0)
          {
            slopes.push_back((points[j].second - points[i].second) / dx);
          }
        }
      }
      return slopes.empty() ? 0.0 : median_(std::move(slopes));
    }

    /// L1 isotonic fit of anchors sorted by (library RT, run RT): one knot per block (median x, block median y).
    RTMap isotonicKnots_(const std::vector<std::pair<double, double>>& points)
    {
      struct Block
      {
        std::vector<double> x, y;
        double level = 0.0;
      };
      std::vector<Block> blocks;
      for (const auto& [x, y] : points)
      {
        blocks.push_back({{x}, {y}, y});
        // Pool adjacent violators: merged blocks take the median of their run RTs (the L1 optimum).
        while (blocks.size() > 1 && blocks[blocks.size() - 2].level > blocks.back().level)
        {
          Block last = std::move(blocks.back());
          blocks.pop_back();
          Block& previous = blocks.back();
          previous.x.insert(previous.x.end(), last.x.begin(), last.x.end());
          previous.y.insert(previous.y.end(), last.y.begin(), last.y.end());
          previous.level = median_(previous.y);
        }
      }
      RTMap map;
      for (const auto& block : blocks)
      {
        const double x = median_(block.x);
        if (!map.library_rt.empty() && x <= map.library_rt.back())
        {
          // Ties in library RT across blocks: keep one knot, at the higher level (still nondecreasing).
          map.run_rt.back() = std::max(map.run_rt.back(), block.level);
          continue;
        }
        map.library_rt.push_back(x);
        map.run_rt.push_back(std::max(block.level, map.run_rt.empty() ? block.level : map.run_rt.back()));
      }
      const std::size_t m = std::max<std::size_t>(3, static_cast<std::size_t>(std::lround(std::sqrt(static_cast<double>(points.size())))));
      const std::size_t k = std::min(m, points.size());
      map.low_slope = std::max(theilSen_({points.begin(), points.begin() + static_cast<std::ptrdiff_t>(k)}), 0.0);
      map.high_slope = std::max(theilSen_({points.end() - static_cast<std::ptrdiff_t>(k), points.end()}), 0.0);
      return map;
    }
  }

  RTMap fitMonotoneRTMap(std::vector<std::pair<double, double>> anchors)
  {
    if (anchors.size() < 20)
    {
      throw std::invalid_argument("RT calibration needs at least 20 anchors");
    }
    for (const auto& anchor : anchors)
    {
      if (!std::isfinite(anchor.first) || !std::isfinite(anchor.second))
      {
        throw std::invalid_argument("RT anchors must be finite");
      }
    }
    std::sort(anchors.begin(), anchors.end());
    constexpr std::size_t FOLDS = 10;
    RTMap map;
    for (int round = 0; round < 4; ++round)
    {
      // Cross-validated residuals: fold f holds out every FOLDS-th anchor of the sorted list.
      std::vector<double> residuals(anchors.size());
      for (std::size_t fold = 0; fold < FOLDS; ++fold)
      {
        std::vector<std::pair<double, double>> train;
        for (std::size_t i = 0; i < anchors.size(); ++i)
        {
          if (i % FOLDS != fold)
          {
            train.push_back(anchors[i]);
          }
        }
        const RTMap fold_map = isotonicKnots_(train);
        for (std::size_t i = fold; i < anchors.size(); i += FOLDS)
        {
          residuals[i] = std::fabs(anchors[i].second - fold_map(anchors[i].first));
        }
      }
      map = isotonicKnots_(anchors);
      map.sigma = 1.4826 * median_(residuals);
      map.anchors = anchors.size();
      // Reject on the held-out residual: an isotonic fit interpolates a lone anchor (in-sample residual 0), so
      // only its neighbours can contradict it. A group of mutually consistent anchors at a real slope change
      // predicts its own held-out members and is kept.
      const double cut = std::max(3.0 * map.sigma, 3.0);
      std::vector<std::pair<double, double>> kept;
      for (std::size_t i = 0; i < anchors.size(); ++i)
      {
        if (residuals[i] <= cut)
        {
          kept.push_back(anchors[i]);
        }
      }
      if (kept.size() == anchors.size())
      {
        break;
      }
      if (kept.size() < 20)
      {
        throw std::invalid_argument("RT calibration rejected too many anchors");
      }
      anchors = std::move(kept);
    }
    return map;
  }
}
