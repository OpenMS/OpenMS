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
      return run_rt.front() + (n > 1 ? end_slope(span, 0) : 0.0) * (x - library_rt.front());
    }
    if (x >= library_rt.back())
    {
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
}
