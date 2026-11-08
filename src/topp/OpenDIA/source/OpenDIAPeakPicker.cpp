// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAPeakPicker.h>

#include <algorithm>
#include <cmath>

namespace OpenMS::OpenDIA
{
  namespace
  {
    // Fractional index where the smoothed trace crosses level between inner (above) and the outward step.
    double crossing(const double* y, std::size_t inner, std::size_t outer, double level)
    {
      const double a = y[inner], b = y[outer];
      const double t = a > b ? std::clamp((a - level) / (a - b), 0.0, 1.0) : 0.0;
      return static_cast<double>(inner) + t * (static_cast<double>(outer) - static_cast<double>(inner));
    }
  }

  PickedPeak pickPeakGroup(const float* data, std::size_t fragments, std::size_t points, std::size_t fragment_stride,
                           std::size_t point_stride, std::size_t hint, const PeakPickerOptions& options,
                           std::vector<double>& scratch)
  {
    PickedPeak peak;
    if (data == nullptr || fragments == 0 || points == 0)
    {
      return peak;
    }
    // Consensus trace: plain sum of the fragments (the brightest fragments dominate, as in a TIC of the group).
    scratch.assign(2 * points, 0.0);
    double* sum = scratch.data();
    double* tmp = scratch.data() + points;
    for (std::size_t f = 0; f < fragments; ++f)
    {
      const float* trace = data + f * fragment_stride;
      for (std::size_t s = 0; s < points; ++s)
      {
        sum[s] += std::max(0.0f, trace[s * point_stride]);
      }
    }
    for (std::size_t pass = 0; pass < options.smoothing_passes && points > 2; ++pass)
    {
      tmp[0] = sum[0];
      tmp[points - 1] = sum[points - 1];
      for (std::size_t s = 1; s + 1 < points; ++s)
      {
        tmp[s] = 0.25 * sum[s - 1] + 0.5 * sum[s] + 0.25 * sum[s + 1];
      }
      std::swap(sum, tmp);
    }
    const double* y = sum;

    hint = std::min(hint, points - 1);
    const std::size_t lo = hint > options.apex_search ? hint - options.apex_search : 0;
    const std::size_t hi = std::min(points - 1, hint + options.apex_search);
    std::size_t apex = hint;
    for (std::size_t s = lo; s <= hi; ++s)
    {
      if (y[s] > y[apex])
      {
        apex = s;
      }
    }
    if (!(y[apex] > 0.0))
    {
      return peak;
    }
    peak.found = true;
    peak.apex = apex;
    peak.apex_height = y[apex];
    const double floor_level = options.boundary_fraction * y[apex];

    // Walk outwards while falling and above the floor.
    std::size_t left = apex;
    peak.left_reason = PeakBoundary::EDGE;
    while (left > 0 && apex - left < options.max_half_span)
    {
      if (y[left - 1] > y[left])
      {
        peak.left_reason = PeakBoundary::VALLEY;
        break;
      }
      --left;
      if (y[left] <= floor_level)
      {
        peak.left_reason = PeakBoundary::FLOOR;
        break;
      }
    }
    std::size_t right = apex;
    peak.right_reason = PeakBoundary::EDGE;
    while (right + 1 < points && right - apex < options.max_half_span)
    {
      if (y[right + 1] > y[right])
      {
        peak.right_reason = PeakBoundary::VALLEY;
        break;
      }
      ++right;
      if (y[right] <= floor_level)
      {
        peak.right_reason = PeakBoundary::FLOOR;
        break;
      }
    }
    peak.left = left;
    peak.right = right;

    // Interpolated half-maximum crossings inside [left, right]; a side that never drops below half stays at its boundary.
    const double half = 0.5 * y[apex];
    std::size_t s = apex;
    while (s > left && y[s - 1] > half)
    {
      --s;
    }
    peak.half_left = s > left ? crossing(y, s, s - 1, half) : static_cast<double>(left);
    s = apex;
    while (s < right && y[s + 1] > half)
    {
      ++s;
    }
    peak.half_right = s < right ? crossing(y, s, s + 1, half) : static_cast<double>(right);

    // Fragment agreement: each fragment's raw maximum inside the peak should sit at the consensus apex.
    for (std::size_t f = 0; f < fragments; ++f)
    {
      const float* trace = data + f * fragment_stride;
      std::size_t best = left;
      float best_value = 0.0f;
      for (std::size_t p = left; p <= right; ++p)
      {
        if (trace[p * point_stride] > best_value)
        {
          best_value = trace[p * point_stride];
          best = p;
        }
      }
      if (best_value > 0.0f)
      {
        ++peak.fragments_with_signal;
        const std::size_t distance = best > apex ? best - apex : apex - best;
        peak.agreeing_fragments += distance <= options.fragment_apex_tolerance;
      }
    }
    return peak;
  }

  double interpolateAt(const std::vector<double>& values, double index)
  {
    if (values.empty())
    {
      return 0.0;
    }
    if (!(index > 0.0))
    {
      return values.front();
    }
    const double last = static_cast<double>(values.size() - 1);
    if (index >= last)
    {
      return values.back();
    }
    const auto i = static_cast<std::size_t>(index);
    const double t = index - static_cast<double>(i);
    return values[i] + t * (values[i + 1] - values[i]);
  }

  const char* peakBoundaryName(PeakBoundary reason)
  {
    switch (reason)
    {
      case PeakBoundary::VALLEY: return "valley";
      case PeakBoundary::FLOOR: return "floor";
      default: return "edge";
    }
  }
}
