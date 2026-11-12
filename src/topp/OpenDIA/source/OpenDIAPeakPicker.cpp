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

  std::vector<XicPeak> pickXicPeaks(const float* data, std::size_t scans, std::size_t fragments, const double* rts,
                                    const XicPickerOptions& options, std::size_t& peaks_found, std::vector<double>& scratch)
  {
    peaks_found = 0;
    std::vector<XicPeak> peaks;
    if (data == nullptr || rts == nullptr || scans == 0 || fragments == 0)
    {
      return peaks;
    }
    const auto value = [data, fragments](std::size_t s, std::size_t k) { return std::max(0.0f, data[s * fragments + k]); };

    // Raw consensus (sum of fragments) and its smoothed copy, which only drives maxima, valleys and boundaries.
    scratch.assign(3 * scans, 0.0);
    double* raw = scratch.data();
    double* y = raw + scans;
    double* tmp = y + scans;
    for (std::size_t s = 0; s < scans; ++s)
    {
      for (std::size_t k = 0; k < fragments; ++k)
      {
        raw[s] += value(s, k);
      }
    }
    std::copy(raw, raw + scans, y);
    for (std::size_t pass = 0; pass < options.smoothing_passes && scans > 2; ++pass)
    {
      tmp[0] = y[0];
      tmp[scans - 1] = y[scans - 1];
      for (std::size_t s = 1; s + 1 < scans; ++s)
      {
        tmp[s] = 0.25 * y[s - 1] + 0.5 * y[s] + 0.25 * y[s + 1];
      }
      std::copy(tmp, tmp + scans, y);
    }

    // Local maxima, then a left-to-right merge: a maximum joins its left neighbour while the valley between
    // them stays above valley_fraction of the lower one, so noise bumps on one elution profile do not split it.
    struct Region
    {
      std::size_t apex;
      std::size_t left;
    };
    std::vector<Region> regions;
    for (std::size_t s = 0; s < scans; ++s)
    {
      const bool rises = s == 0 || y[s] > y[s - 1];
      const bool falls = s + 1 == scans || y[s] >= y[s + 1];
      if (!(y[s] > 0.0) || !rises || !falls)
      {
        continue;
      }
      Region current{s, 0};
      while (!regions.empty())
      {
        const Region& top = regions.back();
        std::size_t valley = top.apex;
        for (std::size_t v = top.apex; v <= current.apex; ++v)
        {
          if (y[v] < y[valley])
          {
            valley = v;
          }
        }
        if (y[valley] > options.valley_fraction * std::min(y[top.apex], y[current.apex]))
        {
          current.apex = y[top.apex] >= y[current.apex] ? top.apex : current.apex;
          current.left = top.left;
          regions.pop_back();
          continue;
        }
        current.left = valley;
        break;
      }
      regions.push_back(current);
    }
    peaks_found = regions.size();

    peaks.reserve(regions.size());
    for (std::size_t r = 0; r < regions.size(); ++r)
    {
      const std::size_t valley_left = regions[r].left;
      const std::size_t valley_right = r + 1 < regions.size() ? regions[r + 1].left : scans - 1;
      const std::size_t top = regions[r].apex;
      const double floor_level = options.boundary_fraction * y[top];
      std::size_t left = top;
      while (left > valley_left && y[left] > floor_level)
      {
        --left;
      }
      std::size_t right = top;
      while (right < valley_right && y[right] > floor_level)
      {
        ++right;
      }

      XicPeak peak;
      peak.left = left;
      peak.right = right;
      peak.apex = static_cast<std::size_t>(std::max_element(raw + left, raw + right + 1) - raw);
      peak.apex_height = raw[peak.apex];
      if (!(peak.apex_height > 0.0))
      {
        continue;
      }
      for (std::size_t s = left; s < right; ++s)
      {
        peak.area += 0.5 * (raw[s] + raw[s + 1]) * (rts[s + 1] - rts[s]);
      }
      for (std::size_t s = left; s <= right; ++s)
      {
        peak.signal_scans += raw[s] > 0.0;
      }

      // FWHM of the raw consensus; a side that stays above half maximum ends at the peak boundary.
      const double half = 0.5 * peak.apex_height;
      std::size_t i = peak.apex;
      while (i > left && raw[i - 1] > half)
      {
        --i;
      }
      const double half_left = i > left ? rts[i - 1] + (half - raw[i - 1]) * (rts[i] - rts[i - 1]) / (raw[i] - raw[i - 1]) : rts[left];
      i = peak.apex;
      while (i < right && raw[i + 1] > half)
      {
        ++i;
      }
      const double half_right = i < right ? rts[i] + (raw[i] - half) * (rts[i + 1] - rts[i]) / (raw[i] - raw[i + 1]) : rts[right];
      peak.fwhm = half_right - half_left;

      // Co-elution over the peak, widened to at least correlation_half_window scans on each side of the apex.
      const std::size_t h = options.correlation_half_window;
      const std::size_t lo = std::min(left, peak.apex > h ? peak.apex - h : 0);
      const std::size_t hi = std::max(right, std::min(scans - 1, peak.apex + h));
      const double length = static_cast<double>(hi - lo + 1);
      for (std::size_t k = 0; k < fragments && hi - lo >= 2; ++k)
      {
        double sx = 0.0, so = 0.0, sxx = 0.0, soo = 0.0, sxo = 0.0;
        for (std::size_t s = lo; s <= hi; ++s)
        {
          const double x = value(s, k);
          const double o = raw[s] - x;
          sx += x;
          so += o;
          sxx += x * x;
          soo += o * o;
          sxo += x * o;
        }
        const double vx = sxx - sx * sx / length;
        const double vo = soo - so * so / length;
        if (vx > 0.0 && vo > 0.0)
        {
          peak.coelution += std::max(0.0, (sxo - sx * so / length) / std::sqrt(vx * vo));
        }
      }
      const std::size_t s_lo = peak.apex > 0 ? peak.apex - 1 : 0;
      const std::size_t s_hi = std::min(scans - 1, peak.apex + 1);
      for (std::size_t k = 0; k < fragments; ++k)
      {
        bool seen = false;
        for (std::size_t s = s_lo; s <= s_hi && !seen; ++s)
        {
          seen = value(s, k) > 0.0f;
        }
        peak.support += seen;
      }
      peaks.push_back(peak);
    }

    std::sort(peaks.begin(), peaks.end(), [](const XicPeak& lhs, const XicPeak& rhs)
    {
      if (lhs.coelution != rhs.coelution) return lhs.coelution > rhs.coelution;
      if (lhs.apex_height != rhs.apex_height) return lhs.apex_height > rhs.apex_height;
      return lhs.apex < rhs.apex;
    });
    if (peaks.size() > options.max_peaks)
    {
      peaks.resize(options.max_peaks);
    }
    return peaks;
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
