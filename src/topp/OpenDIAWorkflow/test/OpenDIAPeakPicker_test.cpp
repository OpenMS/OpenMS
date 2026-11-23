// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAPeakPicker.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace OpenMS::OpenDIA;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
    {
      throw std::runtime_error(message);
    }
  }

  constexpr std::size_t F = 6;

  // points x F interleaved (point-major, like the calibration traces), Gaussian peaks at centre with sigma.
  std::vector<float> gaussianGroup(std::size_t points, double centre, double sigma, double height = 1000.0)
  {
    std::vector<float> data(points * F, 0.0f);
    const std::array<double, F> ratio{1.0, 0.8, 0.6, 0.5, 0.3, 0.2};
    for (std::size_t s = 0; s < points; ++s)
    {
      const double z = (static_cast<double>(s) - centre) / sigma;
      for (std::size_t f = 0; f < F; ++f)
      {
        data[s * F + f] += static_cast<float>(height * ratio[f] * std::exp(-0.5 * z * z));
      }
    }
    return data;
  }

  PickedPeak pick(const std::vector<float>& data, std::size_t points, std::size_t hint, PeakPickerOptions options = {})
  {
    std::vector<double> scratch;
    return pickPeakGroup(data.data(), F, points, 1, F, hint, options, scratch);
  }
}

int main()
{
  try
  {
    // Isolated Gaussian: apex at the centre, FWHM 2.3548 sigma, symmetric FLOOR boundaries, all fragments agree.
    {
      const double sigma = 3.0;
      const auto data = gaussianGroup(61, 30.0, sigma);
      const auto peak = pick(data, 61, 28);
      require(peak.found && peak.apex == 30, "apex not found at the Gaussian centre");
      const double fwhm = peak.half_right - peak.half_left;
      // One [1 2 1]/4 pass adds variance 0.5 point^2.
      const double expected = 2.0 * std::sqrt(2.0 * std::log(2.0)) * std::sqrt(sigma * sigma + 0.5);
      require(std::fabs(fwhm - expected) < 0.15, "interpolated FWHM is off");
      require(peak.left_reason == PeakBoundary::FLOOR && peak.right_reason == PeakBoundary::FLOOR, "isolated peak should stop at the floor");
      require(peak.apex - peak.left == peak.right - peak.apex, "isolated peak boundaries are not symmetric");
      // 5% of the apex lies at sqrt(2 ln 20) ~ 2.45 sigma.
      require(peak.right - peak.apex >= 7 && peak.right - peak.apex <= 8, "floor boundary at the wrong distance");
      require(peak.fragments_with_signal == F && peak.agreeing_fragments == F, "all fragments should agree");
    }
    // Two peaks: boundaries stop at the valley between them and do not swallow the neighbour.
    {
      auto data = gaussianGroup(61, 20.0, 2.0);
      const auto second = gaussianGroup(61, 28.0, 2.0, 600.0);
      for (std::size_t i = 0; i < data.size(); ++i)
      {
        data[i] += second[i];
      }
      const auto peak = pick(data, 61, 20);
      require(peak.apex == 20 && peak.right_reason == PeakBoundary::VALLEY, "neighbouring peak not split at the valley");
      require(peak.right >= 23 && peak.right <= 25, "valley boundary at the wrong position");
      const auto other = pick(data, 61, 28);
      require(other.apex == 28 && other.left_reason == PeakBoundary::VALLEY && other.left == peak.right, "second peak not split at the same valley");
    }
    // Interference on one fragment at another RT: the consensus apex holds and that fragment does not agree.
    {
      auto data = gaussianGroup(61, 30.0, 3.0);
      // Small enough that the consensus keeps falling, so the interfered points stay inside the peak.
      for (std::size_t s = 34; s < 36; ++s)
      {
        data[s * F + 5] += 400.0f;
      }
      const auto peak = pick(data, 61, 30);
      require(peak.found && peak.apex == 30 && peak.right > 35, "interference changed the consensus peak");
      require(peak.agreeing_fragments == F - 1, "interfered fragment should disagree with the apex");
    }
    // Edge and empty traces.
    {
      const auto data = gaussianGroup(20, 1.0, 3.0);
      const auto peak = pick(data, 20, 1);
      require(peak.found && peak.left == 0 && peak.left_reason == PeakBoundary::EDGE, "edge boundary not reported");
      const std::vector<float> zero(20 * F, 0.0f);
      require(!pick(zero, 20, 10).found, "all-zero trace reported a peak");
      std::vector<double> scratch;
      require(!pickPeakGroup(nullptr, F, 20, 1, F, 0, {}, scratch).found, "null data reported a peak");
    }
    // Fragment-major layout (like the 32-bin tensor) gives the same result as point-major.
    {
      const auto data = gaussianGroup(32, 15.0, 1.5);
      std::vector<float> fragment_major(32 * F);
      for (std::size_t s = 0; s < 32; ++s)
      {
        for (std::size_t f = 0; f < F; ++f)
        {
          fragment_major[f * 32 + s] = data[s * F + f];
        }
      }
      std::vector<double> scratch;
      const auto a = pickPeakGroup(data.data(), F, 32, 1, F, 15, {}, scratch);
      const auto b = pickPeakGroup(fragment_major.data(), F, 32, 32, 1, 15, {}, scratch);
      require(a.left == b.left && a.right == b.right && a.half_left == b.half_left && a.half_right == b.half_right, "layouts disagree");
    }
    // interpolateAt.
    {
      const std::vector<double> rt{10.0, 11.0, 13.0};
      require(interpolateAt(rt, 0.5) == 10.5 && interpolateAt(rt, 1.5) == 12.0 && interpolateAt(rt, -1.0) == 10.0 && interpolateAt(rt, 9.0) == 13.0,
              "interpolateAt is wrong");
    }
    // pickXicPeaks: two co-eluting peaks plus a single-fragment interference spike, on a non-uniform RT grid.
    {
      const std::size_t n = 120;
      std::vector<double> rts(n);
      for (std::size_t s = 0; s < n; ++s)
      {
        rts[s] = 100.0 + 1.6 * static_cast<double>(s) + (s % 3 == 0 ? 0.1 : 0.0);
      }
      auto data = gaussianGroup(n, 30.0, 1.5, 1000.0); // strong co-eluting peak, FWHM ~ 2.355 * 1.5 * 1.6 s
      const auto second = gaussianGroup(n, 80.0, 2.0, 300.0);
      for (std::size_t i = 0; i < data.size(); ++i)
      {
        data[i] += second[i];
      }
      data[55 * F + 0] += 5000.0f; // one fragment, one scan: not co-eluting
      // Small noise bump on the first peak's tail must not split it.
      data[33 * F + 1] += 60.0f;
      XicPickerOptions options;
      std::size_t found = 0;
      std::vector<double> scratch;
      const auto peaks = pickXicPeaks(data.data(), n, F, rts.data(), options, found, scratch);
      require(found == 3, "pickXicPeaks: expected exactly three separated peaks");
      require(peaks.size() == 3, "pickXicPeaks: all peaks should be reported below max_peaks");
      require(peaks[0].apex == 30 || peaks[1].apex == 30, "pickXicPeaks: the strong co-eluting peak is not in the top two");
      require(peaks[2].apex == 55, "pickXicPeaks: the single-fragment spike should rank last by co-elution");
      const XicPeak& strong = peaks[0].apex == 30 ? peaks[0] : peaks[1];
      const double expected_fwhm = 2.0 * std::sqrt(2.0 * std::log(2.0)) * 1.5 * 1.6;
      require(std::abs(strong.fwhm - expected_fwhm) < 0.6, "pickXicPeaks: FWHM of the strong peak is off");
      require(strong.support == F, "pickXicPeaks: every fragment supports the strong peak");
      require(strong.coelution > 5.5, "pickXicPeaks: co-eluting fragments should correlate");
      require(strong.left > 20 && strong.right < 40, "pickXicPeaks: boundaries should follow the trace, not a constant");
      require(peaks[2].coelution < 0.5, "pickXicPeaks: an isolated spike has no co-elution");

      options.max_peaks = 1;
      const auto best = pickXicPeaks(data.data(), n, F, rts.data(), options, found, scratch);
      require(best.size() == 1 && found == 3, "pickXicPeaks: max_peaks must cut after ranking");

      const std::vector<float> zeros(10 * F, 0.0f);
      const auto none = pickXicPeaks(zeros.data(), 10, F, rts.data(), options, found, scratch);
      require(none.empty() && found == 0, "pickXicPeaks: an empty trace has no peak");
    }

    // min_profile_scans: two fragments nonzero on one common scan correlate with r = 1 over any window. With the
    // default (2 scans) that coincidence scores 0 and a real but imperfect 2-fragment elution ranks first.
    {
      const std::size_t n = 60;
      std::vector<double> rts(n);
      for (std::size_t s = 0; s < n; ++s)
      {
        rts[s] = 1.6 * static_cast<double>(s);
      }
      std::vector<float> data(n * F, 0.0f);
      data[10 * F + 0] = 900.0f;
      data[10 * F + 1] = 700.0f; // coincidence: two fragments, one scan
      // A real but noisy two-fragment elution: each fragment correlates with the other at r < 1, so its
      // co-elution is just below the coincidence's 2.
      for (std::size_t s = 35; s <= 45; ++s)
      {
        const double z = (static_cast<double>(s) - 40.0) / 1.5;
        for (std::size_t f = 0; f < 2; ++f)
        {
          data[s * F + f] = static_cast<float>(300.0 * (1.0 - 0.2 * f) * std::exp(-0.5 * z * z) + ((s + f) % 2 ? 40.0 : 0.0));
        }
      }
      std::vector<double> scratch;
      std::size_t found = 0;
      XicPickerOptions old_rule;
      old_rule.min_profile_scans = 1;
      const auto before = pickXicPeaks(data.data(), n, F, rts.data(), old_rule, found, scratch);
      require(!before.empty() && before[0].apex == 10 && std::abs(before[0].coelution - 2.0) < 1e-9,
              "min_profile_scans=1: the one-scan coincidence should score exactly 2 and rank first");
      const auto after = pickXicPeaks(data.data(), n, F, rts.data(), XicPickerOptions{}, found, scratch);
      require(!after.empty() && after[0].apex == 40, "min_profile_scans=2: the real elution should rank first");
      for (const auto& p : after)
      {
        require(p.apex != 10 || p.coelution == 0.0, "min_profile_scans=2: the coincidence should have no co-elution");
      }
    }

    // boundary_half_widths: per-side cap in that side's own half widths. On an exponentially tailing peak the right
    // side keeps more scans than the left, and the cap scales with the peak's width (same k, 2x wider peak).
    {
      const auto tailing = [](std::size_t n, double centre, double sigma, double tau)
      {
        std::vector<float> data(n * F, 0.0f);
        const std::array<double, F> ratio{1.0, 0.8, 0.6, 0.5, 0.3, 0.2};
        for (std::size_t s = 0; s < n; ++s)
        {
          const double x = static_cast<double>(s) - centre;
          // Gaussian rise, exponential decay after the centre (a simple tailing shape).
          const double y = x < 0.0 ? std::exp(-0.5 * x * x / (sigma * sigma)) : std::exp(-x / tau);
          for (std::size_t f = 0; f < F; ++f)
          {
            data[s * F + f] = static_cast<float>(1000.0 * ratio[f] * y);
          }
        }
        return data;
      };
      const std::size_t n = 200;
      std::vector<double> rts(n);
      for (std::size_t s = 0; s < n; ++s)
      {
        rts[s] = 1.6 * static_cast<double>(s);
      }
      std::vector<double> scratch;
      std::size_t found = 0;
      XicPickerOptions capped;
      capped.boundary_half_widths = 2.5;
      for (double scale : {1.0, 2.0})
      {
        const auto data = tailing(n, 60.0, 1.5 * scale, 4.0 * scale);
        const auto free_peaks = pickXicPeaks(data.data(), n, F, rts.data(), XicPickerOptions{}, found, scratch);
        const auto cap_peaks = pickXicPeaks(data.data(), n, F, rts.data(), capped, found, scratch);
        require(free_peaks.size() == 1 && cap_peaks.size() == 1, "boundary cap: one peak expected");
        const XicPeak& f = free_peaks[0];
        const XicPeak& c = cap_peaks[0];
        require(c.apex == 60 && c.left >= f.left && c.right <= f.right, "boundary cap: the cap only narrows the peak");
        const double hl = c.half_width_left, hr = c.half_width_right;
        require(std::abs(hl - 1.1774 * 1.5 * scale) < 0.35 && std::abs(hr - std::log(2.0) * 4.0 * scale) < 0.35,
                "boundary cap: half widths at half maximum are off");
        // Each side ends at ceil(k x its own half width), or earlier where the uncapped boundary already stopped.
        const auto reach = [](double h) { return static_cast<std::size_t>(std::ceil(2.5 * std::max(h, 0.5))); };
        require(c.right == std::min(f.right, c.apex + reach(hr)) && c.left == std::max(f.left, c.apex - reach(hl)),
                "boundary cap: each side should end at ceil(k x its own half width)");
        require(c.right < f.right, "boundary cap: the tail side should be cut");
        require(c.right - c.apex > c.apex - c.left, "boundary cap: a tailing peak keeps a longer right side");
        require(c.area < f.area && c.area > 0.9 * f.area, "boundary cap: the cap cuts only the far tail");
        require(c.coelution == f.coelution, "boundary cap: co-elution (peak choice) must not depend on the cap");
      }
    }

    // estimateBoundaryHalfWidths: peaks whose fragment pattern holds out to a known distance and becomes random
    // beyond it. The estimate should land near that distance and follow it when the distance changes.
    {
      const auto make = [](double signal_extent, unsigned seed)
      {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<float> uni(0.0f, 1.0f);
        std::vector<BoundaryTrace> traces;
        const std::array<double, F> ratio{1.0, 0.8, 0.6, 0.5, 0.3, 0.2};
        for (int t = 0; t < 200; ++t)
        {
          BoundaryTrace trace;
          trace.fragments = F;
          const std::size_t n = 120;
          const double sigma = 1.0 + 0.5 * (t % 4); // widths 1 to 2.5 scans (FWHM 2.4 to 5.9)
          const double hw = 1.1774 * sigma;
          for (std::size_t s = 0; s < n; ++s)
          {
            trace.rts.push_back(1.6 * static_cast<double>(s));
          }
          trace.hint = 60;
          trace.data.assign(n * F, 0.0f);
          for (std::size_t s = 0; s < n; ++s)
          {
            const double x = static_cast<double>(s) - 60.0;
            const double u = std::abs(x) / hw;
            const double y = 1000.0 * std::exp(-0.5 * x * x / (sigma * sigma));
            for (std::size_t f = 0; f < F; ++f)
            {
              // Within signal_extent half widths: library pattern; beyond: a random pattern of the same height
              // (an interfering signal), plus low random background everywhere.
              const double pattern = u <= signal_extent ? ratio[f] : uni(rng);
              trace.data[s * F + f] = static_cast<float>(y * pattern + 2.0 * uni(rng));
            }
          }
          traces.push_back(std::move(trace));
        }
        return traces;
      };
      const auto narrow = estimateBoundaryHalfWidths(make(2.0, 3), XicPickerOptions{});
      const auto wide = estimateBoundaryHalfWidths(make(3.5, 4), XicPickerOptions{});
      require(narrow.peaks == 200 && wide.peaks == 200, "estimateBoundaryHalfWidths: every trace should be profiled");
      require(narrow.half_widths > 1.5 && narrow.half_widths < 2.75, "estimateBoundaryHalfWidths: estimate off for extent 2");
      require(wide.half_widths > 3.0 && wide.half_widths < 4.25, "estimateBoundaryHalfWidths: estimate off for extent 3.5");
      require(wide.half_widths > narrow.half_widths + 0.75, "estimateBoundaryHalfWidths: the estimate should follow the extent");
      require(estimateBoundaryHalfWidths({}, XicPickerOptions{}).half_widths == 0.0, "estimateBoundaryHalfWidths: no traces, no estimate");
      std::printf("estimateBoundaryHalfWidths: extent 2.0 -> %.3f, extent 3.5 -> %.3f\n", narrow.half_widths, wide.half_widths);
    }

    // Speed on noisy 32-point windows (informational; the check is only that it is far below the CNN cost).
    {
      std::mt19937 rng(7);
      std::normal_distribution<float> noise(0.0f, 30.0f);
      auto data = gaussianGroup(32, 16.0, 1.5);
      for (auto& value : data)
      {
        value = std::max(0.0f, value + noise(rng));
      }
      std::vector<double> scratch;
      const std::size_t n = 1'000'000;
      std::size_t found = 0;
      const auto start = std::chrono::steady_clock::now();
      for (std::size_t i = 0; i < n; ++i)
      {
        found += pickPeakGroup(data.data(), F, 32, 1, F, 16, {}, scratch).found;
      }
      const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
      std::printf("pickPeakGroup: %.0f ns per 6x32 window (%zu found)\n", 1e9 * seconds / static_cast<double>(n), found);
      require(found == n && seconds < 20.0, "picker is unexpectedly slow");
    }
  }
  catch (const std::exception& error)
  {
    std::cerr << "OpenDIAPeakPicker_test failed: " << error.what() << '\n';
    return 1;
  }
  std::cout << "OpenDIAPeakPicker_test passed\n";
  return 0;
}
