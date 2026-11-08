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
