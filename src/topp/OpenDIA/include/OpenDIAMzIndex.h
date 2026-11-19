// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Coarse m/z bucket index of one m/z-sorted spectrum, so the first point at or above an m/z is found by a
    lower bound inside one bucket instead of over the whole spectrum.

    lowerBound(mz) returns exactly spectrum.MZBegin(mz) (the first point with m/z >= mz): offsets_[b] is the first point
    at or above the start of bucket b, so the answer for an m/z inside bucket b lies in [offsets_[b], offsets_[b + 1]].
    Build it after the spectrum is sorted and use it only while the spectrum is unchanged.
  */
  class SpectrumMzIndex
  {
  public:
    static constexpr double BUCKET_TH = 0.25;

    void build(const MSSpectrum& spectrum)
    {
      spectrum_ = &spectrum;
      offsets_.clear();
      if (spectrum.empty()) return;
      origin_ = std::floor(spectrum.front().getMZ() / BUCKET_TH) * BUCKET_TH;
      const double top = spectrum.back().getMZ();
      if (!std::isfinite(origin_) || !std::isfinite(top) || top < origin_) return;
      const auto buckets = static_cast<std::size_t>((top - origin_) / BUCKET_TH) + 2;
      offsets_.resize(buckets + 1);
      std::size_t i = 0;
      for (std::size_t b = 0; b <= buckets; ++b)
      {
        const double start = origin_ + static_cast<double>(b) * BUCKET_TH;
        while (i < spectrum.size() && spectrum[i].getMZ() < start) ++i;
        offsets_[b] = static_cast<std::uint32_t>(i);
      }
    }

    MSSpectrum::ConstIterator lowerBound(double mz) const
    {
      const MSSpectrum& s = *spectrum_;
      if (offsets_.empty() || !(mz >= origin_)) return s.MZBegin(mz);
      const auto b = static_cast<std::size_t>((mz - origin_) / BUCKET_TH);
      if (b + 1 >= offsets_.size()) return s.end();
      // Guard against rounding at a bucket edge: widen by one bucket on each side.
      const std::size_t lo = offsets_[b == 0 ? 0 : b - 1];
      const std::size_t hi = offsets_[std::min(offsets_.size() - 1, b + 2)];
      return std::lower_bound(s.begin() + static_cast<std::ptrdiff_t>(lo), s.begin() + static_cast<std::ptrdiff_t>(hi), mz,
                              [](const Peak1D& p, double v) { return p.getMZ() < v; });
    }

  private:
    const MSSpectrum* spectrum_ = nullptr;
    double origin_ = 0.0;
    std::vector<std::uint32_t> offsets_;
  };
}
