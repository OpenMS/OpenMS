// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAMzIndex.h"

#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace OpenMS;

namespace
{
  void require(bool value, const char* message)
  {
    if (!value) throw std::runtime_error(message);
  }

  void check(const MSSpectrum& s, const OpenDIA::SpectrumMzIndex& index, double mz)
  {
    require(index.lowerBound(mz) == s.MZBegin(mz), "lowerBound differs from MZBegin");
  }
}

int main()
{
  std::mt19937_64 rng(7);
  for (int trial = 0; trial < 200; ++trial)
  {
    MSSpectrum s;
    std::uniform_real_distribution<double> mz(100.0, 2000.0);
    const int n = static_cast<int>(rng() % 3000);
    for (int i = 0; i < n; ++i) s.emplace_back(mz(rng), 1.0f);
    // Exact bucket edges and duplicates.
    if (trial % 3 == 0)
    {
      for (double e = 100.0; e < 2000.0; e += OpenDIA::SpectrumMzIndex::BUCKET_TH * 7) s.emplace_back(e, 2.0f);
      s.emplace_back(500.0, 1.0f);
      s.emplace_back(500.0, 3.0f);
    }
    s.sortByPosition();
    OpenDIA::SpectrumMzIndex index;
    index.build(s);
    for (int q = 0; q < 2000; ++q) check(s, index, std::uniform_real_distribution<double>(0.0, 2500.0)(rng));
    for (const auto& p : s)
    {
      check(s, index, p.getMZ());
      check(s, index, std::nextafter(p.getMZ(), 0.0));
      check(s, index, std::nextafter(p.getMZ(), 1e9));
    }
    for (double e = 0.0; e < 2500.0; e += OpenDIA::SpectrumMzIndex::BUCKET_TH)
    {
      check(s, index, e);
      check(s, index, std::nextafter(e, 0.0));
    }
    check(s, index, -1.0);
    check(s, index, 1e9);
  }
  MSSpectrum empty;
  OpenDIA::SpectrumMzIndex index;
  index.build(empty);
  require(index.lowerBound(500.0) == empty.MZBegin(500.0), "empty spectrum");
  std::cout << "OpenDIAMzIndex_test passed" << std::endl;
  return 0;
}
