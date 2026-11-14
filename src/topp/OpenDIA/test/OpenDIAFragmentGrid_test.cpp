// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAFragmentGrid.h>

#include <OpenMS/METADATA/Precursor.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

using namespace OpenMS;
using namespace OpenMS::OpenDIA;

namespace
{
  int failures = 0;
  void check(bool ok, const std::string& what)
  {
    if (!ok)
    {
      ++failures;
      std::cerr << "FAIL: " << what << '\n';
    }
  }

  constexpr std::size_t K = FragmentGrid::FRAGMENTS;
  constexpr std::size_t B = FragmentGrid::BINS;

  /// Grids, library intensities, search flags and the 9 reference values (6 evidence + d_max, d6_max, d_lib with noise
  /// 250) computed by the reference numpy formulas (features12.py, deficit_features.py).
  struct Case
  {
    std::array<double, K * B> grid;
    std::array<double, K> library;
    std::array<bool, K> search;
    std::array<double, 9> expected;
  };
  const Case cases[] = {
#include "OpenDIAFragmentGrid_cases.inc"
  };

  MSSpectrum spectrum(double rt, double lo, double hi, const std::vector<std::pair<double, float>>& peaks)
  {
    MSSpectrum s;
    s.setMSLevel(2);
    s.setRT(rt);
    OpenMS::Precursor p;
    p.setMZ(0.5 * (lo + hi));
    p.setIsolationWindowLowerOffset(0.5 * (hi - lo));
    p.setIsolationWindowUpperOffset(0.5 * (hi - lo));
    s.setPrecursors({p});
    for (const auto& [mz, intensity] : peaks)
    {
      Peak1D peak;
      peak.setMZ(mz);
      peak.setIntensity(intensity);
      s.push_back(peak);
    }
    s.sortByPosition();
    return s;
  }
}

int main()
{
  // 1. Feature formulas against the numpy reference.
  for (std::size_t c = 0; c < std::size(cases); ++c)
  {
    const auto& k = cases[c];
    const auto e = FragmentGrid::evidence(k.grid.data(), k.library, k.search);
    std::array<double, K> area{};
    for (std::size_t f = 0; f < K; ++f)
      for (std::size_t b = FragmentGrid::CORE_BEGIN; b < FragmentGrid::CORE_END; ++b) area[f] += k.grid[f * B + b];
    const auto d = FragmentGrid::deficit(area, k.library, k.search, 250.0);
    for (std::size_t j = 0; j < 6; ++j) check(std::abs(e[j] - k.expected[j]) < 1e-5, "case " + std::to_string(c) + " evidence " + FragmentGrid::evidenceNames()[j]);
    for (std::size_t j = 0; j < 3; ++j) check(std::abs(d[j] - k.expected[6 + j]) < 1e-5, "case " + std::to_string(c) + " deficit " + FragmentGrid::deficitNames()[j]);
  }

  // 2. Collector against a brute-force fill: random spectra, two lanes, windows that open and close over the run.
  {
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<FragmentGridCollector::Precursor> precursors(6);
    for (std::size_t p = 0; p < precursors.size(); ++p)
    {
      precursors[p].mz = p < 3 ? 450.0 + p : 700.0 + p;
      for (std::size_t s = 0; s < K; ++s)
      {
        precursors[p].extraction_mz[s] = s < 10 ? 300.0 + 50.0 * s + 3.0 * p : 0.0;
        precursors[p].library[s] = s < 10 ? 1.0 / (1.0 + s) : 0.0;
        precursors[p].search[s] = s < 6;
      }
    }
    std::vector<InterferenceResolver::Window> windows;
    for (std::uint32_t p = 0; p < precursors.size(); ++p)
    {
      for (double apex : {40.0 + 7.0 * p, 95.0 + 3.0 * p})
      {
        InterferenceResolver::Window w;
        w.precursor = p;
        w.rt_apex = apex;
        windows.push_back(w);
      }
    }
    std::vector<MSSpectrum> run;
    for (double rt = 0.0; rt < 140.0; rt += 0.7)
    {
      for (const auto& [lo, hi] : {std::pair{440.0, 460.0}, std::pair{690.0, 710.0}})
      {
        std::vector<std::pair<double, float>> peaks;
        for (std::size_t p = 0; p < precursors.size(); ++p)
          for (std::size_t s = 0; s < 10; ++s)
            if (u(rng) < 0.6) peaks.emplace_back(precursors[p].extraction_mz[s] * (1.0 + (u(rng) - 0.5) * 30e-6), static_cast<float>(1000.0 * u(rng)));
        run.push_back(spectrum(rt, lo, hi, peaks));
      }
    }
    // Apex scans: one spectrum exactly at the first apex of precursor 0.
    const double ppm = 10.0, half = 15.0;
    std::vector<double> brute(windows.size() * K * B, 0.0);
    for (const auto& s : run)
    {
      const double lo = s.getPrecursors()[0].getMZ() - s.getPrecursors()[0].getIsolationWindowLowerOffset();
      const double hi = s.getPrecursors()[0].getMZ() + s.getPrecursors()[0].getIsolationWindowUpperOffset();
      for (std::size_t w = 0; w < windows.size(); ++w)
      {
        const auto& p = precursors[windows[w].precursor];
        if (p.mz < lo || p.mz > hi) continue;
        double d = s.getRT() - windows[w].rt_apex;
        if (!(std::abs(d) < half)) continue;
        if (std::abs(d) < 1e-3) d = 0.0;
        const auto bin = static_cast<std::size_t>(std::clamp(static_cast<long>(std::floor((d + half) / (2 * half) * B)), 0L, 31L));
        for (std::size_t f = 0; f < K; ++f)
        {
          const double t = p.extraction_mz[f];
          if (!(t > 0.0)) continue;
          for (const auto& pk : s)
            if (std::abs((pk.getMZ() - t) / t * 1e6) < ppm) brute[(w * K + f) * B + bin] += pk.getIntensity();
        }
      }
    }
    FragmentGridCollector collector(precursors, windows, half, ppm);
    for (auto s : run) collector.consumeSpectrum(s);
    collector.finish();
    bool profiles = true, areas = true, evidence = true;
    for (std::size_t w = 0; w < windows.size(); ++w)
    {
      std::array<double, K * B> g{};
      for (std::size_t c = 0; c < K * B; ++c) g[c] = static_cast<double>(static_cast<float>(brute[w * K * B + c]));
      for (std::size_t b = 0; b < B; ++b)
      {
        double total = 0.0;
        for (std::size_t f = 0; f < K; ++f) total += g[f * B + b];
        profiles = profiles && std::abs(total - windows[w].profile[b]) <= 1e-6 * std::max(1.0, total);
      }
      for (std::size_t f = 0; f < K; ++f)
      {
        double a = 0.0;
        for (std::size_t b = FragmentGrid::CORE_BEGIN; b < FragmentGrid::CORE_END; ++b) a += g[f * B + b];
        areas = areas && std::abs(a - windows[w].core_area[f]) <= 1e-6 * std::max(1.0, a);
      }
      const auto& p = precursors[windows[w].precursor];
      const auto e = FragmentGrid::evidence(g.data(), p.library, p.search);
      for (std::size_t j = 0; j < FragmentGrid::EVIDENCE; ++j) evidence = evidence && std::abs(e[j] - collector.evidence()[w][j]) < 1e-5;
    }
    check(profiles, "collector profiles equal brute force");
    check(areas, "collector central areas equal brute force");
    check(evidence, "collector evidence equals brute force");
    check(collector.peakOpenGrids() < windows.size(), "grids are released after their window closes");
    const auto iso = collector.isolationIndices();
    check(iso[0] == 0 && iso[3] == 1, "isolation window indices");
    // Out-of-order spectra are rejected.
    FragmentGridCollector strict(precursors, windows, half, ppm);
    auto late = run[10], early = run[2];
    strict.consumeSpectrum(late);
    bool threw = false;
    try { strict.consumeSpectrum(early); } catch (const std::invalid_argument&) { threw = true; }
    check(threw, "out-of-order spectrum throws");
  }

  if (failures) return EXIT_FAILURE;
  std::cerr << "OpenDIAFragmentGrid_test passed\n";
  return EXIT_SUCCESS;
}
