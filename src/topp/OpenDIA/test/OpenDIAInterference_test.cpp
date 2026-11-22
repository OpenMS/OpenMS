// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAInterference.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::OpenDIA;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition) throw std::runtime_error(message);
  }
  bool near(double a, double b, double tol = 1e-5) { return std::abs(a - b) <= tol * std::max(1.0, std::abs(b)); }

  constexpr std::size_t K = InterferenceResolver::FRAGMENTS;
  constexpr std::size_t NB = InterferenceResolver::BINS;

  std::array<double, NB> gaussian(double centre_bin, double sigma)
  {
    std::array<double, NB> p{};
    for (std::size_t b = 0; b < NB; ++b) p[b] = std::exp(-0.5 * std::pow((b - centre_bin) / sigma, 2));
    return p;
  }

  /// Window with the profile prof x abundance x lib per fragment.
  InterferenceResolver::Window window(std::uint32_t precursor, const InterferenceResolver::Precursor& p, double abundance,
                                      const std::array<double, NB>& prof, double apex = 100.0)
  {
    InterferenceResolver::Window w;
    w.precursor = precursor;
    w.rt_apex = apex;
    w.rt_start = apex - 3.0;
    w.rt_end = apex + 3.0;
    for (std::size_t b = 0; b < NB; ++b)
    {
      for (std::size_t s = 0; s < K; ++s)
      {
        const double v = abundance * p.library[s] * prof[b];
        w.profile[b] += v;
        if (b >= InterferenceResolver::CORE_BEGIN && b < InterferenceResolver::CORE_END) w.core_area[s] += v;
      }
    }
    return w;
  }

  InterferenceResolver::Precursor precursor(double base, std::int32_t iso)
  {
    InterferenceResolver::Precursor p;
    p.isolation = iso;
    for (std::size_t s = 0; s < K; ++s)
    {
      p.fragment_mz[s] = base + 17.3 * s;
      p.library[s] = 1.0 - 0.06 * s;
    }
    return p;
  }

  /// The analytic case of features_I.py's synthetic test: a carries only b's signal on its shared fragment 0.
  void testExplainedShare()
  {
    auto pb = precursor(300.0, 0), pa = precursor(1300.0, 0), pc = precursor(2300.0, 1);
    pb.fragment_mz[0] = pa.fragment_mz[0] = 500.0;
    const auto prof = gaussian(15.5, 2.0);
    auto b = window(0, pb, 100.0, prof), a = window(1, pa, 10.0, prof), c = window(2, pc, 10.0, prof);
    // a's fragment 0 holds only b's signal: abundance 100 x lib_b[0].
    a.profile.fill(0.0);
    a.core_area.fill(0.0);
    for (std::size_t bin = 0; bin < NB; ++bin)
    {
      for (std::size_t s = 0; s < K; ++s)
      {
        const double v = (s == 0 ? 100.0 * pb.library[0] : 10.0 * pa.library[s]) * prof[bin];
        a.profile[bin] += v;
        if (bin >= InterferenceResolver::CORE_BEGIN && bin < InterferenceResolver::CORE_END) a.core_area[s] += v;
      }
    }
    const std::vector<InterferenceResolver::Precursor> precursors{pb, pa, pc};
    const std::vector<InterferenceResolver::Window> windows{b, a, c};
    InterferenceResolver::Statistics stats;
    const auto f = InterferenceResolver::compute(precursors, windows, 7.0, 30.0 / 32, 1, &stats);
    double area = 0.0;
    for (double v : a.core_area) area += v;
    // b's scale from its unshared fragments is 100 (exact), frac = 1 (same apex), so a's fragment 0 is fully explained.
    const double expected_share = a.core_area[0] / area;
    require(near(f[4 + 0], expected_share), "explained share of a");
    require(f[4 + 2] > 0.0f, "residual library cosine of a rises");
    require(near(f[4 + 3], std::log1p(1.0)), "one claimant for a");
    // b is explained by a at a's abundance (10), a smaller share.
    require(f[0] > 0.0f && f[0] < f[4], "claimant b gets a smaller explained share");
    // c is alone in its isolation window.
    require(f[8 + 0] == 0.0f && f[8 + 3] == 0.0f && f[8 + 2] == 0.0f, "other isolation window untouched");
    require(stats.links == 2 && stats.used_links == 2 && stats.windows_with_claimant == 2, "link statistics");
  }

  void testNoOverlapNoSharedNoUnshared()
  {
    auto pb = precursor(300.0, 0), pa = precursor(1300.0, 0);
    pb.fragment_mz[0] = pa.fragment_mz[0] = 500.0;
    const auto prof = gaussian(15.5, 2.0);
    auto b = window(0, pb, 100.0, prof, 200.0), a = window(1, pa, 10.0, prof, 100.0);
    auto f = InterferenceResolver::compute({pb, pa}, {b, a}, 7.0, 30.0 / 32);
    require(f[4] == 0.0f && f[7] == 0.0f, "no RT overlap: no link");
    // Same RT, no shared fragment m/z within 7 ppm.
    pb.fragment_mz[0] = 500.01;
    b = window(0, pb, 100.0, prof, 100.0);
    f = InterferenceResolver::compute({pb, pa}, {b, a}, 7.0, 30.0 / 32);
    require(f[4] == 0.0f && f[7] == 0.0f, "no shared fragment: no link");
    // Shared, but the claimant has no unshared fragment with signal: it explains nothing.
    pb.fragment_mz[0] = 500.0;
    for (std::size_t s = 1; s < K; ++s) pb.fragment_mz[s] = 0.0, pb.library[s] = 0.0;
    b = window(0, pb, 100.0, prof, 100.0);
    f = InterferenceResolver::compute({pb, pa}, {b, a}, 7.0, 30.0 / 32);
    require(f[4] == 0.0f && f[7] == 0.0f, "claimant without unshared signal explains nothing");
    // Same precursor's own windows never claim each other.
    auto a2 = window(1, pa, 10.0, prof, 101.0);
    f = InterferenceResolver::compute({pa}, {window(0, pa, 10.0, prof), [&] { auto w = a2; w.precursor = 0; return w; }()}, 7.0, 30.0 / 32);
    require(f[0] == 0.0f && f[3] == 0.0f, "own windows are not claimants");
    // Mobility: the same pair links only when the 1/K0 difference is within the link width (or a 1/K0 is unknown).
    {
      auto mb = precursor(300.0, 0), ma = precursor(1300.0, 0);
      mb.fragment_mz[0] = ma.fragment_mz[0] = 500.0;
      const auto bw = window(0, mb, 100.0, prof), aw = window(1, ma, 10.0, prof);
      mb.im = 0.80;
      ma.im = 1.00;
      const auto apart = InterferenceResolver::compute({mb, ma}, {bw, aw}, 7.0, 30.0 / 32, 1, nullptr, 0.06);
      require(apart[4] == 0.0f && apart[7] == 0.0f, "different 1/K0: no link");
      const auto m_only = InterferenceResolver::compute({mb, ma}, {bw, aw}, 7.0, 30.0 / 32, 1, nullptr, 0.0);
      require(m_only[4] > 0.0f, "link width 0: m/z-only link as before");
      ma.im = 0.83;
      const auto close = InterferenceResolver::compute({mb, ma}, {bw, aw}, 7.0, 30.0 / 32, 1, nullptr, 0.06);
      require(close[4] == m_only[4], "1/K0 within the link width: same link");
      ma.im = 0.0;
      const auto unknown = InterferenceResolver::compute({mb, ma}, {bw, aw}, 7.0, 30.0 / 32, 1, nullptr, 0.06);
      require(unknown[4] == m_only[4], "unknown 1/K0: m/z-only link");
    }
  }

  /// Brute-force reference (all pairs) against the sweep, on random data; also thread independence.
  void testRandomAgainstBruteForce()
  {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    const std::size_t np = 300;
    std::vector<InterferenceResolver::Precursor> precursors(np);
    std::vector<double> pool;
    for (int i = 0; i < 60; ++i) pool.push_back(200.0 + 1000.0 * u(rng));
    for (auto& p : precursors)
    {
      p.isolation = static_cast<std::int32_t>(u(rng) * 3);
      for (std::size_t s = 0; s < K; ++s)
      {
        if (u(rng) < 0.1) continue;
        p.fragment_mz[s] = u(rng) < 0.3 ? pool[static_cast<std::size_t>(u(rng) * pool.size())] * (1.0 + (u(rng) - 0.5) * 8e-6)
                                        : 200.0 + 1500.0 * u(rng);
        p.library[s] = 0.05 + u(rng);
      }
    }
    std::vector<InterferenceResolver::Window> windows;
    for (std::uint32_t p = 0; p < np; ++p)
    {
      for (int k = 0; k < 3; ++k)
      {
        InterferenceResolver::Window w;
        w.precursor = p;
        w.rt_apex = 100.0 + 200.0 * u(rng);
        w.rt_start = w.rt_apex - 1.0 - 4.0 * u(rng);
        w.rt_end = w.rt_apex + 1.0 + 4.0 * u(rng);
        for (auto& v : w.profile) v = u(rng) < 0.3 ? 0.0 : 1000.0 * u(rng);
        for (std::size_t s = 0; s < K; ++s) w.core_area[s] = precursors[p].fragment_mz[s] > 0.0 && u(rng) > 0.2 ? 1000.0 * u(rng) : 0.0;
        windows.push_back(w);
      }
    }
    const auto one = InterferenceResolver::compute(precursors, windows, 7.0, 30.0 / 32, 1);
    const auto four = InterferenceResolver::compute(precursors, windows, 7.0, 30.0 / 32, 4);
    require(one == four, "result independent of thread count");
    // Brute force: every pair, without the sweep.
    const double tol = 7e-6;
    for (std::size_t ai = 0; ai < windows.size(); ++ai)
    {
      const auto& a = windows[ai];
      const auto& pa = precursors[a.precursor];
      std::array<double, K> contribution{};
      std::size_t claimants = 0;
      // nearest claimant window per precursor
      std::vector<long> best(np, -1);
      for (std::size_t bi = 0; bi < windows.size(); ++bi)
      {
        const auto& b = windows[bi];
        if (b.precursor == a.precursor || precursors[b.precursor].isolation != pa.isolation) continue;
        if (b.rt_start > a.rt_end || b.rt_end < a.rt_start) continue;
        bool shared = false;
        for (std::size_t i = 0; i < K && !shared; ++i)
          for (std::size_t j = 0; j < K && !shared; ++j)
            shared = pa.fragment_mz[i] > 0.0 && precursors[b.precursor].fragment_mz[j] > 0.0 &&
                     std::abs(pa.fragment_mz[i] - precursors[b.precursor].fragment_mz[j]) <= tol * pa.fragment_mz[i];
        if (!shared) continue;
        auto& cur = best[b.precursor];
        if (cur < 0 || std::abs(a.rt_apex - b.rt_apex) < std::abs(a.rt_apex - windows[cur].rt_apex)) cur = static_cast<long>(bi);
      }
      for (std::size_t q = 0; q < np; ++q)
      {
        if (best[q] < 0) continue;
        const auto& b = windows[best[q]];
        const auto& pb = precursors[q];
        std::array<int, K> a_to_b, b_to_a;
        a_to_b.fill(-1);
        b_to_a.fill(-1);
        for (std::size_t i = 0; i < K; ++i)
          for (std::size_t j = 0; j < K; ++j)
            if (pa.fragment_mz[i] > 0.0 && pb.fragment_mz[j] > 0.0 && std::abs(pa.fragment_mz[i] - pb.fragment_mz[j]) <= tol * pa.fragment_mz[i])
            {
              if (a_to_b[i] < 0) a_to_b[i] = static_cast<int>(j);
              if (b_to_a[j] < 0) b_to_a[j] = static_cast<int>(i);
            }
        std::vector<double> r;
        for (std::size_t g = 0; g < K; ++g)
          if (b_to_a[g] < 0 && pb.library[g] > 0.0 && b.core_area[g] > 0.0) r.push_back(b.core_area[g] / pb.library[g]);
        if (r.empty()) continue;
        std::sort(r.begin(), r.end());
        const double scale = r.size() % 2 ? r[r.size() / 2] : 0.5 * (r[r.size() / 2 - 1] + r[r.size() / 2]);
        const long shift = static_cast<long>(std::nearbyint((a.rt_apex - b.rt_apex) / (30.0 / 32)));
        const long lo = std::clamp(12L + shift, 0L, 32L), hi = std::clamp(20L + shift, 0L, 32L);
        double own = 0.0, in = 0.0;
        for (long k = 12; k < 20; ++k) own += b.profile[k];
        for (long k = lo; k < hi; ++k) in += b.profile[k];
        const double frac = own > 0.0 ? in / own : 0.0;
        double total = 0.0;
        for (std::size_t f = 0; f < K; ++f)
          if (a_to_b[f] >= 0)
          {
            contribution[f] += scale * pb.library[a_to_b[f]] * frac;
            total += scale * pb.library[a_to_b[f]] * frac;
          }
        if (total > 0.0) ++claimants;
      }
      double observed = 0.0, explained = 0.0;
      for (std::size_t f = 0; f < K; ++f)
      {
        observed += a.core_area[f];
        explained += std::min(a.core_area[f], contribution[f]);
      }
      const double share = observed > 0.0 ? explained / observed : 0.0;
      require(near(one[ai * 4 + 0], share, 1e-4), "explained share equals brute force");
      require(near(one[ai * 4 + 3], std::log1p(static_cast<double>(claimants)), 1e-5), "claimant count equals brute force");
    }
  }
}

int main()
{
  try
  {
    testExplainedShare();
    testNoOverlapNoSharedNoUnshared();
    testRandomAgainstBruteForce();
  }
  catch (const std::exception& error)
  {
    std::cerr << "OpenDIAInterference_test failed: " << error.what() << '\n';
    return 1;
  }
  std::cout << "OpenDIAInterference_test passed\n";
  return 0;
}
