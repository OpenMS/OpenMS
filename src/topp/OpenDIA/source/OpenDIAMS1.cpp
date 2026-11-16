// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAMS1.h>

#include <algorithm>
#include <cmath>

namespace OpenMS::OpenDIA
{
  namespace
  {
    using Composition = std::array<int, 5>; // C, H, N, O, S

    bool residue(char aa, Composition& c)
    {
      static const auto table = []
      {
        std::array<Composition, 26> t{};
        for (auto& r : t) r = {-1, 0, 0, 0, 0};
        const auto set = [&](char a, Composition v) { t[static_cast<std::size_t>(a - 'A')] = v; };
        set('G', {2, 3, 1, 1, 0}); set('A', {3, 5, 1, 1, 0}); set('S', {3, 5, 1, 2, 0}); set('P', {5, 7, 1, 1, 0});
        set('V', {5, 9, 1, 1, 0}); set('T', {4, 7, 1, 2, 0}); set('C', {3, 5, 1, 1, 1}); set('L', {6, 11, 1, 1, 0});
        set('I', {6, 11, 1, 1, 0}); set('N', {4, 6, 2, 2, 0}); set('D', {4, 5, 1, 3, 0}); set('Q', {5, 8, 2, 2, 0});
        set('K', {6, 12, 2, 1, 0}); set('E', {5, 7, 1, 3, 0}); set('M', {5, 9, 1, 1, 1}); set('H', {6, 7, 3, 1, 0});
        set('F', {9, 9, 1, 1, 0}); set('R', {6, 12, 4, 1, 0}); set('Y', {9, 9, 1, 2, 0}); set('W', {11, 10, 2, 1, 0});
        return t;
      }();
      if (aa < 'A' || aa > 'Z') return false;
      const Composition& r = table[static_cast<std::size_t>(aa - 'A')];
      if (r[0] < 0) return false;
      for (std::size_t e = 0; e < 5; ++e) c[e] += r[e];
      return true;
    }

    double pearson(const double* a, const double* b, std::size_t n)
    {
      double ma = 0.0, mb = 0.0;
      for (std::size_t i = 0; i < n; ++i) { ma += a[i]; mb += b[i]; }
      ma /= static_cast<double>(n);
      mb /= static_cast<double>(n);
      double sab = 0.0, saa = 0.0, sbb = 0.0;
      for (std::size_t i = 0; i < n; ++i)
      {
        sab += (a[i] - ma) * (b[i] - mb);
        saa += (a[i] - ma) * (a[i] - ma);
        sbb += (b[i] - mb) * (b[i] - mb);
      }
      return sab / std::sqrt(saa * sbb + 1e-12);
    }
  }

  const std::array<std::string, MS1Evidence::FEATURES>& MS1Evidence::names()
  {
    static const std::array<std::string, FEATURES> n = {"m1_corr_xic", "m1_iso_cos", "m1_iso_corr", "m1_apex_dt", "m1_log_m0"};
    return n;
  }

  std::array<double, 4> MS1Evidence::envelope(const std::string& sequence)
  {
    Composition c = {0, 2, 0, 1, 0};
    for (std::size_t i = 0; i < sequence.size();)
    {
      if (sequence[i] == '(')
      {
        const auto close = sequence.find(')', i);
        if (close == std::string::npos) return {};
        const std::string mod = sequence.substr(i + 1, close - i - 1);
        if (mod == "UniMod:4") { c[0] += 2; c[1] += 3; c[2] += 1; c[3] += 1; }
        else if (mod == "UniMod:35") { c[3] += 1; }
        else return {};
        i = close + 1;
        continue;
      }
      if (!residue(sequence[i], c)) return {};
      ++i;
    }
    // Per element: probabilities of +1 and +2 neutrons (C, H, N, O, S).
    static constexpr double heavy[5][2] = {{0.0107, 0.0}, {0.000115, 0.0}, {0.00364, 0.0}, {0.00038, 0.00205}, {0.0075, 0.0425}};
    std::array<double, 4> p = {1.0, 0.0, 0.0, 0.0};
    for (std::size_t e = 0; e < 5; ++e)
    {
      const std::array<double, 3> atom = {1.0 - heavy[e][0] - heavy[e][1], heavy[e][0], heavy[e][1]};
      for (int k = 0; k < c[e]; ++k)
      {
        std::array<double, 4> q{};
        for (std::size_t a = 0; a < 4; ++a)
          for (std::size_t b = 0; b < 3 && a + b < 4; ++b) q[a + b] += p[a] * atom[b];
        p = q;
      }
    }
    const double s = p[0] + p[1] + p[2] + p[3];
    for (auto& v : p) v /= s;
    return p;
  }

  std::array<double, MS1Evidence::FEATURES> MS1Evidence::features(const IsotopeTraces& iso, const std::array<double, BINS>& fragment_sum,
                                                                  const std::array<double, 4>& envelope)
  {
    constexpr std::size_t n = TRACE_END - TRACE_BEGIN;
    std::array<std::array<double, n>, 4> m{};
    std::array<double, n> x{};
    double peak = 0.0;
    for (std::size_t b = 0; b < n; ++b)
    {
      for (std::size_t i = 0; i < 4; ++i) m[i][b] = iso[i][TRACE_BEGIN + b];
      x[b] = fragment_sum[TRACE_BEGIN + b];
      peak = std::max(peak, m[0][b]);
    }
    std::array<double, 4> area{};
    for (std::size_t i = 0; i < 4; ++i)
      for (std::size_t b = CORE_BEGIN; b < CORE_END; ++b) area[i] += iso[i][b];
    std::array<double, FEATURES> f{};
    f[4] = std::log1p(area[0]);
    if (!(peak > 0.0))
    {
      f[3] = 20.0;
      return f;
    }
    f[0] = pearson(m[0].data(), x.data(), n);
    const double total = area[0] + area[1] + area[2] + area[3];
    double dot = 0.0, no = 0.0, ne = 0.0;
    for (std::size_t i = 0; i < 4; ++i)
    {
      const double o = area[i] / std::max(total, 1e-12);
      dot += o * envelope[i];
      no += o * o;
      ne += envelope[i] * envelope[i];
    }
    f[1] = dot / (std::sqrt(no) * std::sqrt(ne) + 1e-12);
    f[2] = 0.5 * (pearson(m[0].data(), m[1].data(), n) + pearson(m[0].data(), m[2].data(), n));
    const auto am = static_cast<long>(std::max_element(m[0].begin(), m[0].end()) - m[0].begin());
    const auto ax = static_cast<long>(std::max_element(x.begin(), x.end()) - x.begin());
    f[3] = static_cast<double>(std::labs(am - ax));
    return f;
  }
}
