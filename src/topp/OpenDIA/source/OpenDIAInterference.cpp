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
#include <limits>
#include <numeric>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    using Resolver = InterferenceResolver;
    constexpr std::size_t K = Resolver::FRAGMENTS;
    constexpr std::size_t C0 = Resolver::CORE_BEGIN;
    constexpr std::size_t C1 = Resolver::CORE_END;
    constexpr long NB = static_cast<long>(Resolver::BINS);
    /// Scans closer than this to a window apex are the apex scan (the apex is written with ~1e-5 s precision).
    constexpr double APEX_SNAP_S = 1e-3;

    /// Fragments of one precursor sorted by m/z, with their slots, for the binary-search shared-fragment test.
    struct SortedFragments
    {
      std::array<double, K> mz{};
      std::array<std::uint8_t, K> slot{};
      std::uint8_t count = 0;
    };

    SortedFragments sortFragments(const Resolver::Precursor& p)
    {
      SortedFragments s;
      std::array<std::uint8_t, K> order{};
      for (std::size_t k = 0; k < K; ++k)
      {
        if (p.fragment_mz[k] > 0.0) order[s.count++] = static_cast<std::uint8_t>(k);
      }
      std::stable_sort(order.begin(), order.begin() + s.count,
                       [&](std::uint8_t x, std::uint8_t y) { return p.fragment_mz[x] < p.fragment_mz[y]; });
      for (std::size_t i = 0; i < s.count; ++i)
      {
        s.mz[i] = p.fragment_mz[order[i]];
        s.slot[i] = order[i];
      }
      return s;
    }

    /**
      Shared-fragment maps of the ordered pair (a, b), as features_I.py: a_to_b[i] is the lowest slot of b within
      tol x mz_a[i] of a's fragment i, b_to_a[j] the lowest slot of a whose window contains b's fragment j (-1 none).
      The tolerance is always relative to a's fragment m/z. Returns whether any fragment is shared.
    */
    bool sharedFragments(const Resolver::Precursor& a, const SortedFragments& sb, double tol, std::array<std::int8_t, K>& a_to_b,
                         std::array<std::int8_t, K>& b_to_a)
    {
      a_to_b.fill(-1);
      b_to_a.fill(-1);
      bool any = false;
      const double* begin = sb.mz.data();
      const double* end = begin + sb.count;
      for (std::size_t i = 0; i < K; ++i)
      {
        const double m = a.fragment_mz[i];
        if (!(m > 0.0)) continue;
        const double half = tol * m;
        for (const double* p = std::lower_bound(begin, end, m - half); p != end && *p <= m + half; ++p)
        {
          if (!(std::abs(m - *p) <= half)) continue;  // exact form of the reference test
          const auto j = static_cast<std::int8_t>(sb.slot[p - begin]);
          any = true;
          if (a_to_b[i] < 0 || j < a_to_b[i]) a_to_b[i] = j;
          if (b_to_a[j] < 0 || static_cast<std::int8_t>(i) < b_to_a[j]) b_to_a[j] = static_cast<std::int8_t>(i);
        }
      }
      return any;
    }

    /// Median of at most K values (insertion sort; n is small).
    double medianInPlace(double* v, std::size_t n)
    {
      for (std::size_t i = 1; i < n; ++i)
      {
        const double x = v[i];
        std::size_t j = i;
        for (; j > 0 && v[j - 1] > x; --j) v[j] = v[j - 1];
        v[j] = x;
      }
      return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    }

    double sqrtCosine(const std::array<double, K>& x, const std::array<double, K>& lib)
    {
      double dot = 0.0, xs = 0.0, ls = 0.0;
      for (std::size_t k = 0; k < K; ++k)
      {
        dot += std::sqrt(x[k]) * std::sqrt(lib[k]);
        xs += x[k];
        ls += lib[k];
      }
      return dot / std::sqrt(xs * ls + 1e-12);
    }

    struct Claim
    {
      std::uint32_t precursor;
      double distance;
      std::uint32_t window;
      std::array<std::int8_t, K> a_to_b;
      std::array<std::int8_t, K> b_to_a;
    };
  }

  const std::array<const char*, InterferenceResolver::FEATURES>& InterferenceResolver::featureNames()
  {
    static const std::array<const char*, FEATURES> names = {"i_explained", "i_libcos_resid", "i_dlibcos", "i_nclaim"};
    return names;
  }

  std::vector<float> InterferenceResolver::compute(const std::vector<Precursor>& precursors, const std::vector<Window>& windows,
                                                   double fragment_tolerance_ppm, double rt_bin_seconds, int threads, Statistics* stats,
                                                   double im_link_width)
  {
    if (!std::isfinite(fragment_tolerance_ppm) || !(fragment_tolerance_ppm > 0.0) || !std::isfinite(rt_bin_seconds) ||
        !(rt_bin_seconds > 0.0) || threads < 1 || !std::isfinite(im_link_width) || im_link_width < 0.0)
    {
      throw std::invalid_argument("InterferenceResolver: the fragment tolerance, bin width and thread count must be positive.");
    }
    const double tol = fragment_tolerance_ppm * 1e-6;
    const std::size_t n = windows.size();
    for (const auto& w : windows)
    {
      if (w.precursor >= precursors.size()) throw std::invalid_argument("InterferenceResolver: window precursor index out of range.");
    }
    std::vector<SortedFragments> sorted(precursors.size());
    for (std::size_t p = 0; p < precursors.size(); ++p) sorted[p] = sortFragments(precursors[p]);

    // Per window: cumulative profile, observed sqrt cosine and total central area.
    std::vector<std::array<double, BINS + 1>> cumulative(n);
    for (std::size_t i = 0; i < n; ++i)
    {
      cumulative[i][0] = 0.0;
      for (std::size_t b = 0; b < BINS; ++b) cumulative[i][b + 1] = cumulative[i][b] + windows[i].profile[b];
    }

    // Windows of linked isolation windows sorted by rt_start (groups of one isolation window). Each window computes its
    // own features from its partners, so results are written by exactly one iteration and do not depend on the thread
    // count or schedule.
    std::vector<std::uint32_t> order;
    order.reserve(n);
    for (std::uint32_t i = 0; i < n; ++i)
    {
      if (precursors[windows[i].precursor].isolation >= 0) order.push_back(i);
    }
    std::sort(order.begin(), order.end(), [&](std::uint32_t x, std::uint32_t y)
    {
      const auto ix = precursors[windows[x].precursor].isolation, iy = precursors[windows[y].precursor].isolation;
      if (ix != iy) return ix < iy;
      if (windows[x].rt_start != windows[y].rt_start) return windows[x].rt_start < windows[y].rt_start;
      return x < y;
    });
    const std::size_t m = order.size();
    std::vector<std::uint32_t> group_end(m);
    for (std::size_t g = 0; g < m;)
    {
      const auto iso = precursors[windows[order[g]].precursor].isolation;
      std::size_t h = g;
      while (h < m && precursors[windows[order[h]].precursor].isolation == iso) ++h;
      std::fill(group_end.begin() + g, group_end.begin() + h, static_cast<std::uint32_t>(h));
      g = h;
    }

    std::vector<float> out(n * FEATURES, 0.0f);
    // Windows without an isolation window have no claimant.
    for (std::size_t i = 0; i < n; ++i)
    {
      out[i * FEATURES + 1] = static_cast<float>(sqrtCosine(windows[i].core_area, precursors[windows[i].precursor].library));
    }
    std::size_t pairs_total = 0, links_total = 0, used_total = 0, with_claimant = 0;

    // Partners are looked up through a fragment index of one isolation window at a time: every window is entered once
    // per fragment under the log(m/z) bin of width log(1 + tol), sorted by (bin, rt_start). A partner b of a must share a
    // fragment (|m_a - m_b| <= tol x m_a, so ln m_b lies within ln m_a + [ln(1 - tol), ln(1 + tol)]), and its range
    // must overlap a's (b.rt_start in [a.rt_start - widest range of the group, a.rt_end]). The lookup returns a superset
    // of the linked windows (one spare bin on each side), and consider() applies exactly the tests of the former scan
    // over all overlapping windows, so the claims and the features are unchanged. Only pairs_tested counts fewer pairs.
    const double bin_width = std::log1p(tol);
    const double lo_shift = std::log1p(-tol), hi_shift = std::log1p(tol);
    struct Entry
    {
      std::int64_t bin;
      double rt_start;
      std::uint32_t pos;
      std::uint32_t precursor;
      double rt_end;
      double im;
    };
    std::vector<Entry> index;
    // Distinct bins of the index and the start of each bin's entries (bin_start[k + 1] ends bin k).
    std::vector<std::int64_t> bins;
    std::vector<std::uint32_t> bin_start;
    for (std::size_t g = 0; g < m;)
    {
      const std::size_t h = group_end[g];
      index.clear();
      double widest = 0.0;
      for (std::size_t q = g; q < h; ++q)
      {
        const Window& b = windows[order[q]];
        widest = std::max(widest, b.rt_end - b.rt_start);
        const SortedFragments& sb = sorted[b.precursor];
        for (std::size_t k = 0; k < sb.count; ++k)
          index.push_back({static_cast<std::int64_t>(std::floor(std::log(sb.mz[k]) / bin_width)), b.rt_start, static_cast<std::uint32_t>(q),
                           b.precursor, b.rt_end, precursors[b.precursor].im});
      }
      std::sort(index.begin(), index.end(), [](const Entry& x, const Entry& y)
      {
        if (x.bin != y.bin) return x.bin < y.bin;
        if (x.rt_start != y.rt_start) return x.rt_start < y.rt_start;
        return x.pos < y.pos;
      });
      bins.clear();
      bin_start.clear();
      for (std::size_t e = 0; e < index.size(); ++e)
      {
        if (e == 0 || index[e].bin != index[e - 1].bin)
        {
          bins.push_back(index[e].bin);
          bin_start.push_back(static_cast<std::uint32_t>(e));
        }
      }
      bin_start.push_back(static_cast<std::uint32_t>(index.size()));
#pragma omp parallel num_threads(threads) reduction(+ : pairs_total, links_total, used_total, with_claimant)
    {
      std::vector<Claim> claims;
      std::vector<std::uint32_t> partners;
#pragma omp for schedule(dynamic, 512)
      for (std::ptrdiff_t pos = static_cast<std::ptrdiff_t>(g); pos < static_cast<std::ptrdiff_t>(h); ++pos)
      {
        const std::uint32_t ai = order[pos];
        const Window& a = windows[ai];
        const Precursor& pa = precursors[a.precursor];
        claims.clear();
        const auto consider = [&](std::size_t q)
        {
          const std::uint32_t bi = order[q];
          const Window& b = windows[bi];
          if (b.precursor == a.precursor || b.rt_start > a.rt_end || b.rt_end < a.rt_start) return;
          // With mobility, precursors whose extraction windows in 1/K0 do not overlap cannot share signal.
          const Precursor& pb_im = precursors[b.precursor];
          if (im_link_width > 0.0 && pa.im > 0.0 && pb_im.im > 0.0 && std::abs(pa.im - pb_im.im) > im_link_width) return;
          ++pairs_total;
          Claim c;
          if (!sharedFragments(pa, sorted[b.precursor], tol, c.a_to_b, c.b_to_a)) return;
          ++links_total;
          c.precursor = b.precursor;
          c.distance = std::abs(a.rt_apex - b.rt_apex);
          c.window = bi;
          claims.push_back(c);
        };
        partners.clear();
        const double rt_from = a.rt_start - widest;
        const auto before = [](const Entry& x, double rt) { return x.rt_start < rt; };
        // The entry carries the fields of the cheap rejection tests of consider() (same precursor, no RT overlap,
        // 1/K0 too far), so most index hits are rejected without touching the windows; consider() repeats every test.
        const bool test_im = im_link_width > 0.0 && pa.im > 0.0;
        for (std::size_t f = 0; f < K; ++f)
        {
          const double mz = pa.fragment_mz[f];
          if (!(mz > 0.0)) continue;
          const double l = std::log(mz);
          const auto first_bin = static_cast<std::int64_t>(std::floor((l + lo_shift) / bin_width)) - 1;
          const auto last_bin = static_cast<std::int64_t>(std::floor((l + hi_shift) / bin_width)) + 1;
          for (std::size_t k = static_cast<std::size_t>(std::lower_bound(bins.begin(), bins.end(), first_bin) - bins.begin());
               k < bins.size() && bins[k] <= last_bin; ++k)
          {
            const auto bin_end = index.begin() + bin_start[k + 1];
            for (auto it = std::lower_bound(index.begin() + bin_start[k], bin_end, rt_from, before);
                 it != bin_end && it->rt_start <= a.rt_end; ++it)
            {
              if (it->pos == static_cast<std::uint32_t>(pos) || it->precursor == a.precursor || it->rt_end < a.rt_start) continue;
              if (test_im && it->im > 0.0 && std::abs(pa.im - it->im) > im_link_width) continue;
              partners.push_back(it->pos);
            }
          }
        }
        std::sort(partners.begin(), partners.end());
        partners.erase(std::unique(partners.begin(), partners.end()), partners.end());
        for (const std::uint32_t q : partners) consider(q);
        if (claims.empty()) continue;
        ++with_claimant;
        // One claimant window per precursor: nearest apex, ties by window index.
        std::sort(claims.begin(), claims.end(), [](const Claim& x, const Claim& y)
        {
          if (x.precursor != y.precursor) return x.precursor < y.precursor;
          if (x.distance != y.distance) return x.distance < y.distance;
          return x.window < y.window;
        });
        std::array<double, K> contribution{};
        std::size_t claimants = 0;
        double ratios[K];
        for (std::size_t c = 0; c < claims.size(); ++c)
        {
          if (c > 0 && claims[c].precursor == claims[c - 1].precursor) continue;
          ++used_total;
          const Claim& claim = claims[c];
          const Window& b = windows[claim.window];
          const Precursor& pb = precursors[b.precursor];
          std::size_t count = 0;
          for (std::size_t g = 0; g < K; ++g)
          {
            if (claim.b_to_a[g] < 0 && pb.library[g] > 0.0 && b.core_area[g] > 0.0) ratios[count++] = b.core_area[g] / std::max(pb.library[g], 1e-12);
          }
          if (count == 0) continue;
          const double scale = medianInPlace(ratios, count);
          const long shift = static_cast<long>(std::nearbyint((a.rt_apex - b.rt_apex) / rt_bin_seconds));
          const long lo = std::clamp(static_cast<long>(C0) + shift, 0L, NB);
          const long hi = std::clamp(static_cast<long>(C1) + shift, 0L, NB);
          const auto& cb = cumulative[claim.window];
          const double own = cb[C1] - cb[C0];
          const double frac = own > 0.0 ? (cb[hi] - cb[lo]) / std::max(own, 1e-12) : 0.0;
          double total = 0.0;
          for (std::size_t f = 0; f < K; ++f)
          {
            if (claim.a_to_b[f] < 0) continue;
            const double e = scale * pb.library[static_cast<std::size_t>(claim.a_to_b[f])] * frac;
            contribution[f] += e;
            total += e;
          }
          if (total > 0.0) ++claimants;
        }
        std::array<double, K> residual{};
        double observed = 0.0, explained = 0.0;
        for (std::size_t f = 0; f < K; ++f)
        {
          const double e = std::min(a.core_area[f], contribution[f]);
          residual[f] = a.core_area[f] - e;
          observed += a.core_area[f];
          explained += e;
        }
        const double c_res = sqrtCosine(residual, pa.library);
        const double c_obs = sqrtCosine(a.core_area, pa.library);
        float* o = out.data() + static_cast<std::size_t>(ai) * FEATURES;
        o[0] = static_cast<float>(observed > 0.0 ? explained / std::max(observed, 1e-12) : 0.0);
        o[1] = static_cast<float>(c_res);
        o[2] = static_cast<float>(c_res - c_obs);
        o[3] = static_cast<float>(std::log1p(static_cast<double>(claimants)));
      }
    }
      g = h;
    }
    if (stats)
    {
      stats->pairs_tested = pairs_total;
      stats->links = links_total;
      stats->used_links = used_total;
      stats->windows_with_claimant = with_claimant;
    }
    return out;
  }
}
