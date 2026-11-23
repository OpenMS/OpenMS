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
      Both fragment lists are sorted by m/z, so one merge visits the same candidates as a binary search per fragment of a:
      the window start m x (1 - tol) rises with m, so b's first candidate only moves forward. The maps take minima, so
      the visiting order does not change them.
    */
    bool sharedFragments(const SortedFragments& sa, const SortedFragments& sb, double tol, std::array<std::int8_t, K>& a_to_b,
                         std::array<std::int8_t, K>& b_to_a)
    {
      a_to_b.fill(-1);
      b_to_a.fill(-1);
      bool any = false;
      std::size_t start = 0;
      for (std::size_t ia = 0; ia < sa.count; ++ia)
      {
        const double m = sa.mz[ia];
        const double half = tol * m;
        while (start < sb.count && sb.mz[start] < m - half) ++start;
        const auto i = static_cast<std::int8_t>(sa.slot[ia]);
        for (std::size_t p = start; p < sb.count && sb.mz[p] <= m + half; ++p)
        {
          if (!(std::abs(m - sb.mz[p]) <= half)) continue;  // exact form of the reference test
          const auto j = static_cast<std::int8_t>(sb.slot[p]);
          any = true;
          if (a_to_b[i] < 0 || j < a_to_b[i]) a_to_b[i] = j;
          if (b_to_a[j] < 0 || i < b_to_a[j]) b_to_a[j] = i;
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
      std::uint32_t pos;  ///< Position of the claimant window in the isolation group (index into the packed data).
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
    // of the linked windows (one spare bin on each side). The window tests (other precursor, RT overlap, 1/K0) are
    // applied to the index entries; the fragment test depends only on the two precursors, so it runs once per claimant
    // precursor, and each claimant precursor's nearest-apex window is its claimant window (ties by window index).
    // pairs_tested and links count windows, as when every window pair was tested.
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
      double rt_apex;
    };
    std::vector<Entry> index;
    // What a claimant window contributes, copied in group order (by rt_start) so that the claimants of a window, which
    // all overlap it in RT, are read from a few nearby megabytes instead of the whole run's window and precursor arrays.
    struct Packed
    {
      std::array<double, BINS + 1> cumulative;  ///< Cumulative summed profile.
      std::array<double, K> core_area;
      std::array<double, K> library;            ///< Library intensities of the window's precursor.
      SortedFragments fragments;                ///< Fragments of the window's precursor, sorted by m/z.
      double rt_apex;
    };
    std::vector<Packed> packed;
    std::vector<std::int32_t> pid_map(precursors.size(), -1);
    std::vector<std::uint32_t> local_pid, local_to_global;
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
                           b.precursor, b.rt_end, precursors[b.precursor].im, b.rt_apex});
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
      packed.resize(h - g);
      // Group-local precursor ids of the group's windows (dense, first-seen order), for the per-precursor claimant slots.
      local_pid.resize(h - g);
      local_to_global.clear();
      for (std::size_t q = g; q < h; ++q)
      {
        const std::uint32_t p = windows[order[q]].precursor;
        if (pid_map[p] < 0)
        {
          pid_map[p] = static_cast<std::int32_t>(local_to_global.size());
          local_to_global.push_back(p);
        }
        local_pid[q - g] = static_cast<std::uint32_t>(pid_map[p]);
      }
#pragma omp parallel for num_threads(threads) schedule(static)
      for (std::ptrdiff_t q = static_cast<std::ptrdiff_t>(g); q < static_cast<std::ptrdiff_t>(h); ++q)
      {
        const Window& b = windows[order[q]];
        Packed& pk = packed[static_cast<std::size_t>(q) - g];
        pk.cumulative[0] = 0.0;
        for (std::size_t k = 0; k < BINS; ++k) pk.cumulative[k + 1] = pk.cumulative[k] + b.profile[k];
        pk.core_area = b.core_area;
        pk.library = precursors[b.precursor].library;
        pk.fragments = sorted[b.precursor];
        pk.rt_apex = b.rt_apex;
      }
#pragma omp parallel num_threads(threads) reduction(+ : pairs_total, links_total, used_total, with_claimant)
    {
      // Partner windows of a, reduced on the fly to the nearest window per claimant precursor (ties by window index),
      // which is what sorting all partner windows by (precursor, distance, window) and taking each precursor's first
      // gave. slot[local precursor] holds the index into `best` (-1 none); pairs counts distinct partner windows.
      struct Candidate
      {
        std::uint32_t precursor;
        std::uint32_t window;
        double distance;
        std::uint32_t pos;
        std::uint32_t windows;  ///< Distinct partner windows of this precursor.
      };
      std::vector<Candidate> best;
      std::vector<std::int32_t> slot(local_to_global.size(), -1);
      std::vector<std::uint32_t> seen_window;  ///< Windows (group positions) already counted for a, to skip repeats.
      std::vector<std::uint8_t> seen(h - g, 0);
      std::vector<Claim> claims;
#pragma omp for schedule(dynamic, 512)
      for (std::ptrdiff_t pos = static_cast<std::ptrdiff_t>(g); pos < static_cast<std::ptrdiff_t>(h); ++pos)
      {
        const std::uint32_t ai = order[pos];
        const Window& a = windows[ai];
        const Precursor& pa = precursors[a.precursor];
        claims.clear();
        best.clear();
        seen_window.clear();
        const double rt_from = a.rt_start - widest;
        const auto before = [](const Entry& x, double rt) { return x.rt_start < rt; };
        // An entry passes when its window is of another precursor, its peak range overlaps a's and (with mobility)
        // its 1/K0 is within the link width: exactly the window tests of a link.
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
              if (it->precursor == a.precursor || it->rt_end < a.rt_start) continue;
              if (test_im && it->im > 0.0 && std::abs(pa.im - it->im) > im_link_width) continue;
              const std::uint32_t q = it->pos - static_cast<std::uint32_t>(g);
              if (seen[q]) continue;  // the same window found through another fragment
              seen[q] = 1;
              seen_window.push_back(q);
              const std::uint32_t lp = local_pid[q];
              const Candidate cand{it->precursor, order[it->pos], std::abs(a.rt_apex - it->rt_apex), q, 1};
              if (slot[lp] < 0)
              {
                slot[lp] = static_cast<std::int32_t>(best.size());
                best.push_back(cand);
                continue;
              }
              Candidate& b = best[static_cast<std::size_t>(slot[lp])];
              ++b.windows;
              if (cand.distance < b.distance || (cand.distance == b.distance && cand.window < b.window))
              {
                b.window = cand.window;
                b.distance = cand.distance;
                b.pos = cand.pos;
              }
            }
          }
        }
        for (const std::uint32_t q : seen_window) seen[q] = 0;
        for (const auto& b : best) slot[local_pid[b.pos]] = -1;
        // Claimant precursors in ascending id, the order of the former sort (it fixes the summation order below).
        std::sort(best.begin(), best.end(), [](const Candidate& x, const Candidate& y) { return x.precursor < y.precursor; });
        // Fragment sharing depends only on the two precursors, so it is tested once per claimant precursor; that
        // precursor's nearest window is its claimant window.
        for (const Candidate& b : best)
        {
          pairs_total += b.windows;
          Claim claim;
          if (sharedFragments(sorted[a.precursor], packed[b.pos].fragments, tol, claim.a_to_b, claim.b_to_a))
          {
            links_total += b.windows;
            claim.precursor = b.precursor;
            claim.distance = b.distance;
            claim.window = b.window;
            claim.pos = b.pos;
            claims.push_back(claim);
          }
        }
        if (claims.empty()) continue;
        ++with_claimant;
        std::array<double, K> contribution{};
        std::size_t claimants = 0;
        double ratios[K];
        for (std::size_t c = 0; c < claims.size(); ++c)
        {
          ++used_total;
          const Claim& claim = claims[c];
          const Packed& b = packed[claim.pos];
          std::size_t count = 0;
          for (std::size_t g = 0; g < K; ++g)
          {
            if (claim.b_to_a[g] < 0 && b.library[g] > 0.0 && b.core_area[g] > 0.0) ratios[count++] = b.core_area[g] / std::max(b.library[g], 1e-12);
          }
          if (count == 0) continue;
          const double scale = medianInPlace(ratios, count);
          const long shift = static_cast<long>(std::nearbyint((a.rt_apex - b.rt_apex) / rt_bin_seconds));
          const long lo = std::clamp(static_cast<long>(C0) + shift, 0L, NB);
          const long hi = std::clamp(static_cast<long>(C1) + shift, 0L, NB);
          const auto& cb = b.cumulative;
          const double own = cb[C1] - cb[C0];
          const double frac = own > 0.0 ? (cb[hi] - cb[lo]) / std::max(own, 1e-12) : 0.0;
          double total = 0.0;
          for (std::size_t f = 0; f < K; ++f)
          {
            if (claim.a_to_b[f] < 0) continue;
            const double e = scale * b.library[static_cast<std::size_t>(claim.a_to_b[f])] * frac;
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
      for (const std::uint32_t p : local_to_global) pid_map[p] = -1;
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
