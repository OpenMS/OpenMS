// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIAInterference.h"

#include <OpenMS/METADATA/Precursor.h>

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

    // Windows of linked isolation windows sorted by rt_start. For a window a the overlapping partners are the later
    // windows with rt_start <= a.rt_end and the earlier windows with rt_end >= a.rt_start; a running maximum of rt_end
    // bounds the backward scan. Each window computes its own features from its partners, so results are written by
    // exactly one iteration and do not depend on the thread count or schedule.
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
    std::vector<std::uint32_t> group_begin(m), group_end(m);
    std::vector<double> running_end(m);
    for (std::size_t g = 0; g < m;)
    {
      const auto iso = precursors[windows[order[g]].precursor].isolation;
      std::size_t h = g;
      double running = -std::numeric_limits<double>::infinity();
      while (h < m && precursors[windows[order[h]].precursor].isolation == iso)
      {
        running = std::max(running, windows[order[h]].rt_end);
        running_end[h] = running;
        ++h;
      }
      std::fill(group_begin.begin() + g, group_begin.begin() + h, static_cast<std::uint32_t>(g));
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

#pragma omp parallel num_threads(threads) reduction(+ : pairs_total, links_total, used_total, with_claimant)
    {
      std::vector<Claim> claims;
#pragma omp for schedule(dynamic, 512)
      for (std::ptrdiff_t pos = 0; pos < static_cast<std::ptrdiff_t>(m); ++pos)
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
        for (std::size_t q = static_cast<std::size_t>(pos); q-- > group_begin[pos];)
        {
          if (running_end[q] < a.rt_start) break;
          consider(q);
        }
        for (std::size_t q = static_cast<std::size_t>(pos) + 1; q < group_end[pos] && windows[order[q]].rt_start <= a.rt_end; ++q) consider(q);
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
    if (stats)
    {
      stats->pairs_tested = pairs_total;
      stats->links = links_total;
      stats->used_links = used_total;
      stats->windows_with_claimant = with_claimant;
    }
    return out;
  }

  // ------------------------------------------------------------------------------------------------ collector

  InterferenceCollector::InterferenceCollector(std::vector<double> precursor_mz,
                                               std::vector<std::array<double, InterferenceResolver::FRAGMENTS>> extraction_mz,
                                               std::vector<InterferenceResolver::Window>& windows, double rt_half_width, double extraction_ppm) :
    precursor_mz_(std::move(precursor_mz)),
    extraction_mz_(std::move(extraction_mz)),
    windows_(windows),
    half_(rt_half_width),
    ppm_(extraction_ppm)
  {
    const std::size_t np = precursor_mz_.size();
    if (extraction_mz_.size() != np || !(half_ > 0.0) || !(ppm_ > 0.0))
    {
      throw std::invalid_argument("InterferenceCollector: inconsistent precursor arrays or non-positive half width / tolerance.");
    }
    window_offsets_.assign(np + 1, 0);
    for (const auto& w : windows_)
    {
      if (w.precursor >= np) throw std::invalid_argument("InterferenceCollector: window precursor index out of range.");
      ++window_offsets_[w.precursor + 1];
    }
    for (std::size_t p = 0; p < np; ++p) window_offsets_[p + 1] += window_offsets_[p];
    window_index_.assign(windows_.size(), 0);
    std::vector<std::uint32_t> fill(window_offsets_.begin(), window_offsets_.end() - 1);
    rt_low_.assign(np, std::numeric_limits<double>::infinity());
    rt_high_.assign(np, -std::numeric_limits<double>::infinity());
    for (std::uint32_t i = 0; i < windows_.size(); ++i)
    {
      const auto p = windows_[i].precursor;
      window_index_[fill[p]++] = i;
      rt_low_[p] = std::min(rt_low_[p], windows_[i].rt_apex - half_);
      rt_high_[p] = std::max(rt_high_[p], windows_[i].rt_apex + half_);
    }
  }

  InterferenceCollector::Lane& InterferenceCollector::lane_(double lo, double hi)
  {
    const auto key = std::make_pair(lo, hi);
    auto it = lanes_.find(key);
    if (it != lanes_.end()) return it->second;
    Lane lane;
    for (std::uint32_t p = 0; p < precursor_mz_.size(); ++p)
    {
      if (window_offsets_[p + 1] > window_offsets_[p] && precursor_mz_[p] >= lo && precursor_mz_[p] <= hi) lane.precursors.push_back(p);
    }
    std::sort(lane.precursors.begin(), lane.precursors.end(), [&](std::uint32_t x, std::uint32_t y)
    {
      return rt_low_[x] != rt_low_[y] ? rt_low_[x] < rt_low_[y] : x < y;
    });
    return lanes_.emplace(key, std::move(lane)).first->second;
  }

  void InterferenceCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() != 2 || spectrum.empty() || spectrum.getPrecursors().empty()) return;
    const double rt = spectrum.getRT();
    if (!std::isfinite(rt)) return;
    if (rt < last_rt_) ++out_of_order_;
    last_rt_ = std::max(last_rt_, rt);
    if (!spectrum.isSorted()) spectrum.sortByPosition();
    bool used = false;
    for (const OpenMS::Precursor& precursor : spectrum.getPrecursors())
    {
      const double centre = precursor.getMZ();
      const double lower = precursor.getIsolationWindowLowerOffset();
      const double upper = precursor.getIsolationWindowUpperOffset();
      if (!std::isfinite(centre) || !std::isfinite(lower) || !std::isfinite(upper) || lower < 0.0 || upper < 0.0 || !(lower + upper > 0.0)) continue;
      Lane& lane = lane_(centre - lower, centre + upper);
      // Activate precursors whose first window has started; drop those whose last window has ended.
      while (lane.cursor < lane.precursors.size() && rt_low_[lane.precursors[lane.cursor]] < rt) lane.active.push_back(lane.precursors[lane.cursor++]);
      std::size_t kept = 0;
      for (const std::uint32_t p : lane.active)
      {
        if (rt_high_[p] > rt) lane.active[kept++] = p;
      }
      lane.active.resize(kept);
      for (const std::uint32_t p : lane.active)
      {
        // Only precursors with a window whose grid contains this RT (|rt - apex| < half).
        bool open = false;
        for (std::uint32_t k = window_offsets_[p]; k < window_offsets_[p + 1] && !open; ++k)
        {
          open = std::abs(rt - windows_[window_index_[k]].rt_apex) < half_;
        }
        if (!open) continue;
        // Fragment sums of this spectrum, once per precursor.
        std::array<double, InterferenceResolver::FRAGMENTS> sum{};
        bool any = false;
        for (std::size_t s = 0; s < InterferenceResolver::FRAGMENTS; ++s)
        {
          const double target = extraction_mz_[p][s];
          if (!(target > 0.0)) continue;
          const double half = target * ppm_ * 1e-6;
          double total = 0.0;
          for (auto point = spectrum.MZBegin(target - half); point != spectrum.end() && point->getMZ() <= target + half; ++point)
          {
            // |ppm| < tolerance, open interval, as build12.py.
            if (!(std::abs((point->getMZ() - target) / target * 1e6) < ppm_)) continue;
            total += point->getIntensity();
          }
          sum[s] = total;
          any = any || total != 0.0;
        }
        if (!any) continue;
        used = true;
        for (std::uint32_t k = window_offsets_[p]; k < window_offsets_[p + 1]; ++k)
        {
          auto& w = windows_[window_index_[k]];
          double d = rt - w.rt_apex;
          if (!(std::abs(d) < half_)) continue;
          // The apex is itself a scan RT and sits on the edge between the two central bins. Snap that scan onto the edge
          // (the apex text in -candidates is rounded), then floor: an edge point goes to the upper bin, as build12.py does
          // for exact values, so the result does not depend on float noise.
          if (std::abs(d) < APEX_SNAP_S) d = 0.0;
          const double position = (d + half_) / (2.0 * half_) * InterferenceResolver::BINS;
          const long bin = std::clamp(static_cast<long>(std::floor(position)), 0L, static_cast<long>(InterferenceResolver::BINS) - 1);
          const bool core = bin >= static_cast<long>(InterferenceResolver::CORE_BEGIN) && bin < static_cast<long>(InterferenceResolver::CORE_END);
          double total = 0.0;
          for (std::size_t s = 0; s < InterferenceResolver::FRAGMENTS; ++s)
          {
            total += sum[s];
            if (core) w.core_area[s] += sum[s];
          }
          w.profile[static_cast<std::size_t>(bin)] += total;
        }
      }
    }
    if (used) ++spectra_used_;
  }

  std::vector<std::pair<double, double>> InterferenceCollector::isolationWindows() const
  {
    std::vector<std::pair<double, double>> out;
    for (const auto& [key, lane] : lanes_) out.push_back(key);
    return out;
  }

  std::vector<std::int32_t> InterferenceCollector::isolationIndices() const
  {
    const auto iso = isolationWindows();
    std::vector<std::int32_t> out(precursor_mz_.size(), -1);
    for (std::size_t p = 0; p < precursor_mz_.size(); ++p)
    {
      const double mz = precursor_mz_[p];
      // Window with the largest lower edge <= mz; valid only if mz < its upper edge.
      auto it = std::upper_bound(iso.begin(), iso.end(), mz, [](double v, const std::pair<double, double>& w) { return v < w.first; });
      if (it == iso.begin()) continue;
      --it;
      if (mz < it->second) out[p] = static_cast<std::int32_t>(it - iso.begin());
    }
    return out;
  }
}
