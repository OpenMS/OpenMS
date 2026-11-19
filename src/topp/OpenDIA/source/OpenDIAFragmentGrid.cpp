// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIAFragmentGrid.h>

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
    constexpr std::size_t K = FragmentGrid::FRAGMENTS;
    constexpr std::size_t B = FragmentGrid::BINS;
    constexpr std::size_t C0 = FragmentGrid::CORE_BEGIN;
    constexpr std::size_t C1 = FragmentGrid::CORE_END;
    constexpr std::size_t CELLS = K * B;
    /// Scans closer than this to a window apex are the apex scan (the apex is itself a scan RT).
    constexpr double APEX_SNAP_S = 1e-3;

    /// NaN to 0 and infinities to the float range, as numpy.nan_to_num followed by a float32 cast.
    float toFloat(double v)
    {
      if (std::isnan(v)) return 0.0f;
      constexpr double big = std::numeric_limits<float>::max();
      return static_cast<float>(std::clamp(v, -big, big));
    }

    double medianInPlace(double* v, std::size_t n)
    {
      std::sort(v, v + n);
      return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    }
  }

  const std::array<const char*, FragmentGrid::EVIDENCE>& FragmentGrid::evidenceNames()
  {
    static const std::array<const char*, EVIDENCE> names = {"f12_support", "f12_coelution", "f12_libcos",
                                                            "f12_extra_support", "f12_extra_coelution", "f12_extra_share"};
    return names;
  }

  const std::array<const char*, FragmentGrid::DEFICIT>& FragmentGrid::deficitNames()
  {
    static const std::array<const char*, DEFICIT> names = {"d_max", "d6_max", "d_lib"};
    return names;
  }

  std::array<float, FragmentGrid::EVIDENCE> FragmentGrid::evidence(const double* x, const Values& li, const Flags& search)
  {
    std::array<double, B> total{};
    for (std::size_t f = 0; f < K; ++f)
      for (std::size_t b = 0; b < B; ++b) total[b] += x[f * B + b];
    double support = 0.0, r_sum = 0.0, extra_support = 0.0, extra_r = 0.0, dot = 0.0, aa = 0.0, ll = 0.0, extra_area = 0.0, all_area = 0.0;
    std::size_t n_present = 0, n_extra = 0;
    for (std::size_t f = 0; f < K; ++f)
    {
      const double* xf = x + f * B;
      double area = 0.0, core_max = 0.0;
      for (std::size_t b = C0; b < C1; ++b)
      {
        area += xf[b];
        core_max = std::max(core_max, xf[b]);
      }
      const double a = std::sqrt(area), l = std::sqrt(li[f]);
      dot += a * l;
      aa += a * a;
      ll += l * l;
      all_area += area;
      if (!(li[f] > 0.0)) continue;
      const bool extra = !search[f];
      ++n_present;
      n_extra += extra;
      // Pearson correlation of the fragment with the sum of the other fragments over the grid.
      double mx = 0.0, mo = 0.0;
      for (std::size_t b = 0; b < B; ++b)
      {
        mx += xf[b];
        mo += total[b] - xf[b];
      }
      mx /= B;
      mo /= B;
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t b = 0; b < B; ++b)
      {
        const double dx = xf[b] - mx, dy = total[b] - xf[b] - mo;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
      }
      const double r = sxy / std::sqrt(sxx * syy + 1e-12);
      const bool ok = core_max > 0.0 && r > 0.5;
      support += ok;
      r_sum += r;
      if (extra)
      {
        extra_support += ok;
        extra_r += r;
        extra_area += area;
      }
    }
    return {static_cast<float>(support), toFloat(r_sum / static_cast<double>(std::max<std::size_t>(n_present, 1))),
            toFloat(dot / (std::sqrt(aa) * std::sqrt(ll) + 1e-12)), static_cast<float>(extra_support),
            toFloat(extra_r / static_cast<double>(std::max<std::size_t>(n_extra, 1))), toFloat(std::log10((extra_area + 1.0) / (all_area + 1.0)))};
  }

  std::array<float, FragmentGrid::DEFICIT> FragmentGrid::deficit(const Values& area, const Values& li, const Flags& search, double noise)
  {
    // One pass per fragment mask: all library fragments, then the search-library fragments.
    const auto run = [&](bool search_only, double& d_max, double& d_lib)
    {
      double ratios[K];
      std::size_t n_sig = 0;
      bool any = false;
      for (std::size_t f = 0; f < K; ++f)
      {
        if (!(li[f] > 0.0) || (search_only && !search[f])) continue;
        any = true;
        if (area[f] > 0.0) ratios[n_sig++] = area[f] / std::max(li[f], 1e-12);
      }
      const double scale = n_sig ? medianInPlace(ratios, n_sig) : 0.0;
      d_max = any ? -std::numeric_limits<double>::infinity() : 0.0;
      double big = 0.0, all = 0.0;
      for (std::size_t f = 0; f < K; ++f)
      {
        if (!(li[f] > 0.0) || (search_only && !search[f])) continue;
        const double d = std::log10(scale * li[f] + noise) - std::log10(area[f] + noise);
        d_max = std::max(d_max, d);
        all += li[f];
        if (d > 1.0) big += li[f];
      }
      d_lib = big / std::max(all, 1e-12);
    };
    double d_all = 0.0, lib_all = 0.0, d_search = 0.0, lib_search = 0.0;
    run(false, d_all, lib_all);
    run(true, d_search, lib_search);
    return {toFloat(d_all), toFloat(d_search), toFloat(lib_all)};
  }

  double FragmentGrid::noiseLevel(const std::vector<InterferenceResolver::Window>& windows, const std::vector<InterferenceResolver::Precursor>& precursors)
  {
    std::vector<double> v;
    for (const auto& w : windows)
      for (std::size_t f = 0; f < K; ++f)
        if (w.core_area[f] > 0.0 && precursors[w.precursor].library[f] > 0.0) v.push_back(w.core_area[f]);
    if (v.empty()) return 1.0;
    const std::size_t h = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + h, v.end());
    if (v.size() % 2) return v[h];
    return 0.5 * (v[h] + *std::max_element(v.begin(), v.begin() + h));
  }

  // ------------------------------------------------------------------------------------------------ collector

  FragmentGridCollector::FragmentGridCollector(const std::vector<Precursor>& precursors, std::vector<InterferenceResolver::Window>& windows,
                                               double rt_half_width, double extraction_ppm) :
    precursors_(precursors), windows_(windows), half_(rt_half_width), ppm_(extraction_ppm)
  {
    const std::size_t np = precursors_.size();
    if (!(half_ > 0.0) || !(ppm_ > 0.0)) throw std::invalid_argument("FragmentGridCollector: half width and tolerance must be positive.");
    window_offsets_.assign(np + 1, 0);
    for (const auto& w : windows_)
    {
      if (w.precursor >= np) throw std::invalid_argument("FragmentGridCollector: window precursor index out of range.");
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
      windows_[i].core_area.fill(0.0);
      windows_[i].profile.fill(0.0);
    }
    by_apex_.resize(windows_.size());
    std::iota(by_apex_.begin(), by_apex_.end(), 0u);
    std::sort(by_apex_.begin(), by_apex_.end(), [&](std::uint32_t a, std::uint32_t b)
    {
      return windows_[a].rt_apex != windows_[b].rt_apex ? windows_[a].rt_apex < windows_[b].rt_apex : a < b;
    });
    slot_.assign(windows_.size(), -1);
    evidence_.assign(windows_.size(), {});
  }

  void FragmentGridCollector::enableIonMobility(double half_width, bool array_as_vssc)
  {
    if (!(half_width > 0.0) || !std::isfinite(half_width)) throw std::invalid_argument("FragmentGridCollector: mobility half width must be positive.");
    if (last_rt_ > -1e300) throw std::logic_error("FragmentGridCollector: enableIonMobility() after the first spectrum.");
    im_half_ = half_width;
    im_array_as_vssc_ = array_as_vssc;
  }

  FragmentGridCollector::Lane& FragmentGridCollector::lane_(double lo, double hi)
  {
    const auto key = std::make_pair(lo, hi);
    auto it = lanes_.find(key);
    if (it != lanes_.end()) return it->second;
    Lane lane;
    for (std::uint32_t p = 0; p < precursors_.size(); ++p)
    {
      if (window_offsets_[p + 1] > window_offsets_[p] && precursors_[p].mz >= lo && precursors_[p].mz <= hi) lane.precursors.push_back(p);
    }
    std::sort(lane.precursors.begin(), lane.precursors.end(), [&](std::uint32_t x, std::uint32_t y)
    {
      return rt_low_[x] != rt_low_[y] ? rt_low_[x] < rt_low_[y] : x < y;
    });
    return lanes_.emplace(key, std::move(lane)).first->second;
  }

  double* FragmentGridCollector::grid_(std::uint32_t w)
  {
    if (slot_[w] < 0)
    {
      if (free_.empty())
      {
        free_.push_back(static_cast<std::int32_t>(pool_.size() / CELLS));
        pool_.resize(pool_.size() + CELLS, 0.0);
      }
      slot_[w] = free_.back();
      free_.pop_back();
      peak_open_ = std::max(peak_open_, ++open_);
    }
    return pool_.data() + static_cast<std::size_t>(slot_[w]) * CELLS;
  }

  void FragmentGridCollector::close_(std::uint32_t w)
  {
    auto& window = windows_[w];
    const auto& p = precursors_[window.precursor];
    if (slot_[w] < 0)
    {
      static const std::array<double, CELLS> zero{};
      evidence_[w] = FragmentGrid::evidence(zero.data(), p.library, p.search);
      return;
    }
    double* g = pool_.data() + static_cast<std::size_t>(slot_[w]) * CELLS;
    // Bins are stored as float32, as the 12-fragment arrays the features were validated on.
    for (std::size_t c = 0; c < CELLS; ++c) g[c] = static_cast<double>(static_cast<float>(g[c]));
    for (std::size_t f = 0; f < K; ++f)
    {
      for (std::size_t b = 0; b < B; ++b)
      {
        window.profile[b] += g[f * B + b];
        if (b >= C0 && b < C1) window.core_area[f] += g[f * B + b];
      }
    }
    evidence_[w] = FragmentGrid::evidence(g, p.library, p.search);
    std::fill(g, g + CELLS, 0.0);
    free_.push_back(slot_[w]);
    slot_[w] = -1;
    --open_;
  }

  void FragmentGridCollector::closeUntil_(double rt)
  {
    // A spectrum at rt contributes only to windows with |rt - apex| < half, so every window with apex + half <= rt is done.
    while (close_cursor_ < by_apex_.size() && windows_[by_apex_[close_cursor_]].rt_apex + half_ <= rt) close_(by_apex_[close_cursor_++]);
  }

  void FragmentGridCollector::consumeSpectrum(MSSpectrum& spectrum)
  {
    if (spectrum.getMSLevel() != 2 || spectrum.empty() || spectrum.getPrecursors().empty()) return;
    const double rt = spectrum.getRT();
    if (!std::isfinite(rt)) return;
    if (finished_) throw std::logic_error("FragmentGridCollector: spectrum after finish().");
    if (rt < last_rt_) throw std::invalid_argument("FragmentGridCollector: MS2 spectra must be in RT order.");
    last_rt_ = rt;
    closeUntil_(rt);
    if (!spectrum.isSorted()) spectrum.sortByPosition();
    mz_index_.build(spectrum);
    // Per-point (or scalar) 1/K0 of this spectrum, when mobility filtering is on and the data carry VSSC mobility.
    const DataArrays::FloatDataArray* mobilities = nullptr;
    double scalar_im = 0.0;
    if (im_half_ > 0.0)
    {
      if (spectrum.containsIMData())
      {
        const auto [index, unit] = spectrum.getIMData();
        if (spectrum.getFloatDataArrays()[index].size() == spectrum.size() &&
            (unit == DriftTimeUnit::VSSC || (im_array_as_vssc_ && unit != DriftTimeUnit::CCS)))
        {
          mobilities = &spectrum.getFloatDataArrays()[index];
        }
      }
      if (mobilities == nullptr && spectrum.getDriftTimeUnit() == DriftTimeUnit::VSSC && std::isfinite(spectrum.getDriftTime()) &&
          spectrum.getDriftTime() > 0.0)
      {
        scalar_im = spectrum.getDriftTime();
      }
    }
    const bool spectrum_has_im = mobilities != nullptr || scalar_im > 0.0;
    bool used = false;
    // Serially: the precursors this spectrum reaches, in the former visiting order. In parallel: their fragment sums
    // (independent, read-only on the spectrum). Serially again: the grid additions in that same order, so every grid
    // cell receives the same sums in the same order as before (thread-count invariant, byte-identical).
    reached_.clear();
    for (const OpenMS::Precursor& precursor : spectrum.getPrecursors())
    {
      const double centre = precursor.getMZ();
      const double lower = precursor.getIsolationWindowLowerOffset();
      const double upper = precursor.getIsolationWindowUpperOffset();
      if (!std::isfinite(centre) || !std::isfinite(lower) || !std::isfinite(upper) || lower < 0.0 || upper < 0.0 || !(lower + upper > 0.0)) continue;
      Lane& lane = lane_(centre - lower, centre + upper);
      while (lane.cursor < lane.precursors.size() && rt_low_[lane.precursors[lane.cursor]] < rt) lane.active.push_back(lane.precursors[lane.cursor++]);
      std::size_t kept = 0;
      for (const std::uint32_t p : lane.active)
      {
        if (rt_high_[p] > rt) lane.active[kept++] = p;
      }
      lane.active.resize(kept);
      for (const std::uint32_t p : lane.active)
      {
        bool open = false;
        for (std::uint32_t k = window_offsets_[p]; k < window_offsets_[p + 1] && !open; ++k)
        {
          open = std::abs(rt - windows_[window_index_[k]].rt_apex) < half_;
        }
        if (!open) continue;
        const double precursor_im = precursors_[p].im;
        const bool filter_im = spectrum_has_im && precursor_im > 0.0;
        if (filter_im && mobilities == nullptr && std::abs(scalar_im - precursor_im) > im_half_) continue;
        reached_.push_back(p);
      }
    }
    sums_.resize(reached_.size());
    any_.assign(reached_.size(), 0);
#pragma omp parallel for schedule(dynamic, 64) num_threads(threads_) if (threads_ > 1 && reached_.size() >= 256)
    for (std::ptrdiff_t r = 0; r < static_cast<std::ptrdiff_t>(reached_.size()); ++r)
    {
      const std::uint32_t p = reached_[static_cast<std::size_t>(r)];
      // Fragment sums of this spectrum, once per precursor: |ppm| < tolerance (open interval), and, with mobility
      // filtering, |1/K0 - precursor 1/K0| <= im half width.
      const double precursor_im = precursors_[p].im;
      const bool filter_im = spectrum_has_im && precursor_im > 0.0;
      std::array<double, K>& sum = sums_[static_cast<std::size_t>(r)];
      sum.fill(0.0);
      bool any = false;
      for (std::size_t s = 0; s < K; ++s)
      {
        const double target = precursors_[p].extraction_mz[s];
        if (!(target > 0.0)) continue;
        const double width = target * ppm_ * 1e-6;
        double total = 0.0;
        for (auto point = mz_index_.lowerBound(target - width); point != spectrum.end() && point->getMZ() <= target + width; ++point)
        {
          if (!(std::abs((point->getMZ() - target) / target * 1e6) < ppm_)) continue;
          if (filter_im && mobilities != nullptr && !(std::abs((*mobilities)[point - spectrum.begin()] - precursor_im) <= im_half_)) continue;
          total += point->getIntensity();
        }
        sum[s] = total;
        any = any || total != 0.0;
      }
      any_[static_cast<std::size_t>(r)] = any ? 1 : 0;
    }
    for (std::size_t r = 0; r < reached_.size(); ++r)
    {
      if (!any_[r]) continue;
      const std::uint32_t p = reached_[r];
      const std::array<double, K>& sum = sums_[r];
      used = true;
      {
        for (std::uint32_t k = window_offsets_[p]; k < window_offsets_[p + 1]; ++k)
        {
          const std::uint32_t w = window_index_[k];
          double d = rt - windows_[w].rt_apex;
          if (!(std::abs(d) < half_)) continue;
          // The apex sits on the edge of the two central bins; snap its own scan onto the edge and floor (upper bin).
          if (std::abs(d) < APEX_SNAP_S) d = 0.0;
          const long bin = std::clamp(static_cast<long>(std::floor((d + half_) / (2.0 * half_) * B)), 0L, static_cast<long>(B) - 1);
          double* g = grid_(w);
          for (std::size_t s = 0; s < K; ++s) g[s * B + static_cast<std::size_t>(bin)] += sum[s];
        }
      }
    }
    if (used) ++spectra_used_;
  }

  void FragmentGridCollector::finish()
  {
    if (finished_) return;
    closeUntil_(std::numeric_limits<double>::infinity());
    finished_ = true;
    std::vector<double>().swap(pool_);
  }

  std::vector<std::int32_t> FragmentGridCollector::isolationIndices() const
  {
    std::vector<std::pair<double, double>> iso;
    for (const auto& [key, lane] : lanes_) iso.push_back(key);
    std::vector<std::int32_t> out(precursors_.size(), -1);
    for (std::size_t p = 0; p < precursors_.size(); ++p)
    {
      const double mz = precursors_[p].mz;
      auto it = std::upper_bound(iso.begin(), iso.end(), mz, [](double v, const std::pair<double, double>& w) { return v < w.first; });
      if (it == iso.begin()) continue;
      --it;
      if (mz < it->second) out[p] = static_cast<std::int32_t>(it - iso.begin());
    }
    return out;
  }
}
