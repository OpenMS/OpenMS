// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include "NuXLConstants.h"
#include "NuXLFragmentAdductDefinition.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Search-time lookup structures for the OpenNuXL database search.

    Everything in here only changes how values are looked up, never what is computed: the lookups return the
    same results as MSSpectrum::findNearest or std::lower_bound / std::upper_bound on the full containers.
  */

  /// Per-spectrum data that is computed once after preprocessing and reused for every candidate peptide.
  struct NuXLSpectrumIndex
  {
    const MSSpectrum* spec = nullptr;
    std::vector<double> mz;        ///< peak m/z (copy, contiguous)
    std::vector<float> intensity;  ///< peak intensities (copy, contiguous)
    std::vector<int> charge;       ///< peak charges (IntegerDataArray NuXLConstants::IA_CHARGE_INDEX)
    static constexpr double BIN = 2.0;
    static constexpr double INV_BIN = 1.0 / BIN;
    std::vector<uint16_t> lb;      ///< lb[k] = index of the first peak with m/z >= k * BIN
    float marker_ions_sub_score = 0; ///< MorpheusScore of the spectrum against all marker ions (independent of the candidate)
    std::vector<uint8_t> has_tag_for_fa_id; ///< fragment adduct name (by registry id) is listed in StringDataArray 0 (nucleotide tags)

    void build(const MSSpectrum& s)
    {
      spec = &s;
      const size_t n = s.size();
      mz.resize(n); intensity.resize(n); charge.resize(n);
      const auto& z = s.getIntegerDataArrays()[NuXLConstants::IA_CHARGE_INDEX];
      for (size_t i = 0; i != n; ++i) { mz[i] = s[i].getMZ(); intensity[i] = s[i].getIntensity(); charge[i] = z[i]; }
      if (n == 0 || n >= std::numeric_limits<uint16_t>::max()) { lb.clear(); return; } // fall back to binary search
      const size_t nb = (size_t)(mz.back() * INV_BIN) + 2;
      lb.resize(nb + 1);
      size_t p = 0;
      for (size_t k = 0; k <= nb; ++k)
      {
        const double edge = k * BIN;
        while (p < n && mz[p] < edge) ++p;
        lb[k] = (uint16_t)p;
      }
    }

    /// same result as std::lower_bound on the peak m/z
    inline size_t lowerBound(double x) const
    {
      if (lb.empty()) return std::lower_bound(mz.begin(), mz.end(), x) - mz.begin();
      const long nb = (long)lb.size() - 1;
      const long k = (long)(x * INV_BIN);
      size_t lo = lb[std::clamp(k - 1, 0L, nb)];
      const size_t hi = lb[std::clamp(k + 2, 0L, nb)];
      while (lo < hi && mz[lo] < x) ++lo;
      return lo;
    }

    /// same result as MSSpectrum::findNearest(x)
    inline size_t nearest(double x) const
    {
      const size_t n = mz.size();
      const size_t it = lowerBound(x);
      if (it == 0) return 0;
      if (it == n) return n - 1;
      if (std::fabs(mz[it] - x) < std::fabs(mz[it - 1] - x)) return it;
      return it - 1;
    }
  };

  /// Dense ids for the fragment adduct definitions used in the search, a rank consistent with
  /// NuXLFragmentAdductDefinition::operator< (avoids copies and toString() calls when sorting) and
  /// the hoisted lookups into the ambiguous-mass lists.
  struct NuXLFragmentAdductRegistry
  {
    std::vector<const NuXLFragmentAdductDefinition*> by_id;
    std::unordered_map<std::string, int> name_id; ///< has_tag only depends on the name
    std::vector<int> id_to_name_id;
    std::vector<int> rank;
    std::vector<const std::vector<std::vector<double>>*> block_lists_by_name_id;

    int add(const NuXLFragmentAdductDefinition* f)
    {
      by_id.push_back(f);
      auto it = name_id.emplace(f->name, (int)name_id.size()).first;
      id_to_name_id.push_back(it->second);
      return (int)by_id.size() - 1;
    }

    void finalize(const std::map<std::string, std::vector<std::vector<double>>>& block)
    {
      std::vector<int> order(by_id.size());
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](int a, int b) { return *by_id[a] < *by_id[b]; });
      rank.assign(by_id.size(), 0);
      int r = 0;
      for (size_t i = 0; i != order.size(); ++i)
      {
        if (i > 0 && *by_id[order[i - 1]] < *by_id[order[i]]) ++r;
        rank[order[i]] = r;
      }
      block_lists_by_name_id.assign(name_id.size(), nullptr);
      for (const auto& [name, id] : name_id)
      {
        auto it = block.find(name);
        if (it != block.end()) block_lists_by_name_id[id] = &it->second;
      }
    }

    void annotateSpectrum(NuXLSpectrumIndex& si) const
    {
      si.has_tag_for_fa_id.assign(name_id.size(), 0);
      if (si.spec->getStringDataArrays().empty()) return;
      const auto& tags = si.spec->getStringDataArrays()[0];
      for (const auto& [name, id] : name_id)
      {
        si.has_tag_for_fa_id[id] = std::find(tags.begin(), tags.end(), name) != tags.end();
      }
    }
  };

  /// Fragment adducts of one cross-linked nucleotide with their registry ids and ranks.
  struct NuXLIndexedAdductList
  {
    char nucleotide = '?';
    const std::vector<NuXLFragmentAdductDefinition>* fas = nullptr;
    std::vector<int> name_ids;
    std::vector<int> ranks;
  };

  /// Sorted precursor masses (same order as the multimap they are created from) with a bucket index.
  struct NuXLPrecursorIndex
  {
    std::vector<double> mass;
    std::vector<std::pair<Size, int>> scan_and_isotope;
    double lo = 0;
    static constexpr double BIN = 0.05;
    std::vector<uint32_t> bucket; ///< bucket[k] = index of the first precursor with mass >= lo + k * BIN

    template <class MultiMap>
    void build(const MultiMap& mm)
    {
      for (const auto& e : mm) { mass.push_back(e.first); scan_and_isotope.push_back(e.second); }
      lo = mass.empty() ? 0 : mass.front();
      const size_t nb = mass.empty() ? 1 : (size_t)((mass.back() - lo) / BIN) + 2;
      bucket.resize(nb + 1);
      size_t q = 0;
      for (size_t k = 0; k <= nb; ++k)
      {
        const double edge = lo + k * BIN;
        while (q < mass.size() && mass[q] < edge) ++q;
        bucket[k] = (uint32_t)q;
      }
    }
    std::pair<size_t, size_t> range_(double x) const
    {
      const long nb = (long)bucket.size() - 1;
      const long k = (long)std::floor((x - lo) / BIN);
      return {bucket[std::clamp(k - 1, 0L, nb)], bucket[std::clamp(k + 2, 0L, nb)]};
    }
    /// == multimap::lower_bound / upper_bound (as index)
    size_t lowerBound(double x) const { auto [a, b] = range_(x); return std::lower_bound(mass.begin() + a, mass.begin() + b, x) - mass.begin(); }
    size_t upperBound(double x) const { auto [a, b] = range_(x); return std::upper_bound(mass.begin() + a, mass.begin() + b, x) - mass.begin(); }
  };

  /// logfactorial with the exact same summation sequence as the loop it replaces
  struct NuXLLogFactorialTable
  {
    std::vector<double> t;
    NuXLLogFactorialTable()
    {
      t.resize(4096);
      t[0] = 0; t[1] = 0;
      double z(0);
      for (size_t x = 2; x < t.size(); ++x) { z += std::log((double)x); t[x] = z; }
    }
    double operator()(UInt x) const
    {
      if (x < t.size()) return t[x];
      const double dx = static_cast<double>(x);
      double z(0);
      for (double y = 2; y <= dx; ++y) { z += std::log(y); }
      return z;
    }
  };

  /// per-thread buffers (no heap allocations per candidate) and memoized match odds
  struct NuXLScoringScratch
  {
    std::vector<double> intensity_sum, b, y, a, intensity_xl, b_xl, y_xl, a_xl;
    std::vector<uint8_t> peak_matched, peak_matched_xl;
    static constexpr int N_ION_TYPES = 6; ///< shifted b, a, y (pre-selection) and b-H2O, y-H2O, y-NH3 (neutral loss peaks)
    struct FirstPassEntry { double mass; uint32_t cnt[3][2]; size_t off; };
    std::vector<FirstPassEntry> fa_entries; ///< first pass of the shifted ion scoring, per fragment adduct mass
    std::vector<int32_t> fa_idx;            ///< matched peak index (or -1) per [entry][ion type][charge][ion]
    std::vector<std::tuple<size_t, size_t, const NuXLFragmentAdductDefinition*, int>> matches_z_fa;
    std::vector<std::vector<double>> match_odds; ///< [N][matches], NaN = not yet computed
    void resetCandidate() { fa_entries.clear(); fa_idx.clear(); }

    template <class F>
    double matchOdds(size_t N, size_t matches, double p, F f)
    {
      if (N >= match_odds.size()) match_odds.resize(N + 1);
      auto& row = match_odds[N];
      if (row.empty()) row.assign(N + 1, std::numeric_limits<double>::quiet_NaN());
      if (matches > N) return f(N, matches, p);
      double& v = row[matches];
      if (std::isnan(v)) v = f(N, matches, p);
      return v;
    }
  };
} // namespace OpenMS
