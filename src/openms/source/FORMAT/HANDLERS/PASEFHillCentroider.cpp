// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/PASEFHillCentroider.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <unordered_set>

namespace OpenMS::Internal::PASEFHillCentroider
{
  namespace
  {
    // A hill is a trace of one ion across consecutive IM scans within a frame.
    // Arrays are parallel and ordered by IM scan (i.e. by time/scan index).
    struct HillInternal
    {
      std::vector<double> mz_values;
      std::vector<double> intensities;
      std::vector<double> im_values;
      double intensity_sum  = 0.0;
      double intensity_apex = 0.0;
      Size   apex_idx       = 0; // index into the parallel arrays
    };

    // 3-tap symmetric mean filter with implicit zero padding outside the
    // signal. Kernel length is 3 (window=1 in Biosaur2's meanFilter_).
    std::vector<double> meanFilter3(const std::vector<double>& data)
    {
      std::vector<double> out(data.size(), 0.0);
      const Size n = data.size();
      for (Size i = 0; i < n; ++i)
      {
        double sum = data[i];
        if (i > 0)     sum += data[i - 1];
        if (i + 1 < n) sum += data[i + 1];
        out[i] = sum / 3.0;
      }
      return out;
    }

    // Faithful port of Biosaur2's splitHills_ (Biosaur2Algorithm.cpp:1789),
    // axis-agnostic: operates on a hill's intensity sequence regardless of
    // whether the underlying axis is RT or IM-scan.
    std::vector<HillInternal> splitHills(const std::vector<HillInternal>& hills,
                                         double hvf, Size min_length)
    {
      std::vector<HillInternal> out;
      out.reserve(hills.size());

      for (const auto& hill : hills)
      {
        const Size hill_len = hill.intensities.size();
        if (hill_len < 2 * min_length)
        {
          out.push_back(hill);
          continue;
        }

        const std::vector<double> smoothed = meanFilter3(hill.intensities);

        std::vector<Size> min_idx_list;
        std::vector<Size> recheck_positions;

        const int min_len = static_cast<int>(min_length);
        const int c_len = static_cast<int>(hill_len) - min_len;
        int idx = min_len - 1;
        Size l_idx = 0;
        double min_val = 0.0;

        while (idx <= c_len)
        {
          if (!min_idx_list.empty() &&
              static_cast<Size>(idx) >= min_idx_list.back() + min_length)
          {
            l_idx = min_idx_list.back();
          }
          const double valley_intensity = smoothed[static_cast<Size>(idx)];
          if (valley_intensity <= 0.0) { ++idx; continue; }

          double left_max = 0.0;
          for (Size j = l_idx; j < static_cast<Size>(idx); ++j)
            left_max = std::max(left_max, smoothed[j]);
          if (left_max == 0.0) { ++idx; continue; }

          const double l_r = left_max / valley_intensity;
          if (l_r >= hvf)
          {
            double right_max = 0.0;
            for (Size j = static_cast<Size>(idx) + 1; j < hill_len; ++j)
              right_max = std::max(right_max, smoothed[j]);
            if (right_max > 0.0)
            {
              const double r_r = right_max / valley_intensity;
              if (r_r >= hvf)
              {
                const double mult_val = l_r * r_r;
                const int include_factor = (l_r > r_r) ? 1 : 0;
                const int candidate_pos_int = idx + include_factor;
                if (min_len <= candidate_pos_int && candidate_pos_int <= c_len)
                {
                  const Size candidate_pos = static_cast<Size>(candidate_pos_int);
                  if (min_idx_list.empty() ||
                      candidate_pos >= min_idx_list.back() + min_length)
                  {
                    min_idx_list.push_back(candidate_pos);
                    recheck_positions.push_back(static_cast<Size>(idx));
                    min_val = mult_val;
                  }
                  else if (mult_val > min_val)
                  {
                    min_idx_list.back()      = candidate_pos;
                    recheck_positions.back() = static_cast<Size>(idx);
                    min_val = mult_val;
                  }
                }
              }
            }
          }
          ++idx;
        }

        // Second pass: confirm right-side max above hvf for each candidate.
        std::vector<Size> final_splits;
        for (Size k = 0; k < min_idx_list.size(); ++k)
        {
          const Size recheck_idx = recheck_positions[k];
          const Size end_idx = (k + 1 < min_idx_list.size()) ? min_idx_list[k + 1]
                                                             : hill_len;
          if (recheck_idx + 1 >= end_idx || recheck_idx >= hill_len) continue;
          const double recheck_val = smoothed[recheck_idx];
          if (recheck_val <= 0.0) continue;
          double right_max = 0.0;
          for (Size j = recheck_idx + 1; j < end_idx; ++j)
            right_max = std::max(right_max, smoothed[j]);
          if (right_max / recheck_val >= hvf)
            final_splits.push_back(min_idx_list[k]);
        }

        if (final_splits.empty())
        {
          out.push_back(hill);
          continue;
        }

        std::vector<Size> boundaries;
        boundaries.reserve(final_splits.size() + 2);
        boundaries.push_back(0);
        for (Size pos : final_splits) boundaries.push_back(pos);
        boundaries.push_back(hill_len);

        for (Size s = 0; s + 1 < boundaries.size(); ++s)
        {
          const Size seg_start = boundaries[s];
          const Size seg_end   = boundaries[s + 1];
          if (seg_end <= seg_start) continue;
          const Size seg_len = seg_end - seg_start;
          if (seg_len < min_length) continue;

          HillInternal seg;
          seg.mz_values.assign(hill.mz_values.begin()   + seg_start,
                               hill.mz_values.begin()   + seg_end);
          seg.intensities.assign(hill.intensities.begin() + seg_start,
                                 hill.intensities.begin() + seg_end);
          seg.im_values.assign(hill.im_values.begin()   + seg_start,
                               hill.im_values.begin()   + seg_end);

          double isum = 0.0, iapex = 0.0;
          Size aidx = 0;
          for (Size k = 0; k < seg.intensities.size(); ++k)
          {
            isum += seg.intensities[k];
            if (seg.intensities[k] > iapex) { iapex = seg.intensities[k]; aidx = k; }
          }
          seg.intensity_sum  = isum;
          seg.intensity_apex = iapex;
          seg.apex_idx       = aidx;
          out.push_back(std::move(seg));
        }
      }
      return out;
    }

    // One tail per active hill: the most recent peak of that hill and how
    // many IM scans have passed since it was last extended (age 0 means the
    // hill was extended in the most recent scan; age N means we've skipped N
    // scans where no peak at this m/z appeared).
    struct HillTail
    {
      Size   hill_id;
      double mz;
      double intensity;
      int    fm;
      Size   age;  // 0 = just extended; bumped after each scan that doesn't extend it
    };

    // Link peaks of one IM scan against active hill tails. Tails older than
    // `max_scan_gap` have already been pruned by the caller. Updates `tails`
    // in place: extended hills get age=0 and (mz, intensity, fm) of the new
    // tail; unextended tails are aged by 1 after the scan.
    void linkScan(
        const std::vector<double>& cur_mz,
        const std::vector<double>& cur_int,
        const std::vector<int>&    cur_fm,
        double                     scan_im,
        double                     mz_tol_ppm,
        std::vector<HillInternal>& hills,
        std::vector<HillTail>&     tails,
        std::map<int, std::vector<int>>& tails_by_fm)
    {
      const Size m = cur_mz.size();
      auto append_to_hill = [&](Size hid, double mz, double inten, double im) {
        auto& h = hills[hid];
        h.mz_values.push_back(mz);
        h.intensities.push_back(inten);
        h.im_values.push_back(im);
        h.intensity_sum += inten;
        if (inten > h.intensity_apex)
        {
          h.intensity_apex = inten;
          h.apex_idx = h.intensities.size() - 1;
        }
      };

      if (m == 0)
      {
        // Empty scan: every active tail gets one scan older. Pruning past
        // max_scan_gap happens in the caller.
        for (auto& t : tails) ++t.age;
        return;
      }

      // Greedy assignment: process current-scan peaks in descending intensity.
      std::vector<Size> order(m);
      std::iota(order.begin(), order.end(), Size(0));
      std::sort(order.begin(), order.end(),
                [&](Size a, Size b) { return cur_int[a] > cur_int[b]; });

      // Which existing tails were claimed by a current-scan peak this round.
      std::unordered_set<int> claimed_tail_idx;
      // For each current-scan peak:
      //   v >= 0           : extended tails[v]
      //   in [INT_MIN+1, -1]: started a new hill, encoded as -hill_id - 1
      //   INT_MIN          : sentinel for "skipped" (e.g. pathological mz_cur <= 0)
      constexpr int k_unassigned = std::numeric_limits<int>::min();
      std::vector<int> cur_peak_to_tail(m, k_unassigned);

      for (Size pos = 0; pos < m; ++pos)
      {
        const Size i = order[pos];
        const int fm = cur_fm[i];
        const double mz_cur = cur_mz[i];
        if (mz_cur <= 0.0) continue;  // pathological input; skip rather than div-by-zero

        double best_intensity = 0.0;
        int best_tail_idx = -1;
        for (int fb : {fm, fm - 1, fm + 1})
        {
          auto it = tails_by_fm.find(fb);
          if (it == tails_by_fm.end()) continue;
          for (int tidx : it->second)
          {
            if (tidx < 0 || static_cast<Size>(tidx) >= tails.size()) continue;
            if (claimed_tail_idx.count(tidx)) continue;
            const auto& t = tails[static_cast<Size>(tidx)];
            const double dppm = std::abs(mz_cur - t.mz) / mz_cur * 1e6;
            if (dppm > mz_tol_ppm) continue;
            // Biosaur2 tie-break: '>=' prefers later (more recently created)
            // tails with equal intensity.
            if (t.intensity >= best_intensity)
            {
              best_intensity = t.intensity;
              best_tail_idx  = tidx;
            }
          }
        }

        if (best_tail_idx >= 0)
        {
          claimed_tail_idx.insert(best_tail_idx);
          cur_peak_to_tail[i] = best_tail_idx;
          const Size hid = tails[static_cast<Size>(best_tail_idx)].hill_id;
          append_to_hill(hid, mz_cur, cur_int[i], scan_im);
        }
        else
        {
          // Start a new hill.
          const Size new_hid = hills.size();
          hills.emplace_back();
          append_to_hill(new_hid, mz_cur, cur_int[i], scan_im);
          // Encode "new hill" as a sentinel; we'll create the tail below.
          cur_peak_to_tail[i] = -static_cast<int>(new_hid) - 1;
        }
      }

      // Age every tail by 1; replace claimed tails' (mz, intensity, fm)
      // with the current scan's peak data and reset their age to 0.
      // Tails not claimed by any current peak keep their old m/z/fm.
      for (auto& t : tails) ++t.age;

      // Apply current-scan extensions.
      for (Size i = 0; i < m; ++i)
      {
        const int v = cur_peak_to_tail[i];
        if (v == k_unassigned) continue;  // peak was skipped (pathological input)
        if (v >= 0)
        {
          // Extended existing tail.
          auto& t = tails[static_cast<Size>(v)];
          t.mz        = cur_mz[i];
          t.intensity = cur_int[i];
          t.fm        = cur_fm[i];
          t.age       = 0;
        }
        else
        {
          // Started a new hill: register a fresh tail. Decoding -v - 1 is
          // safe here because v is guaranteed >= INT_MIN+1 (the sentinel
          // was filtered out above).
          const Size new_hid = static_cast<Size>(-v - 1);
          HillTail t;
          t.hill_id   = new_hid;
          t.mz        = cur_mz[i];
          t.intensity = cur_int[i];
          t.fm        = cur_fm[i];
          t.age       = 0;
          tails.push_back(t);
        }
      }
    }

    // Drop tails older than max_scan_gap and rebuild the fm dict over the
    // surviving tails. Called once per scan after linkScan().
    //
    // Additional rule: a length-1 hill (a brand-new "starter") can never
    // bridge a gap. Such a tail must be extended in the very next scan, or
    // it is dropped. This prevents single-scan noise blips from chaining
    // across an empty scan into false length-2 hills. Once a hill has been
    // paired (length >= 2) further gaps up to max_scan_gap are tolerated.
    void pruneAndIndexTails(std::vector<HillTail>& tails,
                            const std::vector<HillInternal>& hills,
                            std::map<int, std::vector<int>>& tails_by_fm,
                            Size max_scan_gap)
    {
      tails.erase(std::remove_if(tails.begin(), tails.end(),
                                 [&hills, max_scan_gap](const HillTail& t)
                                 {
                                   const Size length =
                                       hills[t.hill_id].intensities.size();
                                   if (length < 2 && t.age > 0) return true;
                                   return t.age > max_scan_gap;
                                 }),
                  tails.end());
      tails_by_fm.clear();
      for (Size i = 0; i < tails.size(); ++i)
        tails_by_fm[tails[i].fm].push_back(static_cast<int>(i));
    }

  } // anon

  std::vector<Centroid>
  centroidFrame(const double* mz_values,
                const double* intensities,
                const double* im_values,
                std::size_t n,
                const Params& p_in,
                const std::uint32_t* scan_indices)
  {
    std::vector<Centroid> result;
    if (n == 0) return result;

    // splitHills indexes from `min_length - 1`, so a 0 here would underflow
    // a signed counter to -1 and read smoothed[SIZE_MAX]. A hill must have
    // at least one peak by definition, so clamp silently.
    Params p = p_in;
    if (p.min_hill_length == 0) p.min_hill_length = 1;

    // Filter by min intensity, find max m/z for the bin width.
    std::vector<Size> kept;
    kept.reserve(n);
    double max_mz = 0.0;
    for (Size i = 0; i < n; ++i)
    {
      if (intensities[i] < p.min_intensity) continue;
      if (mz_values[i]   > max_mz) max_mz = mz_values[i];
      kept.push_back(i);
    }
    if (kept.empty()) return result;

    // m/z bin width: htol_ppm * max_mz, mirroring Biosaur2 computeHillMzStep_.
    const double mz_step = (p.mz_tol_ppm * 1e-6) * max_mz;
    if (mz_step <= 0.0) return result;

    std::vector<HillInternal> hills;
    hills.reserve(kept.size() / 4 + 1);

    // Linking state across scans: active hill tails, possibly older than
    // the most recent scan (up to p.max_scan_gap scans ago).
    std::vector<HillTail>           tails;
    std::map<int, std::vector<int>> tails_by_fm;

    if (scan_indices != nullptr)
    {
      // Sort by (scan_index ASC, mz ASC) so peaks within one IM scan are
      // contiguous and m/z-sorted.
      std::sort(kept.begin(), kept.end(),
                [&](Size a, Size b) {
                  if (scan_indices[a] != scan_indices[b])
                    return scan_indices[a] < scan_indices[b];
                  return mz_values[a] < mz_values[b];
                });

      const std::uint32_t scan_min = scan_indices[kept.front()];
      const std::uint32_t scan_max = scan_indices[kept.back()];

      // Walk every scan index from min to max, calling linkScan with the
      // appropriate (possibly empty) batch so empty scans age tails.
      Size i0 = 0;
      for (std::uint32_t s = scan_min; s <= scan_max; ++s)
      {
        Size i1 = i0;
        while (i1 < kept.size() && scan_indices[kept[i1]] == s) ++i1;
        const Size m = i1 - i0;

        std::vector<double> cur_mz(m), cur_int(m);
        std::vector<int>    cur_fm(m);
        double scan_im = 0.0;
        for (Size k = 0; k < m; ++k)
        {
          const Size oi = kept[i0 + k];
          cur_mz[k]  = mz_values[oi];
          cur_int[k] = intensities[oi];
          cur_fm[k]  = static_cast<int>(mz_values[oi] / mz_step);
          scan_im    = im_values[oi];  // all peaks in one scan share IM
        }

        linkScan(cur_mz, cur_int, cur_fm, scan_im, p.mz_tol_ppm, hills,
                 tails, tails_by_fm);
        pruneAndIndexTails(tails, hills, tails_by_fm, p.max_scan_gap);

        i0 = i1;
      }
    }
    else
    {
      // No scan-index info: group by IM value with im_group_tolerance and
      // iterate only non-empty groups. Without knowing the actual scan-id
      // spacing we cannot tell two adjacent groups apart from two groups
      // separated by many empty scans, so max_scan_gap is forced to 0 in
      // this branch to avoid silently merging features across unknown gaps.
      // Callers wanting gap bridging must supply scan_indices.
      std::sort(kept.begin(), kept.end(),
                [&](Size a, Size b) {
                  if (im_values[a] != im_values[b])
                    return im_values[a] < im_values[b];
                  return mz_values[a] < mz_values[b];
                });

      Size i0 = 0;
      while (i0 < kept.size())
      {
        const double scan_im = im_values[kept[i0]];
        Size i1 = i0;
        while (i1 < kept.size() &&
               std::abs(im_values[kept[i1]] - scan_im) <= p.im_group_tolerance)
        {
          ++i1;
        }

        const Size m = i1 - i0;
        std::vector<double> cur_mz(m), cur_int(m);
        std::vector<int>    cur_fm(m);
        for (Size k = 0; k < m; ++k)
        {
          const Size oi = kept[i0 + k];
          cur_mz[k]  = mz_values[oi];
          cur_int[k] = intensities[oi];
          cur_fm[k]  = static_cast<int>(mz_values[oi] / mz_step);
        }

        linkScan(cur_mz, cur_int, cur_fm, scan_im, p.mz_tol_ppm, hills,
                 tails, tails_by_fm);
        // Force max_scan_gap=0 here: see comment above the loop.
        pruneAndIndexTails(tails, hills, tails_by_fm, 0);

        i0 = i1;
      }
    }

    // Split hills at intensity valleys.
    auto split = splitHills(hills, p.valley_factor, p.min_hill_length);

    // Drop short hills and collapse the rest to centroids.
    result.reserve(split.size());
    for (const auto& h : split)
    {
      const Size len = h.intensities.size();
      if (len < p.min_hill_length) continue;

      double w_sum = 0.0;
      double mz_w  = 0.0;
      double im_lower = h.im_values.front();
      double im_upper = h.im_values.front();
      double mz_lower = h.mz_values.front();
      double mz_upper = h.mz_values.front();
      for (Size k = 0; k < len; ++k)
      {
        w_sum += h.intensities[k];
        mz_w  += h.mz_values[k] * h.intensities[k];
        if (h.im_values[k] < im_lower) im_lower = h.im_values[k];
        if (h.im_values[k] > im_upper) im_upper = h.im_values[k];
        if (h.mz_values[k] < mz_lower) mz_lower = h.mz_values[k];
        if (h.mz_values[k] > mz_upper) mz_upper = h.mz_values[k];
      }
      if (w_sum <= 0.0) continue;

      Centroid c;
      c.mz        = mz_w / w_sum;
      c.intensity = h.intensity_sum;
      c.im_apex   = h.im_values[h.apex_idx];
      c.im_lower  = im_lower;
      c.im_upper  = im_upper;
      c.mz_lower  = mz_lower;
      c.mz_upper  = mz_upper;
      c.length    = len;
      result.push_back(c);
    }

    std::sort(result.begin(), result.end(),
              [](const Centroid& a, const Centroid& b) { return a.mz < b.mz; });

    return result;
  }

} // namespace OpenMS::Internal::PASEFHillCentroider
