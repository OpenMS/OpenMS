// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
#include <OpenMS/FORMAT/HANDLERS/PASEFHillCentroider.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <opentims++/opentims.h>
#include <opentims++/tof2mz_converter.h>
#include <opentims++/scan2inv_ion_mobility_converter.h>
#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
#include <SQLiteCpp/SQLiteCpp.h>

#include <memory>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <optional>
#include <stdexcept>

// Open-addressing flat hash map (Boost 1.81+). Used in FrameAggregator's
// hot path because it eliminates per-node allocations and pointer chasing,
// giving 2-3x faster lookup/insert vs std::unordered_map for small values.
#include <boost/unordered/unordered_flat_map.hpp>

namespace OpenMS
{

  // --------------------------------------------------------------------------
  // TOF-domain spectrum processing (ported from timsrust)
  //
  // Pipeline: raw peaks → group_and_sum → smooth → centroid
  // All processing in sparse TOF-index space (integers), before m/z conversion.
  // --------------------------------------------------------------------------

  /// Sort TOF indices and sum intensities at identical TOF bins.
  /// Merges duplicate detections across scans/frames into single entries.
  /// Input arrays are modified in-place; returns the new size.
  static size_t groupAndSum(std::vector<uint32_t>& tofs, std::vector<uint64_t>& intensities)
  {
    if (tofs.empty()) return 0;
    // Sort by TOF index (co-sort intensities)
    std::vector<size_t> order(tofs.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(),
              [&tofs](size_t a, size_t b) { return tofs[a] < tofs[b]; });

    std::vector<uint32_t> sorted_tofs(tofs.size());
    std::vector<uint64_t> sorted_int(intensities.size());
    for (size_t i = 0; i < order.size(); ++i)
    {
      sorted_tofs[i] = tofs[order[i]];
      sorted_int[i] = intensities[order[i]];
    }

    // Merge consecutive entries with same TOF index
    size_t write = 0;
    for (size_t read = 0; read < sorted_tofs.size(); ++read)
    {
      if (write > 0 && sorted_tofs[read] == tofs[write - 1])
      {
        intensities[write - 1] += sorted_int[read];
      }
      else
      {
        tofs[write] = sorted_tofs[read];
        intensities[write] = sorted_int[read];
        ++write;
      }
    }
    tofs.resize(write);
    intensities.resize(write);
    return write;
  }

  /// Sparse symmetric neighbor intensity sharing.
  /// For each pair within `window` TOF units, add each other's original intensity.
  /// Reads from original intensities, writes to output (no cascading).
  static void smoothTofSpectrum(const std::vector<uint32_t>& tofs,
                                std::vector<uint64_t>& intensities,
                                uint32_t window)
  {
    if (window == 0 || tofs.size() <= 1) return;
    std::vector<uint64_t> smoothed(intensities); // copy original
    for (size_t i = 0; i < tofs.size(); ++i)
    {
      for (size_t j = i + 1; j < tofs.size(); ++j)
      {
        if (tofs[j] - tofs[i] <= window)
        {
          smoothed[i] += intensities[j];
          smoothed[j] += intensities[i];
        }
        else
        {
          break; // sorted, no more neighbors
        }
      }
    }
    intensities = std::move(smoothed);
  }

  /// Sparse local maximum apex picking.
  /// For each pair within `window` TOF units, suppress the lower-intensity one.
  /// Returns a boolean mask (true = keep).
  static std::vector<bool> findSparseLocalMaxima(const std::vector<uint32_t>& tofs,
                                                  const std::vector<uint64_t>& intensities,
                                                  uint32_t window)
  {
    std::vector<bool> keep(tofs.size(), true);
    if (window == 0) return keep;
    for (size_t i = 0; i < tofs.size(); ++i)
    {
      for (size_t j = i + 1; j < tofs.size(); ++j)
      {
        if (tofs[j] - tofs[i] <= window)
        {
          if (intensities[i] < intensities[j])
            keep[i] = false;
          else
            keep[j] = false;
        }
        else
        {
          break;
        }
      }
    }
    return keep;
  }

  // --------------------------------------------------------------------------

  // IM-aware peak for centroiding. Adapted from Sage's ImsPeak.
  struct ImsPeak
  {
    float mz;
    float intensity;
    float im;
  };

  // Per-frame parallel buffers that travel together through the
  // centroiding pipeline: m/z, intensity, IM (1/K0), and TIMS scan index.
  // All four vectors stay in lockstep (same size, same ordering).
  // Downstream centroiders consume them via .data() pointers; the SoA
  // layout is what PASEFHillCentroider / FrameCentroider / PeakPickerIM
  // already expect.
  struct FrameBuffers
  {
    std::vector<double>   mz;
    std::vector<double>   intensity;
    std::vector<double>   im;
    std::vector<uint32_t> scan_ids;

    size_t size() const { return mz.size(); }
    bool   empty() const { return mz.empty(); }

    void clear()
    {
      mz.clear(); intensity.clear(); im.clear(); scan_ids.clear();
    }
    void resize(size_t n)
    {
      mz.resize(n); intensity.resize(n); im.resize(n); scan_ids.resize(n);
    }
    void reserve(size_t n)
    {
      mz.reserve(n); intensity.reserve(n); im.reserve(n); scan_ids.reserve(n);
    }

    // Compact in place to the entries where keep(i) returns true. The
    // ordering of survivors is preserved.
    template <class Keep>
    void filter(Keep keep)
    {
      const size_t n = mz.size();
      size_t w = 0;
      for (size_t r = 0; r < n; ++r)
      {
        if (!keep(r)) continue;
        if (w != r)
        {
          mz[w]        = mz[r];
          intensity[w] = intensity[r];
          im[w]        = im[r];
          scan_ids[w]  = scan_ids[r];
        }
        ++w;
      }
      resize(w);
    }
  };

  // Sort the four buffers in lockstep by m/z ascending. The downstream
  // centroiders all re-sort internally if they need a specific order,
  // so the only caller that strictly needs a sort is the isotopic
  // prefilter (which binary-searches m/z windows).
  static void sortFrameBuffersByMz(FrameBuffers& buf)
  {
    const size_t n = buf.size();
    if (n <= 1) return;
    std::vector<size_t> idx(n);
    std::iota(idx.begin(), idx.end(), size_t(0));
    std::sort(idx.begin(), idx.end(),
              [&buf](size_t a, size_t b) { return buf.mz[a] < buf.mz[b]; });
    // Build a permuted copy then move back; cheaper than in-place
    // permutation for typical PASEF frame sizes (10^4–10^5 peaks).
    FrameBuffers tmp;
    tmp.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
      tmp.mz[i]        = buf.mz[idx[i]];
      tmp.intensity[i] = buf.intensity[idx[i]];
      tmp.im[i]        = buf.im[idx[i]];
      tmp.scan_ids[i]  = buf.scan_ids[idx[i]];
    }
    buf = std::move(tmp);
  }

  // Isotopic-partner prefilter. Drops peaks lacking at least one
  // isotopic partner at m/z ± C13C12_MASSDIFF_U/q (for q in {1..5})
  // within ±tol_ppm (mass-relative) AND |Δscan_id| <= 1. m/z lookup uses
  // binary search, so the buffer is sorted by m/z first (in lockstep).
  // Returns the number of dropped peaks. Output buffers are m/z-sorted.
  //
  // Charges 1..5 cover essentially all peptide MS1 / DIA-MS2 precursors;
  // tol_ppm is intentionally broad (default 50 ppm from the Config) to
  // tolerate per-scan m/z calibration jitter at the prefilter stage.
  // No intensity-ratio check is applied — pure existence test.
  static size_t isotopicPrefilter(FrameBuffers& buf, double tol_ppm)
  {
    if (buf.empty()) return 0;

    sortFrameBuffersByMz(buf);

    const size_t n = buf.size();
    std::vector<bool> keep(n, false);
    static constexpr std::array<int, 5> CHARGES = {1, 2, 3, 4, 5};

    for (size_t i = 0; i < n; ++i)
    {
      if (keep[i]) continue;  // already validated as partner of another
      bool found = false;
      for (int q : CHARGES)
      {
        const double dmz = Constants::C13C12_MASSDIFF_U / static_cast<double>(q);
        for (double sign : {-1.0, 1.0})
        {
          const double target = buf.mz[i] + sign * dmz;
          // Convert ppm tolerance to absolute Da at the target m/z.
          const double tol_da = target * tol_ppm * 1e-6;
          auto it_lo = std::lower_bound(buf.mz.begin(), buf.mz.end(), target - tol_da);
          auto it_hi = std::upper_bound(buf.mz.begin(), buf.mz.end(), target + tol_da);
          for (auto it = it_lo; it != it_hi; ++it)
          {
            const size_t j = static_cast<size_t>(it - buf.mz.begin());
            if (j == i) continue;
            const int64_t ds = static_cast<int64_t>(buf.scan_ids[j])
                             - static_cast<int64_t>(buf.scan_ids[i]);
            if (std::abs(ds) <= 1)
            {
              keep[i] = true;
              keep[j] = true;   // mutual confirmation; saves later work
              found = true;
              break;
            }
          }
          if (found) break;
        }
        if (found) break;
      }
    }

    const size_t before = n;
    buf.filter([&keep](size_t i) { return keep[i]; });
    return before - buf.size();
  }

  // Default cap on the number of centroided peaks per frame. Overridable via
  // BrukerTimsFile::Config::ms1_centroid_max_peaks. For per-frame MS1 (~200k
  // peaks) the 10k default covers the top-intensity peaks well; for aggregated
  // MS1 where the frame-stacked grid can exceed 600k peaks, users should raise
  // this (100k is the new CLI default).
  static constexpr size_t DEFAULT_CENTROID_MAX_PEAKS = 10000;

  // Reusable buffer for IM-dimension centroiding of a single frame.
  // Translated from Sage's PeakBuffer (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
  // Buffers persist across frames within a single load() call for memory reuse.
  struct FrameCentroider
  {
    std::vector<ImsPeak> peaks;
    std::vector<size_t> order;       // indices sorted by descending intensity
    std::vector<ImsPeak> agg_buff;   // centroided output

    // Cap-hit tracking across frames within one load() call. Emitted as a
    // single summary line at the end of the load via reportCapSummary().
    size_t cap_hits_ = 0;
    size_t total_calls_ = 0;
    float  max_dropped_intensity_ = 0.0f;

    void clear()
    {
      peaks.clear();
      order.clear();
      agg_buff.clear();
    }

    // Load already-converted frame data into the buffer.
    // mz_values: converted m/z (double, from save_to_buffs)
    // intensities: corrected uint32_t from save_to_buffs (cast to float internally)
    // im_values: converted per-peak inv_ion_mobility (double, from save_to_buffs)
    void loadFrame(const double* mz_values, const uint32_t* intensities,
                   const double* im_values, uint32_t count)
    {
      clear();
      peaks.reserve(count);
      for (uint32_t i = 0; i < count; ++i)
      {
        peaks.push_back({static_cast<float>(mz_values[i]),
                         static_cast<float>(intensities[i]),
                         static_cast<float>(im_values[i])});
      }

      // Sort by m/z for binary-search neighbor finding
      std::sort(peaks.begin(), peaks.end(),
                [](const ImsPeak& a, const ImsPeak& b) { return a.mz < b.mz; });

      // Build intensity-descending index
      order.resize(count);
      std::iota(order.begin(), order.end(), size_t(0));
      std::sort(order.begin(), order.end(),
                [this](size_t a, size_t b)
                {
                  return peaks[b].intensity < peaks[a].intensity; // descending
                });
    }

    // Overload: accepts double* intensities (from FrameAggregator::finalize() output,
    // where summed intensities across aggregated frames can exceed uint32_t range).
    void loadFrame(const double* mz_values, const double* intensities,
                   const double* im_values, uint32_t count)
    {
      clear();
      peaks.reserve(count);
      for (uint32_t i = 0; i < count; ++i)
      {
        peaks.push_back({static_cast<float>(mz_values[i]),
                         static_cast<float>(intensities[i]),
                         static_cast<float>(im_values[i])});
      }

      std::sort(peaks.begin(), peaks.end(),
                [](const ImsPeak& a, const ImsPeak& b) { return a.mz < b.mz; });

      order.resize(count);
      std::iota(order.begin(), order.end(), size_t(0));
      std::sort(order.begin(), order.end(),
                [this](size_t a, size_t b)
                {
                  return peaks[b].intensity < peaks[a].intensity;
                });
    }

    // Centroid the loaded frame by collapsing the IM dimension.
    // Iterates peaks in descending intensity order. For each apex peak, aggregates
    // all unconsumed neighbors within m/z (ppm) and IM (percent) tolerances.
    // Output is capped at max_peaks entries (dropped peaks logged at WARN level
    // when high-intensity, at DEBUG otherwise).
    void centroid(float mz_ppm, float im_pct,
                  std::vector<double>& out_mz,
                  std::vector<double>& out_intensity,
                  std::vector<float>& out_im,
                  size_t max_peaks = DEFAULT_CENTROID_MAX_PEAKS)
    {
      agg_buff.clear();
      agg_buff.reserve(std::min(peaks.size(), max_peaks));
      ++total_calls_;

      const float utol = mz_ppm / 1e6f;
      const float im_tol_frac = im_pct / 100.0f;
      size_t global_consumed = 0;

      for (size_t idx : order)
      {
        if (peaks[idx].intensity <= 0.0f) continue;  // already consumed

        if (agg_buff.size() >= max_peaks)
        {
          // Only count as a "real" cap hit if the next dropped peak is above
          // the noise floor (~200 counts on timsTOF). Dropping noise-floor
          // peaks is the expected behavior of the cap and not worth reporting.
          if (peaks[idx].intensity > 200.0f)
          {
            ++cap_hits_;
            if (peaks[idx].intensity > max_dropped_intensity_)
              max_dropped_intensity_ = peaks[idx].intensity;
          }
          break;
        }

        const float mz = peaks[idx].mz;
        const float im = peaks[idx].im;
        const float da_tol = mz * utol;
        const float left_mz = mz - da_tol;
        const float right_mz = mz + da_tol;

        // Binary search for m/z neighbor range
        auto it_start = std::lower_bound(peaks.begin(), peaks.end(), left_mz,
                          [](const ImsPeak& p, float val) { return p.mz < val; });
        auto it_end = std::upper_bound(peaks.begin(), peaks.end(), right_mz,
                          [](float val, const ImsPeak& p) { return val < p.mz; });

        const float abs_im_tol = im * im_tol_frac;
        const float left_im = im - abs_im_tol;
        const float right_im = im + abs_im_tol;

        float curr_intensity = 0.0f;
        size_t num_consumed = 0;

        for (auto it = it_start; it != it_end; ++it)
        {
          if (it->intensity <= 0.0f) continue;  // already consumed by earlier apex
          if (it->im >= left_im && it->im <= right_im)
          {
            curr_intensity += it->intensity;
            it->intensity = -1.0f;  // mark consumed
            ++num_consumed;
          }
        }

        agg_buff.push_back({mz, curr_intensity, im});
        global_consumed += num_consumed;

        if (global_consumed == peaks.size()) break;  // all peaks assigned
      }

      // Sort output by m/z
      std::sort(agg_buff.begin(), agg_buff.end(),
                [](const ImsPeak& a, const ImsPeak& b) { return a.mz < b.mz; });

      // Write to output vectors
      out_mz.resize(agg_buff.size());
      out_intensity.resize(agg_buff.size());
      out_im.resize(agg_buff.size());
      for (size_t i = 0; i < agg_buff.size(); ++i)
      {
        out_mz[i] = static_cast<double>(agg_buff[i].mz);
        out_intensity[i] = static_cast<double>(agg_buff[i].intensity);
        out_im[i] = agg_buff[i].im;
      }
    }

    // Emit one summary WARN if the cap fired on any frames during this load()
    // call. Called by each top-level load path (loadDIA_, loadDDA_,
    // loadDIAStreaming, loadFrames_) after all centroiding is done.
    void reportCapSummary(size_t max_peaks)
    {
      if (cap_hits_ == 0) return;
      OPENMS_LOG_WARN << "Warning: MS1 centroiding hit ms1_centroid_max_peaks (="
                      << max_peaks << ") on " << cap_hits_ << " of "
                      << total_calls_ << " spectra; highest dropped-peak "
                      << "intensity was " << max_dropped_intensity_
                      << ". Consider raising -bruker:ms1_centroid_max_peaks "
                      << "to retain more peaks." << std::endl;
      // Reset so successive load() calls on the same instance don't carry
      // state forward.
      cap_hits_ = 0;
      total_calls_ = 0;
      max_dropped_intensity_ = 0.0f;
    }
  };

  // Resolve the effective MS1 centroiding algorithm, taking the legacy
  // (mz_ppm, im_pct) field combination as an implicit Greedy2D selection
  // when ms1_centroid_algo == Off. @p emit_warnings controls whether
  // partial-config warnings are logged; pass true at most once per
  // top-level load() call to avoid log spam.
  static BrukerTimsFile::Config::CentroidAlgo
  effectiveMS1Algo(const BrukerTimsFile::Config& config, bool emit_warnings)
  {
    using CA = BrukerTimsFile::Config::CentroidAlgo;
    if (config.ms1_centroid_algo != CA::OFF)
    {
      if (config.ms1_centroid_mz_ppm <= 0.0f)
      {
        if (emit_warnings)
          OPENMS_LOG_WARN << "Warning: ms1_centroid_algo selected but "
                          << "ms1_centroid_mz_ppm <= 0; MS1 centroiding disabled."
                          << std::endl;
        return CA::OFF;
      }
      if (config.ms1_centroid_algo == CA::GREEDY2D &&
          config.ms1_centroid_im_pct <= 0.0f)
      {
        if (emit_warnings)
          OPENMS_LOG_WARN << "Warning: Greedy2D MS1 centroiding requires "
                          << "ms1_centroid_im_pct > 0; centroiding disabled."
                          << std::endl;
        return CA::OFF;
      }
      return config.ms1_centroid_algo;
    }
    // Back-compat: legacy auto-enable when both mz/im fields are set > 0.
    const bool has_mz = config.ms1_centroid_mz_ppm > 0.0f;
    const bool has_im = config.ms1_centroid_im_pct > 0.0f;
    if (has_mz != has_im)
    {
      if (emit_warnings)
        OPENMS_LOG_WARN << "Warning: ms1_centroid_mz_ppm and "
                        << "ms1_centroid_im_pct must both be > 0 to enable "
                        << "MS1 centroiding. Centroiding disabled."
                        << std::endl;
      return CA::OFF;
    }
    return (has_mz && has_im) ? CA::GREEDY2D : CA::OFF;
  }

  // Resolve the effective DIA-MS2 centroiding algorithm. The legacy
  // boolean dia_ms2_centroid is treated as an implicit Greedy2D selection
  // (Gaussian smoothing + local maxima in finalizeCentroided). HillBased
  // requires ms2_centroid_mz_ppm > 0; warns and falls back to Off otherwise.
  static BrukerTimsFile::Config::CentroidAlgo
  effectiveMS2Algo(const BrukerTimsFile::Config& config, bool emit_warnings)
  {
    using CA = BrukerTimsFile::Config::CentroidAlgo;
    if (config.ms2_centroid_algo != CA::OFF)
    {
      if (config.ms2_centroid_algo == CA::HILL_BASED &&
          config.ms2_centroid_mz_ppm <= 0.0f)
      {
        if (emit_warnings)
          OPENMS_LOG_WARN << "Warning: HillBased DIA-MS2 centroiding requires "
                          << "ms2_centroid_mz_ppm > 0; centroiding disabled."
                          << std::endl;
        return CA::OFF;
      }
      return config.ms2_centroid_algo;
    }
    // Back-compat: dia_ms2_centroid=true implies Greedy2D.
    return config.dia_ms2_centroid ? CA::GREEDY2D : CA::OFF;
  }

  // Convenience wrapper preserving the old call-site semantics.
  static bool isCentroidingEnabled(const BrukerTimsFile::Config& config)
  {
    return effectiveMS1Algo(config, /*emit_warnings=*/false) !=
           BrukerTimsFile::Config::CentroidAlgo::OFF;
  }

  // Run the selected MS1 centroiding algorithm on parallel (mz, intensity, im)
  // arrays. Output arrays are sorted by m/z ascending and parallel to each
  // other. The Greedy2D path uses the existing FrameCentroider; the HillBased
  // path uses PASEFHillCentroider. The caller must already have resolved the
  // algorithm via effectiveMS1Algo() and confirmed it is not Off.
  // Bounding box arrays for hill-based centroiding output. When non-null,
  // populated parallel to the (mz, intensity, im) outputs and exposed as
  // extra FloatDataArrays via expose_hill_bounds. Greedy2D is not a hill
  // algorithm, so these are only populated for the HillBased branch.
  struct HillBoundsOut
  {
    std::vector<float>* im_lower;
    std::vector<float>* im_upper;
    std::vector<float>* mz_lower;
    std::vector<float>* mz_upper;
  };

  static void fillHillBounds(
      const std::vector<Internal::PASEFHillCentroider::Centroid>& centroids,
      const HillBoundsOut& b)
  {
    if (!b.im_lower) return;
    b.im_lower->resize(centroids.size());
    b.im_upper->resize(centroids.size());
    b.mz_lower->resize(centroids.size());
    b.mz_upper->resize(centroids.size());
    for (std::size_t i = 0; i < centroids.size(); ++i)
    {
      (*b.im_lower)[i] = static_cast<float>(centroids[i].im_lower);
      (*b.im_upper)[i] = static_cast<float>(centroids[i].im_upper);
      (*b.mz_lower)[i] = static_cast<float>(centroids[i].mz_lower);
      (*b.mz_upper)[i] = static_cast<float>(centroids[i].mz_upper);
    }
  }

  static void runMS1Centroider(
      BrukerTimsFile::Config::CentroidAlgo algo,
      const BrukerTimsFile::Config& config,
      FrameCentroider& greedy,
      const double* in_mz,
      const double* in_intensity,
      const double* in_im,
      uint32_t n,
      std::vector<double>& out_mz,
      std::vector<double>& out_intensity,
      std::vector<float>& out_im,
      const std::uint32_t* in_scan_ids = nullptr,
      const HillBoundsOut& bounds = HillBoundsOut{nullptr,nullptr,nullptr,nullptr})
  {
    using CA = BrukerTimsFile::Config::CentroidAlgo;
    if (algo == CA::HILL_BASED)
    {
      Internal::PASEFHillCentroider::Params p;
      p.mz_tol_ppm        = static_cast<double>(config.ms1_centroid_mz_ppm);
      p.valley_factor     = config.centroid_valley_factor;
      p.min_hill_length   = config.ms1_centroid_min_hill_length;
      p.max_scan_gap      = config.centroid_max_scan_gap;
      auto centroids = Internal::PASEFHillCentroider::centroidFrame(
          in_mz, in_intensity, in_im, static_cast<std::size_t>(n), p,
          in_scan_ids);
      out_mz.resize(centroids.size());
      out_intensity.resize(centroids.size());
      out_im.resize(centroids.size());
      for (std::size_t i = 0; i < centroids.size(); ++i)
      {
        out_mz[i]        = centroids[i].mz;
        out_intensity[i] = centroids[i].intensity;
        out_im[i]        = static_cast<float>(centroids[i].im_apex);
      }
      fillHillBounds(centroids, bounds);
      return;
    }
    // Greedy2D path: matches the legacy behavior bit-for-bit.
    greedy.loadFrame(in_mz, in_intensity, in_im, n);
    greedy.centroid(config.ms1_centroid_mz_ppm, config.ms1_centroid_im_pct,
                    out_mz, out_intensity, out_im,
                    static_cast<std::size_t>(config.ms1_centroid_max_peaks));
  }

  // Run hill-based centroiding on a DIA-MS2 window's aggregated peaks. The
  // Greedy2D MS2 path is handled inline via FrameAggregator::finalizeCentroided
  // (Gaussian smoothing + local maxima); only the HillBased path needs this
  // helper. Returns true if centroids were produced; false if the input was
  // empty or the helper rejected all hills.
  static bool runMS2HillCentroider(
      const BrukerTimsFile::Config& config,
      const double* in_mz,
      const double* in_intensity,
      const double* in_im,
      std::size_t n,
      std::vector<double>& out_mz,
      std::vector<double>& out_intensity,
      std::vector<float>& out_im,
      const std::uint32_t* in_scan_ids = nullptr,
      const HillBoundsOut& bounds = HillBoundsOut{nullptr,nullptr,nullptr,nullptr})
  {
    Internal::PASEFHillCentroider::Params p;
    p.mz_tol_ppm        = static_cast<double>(config.ms2_centroid_mz_ppm);
    p.valley_factor     = config.centroid_valley_factor;
    p.min_hill_length   = config.ms2_centroid_min_hill_length;
    p.max_scan_gap      = config.centroid_max_scan_gap;
    auto centroids = Internal::PASEFHillCentroider::centroidFrame(
        in_mz, in_intensity, in_im, n, p, in_scan_ids);
    out_mz.resize(centroids.size());
    out_intensity.resize(centroids.size());
    out_im.resize(centroids.size());
    for (std::size_t i = 0; i < centroids.size(); ++i)
    {
      out_mz[i]        = centroids[i].mz;
      out_intensity[i] = centroids[i].intensity;
      out_im[i]        = static_cast<float>(centroids[i].im_apex);
    }
    fillHillBounds(centroids, bounds);
    return !centroids.empty();
  }

  // Helper to attach the four hill-bounds FloatDataArrays to a spectrum.
  // Centralized so emission sites remain uniform.
  static void attachHillBoundsArrays(
      MSSpectrum& spec,
      const std::vector<float>& im_lower,
      const std::vector<float>& im_upper,
      const std::vector<float>& mz_lower,
      const std::vector<float>& mz_upper)
  {
    auto make = [](const std::vector<float>& src, const String& name) {
      DataArrays::FloatDataArray a;
      a.setName(name);
      a.assign(src.begin(), src.end());
      return a;
    };
    spec.getFloatDataArrays().push_back(make(im_lower, "im lower bound"));
    spec.getFloatDataArrays().push_back(make(im_upper, "im upper bound"));
    spec.getFloatDataArrays().push_back(make(mz_lower, "m/z lower bound"));
    spec.getFloatDataArrays().push_back(make(mz_upper, "m/z upper bound"));
  }

  // Predicate for the Config::frame_id_min/frame_id_max filter. A default
  // Config (0, UINT32_MAX) trivially matches every frame, so callers can
  // gate every frame iteration on this without branching for the no-filter
  // case.
  static inline bool frameInRange(uint32_t fid, const BrukerTimsFile::Config& config)
  {
    return fid >= config.frame_id_min && fid <= config.frame_id_max;
  }

  // One-shot validation of the frame_id and RT ranges against the file's
  // actual frame bounds + RT resolution. Called from each top-level public
  // entry point. Returns an effective Config in which frame_id_min/max have
  // been narrowed to the intersection of:
  //   - the user-supplied frame_id range
  //   - the frame_id range derived from the user-supplied RT range
  //     (SELECT MIN/MAX Id FROM Frames WHERE Time IN [rt_min, rt_max])
  // so downstream code only needs to honor frame_id_min/frame_id_max.
  // Emits warnings for inverted or fully-out-of-bounds ranges.
  static BrukerTimsFile::Config resolveEffectiveConfig(
      SQLite::Database& db,
      const BrukerTimsFile::Config& config,
      uint32_t file_min, uint32_t file_max)
  {
    BrukerTimsFile::Config eff = config;

    // --- 1) Validate frame_id range ---
    const bool fid_set = (config.frame_id_min != 0) ||
                         (config.frame_id_max != std::numeric_limits<uint32_t>::max());
    if (fid_set)
    {
      if (config.frame_id_min > config.frame_id_max)
      {
        OPENMS_LOG_WARN << "Warning: BrukerTimsFile::Config::frame_id_min ("
                        << config.frame_id_min << ") > frame_id_max ("
                        << config.frame_id_max
                        << "). No frames will be loaded." << std::endl;
      }
      else if (config.frame_id_max < file_min || config.frame_id_min > file_max)
      {
        OPENMS_LOG_WARN << "Warning: BrukerTimsFile::Config frame range ["
                        << config.frame_id_min << ", " << config.frame_id_max
                        << "] does not intersect file frame range ["
                        << file_min << ", " << file_max
                        << "]. No frames will be loaded." << std::endl;
      }
    }

    // --- 2) Validate and resolve RT range ---
    const bool rt_set = (config.rt_min_sec > 0.0) ||
                        std::isfinite(config.rt_max_sec);
    if (!rt_set) return eff;

    if (config.rt_min_sec > config.rt_max_sec)
    {
      OPENMS_LOG_WARN << "Warning: BrukerTimsFile::Config::rt_min_sec ("
                      << config.rt_min_sec << ") > rt_max_sec ("
                      << config.rt_max_sec
                      << "). No frames will be loaded." << std::endl;
      // Force the effective frame_id range empty so downstream loops produce 0 spectra.
      eff.frame_id_min = 1;
      eff.frame_id_max = 0;
      return eff;
    }

    try
    {
      // TIMS frame IDs are monotonic in acquisition time, so the matching
      // (MIN Id, MAX Id) is a contiguous range. Use an inclusive
      // Time filter; if the user's max is +inf, the comparison still works
      // (SQLite stores doubles natively).
      SQLite::Statement q(db,
        "SELECT MIN(Id), MAX(Id) FROM Frames WHERE Time >= ? AND Time <= ?");
      q.bind(1, config.rt_min_sec);
      q.bind(2, std::isfinite(config.rt_max_sec)
                  ? config.rt_max_sec
                  : std::numeric_limits<double>::max());

      if (q.executeStep() && !q.getColumn(0).isNull())
      {
        const uint32_t rt_fid_min = static_cast<uint32_t>(q.getColumn(0).getInt());
        const uint32_t rt_fid_max = static_cast<uint32_t>(q.getColumn(1).getInt());
        // Save pre-intersection values to detect disjoint ranges
        const uint32_t pre_min = eff.frame_id_min;
        const uint32_t pre_max = eff.frame_id_max;
        eff.frame_id_min = std::max(eff.frame_id_min, rt_fid_min);
        eff.frame_id_max = std::min(eff.frame_id_max, rt_fid_max);
        // Check if the intersection resulted in an empty range
        if (pre_min <= pre_max && eff.frame_id_min > eff.frame_id_max)
        {
          OPENMS_LOG_WARN << "Warning: BrukerTimsFile::Config frame_id range ["
                          << pre_min << ", " << pre_max
                          << "] and RT range [" << config.rt_min_sec << ", "
                          << config.rt_max_sec << "] s do not intersect. "
                          << "No frames will be loaded." << std::endl;
        }
      }
      else
      {
        // No frame's Time falls in the user's RT range — produce empty output.
        OPENMS_LOG_WARN << "Warning: BrukerTimsFile::Config RT range ["
                        << config.rt_min_sec << ", " << config.rt_max_sec
                        << "] s does not intersect any frame's acquisition "
                        << "time. No frames will be loaded." << std::endl;
        eff.frame_id_min = 1;
        eff.frame_id_max = 0;
      }
    }
    catch (const SQLite::Exception& ex)
    {
      OPENMS_LOG_WARN << "Warning: failed to resolve RT range via SQL: "
                      << ex.what() << ". RT filter will be ignored." << std::endl;
    }

    return eff;
  }

  // Emit one-shot partial-config warnings for the MS1 aggregation knobs.
  // Called once per top-level load path, after ms1_frame_ids is populated.
  static void warnPartialMS1AggregationConfig(const BrukerTimsFile::Config& config,
                                               size_t ms1_frame_count)
  {
    if (config.ms1_min_support > 0 && config.ms1_n_neighbors == 0)
    {
      OPENMS_LOG_WARN << "Warning: ms1_min_support (=" << config.ms1_min_support
                      << ") ignored: ms1_n_neighbors=0 disables aggregation." << std::endl;
    }
    if (config.ms1_max_rt_distance_sec > 0.0 && config.ms1_n_neighbors == 0)
    {
      OPENMS_LOG_WARN << "Warning: ms1_max_rt_distance_sec (="
                      << config.ms1_max_rt_distance_sec
                      << ") ignored: ms1_n_neighbors=0 disables aggregation." << std::endl;
    }
    if (config.ms1_n_neighbors > 0 && ms1_frame_count > 0 &&
        static_cast<size_t>(2 * config.ms1_n_neighbors + 1) > ms1_frame_count)
    {
      OPENMS_LOG_WARN << "Warning: ms1_n_neighbors (=" << config.ms1_n_neighbors
                      << ") exceeds half the run's MS1 frame count (="
                      << ms1_frame_count << "); effective window will be clamped "
                      << "at run edges." << std::endl;
    }
    if (!config.load_ms1 && config.ms1_n_neighbors > 0)
    {
      OPENMS_LOG_WARN << "Warning: ms1_n_neighbors (=" << config.ms1_n_neighbors
                      << ") ignored: load_ms1=false disables MS1 loading entirely."
                      << std::endl;
    }

    // One-shot validation of the centroiding config: emits warnings for
    // partial/inconsistent inputs but does not change behavior. Per-spectrum
    // dispatch later calls effectiveMS{1,2}Algo with emit_warnings=false.
    // Both DDA-MS2 and DIA-MS2 HillBased work without any frame aggregation
    // (the hill centroider operates within a single PASEF frame across its
    // own IM scans); dia_ms2_n_neighbors only controls cross-RT summing.
    (void)effectiveMS1Algo(config, /*emit_warnings=*/true);
    (void)effectiveMS2Algo(config, /*emit_warnings=*/true);
  }

  // =====================================================================
  // TdfMetadataReader: private helper functions for SQL-based metadata
  // =====================================================================
  namespace
  {

    // Check if the dataset is DIA by looking for SWATH window tables with data
    bool isDIA(SQLite::Database& db)
    {
      // Try DiaFrameMsMsWindows first (newer format)
      try
      {
        SQLite::Statement q(db, "SELECT COUNT(*) FROM DiaFrameMsMsWindows");
        if (q.executeStep() && q.getColumn(0).getInt() > 0)
          return true;
      }
      catch (const SQLite::Exception&)
      {
        // Table does not exist, try fallback
      }

      // Fallback: DiaFrameMsMsInfo (older format)
      try
      {
        SQLite::Statement q(db, "SELECT COUNT(*) FROM DiaFrameMsMsInfo");
        if (q.executeStep() && q.getColumn(0).getInt() > 0)
          return true;
      }
      catch (const SQLite::Exception&)
      {
        // Table does not exist either
      }

      return false;
    }

    /// Grid-based aggregator for stacking peaks from adjacent frames into a sparse
    /// (mz_bin, scan_id) grid, with optional 3x3 spatial denoising. Used by both
    /// DIA MS2 aggregation and MS1 aggregation (with different bin widths).
    class FrameAggregator
    {
    public:
      explicit FrameAggregator(double mz_bin_width = 0.02)
        : mz_bin_width_(mz_bin_width) {}

      struct OutputPeak
      {
        double mz;        // intensity-weighted mean m/z
        double intensity;  // summed intensity
        uint32_t scan_id;  // native scan index (for IM conversion)
      };

      /// Returns the grid key of the neighbor cell at signed (dm, ds) offset.
      /// Returns std::nullopt when the offset would wrap past 0 on either axis.
      static inline std::optional<uint64_t>
      neighborKey(uint32_t mz_bin, uint32_t scan_id, int dm, int ds)
      {
        if (dm < 0 && mz_bin < static_cast<uint32_t>(-dm)) return std::nullopt;
        if (ds < 0 && scan_id < static_cast<uint32_t>(-ds)) return std::nullopt;
        uint32_t new_mz = mz_bin + static_cast<uint32_t>(dm);
        uint32_t new_scan = scan_id + static_cast<uint32_t>(ds);
        return (static_cast<uint64_t>(new_mz) << 32) | new_scan;
      }

      /// Add a peak to the grid. Call for every peak from every neighbor frame.
      void addPeak(double mz, uint32_t intensity, uint32_t scan_id)
      {
        int64_t mz_bin = static_cast<int64_t>(std::round(mz / mz_bin_width_));
        uint64_t key = (static_cast<uint64_t>(static_cast<uint32_t>(mz_bin)) << 32) | scan_id;

        auto& cell = grid_[key];
        double int_d = static_cast<double>(intensity);
        cell.intensity_sum += int_d;
        cell.mz_weighted_sum += mz * int_d;
      }

      /// Apply spatial denoising and return surviving peaks.
      /// min_support: minimum occupied neighbors in 3x3 grid (center excluded).
      /// If skip_denoise is true (e.g., only 1 frame), return all peaks without filtering.
      std::vector<OutputPeak> finalize(int min_support, bool skip_denoise) const
      {
        std::vector<OutputPeak> result;
        result.reserve(grid_.size());

        for (const auto& [key, cell] : grid_)
        {
          if (!skip_denoise)
          {
            int neighbors = 0;
            uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
            uint32_t mz_bin = static_cast<uint32_t>(key >> 32);

            // Check 3x3 neighborhood (center excluded)
            for (int dm = -1; dm <= 1; ++dm)
            {
              for (int ds = -1; ds <= 1; ++ds)
              {
                if (dm == 0 && ds == 0) continue;
                if (auto nkey = neighborKey(mz_bin, scan_id, dm, ds))
                {
                  if (grid_.count(*nkey)) ++neighbors;
                }
              }
            }

            if (neighbors < min_support) continue;
          }

          uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
          double mz = cell.mz_weighted_sum / cell.intensity_sum;

          result.push_back({mz, cell.intensity_sum, scan_id});
        }

        // Sort by m/z for spectrum output
        std::sort(result.begin(), result.end(),
          [](const OutputPeak& a, const OutputPeak& b) { return a.mz < b.mz; });

        return result;
      }

      /// Apply 2D Gaussian smoothing + local maxima peak picking to the denoised grid.
      /// Returns centroided peaks with sub-bin (m/z, scan_id) precision.
      std::vector<OutputPeak> finalizeCentroided(int min_support, bool skip_denoise) const
      {
        // Step 1: Denoise (same as finalize)
        boost::unordered_flat_map<uint64_t, Cell> denoised;
        denoised.reserve(grid_.size());
        for (const auto& [key, cell] : grid_)
        {
          if (!skip_denoise)
          {
            int neighbors = 0;
            uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
            uint32_t mz_bin = static_cast<uint32_t>(key >> 32);
            for (int dm = -1; dm <= 1; ++dm)
            {
              for (int ds = -1; ds <= 1; ++ds)
              {
                if (dm == 0 && ds == 0) continue;
                if (auto nkey = neighborKey(mz_bin, scan_id, dm, ds))
                {
                  if (grid_.count(*nkey)) ++neighbors;
                }
              }
            }
            if (neighbors < min_support) continue;
          }
          denoised[key] = cell;
        }

        // Step 2: Gaussian smooth (sigma=1 bin, 5x5 support)
        static const double kernel[5][5] = {
          {0.003, 0.013, 0.022, 0.013, 0.003},
          {0.013, 0.059, 0.097, 0.059, 0.013},
          {0.022, 0.097, 0.159, 0.097, 0.022},
          {0.013, 0.059, 0.097, 0.059, 0.013},
          {0.003, 0.013, 0.022, 0.013, 0.003}
        };

        boost::unordered_flat_map<uint64_t, double> smoothed;
        smoothed.reserve(denoised.size());
        for (const auto& [key, cell] : denoised)
        {
          uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
          uint32_t mz_bin = static_cast<uint32_t>(key >> 32);
          double weighted_sum = 0.0;
          for (int dm = -2; dm <= 2; ++dm)
          {
            for (int ds = -2; ds <= 2; ++ds)
            {
              if (auto nkey = neighborKey(mz_bin, scan_id, dm, ds))
              {
                auto it = denoised.find(*nkey);
                if (it != denoised.end())
                {
                  weighted_sum += it->second.intensity_sum * kernel[dm + 2][ds + 2];
                }
              }
            }
          }
          smoothed[key] = weighted_sum;
        }

        // Step 3: Find local maxima in smoothed grid
        std::vector<uint64_t> maxima;
        for (const auto& [key, val] : smoothed)
        {
          uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
          uint32_t mz_bin = static_cast<uint32_t>(key >> 32);
          bool is_max = true;
          for (int dm = -1; dm <= 1 && is_max; ++dm)
          {
            for (int ds = -1; ds <= 1 && is_max; ++ds)
            {
              if (dm == 0 && ds == 0) continue;
              if (auto nkey = neighborKey(mz_bin, scan_id, dm, ds))
              {
                auto it = smoothed.find(*nkey);
                if (it != smoothed.end() && it->second > val) is_max = false;
              }
            }
          }
          if (is_max) maxima.push_back(key);
        }

        // Step 4: Centroid each maximum from original denoised cells within ±2 radius
        std::vector<OutputPeak> result;
        result.reserve(maxima.size());
        for (uint64_t max_key : maxima)
        {
          uint32_t center_scan = static_cast<uint32_t>(max_key & 0xFFFFFFFF);
          uint32_t center_mz_bin = static_cast<uint32_t>(max_key >> 32);

          double total_intensity = 0.0;
          double mz_weighted = 0.0;
          double scan_weighted = 0.0;

          for (int dm = -2; dm <= 2; ++dm)
          {
            for (int ds = -2; ds <= 2; ++ds)
            {
              if (auto nkey = neighborKey(center_mz_bin, center_scan, dm, ds))
              {
                auto it = denoised.find(*nkey);
                if (it != denoised.end())
                {
                  double int_val = it->second.intensity_sum;
                  double mz_val = it->second.mz_weighted_sum / it->second.intensity_sum;
                  total_intensity += int_val;
                  mz_weighted += mz_val * int_val;
                  scan_weighted += static_cast<double>(static_cast<uint32_t>(*nkey & 0xFFFFFFFF)) * int_val;
                }
              }
            }
          }

          double centroid_mz = mz_weighted / total_intensity;
          uint32_t centroid_scan = static_cast<uint32_t>(std::round(scan_weighted / total_intensity));
          result.push_back({centroid_mz, total_intensity, centroid_scan});
        }

        std::sort(result.begin(), result.end(),
          [](const OutputPeak& a, const OutputPeak& b) { return a.mz < b.mz; });
        return result;
      }

      /// Clear grid for reuse
      void clear() { grid_.clear(); }

      /// Pre-reserve bucket space for N cells. Call once before a hot loop
      /// to avoid rehash cascades. Safe to over-allocate — clear() preserves
      /// the bucket array across center frames.
      void reserve(size_t n) { grid_.reserve(n); }

    private:
      double mz_bin_width_;
      struct Cell
      {
        double intensity_sum = 0.0;
        double mz_weighted_sum = 0.0;
      };
      boost::unordered_flat_map<uint64_t, Cell> grid_;
    };

    // Fill a FrameBuffers struct from FrameAggregator output peaks. Converts
    // each peak's scan_id to inverse ion mobility (1/K0) via the handle's
    // converter using @p center_frame_id as the calibration reference. Used
    // by the sites that consume aggregator output: loadAggregatedMS1Spectrum
    // and the DIA-MS2 aggregated hill path.
    static void aggregatorPeaksToBuffers(
        const std::vector<FrameAggregator::OutputPeak>& peaks,
        TimsDataHandle& handle,
        uint32_t center_frame_id,
        FrameBuffers& out)
    {
      const size_t n = peaks.size();
      out.resize(n);
      for (size_t i = 0; i < n; ++i)
      {
        out.mz[i]        = peaks[i].mz;
        out.intensity[i] = peaks[i].intensity;
        double im_val = 0.0;
        handle.scan2inv_ion_mobility_converter->convert(
          center_frame_id, &im_val, &peaks[i].scan_id, 1);
        out.im[i]        = im_val;
        out.scan_ids[i]  = peaks[i].scan_id;
      }
    }

    // SWATH window descriptor
    struct DIAWindow
    {
      int window_group;
      double mz_center;
      double mz_width;
      double im_lower;  // 1/K0 lower bound
      double im_upper;  // 1/K0 upper bound
      uint32_t scan_begin;  // raw scan index lower bound (for grid-based IM filtering)
      uint32_t scan_end;    // raw scan index upper bound
    };

    // Read DIA SWATH windows from SQL, converting scan bounds to IM
    std::vector<DIAWindow> readDIAWindows(SQLite::Database& db,
                                          Scan2InvIonMobilityConverter& im_converter)
    {
      std::vector<DIAWindow> windows;

      // Try newer DiaFrameMsMsWindows table first
      bool use_new_table = false;
      try
      {
        SQLite::Statement check(db, "SELECT COUNT(*) FROM DiaFrameMsMsWindows");
        if (check.executeStep() && check.getColumn(0).getInt() > 0)
          use_new_table = true;
      }
      catch (const SQLite::Exception&) {}

      if (use_new_table)
      {
        SQLite::Statement q(db,
          "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth "
          "FROM DiaFrameMsMsWindows ORDER BY WindowGroup, ScanNumBegin");

        while (q.executeStep())
        {
          DIAWindow w;
          w.window_group = q.getColumn(0).getInt();
          uint32_t scan_begin = static_cast<uint32_t>(q.getColumn(1).getInt());
          uint32_t scan_end   = static_cast<uint32_t>(q.getColumn(2).getInt());
          w.scan_begin = scan_begin;
          w.scan_end = scan_end;
          w.mz_center = q.getColumn(3).getDouble();
          w.mz_width  = q.getColumn(4).getDouble();

          // Convert scan bounds to 1/K0 (frame_id=1 is used as reference; linear model is frame-independent)
          double im_begin = 0.0, im_end = 0.0;
          im_converter.convert(1, &im_begin, &scan_begin, 1);
          im_converter.convert(1, &im_end, &scan_end, 1);
          w.im_lower = std::min(im_begin, im_end);
          w.im_upper = std::max(im_begin, im_end);

          windows.push_back(w);
        }
      }
      else
      {
        // Fallback: DiaFrameMsMsInfo (older format, may have per-frame entries)
        // Deduplicate by (IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd)
        SQLite::Statement q(db,
          "SELECT DISTINCT ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth "
          "FROM DiaFrameMsMsInfo ORDER BY IsolationMz, ScanNumBegin");

        int group_id = 0;
        while (q.executeStep())
        {
          DIAWindow w;
          w.window_group = group_id++;
          uint32_t scan_begin = static_cast<uint32_t>(q.getColumn(0).getInt());
          uint32_t scan_end   = static_cast<uint32_t>(q.getColumn(1).getInt());
          w.scan_begin = scan_begin;
          w.scan_end = scan_end;
          w.mz_center = q.getColumn(2).getDouble();
          w.mz_width  = q.getColumn(3).getDouble();

          double im_begin = 0.0, im_end = 0.0;
          im_converter.convert(1, &im_begin, &scan_begin, 1);
          im_converter.convert(1, &im_end, &scan_end, 1);
          w.im_lower = std::min(im_begin, im_end);
          w.im_upper = std::max(im_begin, im_end);

          windows.push_back(w);
        }
      }

      return windows;
    }

    // Build frame_id -> WindowGroup mapping from SQL.
    // Returns map<window_group, vector<frame_id>> sorted by frame_id within each group.
    std::map<int, std::vector<uint32_t>> readFrameToWindowGroupMapping(
      SQLite::Database& db,
      const std::vector<DIAWindow>& windows)
    {
      std::map<int, std::vector<uint32_t>> mapping;

      // Try DiaFrameMsMsInfo with WindowGroup column first (newer format)
      bool has_window_group_column = false;
      try
      {
        SQLite::Statement check(db,
          "SELECT WindowGroup FROM DiaFrameMsMsInfo LIMIT 1");
        if (check.executeStep()) has_window_group_column = true;
      }
      catch (const SQLite::Exception&) {}

      if (has_window_group_column)
      {
        SQLite::Statement q(db,
          "SELECT DISTINCT Frame, WindowGroup FROM DiaFrameMsMsInfo ORDER BY Frame");
        while (q.executeStep())
        {
          uint32_t frame_id = static_cast<uint32_t>(q.getColumn(0).getInt());
          int group = q.getColumn(1).getInt();
          mapping[group].push_back(frame_id);
        }
      }
      else
      {
        // Older format: DiaFrameMsMsInfo has (Frame, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth)
        // but no WindowGroup column. Match each frame's geometry to the synthetic group IDs
        // assigned by readDIAWindows().
        struct WindowKey
        {
          double mz; double width; uint32_t sb; uint32_t se;
          bool operator<(const WindowKey& o) const
          {
            if (mz != o.mz) return mz < o.mz;
            if (width != o.width) return width < o.width;
            if (sb != o.sb) return sb < o.sb;
            return se < o.se;
          }
        };
        std::map<WindowKey, int> key_to_group;
        for (const auto& w : windows)
        {
          key_to_group[{w.mz_center, w.mz_width, w.scan_begin, w.scan_end}] = w.window_group;
        }

        SQLite::Statement q(db,
          "SELECT Frame, IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd "
          "FROM DiaFrameMsMsInfo ORDER BY Frame");
        while (q.executeStep())
        {
          uint32_t frame_id = static_cast<uint32_t>(q.getColumn(0).getInt());
          WindowKey key{
            q.getColumn(1).getDouble(),
            q.getColumn(2).getDouble(),
            static_cast<uint32_t>(q.getColumn(3).getInt()),
            static_cast<uint32_t>(q.getColumn(4).getInt())
          };
          auto it = key_to_group.find(key);
          if (it != key_to_group.end())
          {
            auto& frames = mapping[it->second];
            if (frames.empty() || frames.back() != frame_id)
            {
              frames.push_back(frame_id);
            }
          }
        }
      }

      return mapping;
    }

    // DDA precursor info (joined from Precursors + PasefFrameMsMsInfo)
    struct DDAPrecursorInfo
    {
      int precursor_id;
      double monoisotopic_mz;
      int charge;
      double intensity;      // precursor intensity (may be NaN/NULL)
      uint32_t frame_id;
      uint32_t scan_begin;
      uint32_t scan_end;
      double isolation_mz;
      double isolation_width;
    };

    // Read DDA precursors from SQL
    std::vector<DDAPrecursorInfo> readDDAPrecursors(SQLite::Database& db)
    {
      std::vector<DDAPrecursorInfo> precursors;

      SQLite::Statement q(db,
        "SELECT p.Id, p.MonoisotopicMz, p.Charge, p.Intensity, "
        "       pfm.Frame, pfm.ScanNumBegin, pfm.ScanNumEnd, "
        "       pfm.IsolationMz, pfm.IsolationWidth "
        "FROM Precursors p "
        "JOIN PasefFrameMsMsInfo pfm ON p.Id = pfm.Precursor "
        "ORDER BY p.Id, pfm.Frame, pfm.ScanNumBegin");

      while (q.executeStep())
      {
        DDAPrecursorInfo info;
        info.precursor_id   = q.getColumn(0).getInt();
        info.monoisotopic_mz = q.getColumn(1).getDouble();
        info.charge          = q.getColumn(2).isNull() ? 0 : q.getColumn(2).getInt();
        info.intensity       = q.getColumn(3).isNull() ? std::numeric_limits<double>::quiet_NaN() : q.getColumn(3).getDouble();
        info.frame_id        = static_cast<uint32_t>(q.getColumn(4).getInt());
        info.scan_begin      = static_cast<uint32_t>(q.getColumn(5).getInt());
        info.scan_end        = static_cast<uint32_t>(q.getColumn(6).getInt());
        info.isolation_mz    = q.getColumn(7).getDouble();
        info.isolation_width = q.getColumn(8).getDouble();
        precursors.push_back(info);
      }

      return precursors;
    }

    // Frame count info for progress reporting
    struct FrameCounts
    {
      uint32_t total = 0;
      uint32_t ms1 = 0;
      uint32_t ms2 = 0;
    };

    FrameCounts readFrameCounts(SQLite::Database& db,
                                const BrukerTimsFile::Config& config)
    {
      FrameCounts counts;

      const std::string range_clause =
        " WHERE Id >= " + std::to_string(config.frame_id_min) +
        " AND Id <= " + std::to_string(config.frame_id_max);

      SQLite::Statement q_total(db, "SELECT COUNT(*) FROM Frames" + range_clause);
      if (q_total.executeStep())
        counts.total = static_cast<uint32_t>(q_total.getColumn(0).getInt());

      SQLite::Statement q_ms1(db,
        "SELECT COUNT(*) FROM Frames" + range_clause + " AND MsMsType = 0");
      if (q_ms1.executeStep())
        counts.ms1 = static_cast<uint32_t>(q_ms1.getColumn(0).getInt());

      // MS2 = everything that is not MS1 (MsMsType != 0)
      SQLite::Statement q_ms2(db,
        "SELECT COUNT(*) FROM Frames" + range_clause + " AND MsMsType != 0");
      if (q_ms2.executeStep())
        counts.ms2 = static_cast<uint32_t>(q_ms2.getColumn(0).getInt());

      return counts;
    }

    // Count MS2 spectra that loadDDA_ will emit under the active frame range:
    // one spectrum per distinct precursor that still has at least one
    // PasefFrameMsMsInfo row pointing at an in-range MS2 frame. With the
    // default range [0, UINT32_MAX] this reduces to "every distinct
    // precursor", matching the legacy behavior.
    uint32_t countDDAPrecursors(SQLite::Database& db,
                                const BrukerTimsFile::Config& config)
    {
      try
      {
        SQLite::Statement q(db,
          "SELECT COUNT(DISTINCT Precursor) FROM PasefFrameMsMsInfo"
          " WHERE Frame >= " + std::to_string(config.frame_id_min) +
          " AND Frame <= " + std::to_string(config.frame_id_max));
        if (q.executeStep())
          return static_cast<uint32_t>(q.getColumn(0).getInt());
      }
      catch (const SQLite::Exception&) {}
      return 0;
    }

  } // anonymous namespace

  // Map OpenMS PressureCompensation enum to opentims pressure_compensation_strategy
  static pressure_compensation_strategy mapPressureCompensation(
    BrukerTimsFile::Config::PressureCompensation pc)
  {
    switch (pc)
    {
      case BrukerTimsFile::Config::PressureCompensation::GLOBAL:
        return AnalyisGlobalPressureCompensation;
      case BrukerTimsFile::Config::PressureCompensation::PER_FRAME:
        return PerFramePressureCompensationWithMissingReference;
      default:
        return NoPressureCompensation;
    }
  }

  // =====================================================================
  // Helper: open TimsDataHandle with tiered calibration strategy
  // =====================================================================
  using Config = BrukerTimsFile::Config;

  static std::unique_ptr<TimsDataHandle> openTimsDataHandle(
    const String& path, const Config& config = Config())
  {
    std::string path_string = path;

    // 1. Always create handle with linear converters (safe baseline)
    auto& tof_factory = OpenSourceTof2MzConverterFactory::instance();
    auto& im_factory = OpenSourceScan2ImConverterFactory::instance();

    std::unique_ptr<TimsDataHandle> handle;
    try
    {
      handle = std::make_unique<TimsDataHandle>(
        path_string, NoPressureCompensation, &tof_factory, &im_factory);
    }
    catch (const std::exception& e)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path + " (opentims: " + String(e.what()) + ")");
    }
    catch (...)
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        path + " (opentims: unknown error)");
    }

    using Strategy = Config::TimsCalibrationStrategy;
    auto strategy = config.tims_calibration_strategy;

    // Warn if pressure compensation requested without SDK strategy
    if (config.pressure_compensation != Config::PressureCompensation::NONE
        && strategy != Strategy::BRUKER_SDK
        && strategy != Strategy::AUTO)
    {
      OPENMS_LOG_WARN << "Pressure compensation requires BRUKER_SDK strategy, ignoring"
                      << std::endl;
    }

    // 2. Try Bruker SDK (AUTO or BRUKER_SDK)
    if (strategy == Strategy::AUTO || strategy == Strategy::BRUKER_SDK)
    {
      std::string sdk_path = config.bruker_sdk_path;
      if (sdk_path.empty())
      {
        const char* env = std::getenv("OPENMS_BRUKER_SDK_PATH");
        if (env) sdk_path = env;
      }

      if (!sdk_path.empty())
      {
        try
        {
          auto pcs = mapPressureCompensation(config.pressure_compensation);
          handle->scan2inv_ion_mobility_converter =
            BrukerScan2InvIonMobilityConverterFactory::instance(sdk_path)
              .produce(*handle, pcs);
          OPENMS_LOG_INFO << "TIMS calibration: Bruker SDK"
                          << (pcs != NoPressureCompensation ? " (pressure_comp=on)" : "")
                          << std::endl;
          return handle;
        }
        catch (const std::exception& e)
        {
          if (strategy == Strategy::BRUKER_SDK)
          {
            throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              String("Bruker SDK failed: ") + e.what());
          }
          OPENMS_LOG_DEBUG << "Bruker SDK not available (" << e.what()
                          << "), trying rational" << std::endl;
        }
      }
      else if (strategy == Strategy::BRUKER_SDK)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Bruker SDK path not set (use OPENMS_BRUKER_SDK_PATH or Config::bruker_sdk_path)");
      }
    }

    // 3. Try rational model from TimsCalibration table (AUTO or RATIONAL)
    if (strategy == Strategy::AUTO || strategy == Strategy::RATIONAL)
    {
      auto converter = tryCreateRationalConverter(path_string);
      if (converter)
      {
        handle->scan2inv_ion_mobility_converter = std::move(converter);
        return handle;
      }
      else if (strategy == Strategy::RATIONAL)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "TimsCalibration table not found or unsupported in " + path);
      }
    }

    // 4. Linear fallback (already set from step 1)
    OPENMS_LOG_INFO << "TIMS calibration: linear (GlobalMetadata)" << std::endl;
    return handle;
  }

  // =====================================================================
  // Helper: convert a single opentims frame to MSSpectrum
  // =====================================================================
  static void frameToSpectrum(TimsFrame& frame, MSSpectrum& spec, int ms_level)
  {
    spec.clear(true);
    spec.setRT(frame.time);
    spec.setMSLevel(ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    // TDF peaks are detector-centroided in m/z (peak-list format); label as
    // centroid so downstream tools (e.g. CometAdapter) don't reject as profile.
    spec.setType(SpectrumSettings::SpectrumType::CENTROID);
    spec.setNativeID("frame=" + String(frame.id));

    if (frame.num_peaks == 0) return;

    // Allocate output buffers
    std::vector<uint32_t> scan_ids(frame.num_peaks);
    std::vector<uint32_t> intensities(frame.num_peaks);
    std::vector<double> mzs(frame.num_peaks);
    std::vector<double> inv_ion_mobilities(frame.num_peaks);

    frame.save_to_buffs(nullptr, scan_ids.data(), nullptr, intensities.data(),
                        mzs.data(), inv_ion_mobilities.data(), nullptr);

    // Fill spectrum peaks
    spec.reserve(frame.num_peaks);
    for (uint32_t i = 0; i < frame.num_peaks; ++i)
    {
      Peak1D peak;
      peak.setMZ(mzs[i]);
      peak.setIntensity(static_cast<double>(intensities[i]));
      spec.push_back(peak);
    }

    // Set IM float data array with correct CV term
    DataArrays::FloatDataArray im_array;
    im_array.resize(frame.num_peaks);
    for (uint32_t i = 0; i < frame.num_peaks; ++i)
    {
      im_array[i] = static_cast<float>(inv_ion_mobilities[i]);
    }
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    spec.getFloatDataArrays().push_back(std::move(im_array));
    spec.setIMPeakType(IMPeakType::IM_PROFILE);
  }

  // =====================================================================
  // Helper: centroid a single MS1 frame
  // =====================================================================
  static void centroidMS1Frame(TimsFrame& frame, MSSpectrum& spec,
                               const BrukerTimsFile::Config& config,
                               FrameCentroider& centroider)
  {
    spec.clear(true);
    spec.setRT(frame.time);
    spec.setMSLevel(1);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    spec.setType(SpectrumSettings::SpectrumType::CENTROID);
    spec.setNativeID("frame=" + String(frame.id));

    if (frame.num_peaks == 0) return;

    // Extract frame data into FrameBuffers. save_to_buffs writes uint32
    // intensities; we widen to double for the unified centroider dispatch
    // (PASEFHillCentroider + FrameCentroider double overloads + PeakPickerIM).
    FrameBuffers buf;
    buf.resize(frame.num_peaks);
    std::vector<uint32_t> raw_intensities(frame.num_peaks);
    frame.save_to_buffs(nullptr, buf.scan_ids.data(), nullptr,
                        raw_intensities.data(),
                        buf.mz.data(), buf.im.data(), nullptr);
    for (uint32_t i = 0; i < frame.num_peaks; ++i)
      buf.intensity[i] = static_cast<double>(raw_intensities[i]);

    if (config.isotopic_prefilter)
    {
      const size_t before = buf.size();
      const size_t dropped = isotopicPrefilter(buf, config.isotopic_prefilter_tol_ppm);
      OPENMS_LOG_INFO << "Isotopic prefilter (MS1, frame=" << frame.id
                      << "): dropped " << dropped << " / " << before
                      << " peaks ("
                      << (100.0 * dropped / std::max<size_t>(1, before)) << "%)"
                      << std::endl;
      if (buf.empty()) return;
    }

    const auto algo = effectiveMS1Algo(config, /*emit_warnings=*/false);
    std::vector<double> cent_mz, cent_intensity;
    std::vector<float> cent_im;
    std::vector<float> b_im_lo, b_im_hi, b_mz_lo, b_mz_hi;
    const bool want_bounds = config.expose_hill_bounds &&
        algo == BrukerTimsFile::Config::CentroidAlgo::HILL_BASED;
    HillBoundsOut bounds = want_bounds
        ? HillBoundsOut{&b_im_lo, &b_im_hi, &b_mz_lo, &b_mz_hi}
        : HillBoundsOut{nullptr, nullptr, nullptr, nullptr};
    runMS1Centroider(algo, config, centroider,
                     buf.mz.data(), buf.intensity.data(),
                     buf.im.data(),
                     static_cast<uint32_t>(buf.size()),
                     cent_mz, cent_intensity, cent_im,
                     buf.scan_ids.data(), bounds);

    // Build MSSpectrum from centroided output
    spec.reserve(cent_mz.size());
    for (size_t i = 0; i < cent_mz.size(); ++i)
    {
      Peak1D peak;
      peak.setMZ(cent_mz[i]);
      peak.setIntensity(cent_intensity[i]);
      spec.push_back(peak);
    }

    // Set IM float data array
    DataArrays::FloatDataArray im_array;
    im_array.assign(cent_im.begin(), cent_im.end());
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    spec.getFloatDataArrays().push_back(std::move(im_array));
    if (want_bounds)
      attachHillBoundsArrays(spec, b_im_lo, b_im_hi, b_mz_lo, b_mz_hi);
    spec.setIMPeakType(IMPeakType::IM_CENTROIDED);
  }

  // =====================================================================
  // Helper: build one MS1 spectrum from a frame, with optional IM centroiding
  // =====================================================================
  static void loadMS1Spectrum(TimsFrame& frame, MSSpectrum& spec,
                              const BrukerTimsFile::Config& config,
                              FrameCentroider& centroider)
  {
    if (isCentroidingEnabled(config))
      centroidMS1Frame(frame, spec, config, centroider);
    else
      frameToSpectrum(frame, spec, 1);
  }

  // =====================================================================
  // Helper: build one aggregated MS1 spectrum from a sliding window of
  // adjacent MS1 frames. Composes with within-frame IM centroiding
  // (FrameCentroider) when ms1_centroid_mz_ppm/pct are set.
  //
  // Pre: center_idx is the output-emitting frame; N and RT cap come from
  // config. Post: spec is populated with center-frame metadata and
  // per-peak IM data. IMPeakType is IM_PROFILE or IM_CENTROIDED based on
  // isCentroidingEnabled(config).
  // =====================================================================
  static void loadAggregatedMS1Spectrum(
      TimsDataHandle& handle,
      const std::vector<uint32_t>& ms1_frame_ids,
      size_t center_idx,
      const BrukerTimsFile::Config& config,
      FrameAggregator& aggregator,
      FrameCentroider& centroider,
      MSSpectrum& spec)
  {
    const int N = config.ms1_n_neighbors;
    const size_t size = ms1_frame_ids.size();
    const size_t lo = (static_cast<size_t>(N) > center_idx) ? 0 : center_idx - N;
    const size_t hi = std::min(center_idx + static_cast<size_t>(N), size - 1);

    TimsFrame& center_frame = handle.get_frame(ms1_frame_ids[center_idx]);

    aggregator.clear();
    size_t contributing_frames = 0;

    for (size_t ni = lo; ni <= hi; ++ni)
    {
      TimsFrame& nframe = handle.get_frame(ms1_frame_ids[ni]);

      // Center is always included regardless of RT cap; neighbors are capped.
      const bool is_center = (ni == center_idx);
      if (!is_center && config.ms1_max_rt_distance_sec > 0.0 &&
          std::abs(nframe.time - center_frame.time) > config.ms1_max_rt_distance_sec)
      {
        continue;
      }

      if (nframe.num_peaks == 0) continue;

      std::vector<uint32_t> scan_ids(nframe.num_peaks);
      std::vector<uint32_t> intensities(nframe.num_peaks);
      std::vector<double> mzs(nframe.num_peaks);

      nframe.save_to_buffs(nullptr, scan_ids.data(), nullptr, intensities.data(),
                           mzs.data(), nullptr, nullptr);

      bool frame_contributed = false;
      for (uint32_t p = 0; p < nframe.num_peaks; ++p)
      {
        aggregator.addPeak(mzs[p], intensities[p], scan_ids[p]);
        frame_contributed = true;
      }
      if (frame_contributed) ++contributing_frames;
    }

    const bool skip_denoise = (contributing_frames <= 1) || (config.ms1_min_support <= 0);
    auto peaks = aggregator.finalize(config.ms1_min_support, skip_denoise);

    // Populate center-frame metadata regardless of whether peaks is empty —
    // preserves the per-center-frame cadence invariant for streaming consumers.
    spec.clear(true);
    spec.setRT(center_frame.time);
    spec.setMSLevel(1);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    spec.setNativeID("frame=" + String(center_frame.id));

    DataArrays::FloatDataArray im_array;
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

    if (peaks.empty())
    {
      spec.getFloatDataArrays().push_back(std::move(im_array));
      spec.setIMPeakType(IMPeakType::IM_PROFILE);
      return;
    }

    // Build FrameBuffers once for all algos (including OFF). scan_ids are
    // needed by HillBased (max_scan_gap) and the isotopic prefilter; for
    // OFF they're unused but cheap.
    FrameBuffers buf;
    aggregatorPeaksToBuffers(peaks, handle, center_frame.id, buf);

    if (config.isotopic_prefilter)
    {
      const size_t before = buf.size();
      const size_t dropped = isotopicPrefilter(buf, config.isotopic_prefilter_tol_ppm);
      OPENMS_LOG_INFO << "Isotopic prefilter (MS1 aggregated, frame="
                      << center_frame.id << "): dropped " << dropped << " / "
                      << before << " peaks ("
                      << (100.0 * dropped / std::max<size_t>(1, before)) << "%)"
                      << std::endl;
      if (buf.empty())
      {
        spec.getFloatDataArrays().push_back(std::move(im_array));
        spec.setIMPeakType(IMPeakType::IM_PROFILE);
        return;
      }
    }

    const auto algo = effectiveMS1Algo(config, /*emit_warnings=*/false);
    if (algo != BrukerTimsFile::Config::CentroidAlgo::OFF)
    {
      // Unpack aggregator output into parallel arrays for the centroider,
      // including per-peak scan_id so HillBased can do gap detection.
      std::vector<double> cent_mzs(peaks.size());
      std::vector<double> cent_intensities(peaks.size());
      std::vector<double> cent_ims(peaks.size());
      std::vector<uint32_t> cent_scan_ids(peaks.size());
      for (size_t i = 0; i < peaks.size(); ++i)
      {
        cent_mzs[i] = peaks[i].mz;
        cent_intensities[i] = peaks[i].intensity;
        double im_val = 0.0;
        handle.scan2inv_ion_mobility_converter->convert(
          center_frame.id, &im_val, &peaks[i].scan_id, 1);
        cent_ims[i] = im_val;
        cent_scan_ids[i] = peaks[i].scan_id;
      }

      std::vector<double> out_mz, out_intensity;
      std::vector<float> out_im;
      std::vector<float> b_im_lo, b_im_hi, b_mz_lo, b_mz_hi;
      const bool want_bounds = config.expose_hill_bounds &&
          algo == BrukerTimsFile::Config::CentroidAlgo::HILL_BASED;
      HillBoundsOut bounds = want_bounds
          ? HillBoundsOut{&b_im_lo, &b_im_hi, &b_mz_lo, &b_mz_hi}
          : HillBoundsOut{nullptr, nullptr, nullptr, nullptr};
      runMS1Centroider(algo, config, centroider,
                       buf.mz.data(), buf.intensity.data(),
                       buf.im.data(),
                       static_cast<uint32_t>(buf.size()),
                       out_mz, out_intensity, out_im,
                       buf.scan_ids.data(), bounds);

      spec.reserve(out_mz.size());
      for (size_t i = 0; i < out_mz.size(); ++i)
      {
        Peak1D peak;
        peak.setMZ(out_mz[i]);
        peak.setIntensity(static_cast<float>(out_intensity[i]));
        spec.push_back(peak);
      }
      im_array.assign(out_im.begin(), out_im.end());
      spec.setType(SpectrumSettings::SpectrumType::CENTROID);
      spec.setIMPeakType(IMPeakType::IM_CENTROIDED);
      spec.getFloatDataArrays().push_back(std::move(im_array));
      if (want_bounds)
        attachHillBoundsArrays(spec, b_im_lo, b_im_hi, b_mz_lo, b_mz_hi);
      return;
    }
    // Off path — emit from prefiltered (or raw) buf.
    spec.reserve(buf.size());
    im_array.reserve(buf.size());
    for (size_t i = 0; i < buf.size(); ++i)
    {
      spec.emplace_back(buf.mz[i], static_cast<float>(buf.intensity[i]));
      im_array.push_back(static_cast<float>(buf.im[i]));
    }
    spec.setIMPeakType(IMPeakType::IM_PROFILE);
    spec.getFloatDataArrays().push_back(std::move(im_array));
  }

  // =====================================================================
  // isDIA_: detect DDA vs DIA by querying for SWATH window tables
  // =====================================================================
  bool BrukerTimsFile::isDIA_(const String& tdf_path) const
  {
    try
    {
      SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);
      return isDIA(db);
    }
    catch (const SQLite::Exception& e)
    {
      OPENMS_LOG_WARN << "Warning: could not query DIA status from " << tdf_path
                      << ": " << e.what() << std::endl;
      return false;
    }
  }

  // =====================================================================
  // loadExperimentalSettings_: populate SourceFile metadata
  // =====================================================================
  void BrukerTimsFile::loadExperimentalSettings_(const String& path, ExperimentalSettings& settings)
  {
    SourceFile sf;
    sf.setNameOfFile(File::basename(path));
    sf.setPathToFile(File::path(path));
    sf.setFileType("Bruker TDF");
    sf.setNativeIDType("Bruker TDF nativeID format");
    sf.setNativeIDTypeAccession("MS:1002818");
    settings.getSourceFiles().push_back(sf);
  }

  // =====================================================================
  // readDIAMetadata: SQL-only boundary and count extraction
  // =====================================================================
  BrukerTimsFile::DIAStreamingMetadata BrukerTimsFile::readDIAMetadata(
      const String& path, ExperimentalSettings& exp_settings)
  {
    return readDIAMetadata(path, exp_settings, Config());
  }

  BrukerTimsFile::DIAStreamingMetadata BrukerTimsFile::readDIAMetadata(
      const String& path, ExperimentalSettings& exp_settings, const Config& config)
  {
    auto handle = openTimsDataHandle(path, config);
    String tdf_path = path + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    if (!isDIA(db))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "readDIAMetadata() requires a DIA dataset, but '" + path + "' appears to be DDA.");
    }

    loadExperimentalSettings_(path, exp_settings);

    // Read DIA windows (with IM conversion via handle's calibration)
    auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);

    DIAStreamingMetadata meta;
    if (windows.empty())
    {
      OPENMS_LOG_WARN << "Warning: DIA dataset detected but no SWATH windows found in '" << path << "'" << std::endl;
      return meta;
    }

    // Build SwathMap boundaries from DIAWindow structs
    meta.boundaries.reserve(windows.size());
    for (const auto& w : windows)
    {
      meta.boundaries.emplace_back(
        w.mz_center - w.mz_width / 2.0,  // lower
        w.mz_center + w.mz_width / 2.0,  // upper
        w.mz_center,                       // center
        w.im_lower,                        // imLower
        w.im_upper,                        // imUpper
        false);                            // ms1 = false
    }

    const auto eff = resolveEffectiveConfig(db, config,
                                             handle->min_frame_id(),
                                             handle->max_frame_id());

    // Count MS1 frames
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (!frameInRange(fid, eff)) continue;
      if (eff.load_ms1 && handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ++meta.nr_ms1_spectra;
    }

    // Count MS2 spectra per window (= frames per WindowGroup).
    // Apply the frame range filter so counts match what loadDIAStreaming
    // will actually emit (CachedSwathFileConsumer sizing depends on it).
    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
    for (auto& [group, frame_ids] : group_to_frames)
    {
      frame_ids.erase(
        std::remove_if(frame_ids.begin(), frame_ids.end(),
                       [&eff](uint32_t fid) { return !frameInRange(fid, eff); }),
        frame_ids.end());
    }
    meta.nr_ms2_spectra.reserve(windows.size());
    for (const auto& w : windows)
    {
      auto it = group_to_frames.find(w.window_group);
      meta.nr_ms2_spectra.push_back(
        it != group_to_frames.end() ? static_cast<int>(it->second.size()) : 0);
    }

    return meta;
  }

  // =====================================================================
  // loadDIAStreaming: stream DIA spectra to consumer one-at-a-time
  // =====================================================================
  void BrukerTimsFile::loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer)
  {
    loadDIAStreaming(path, consumer, Config());
  }

  void BrukerTimsFile::loadDIAStreaming(
      const String& path, FullSwathFileConsumer& consumer, const Config& config)
  {
    auto handle = openTimsDataHandle(path, config);
    String tdf_path = handle->get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    const auto eff = resolveEffectiveConfig(db, config,
                                             handle->min_frame_id(),
                                             handle->max_frame_id());

    // --- MS1 frames ---
    FrameCentroider centroider;
    FrameAggregator ms1_aggregator(0.01);

    std::vector<uint32_t> ms1_frame_ids;
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (!frameInRange(fid, eff)) continue;
      if (handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ms1_frame_ids.push_back(fid);
    }
    warnPartialMS1AggregationConfig(eff, ms1_frame_ids.size());
    if (eff.ms1_n_neighbors > 0) ms1_aggregator.reserve(300'000);

    const size_t ms1_to_emit = eff.load_ms1 ? ms1_frame_ids.size() : 0;
    startProgress(0, ms1_to_emit, "Streaming DIA-PASEF MS1 frames");
    for (size_t i = 0; i < ms1_to_emit; ++i)
    {
      MSSpectrum spec;
      if (eff.ms1_n_neighbors > 0)
      {
        loadAggregatedMS1Spectrum(*handle, ms1_frame_ids, i, eff,
                                  ms1_aggregator, centroider, spec);
      }
      else
      {
        TimsFrame& frame = handle->get_frame(ms1_frame_ids[i]);
        loadMS1Spectrum(frame, spec, eff, centroider);
      }
      // Sort peaks by m/z (and the associated IM float data array alongside).
      // TIMS save_to_buffs returns peaks in (scan_id, m/z-within-scan) order,
      // which is NOT globally m/z-sorted. Downstream consumers
      // (OpenSwath's chromatogram extraction) assume sorted m/z for lower_bound /
      // upper_bound range queries; without this sort, queries silently miss
      // peaks, producing empty chromatograms and breaking iRT calibration.
      // The non-streaming load() path gets this sort for free via exp.sortSpectra(true).
      spec.sortByPosition();
      consumer.consumeSpectrum(spec);
      setProgress(i);
    }
    endProgress();
    centroider.reportCapSummary(static_cast<size_t>(eff.ms1_centroid_max_peaks));

    // --- MS2 frames: raw per-WindowGroup iteration (no aggregation) ---
    auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);
    if (windows.empty())
    {
      OPENMS_LOG_WARN << "Warning: DIA dataset detected but no SWATH windows found" << std::endl;
      return;
    }

    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
    for (auto& [group, frame_ids] : group_to_frames)
    {
      frame_ids.erase(
        std::remove_if(frame_ids.begin(), frame_ids.end(),
                       [&eff](uint32_t fid) { return !frameInRange(fid, eff); }),
        frame_ids.end());
    }

    std::map<int, std::vector<const DIAWindow*>> group_to_windows;
    for (const auto& w : windows)
      group_to_windows[w.window_group].push_back(&w);

    Size total_work = 0;
    for (const auto& [group, frames] : group_to_frames)
    {
      auto wit = group_to_windows.find(group);
      if (wit != group_to_windows.end())
        total_work += frames.size() * wit->second.size();
    }

    startProgress(0, total_work, "Streaming DIA-PASEF MS2 frames");
    Size progress_count = 0;

    for (const auto& [group, frame_ids] : group_to_frames)
    {
      auto wit = group_to_windows.find(group);
      if (wit == group_to_windows.end()) continue;
      const auto& dia_windows = wit->second;

      for (const DIAWindow* win : dia_windows)
      {
        for (size_t i = 0; i < frame_ids.size(); ++i)
        {
          setProgress(progress_count++);
          TimsFrame& frame = handle->get_frame(frame_ids[i]);
          if (frame.num_peaks == 0) continue;

          std::vector<uint32_t> scan_ids(frame.num_peaks);
          std::vector<uint32_t> intensities(frame.num_peaks);
          std::vector<double> mzs(frame.num_peaks);
          std::vector<double> inv_ion_mobilities(frame.num_peaks);

          frame.save_to_buffs(nullptr, scan_ids.data(), nullptr, intensities.data(),
                              mzs.data(), inv_ion_mobilities.data(), nullptr);

          MSSpectrum spec;
          spec.setRT(frame.time);
          spec.setMSLevel(2);
          spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
          spec.setType(SpectrumSettings::SpectrumType::CENTROID);
          spec.setNativeID("frame=" + String(frame.id) + " windowGroup=" + String(win->window_group) + " scan=" + String(win->scan_begin));

          Precursor prec;
          prec.setMZ(win->mz_center);
          prec.setIsolationWindowLowerOffset(win->mz_width / 2.0);
          prec.setIsolationWindowUpperOffset(win->mz_width / 2.0);
          spec.setPrecursors({prec});

          spec.setMetaValue("ion mobility lower limit", win->im_lower);
          spec.setMetaValue("ion mobility upper limit", win->im_upper);

          DataArrays::FloatDataArray im_array;
          IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

          for (uint32_t p = 0; p < frame.num_peaks; ++p)
          {
            if (inv_ion_mobilities[p] >= win->im_lower && inv_ion_mobilities[p] <= win->im_upper)
            {
              Peak1D peak;
              peak.setMZ(mzs[p]);
              peak.setIntensity(static_cast<double>(intensities[p]));
              spec.push_back(peak);
              im_array.push_back(static_cast<float>(inv_ion_mobilities[p]));
            }
          }

          if (!spec.empty())
          {
            spec.getFloatDataArrays().push_back(std::move(im_array));
            spec.setIMPeakType(IMPeakType::IM_PROFILE);
            // Sort peaks by m/z (and the IM float data array alongside).
            // See the matching comment in the MS1 loop above — TIMS save_to_buffs
            // returns peaks in (scan_id, m/z-within-scan) order; OpenSwath's
            // chromatogram extraction assumes globally m/z-sorted input.
            spec.sortByPosition();
            consumer.consumeSpectrum(spec);
          }
        }
      }
    }
    endProgress();
  }

  // =====================================================================
  // load() overloads
  // =====================================================================
  void BrukerTimsFile::load(const String& path, MSExperiment& exp)
  {
    load(path, exp, Config());
  }

  void BrukerTimsFile::load(const String& path, MSExperiment& exp, const Config& config)
  {
    exp.clear(true);
    auto handle = openTimsDataHandle(path, config);

    String tdf_path = path + "/analysis.tdf";

    // Resolve RT range (if any) to an effective frame_id range; also
    // validates the user-supplied frame_id/rt ranges and emits warnings.
    Config eff;
    {
      SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);
      eff = resolveEffectiveConfig(db, config,
                                    handle->min_frame_id(),
                                    handle->max_frame_id());
    }

    bool is_dia = isDIA_(tdf_path);
    Config::ExportMode mode = eff.export_mode;

    if (mode == Config::FRAME)
    {
      loadFrames_(*handle, exp, eff);
    }
    else if (mode == Config::SPECTRUM || (mode == Config::AUTO && !is_dia))
    {
      loadDDA_(*handle, exp, eff);
    }
    else // AUTO + DIA
    {
      loadDIA_(*handle, exp, eff);
    }

    loadExperimentalSettings_(path, exp);

    // Sort by RT, interleaved across MS levels
    exp.sortSpectra(true);
    exp.updateRanges();
  }

  // =====================================================================
  // transform() overloads
  // =====================================================================
  void BrukerTimsFile::transform(const String& path, Interfaces::IMSDataConsumer* consumer)
  {
    transform(path, consumer, Config());
  }

  void BrukerTimsFile::transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config)
  {
    auto handle = openTimsDataHandle(path, config);

    String tdf_path = path + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    const auto eff = resolveEffectiveConfig(db, config,
                                             handle->min_frame_id(),
                                             handle->max_frame_id());

    bool is_dia = isDIA(db);
    Config::ExportMode mode = eff.export_mode;

    // Compute expected size for consumer. Must mirror what loadFrames_/loadDIA_/
    // loadDDA_ actually emit under the active frame_id range and load_ms1 knob.
    FrameCounts counts = readFrameCounts(db, eff);
    const size_t ms1_emitted = eff.load_ms1 ? counts.ms1 : 0;
    size_t expected = 0;
    if (mode == Config::FRAME)
    {
      expected = ms1_emitted + counts.ms2;
    }
    else if (is_dia && mode != Config::SPECTRUM)
    {
      // DIA: each MS2 frame belongs to exactly ONE WindowGroup and only
      // produces spectra for that group's windows. So the emit count is
      //   Σ_group (frames_in_group × windows_in_group)
      // — NOT frames × all_windows. Mirrors the total_work computation in
      // loadDIA_.
      auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);
      auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
      for (auto& [group, frame_ids] : group_to_frames)
      {
        frame_ids.erase(
          std::remove_if(frame_ids.begin(), frame_ids.end(),
                         [&eff](uint32_t fid) { return !frameInRange(fid, eff); }),
          frame_ids.end());
      }
      std::map<int, size_t> group_window_count;
      for (const auto& w : windows) ++group_window_count[w.window_group];
      size_t dia_ms2_spectra = 0;
      for (const auto& [group, frames] : group_to_frames)
      {
        auto wit = group_window_count.find(group);
        if (wit != group_window_count.end())
          dia_ms2_spectra += frames.size() * wit->second;
      }
      expected = ms1_emitted + dia_ms2_spectra;
    }
    else
    {
      // DDA: MS1 frames + per-precursor MS2 spectra
      uint32_t num_precursors = countDDAPrecursors(db, eff);
      expected = ms1_emitted + num_precursors;
    }

    consumer->setExpectedSize(expected, 0);

    // Populate source file metadata (same as load())
    ExperimentalSettings settings;
    loadExperimentalSettings_(path, settings);
    consumer->setExperimentalSettings(settings);

    // NOTE: This loads into a temporary experiment then feeds to consumer.
    // Not truly constant-memory -- a future optimization should iterate
    // frame-by-frame and call consumer->consumeSpectrum() inline.
    OPENMS_LOG_INFO << "BrukerTimsFile::transform(): loading full dataset (streaming optimization pending)" << std::endl;
    MSExperiment exp;
    if (mode == Config::FRAME)
      loadFrames_(*handle, exp, eff);
    else if (is_dia && mode != Config::SPECTRUM)
      loadDIA_(*handle, exp, eff);
    else
      loadDDA_(*handle, exp, eff);

    exp.sortSpectra(true);
    for (auto& spec : exp)
    {
      consumer->consumeSpectrum(spec);
    }
  }

  // =====================================================================
  // loadDDA_: MS1 frames (IM_PEAK) + MS2 spectra (scalar IM)
  // =====================================================================
  void BrukerTimsFile::loadDDA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config)
  {
    String tdf_path = handle.get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    // Collect frame IDs by MS level
    std::vector<uint32_t> ms1_frame_ids;
    std::vector<uint32_t> ms2_frame_ids;

    for (uint32_t fid = handle.min_frame_id(); fid <= handle.max_frame_id(); ++fid)
    {
      if (!frameInRange(fid, config)) continue;
      if (!handle.has_frame(fid)) continue;
      TimsFrame& frame = handle.get_frame(fid);
      if (frame.msms_type == 0)
        ms1_frame_ids.push_back(fid);
      else
        ms2_frame_ids.push_back(fid);
    }
    warnPartialMS1AggregationConfig(config, ms1_frame_ids.size());

    // Read DDA precursors from SQL. Drop entries whose MS2 frame is outside
    // the active range so the per-precursor MS2 spectra match what the
    // ms2_frame_ids filter selected above; the existing has_frame() safety
    // checks in the loops downstream silently skip any orphans, but stripping
    // them here keeps progress counts and the precursor_groups loop honest.
    std::vector<DDAPrecursorInfo> precursor_entries = readDDAPrecursors(db);
    const bool range_active =
      (config.frame_id_min != 0) ||
      (config.frame_id_max != std::numeric_limits<uint32_t>::max());
    if (range_active)
    {
      // Tally how many distinct precursors had at least one source MS2 frame
      // dropped by the range filter (= partial reconstruction). One-shot
      // warning at the end so the user knows the reconstruction is partial
      // at the range boundaries.
      std::map<int, std::pair<size_t, size_t>> per_pid_kept_dropped;  // pid -> (kept, dropped)
      for (const auto& e : precursor_entries)
      {
        auto& kd = per_pid_kept_dropped[e.precursor_id];
        if (frameInRange(e.frame_id, config)) ++kd.first;
        else                                   ++kd.second;
      }
      size_t partial = 0;
      size_t fully_dropped = 0;
      for (const auto& [pid, kd] : per_pid_kept_dropped)
      {
        if (kd.first == 0)                      ++fully_dropped;
        else if (kd.second > 0)                 ++partial;
      }
      if (partial > 0)
      {
        OPENMS_LOG_WARN << "Warning: " << partial
                        << " DDA precursor(s) had source MS2 frames outside the active "
                        << "frame range [" << config.frame_id_min << ", "
                        << config.frame_id_max << "]; their MS2 spectra will be "
                        << "reconstructed from the in-range subset only (partial "
                        << "spectra). Widen the range to avoid this." << std::endl;
      }
      if (fully_dropped > 0)
      {
        OPENMS_LOG_INFO << fully_dropped
                        << " DDA precursor(s) fully outside the active frame range "
                        << "were skipped." << std::endl;
      }
    }
    precursor_entries.erase(
      std::remove_if(precursor_entries.begin(), precursor_entries.end(),
                     [&config](const DDAPrecursorInfo& e) {
                       return !frameInRange(e.frame_id, config);
                     }),
      precursor_entries.end());

    // Group by precursor ID
    std::map<int, std::vector<const DDAPrecursorInfo*>> precursor_groups;
    for (const auto& entry : precursor_entries)
    {
      precursor_groups[entry.precursor_id].push_back(&entry);
    }

    uint32_t num_ms2 = static_cast<uint32_t>(precursor_groups.size());
    const size_t ms1_to_emit = config.load_ms1 ? ms1_frame_ids.size() : 0;
    exp.reserveSpaceSpectra(static_cast<unsigned int>(ms1_to_emit) + num_ms2);

    startProgress(0, ms1_to_emit + num_ms2, "Loading DDA-PASEF data");

    // --- MS1 frames ---
    FrameCentroider centroider;
    FrameAggregator ms1_aggregator(0.01);
    if (config.load_ms1 && config.ms1_n_neighbors > 0) ms1_aggregator.reserve(300'000);
    for (size_t i = 0; i < ms1_to_emit; ++i)
    {
      MSSpectrum spec;
      if (config.ms1_n_neighbors > 0)
      {
        loadAggregatedMS1Spectrum(handle, ms1_frame_ids, i, config,
                                  ms1_aggregator, centroider, spec);
      }
      else
      {
        TimsFrame& frame = handle.get_frame(ms1_frame_ids[i]);
        loadMS1Spectrum(frame, spec, config, centroider);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }

    // --- MS2 spectra (one per precursor, reconstructed from frames) ---
    // Pre-extract frame data for all MS2 frames that are referenced by precursors.
    // Cache: frame_id -> (scan_ids, intensities, mzs, inv_ion_mobilities)
    struct FrameData
    {
      std::vector<uint32_t> scan_ids;
      std::vector<uint32_t> tofs;
      std::vector<uint32_t> intensities;
      std::vector<double> mzs;
      std::vector<double> inv_ion_mobilities;
      std::vector<uint32_t> scan_offsets; // scan_offsets[s] = first peak index for scan s; size = num_scans + 1
      double rt = 0.0;

      // Get peak index range [begin, end) for scans in [scan_begin, scan_end)
      std::pair<uint32_t, uint32_t> peakRangeForScans(uint32_t scan_begin, uint32_t scan_end) const
      {
        if (scan_offsets.empty() || scan_begin >= scan_offsets.size() - 1) return {0, 0};
        scan_end = std::min(scan_end, static_cast<uint32_t>(scan_offsets.size() - 1));
        return {scan_offsets[scan_begin], scan_offsets[scan_end]};
      }
    };

    // Determine which frames are needed
    std::set<uint32_t> needed_frames;
    for (const auto& entry : precursor_entries)
    {
      needed_frames.insert(entry.frame_id);
    }

    // Extract data for needed frames (lazy, on demand per precursor group)
    // We cache extracted frame data to avoid re-extraction when multiple
    // precursors share the same frame.
    std::unordered_map<uint32_t, FrameData> frame_cache;

    auto getFrameData = [&](uint32_t frame_id) -> const FrameData&
    {
      auto it = frame_cache.find(frame_id);
      if (it != frame_cache.end()) return it->second;

      TimsFrame& frame = handle.get_frame(frame_id);
      FrameData fd;
      fd.rt = frame.time;

      if (frame.num_peaks > 0)
      {
        fd.scan_ids.resize(frame.num_peaks);
        fd.tofs.resize(frame.num_peaks);
        fd.intensities.resize(frame.num_peaks);
        fd.mzs.resize(frame.num_peaks);
        fd.inv_ion_mobilities.resize(frame.num_peaks);

        frame.save_to_buffs(nullptr, fd.scan_ids.data(), fd.tofs.data(),
                            fd.intensities.data(), fd.mzs.data(),
                            fd.inv_ion_mobilities.data(), nullptr);

        // Build scan_offsets for O(1) scan-range lookups.
        // Peaks are ordered by scan (opentims fills scan_ids monotonically).
        fd.scan_offsets.resize(frame.num_scans + 1, frame.num_peaks);
        uint32_t current_scan = 0;
        fd.scan_offsets[0] = 0;
        for (uint32_t p = 0; p < frame.num_peaks; ++p)
        {
          while (current_scan < fd.scan_ids[p])
          {
            ++current_scan;
            fd.scan_offsets[current_scan] = p;
          }
        }
        // Fill remaining scans that have no peaks
        for (uint32_t s = current_scan + 1; s <= frame.num_scans; ++s)
        {
          fd.scan_offsets[s] = frame.num_peaks;
        }
      }

      auto result = frame_cache.emplace(frame_id, std::move(fd));
      return result.first->second;
    };

    // OLS recalibration for DDA (when config.calibrate == true)
    const double cal_tol = (config.calibration_tolerance > 0.0) ? config.calibration_tolerance : 0.1;
    if (config.calibrate)
    {
      // Collect (tof_index, sqrt(monoisotopic_mz)) pairs for regression
      std::vector<double> tof_vals;
      std::vector<double> sqrt_mz_vals;

      for (const auto& entry : precursor_entries)
      {
        if (!handle.has_frame(entry.frame_id)) continue;
        if (std::isnan(entry.monoisotopic_mz) || entry.monoisotopic_mz <= 0.0) continue;

        const FrameData& fd = getFrameData(entry.frame_id);
        if (fd.mzs.empty()) continue;

        // Find the highest-intensity peak in the scan range [scan_begin, scan_end)
        auto [p_begin, p_end] = fd.peakRangeForScans(entry.scan_begin, entry.scan_end);
        double best_intensity = -1.0;
        uint32_t best_tof = 0;
        double best_mz = 0.0;

        for (uint32_t p = p_begin; p < p_end; ++p)
        {
          if (static_cast<double>(fd.intensities[p]) > best_intensity)
          {
            best_intensity = static_cast<double>(fd.intensities[p]);
            best_tof = fd.tofs[p];
            best_mz = fd.mzs[p];
          }
        }

        if (best_intensity > 0.0 &&
            std::abs(best_mz - entry.monoisotopic_mz) < cal_tol)
        {
          tof_vals.push_back(static_cast<double>(best_tof));
          sqrt_mz_vals.push_back(std::sqrt(entry.monoisotopic_mz));
        }
      }

      if (tof_vals.size() >= 2)
      {
        // Simple OLS regression: sqrt(mz) = intercept + slope * tof
        double n = static_cast<double>(tof_vals.size());
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_xx = 0;
        for (size_t i = 0; i < tof_vals.size(); ++i)
        {
          sum_x += tof_vals[i];
          sum_y += sqrt_mz_vals[i];
          sum_xy += tof_vals[i] * sqrt_mz_vals[i];
          sum_xx += tof_vals[i] * tof_vals[i];
        }
        double denom = n * sum_xx - sum_x * sum_x;
        if (std::abs(denom) > 1e-15)
        {
          double new_slope = (n * sum_xy - sum_x * sum_y) / denom;
          double new_intercept = (sum_y - new_slope * sum_x) / n;

          OPENMS_LOG_INFO << "OLS recalibration: " << tof_vals.size() << " calibration points, "
                          << "new intercept=" << new_intercept << ", slope=" << new_slope << std::endl;

          // Update the converter
          auto* converter = dynamic_cast< ::OpenSourceTof2MzConverter*>(handle.tof2mz_converter.get());
          if (converter)
          {
            converter->updateCalibration(new_intercept, new_slope);
            // Clear the frame cache so frames are re-extracted with new calibration
            frame_cache.clear();
          }
          else
          {
            OPENMS_LOG_WARN << "Warning: OLS recalibration requested but converter type is not OpenSourceTof2MzConverter" << std::endl;
          }
        }
      }
      else
      {
        OPENMS_LOG_WARN << "Warning: OLS recalibration requested but only "
                        << tof_vals.size() << " calibration points found (need >= 2)" << std::endl;
      }
    }

    // Build one MS2 spectrum per precursor
    // TOF-domain smoothing and centroiding windows (matching timsrust defaults)
    constexpr uint32_t SMOOTH_WINDOW = 1;
    constexpr uint32_t CENTROID_WINDOW = 1;

    uint32_t precursor_idx = 0;
    for (const auto& [prec_id, entries] : precursor_groups)
    {
      // Collect raw peaks (TOF space) from all frames for this precursor
      std::vector<uint32_t> all_tofs;
      std::vector<uint64_t> all_intensities;
      std::vector<double> all_im;
      double rt_sum = 0.0;

      for (const DDAPrecursorInfo* entry : entries)
      {
        if (!handle.has_frame(entry->frame_id)) continue;

        const FrameData& fd = getFrameData(entry->frame_id);
        rt_sum += fd.rt;

        // Extract peaks in the scan range [scan_begin, scan_end) via precomputed offsets
        auto [p_begin, p_end] = fd.peakRangeForScans(entry->scan_begin, entry->scan_end);
        for (uint32_t p = p_begin; p < p_end; ++p)
        {
          all_tofs.push_back(fd.tofs[p]);
          all_intensities.push_back(static_cast<uint64_t>(fd.intensities[p]));
          all_im.push_back(fd.inv_ion_mobilities[p]);
        }
      }

      if (all_tofs.empty())
      {
        ++precursor_idx;
        setProgress(ms1_to_emit + precursor_idx);
        continue;
      }

      // Branch: HillBased DDA-MS2 centroiding bypasses the TOF-domain pipeline
      // entirely. It converts each raw peak's TOF to m/z, then runs IM-axis
      // hill detection on the resulting (m/z, intensity, IM) triples. Output
      // schema is unchanged: scalar drift_time per spectrum (no per-peak IM
      // array). For DDA-MS2, peaks come from possibly several PASEF frames
      // covering one precursor; per-frame IM-calibration jitter is absorbed
      // by the linker (it chains across IM-sorted scans by m/z proximity, not
      // IM proximity). Greedy/Off paths fall through to the existing TOF
      // pipeline below.
      const auto dda_ms2_algo = effectiveMS2Algo(config, /*emit_warnings=*/false);
      if (dda_ms2_algo == BrukerTimsFile::Config::CentroidAlgo::HILL_BASED)
      {
        std::vector<double> hb_mz(all_tofs.size());
        std::vector<double> hb_int(all_tofs.size());
        for (std::size_t i = 0; i < all_tofs.size(); ++i)
        {
          handle.tof2mz_converter->convert(0, &hb_mz[i], &all_tofs[i], 1);
          hb_int[i] = static_cast<double>(all_intensities[i]);
        }
        Internal::PASEFHillCentroider::Params p;
        p.mz_tol_ppm        = static_cast<double>(config.ms2_centroid_mz_ppm);
        p.valley_factor     = config.centroid_valley_factor;
        p.min_hill_length   = config.ms2_centroid_min_hill_length;
        auto centroids = Internal::PASEFHillCentroider::centroidFrame(
            hb_mz.data(), hb_int.data(), all_im.data(),
            all_tofs.size(), p);

        if (centroids.empty())
        {
          ++precursor_idx;
          setProgress(ms1_to_emit + precursor_idx);
          continue;
        }

        std::vector<double> all_mz, all_intensity_d;
        std::vector<double> kept_im;
        all_mz.reserve(centroids.size());
        all_intensity_d.reserve(centroids.size());
        kept_im.reserve(centroids.size());
        for (const auto& c : centroids)
        {
          all_mz.push_back(c.mz);
          all_intensity_d.push_back(c.intensity);
          kept_im.push_back(c.im_apex);
        }

        std::vector<size_t> sort_idx(all_mz.size());
        std::iota(sort_idx.begin(), sort_idx.end(), 0u);
        std::sort(sort_idx.begin(), sort_idx.end(),
                  [&all_mz](size_t a, size_t b) { return all_mz[a] < all_mz[b]; });

        double im_weighted_sum = 0.0;
        double intensity_sum = 0.0;
        for (std::size_t i = 0; i < all_mz.size(); ++i)
        {
          im_weighted_sum += kept_im[i] * all_intensity_d[i];
          intensity_sum   += all_intensity_d[i];
        }
        double scalar_im = (intensity_sum > 0.0) ? (im_weighted_sum / intensity_sum) : 0.0;

        const DDAPrecursorInfo* first_entry = entries.front();
        double avg_rt = rt_sum / static_cast<double>(entries.size());

        MSSpectrum spec;
        spec.setRT(avg_rt);
        spec.setMSLevel(2);
        spec.setDriftTime(scalar_im);
        spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
        spec.setType(SpectrumSettings::SpectrumType::CENTROID);
        spec.setNativeID("frame=" + String(first_entry->frame_id)
                         + " scan=" + String(first_entry->scan_begin)
                         + " precursor=" + String(prec_id));
        spec.setMetaValue("bruker_precursor_id", static_cast<int>(prec_id));

        spec.reserve(all_mz.size());
        for (size_t idx : sort_idx)
        {
          Peak1D peak;
          peak.setMZ(all_mz[idx]);
          peak.setIntensity(all_intensity_d[idx]);
          spec.push_back(peak);
        }

        Precursor prec;
        prec.setMZ(first_entry->monoisotopic_mz);
        prec.setCharge(first_entry->charge);
        if (!std::isnan(first_entry->intensity))
          prec.setIntensity(static_cast<float>(first_entry->intensity));
        prec.setDriftTime(scalar_im);
        prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
        double iso_offset_lower = first_entry->isolation_width / 2.0 +
                                  (first_entry->monoisotopic_mz - first_entry->isolation_mz);
        double iso_offset_upper = first_entry->isolation_width / 2.0 -
                                  (first_entry->monoisotopic_mz - first_entry->isolation_mz);
        prec.setIsolationWindowLowerOffset(std::max(0.0, iso_offset_lower));
        prec.setIsolationWindowUpperOffset(std::max(0.0, iso_offset_upper));
        std::vector<Precursor> precursors;
        precursors.push_back(std::move(prec));
        spec.setPrecursors(std::move(precursors));

        exp.addSpectrum(std::move(spec));
        ++precursor_idx;
        setProgress(ms1_to_emit + precursor_idx);
        continue;
      }

      // --- TOF-domain processing pipeline (ported from timsrust) ---
      // Step 1: group_and_sum — sort by TOF, merge duplicates
      // We need to co-sort IM values with the TOF/intensity arrays.
      // Build a sort order, apply it, then group.
      {
        std::vector<size_t> order(all_tofs.size());
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(),
                  [&all_tofs](size_t a, size_t b) { return all_tofs[a] < all_tofs[b]; });

        std::vector<uint32_t> sorted_tofs(all_tofs.size());
        std::vector<uint64_t> sorted_int(all_intensities.size());
        std::vector<double> sorted_im(all_im.size());
        for (size_t i = 0; i < order.size(); ++i)
        {
          sorted_tofs[i] = all_tofs[order[i]];
          sorted_int[i] = all_intensities[order[i]];
          sorted_im[i] = all_im[order[i]];
        }

        // Merge consecutive entries with same TOF index (intensity-weighted IM merge)
        size_t write = 0;
        for (size_t read = 0; read < sorted_tofs.size(); ++read)
        {
          if (write > 0 && sorted_tofs[read] == all_tofs[write - 1])
          {
            // Weighted IM update: new_im = (old_im * old_int + new_im * new_int) / total_int
            double total = static_cast<double>(all_intensities[write - 1] + sorted_int[read]);
            if (total > 0)
              all_im[write - 1] = (all_im[write - 1] * static_cast<double>(all_intensities[write - 1])
                                   + sorted_im[read] * static_cast<double>(sorted_int[read])) / total;
            all_intensities[write - 1] += sorted_int[read];
          }
          else
          {
            all_tofs[write] = sorted_tofs[read];
            all_intensities[write] = sorted_int[read];
            all_im[write] = sorted_im[read];
            ++write;
          }
        }
        all_tofs.resize(write);
        all_intensities.resize(write);
        all_im.resize(write);
      }

      // Step 2: smooth — symmetric neighbor intensity sharing in TOF space
      smoothTofSpectrum(all_tofs, all_intensities, SMOOTH_WINDOW);

      // Step 3: centroid — sparse local maximum apex picking in TOF space
      std::vector<bool> keep = findSparseLocalMaxima(all_tofs, all_intensities, CENTROID_WINDOW);

      // Apply centroid mask and convert TOF to m/z
      std::vector<double> all_mz;
      std::vector<double> all_intensity_d;
      std::vector<double> kept_im;
      all_mz.reserve(all_tofs.size());
      all_intensity_d.reserve(all_tofs.size());
      kept_im.reserve(all_tofs.size());

      for (size_t i = 0; i < all_tofs.size(); ++i)
      {
        if (keep[i])
        {
          // Convert TOF to m/z using the registered converter
          double mz = 0;
          handle.tof2mz_converter->convert(0, &mz, &all_tofs[i], 1);
          all_mz.push_back(mz);
          all_intensity_d.push_back(static_cast<double>(all_intensities[i]));
          kept_im.push_back(all_im[i]);
        }
      }

      if (all_mz.empty())
      {
        ++precursor_idx;
        setProgress(ms1_to_emit + precursor_idx);
        continue;
      }

      // Sort by m/z
      std::vector<size_t> sort_idx(all_mz.size());
      std::iota(sort_idx.begin(), sort_idx.end(), 0u);
      std::sort(sort_idx.begin(), sort_idx.end(),
                [&all_mz](size_t a, size_t b) { return all_mz[a] < all_mz[b]; });

      // Compute intensity-weighted mean IM (scalar drift time for DDA MS2)
      double im_weighted_sum = 0.0;
      double intensity_sum = 0.0;
      for (size_t i = 0; i < all_mz.size(); ++i)
      {
        im_weighted_sum += kept_im[i] * all_intensity_d[i];
        intensity_sum += all_intensity_d[i];
      }
      double scalar_im = (intensity_sum > 0.0) ? (im_weighted_sum / intensity_sum) : 0.0;

      // Use the first entry's metadata for the precursor info
      const DDAPrecursorInfo* first_entry = entries.front();
      double avg_rt = rt_sum / static_cast<double>(entries.size());

      MSSpectrum spec;
      spec.setRT(avg_rt);
      spec.setMSLevel(2);
      spec.setDriftTime(scalar_im);
      spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      spec.setType(SpectrumSettings::SpectrumType::CENTROID);
      // Native ID: "frame=<F> scan=<S> precursor=<P>".
      //
      // This extends the PSI-MS MS:1002818 "Bruker TDF nativeID format"
      // pattern ("frame=<FRAME_ID> scan=<SCAN_ID>") with a trailing
      // "precursor=<Precursors.Id>" token. The extension is REQUIRED for
      // uniqueness under OpenMS's per-precursor aggregation:
      //
      //   - Unlike pwiz/msconvert (which emits one MS2 spectrum per
      //     Bruker PASEF mobility-scan row, where (frame, scan_begin) is
      //     inherently unique), this reader ports timsrust's strategy of
      //     merging all PasefFrameMsMsInfo entries sharing the same
      //     Precursors.Id (potentially across multiple frames) into ONE
      //     aggregated MS2 spectrum — see the loop above and commit
      //     41ce6cfeb (PR #8975). first_entry->(frame_id, scan_begin) is
      //     just the min-ordered row of that group and is NOT guaranteed
      //     unique across distinct precursors (PASEF can co-isolate
      //     multiple precursors in the same (Frame, ScanNumBegin) mobility
      //     window differing only by IsolationMz; multi-cycle precursors
      //     can also collide on their first entry).
      //
      //   - Strict MS:1002818 validators may flag the trailing token, but
      //     the ecosystem convention (pwiz combined mode uses
      //     "merged=N frame=F scanStart=A scanEnd=B" under the SAME
      //     MS:1002818 accession; our own DIA path emits
      //     "frame=F windowGroup=G scan=S") is that vendor-specific
      //     disambiguators are added while keeping MS:1002818. Comet and
      //     Sage dispatch on nativeID content, not the CV accession, so
      //     compatibility is preserved in practice.
      //
      //   - Bruker Precursors.Id is also duplicated as a typed MetaValue
      //     "bruker_precursor_id" for direct programmatic lookup without
      //     parsing the nativeID string.
      spec.setNativeID("frame=" + String(first_entry->frame_id)
                       + " scan=" + String(first_entry->scan_begin)
                       + " precursor=" + String(prec_id));
      spec.setMetaValue("bruker_precursor_id", static_cast<int>(prec_id));

      // Copy sorted peak data
      spec.reserve(all_mz.size());
      for (size_t idx : sort_idx)
      {
        Peak1D peak;
        peak.setMZ(all_mz[idx]);
        peak.setIntensity(all_intensity_d[idx]);
        spec.push_back(peak);
      }

      // Precursor metadata
      Precursor prec;
      prec.setMZ(first_entry->monoisotopic_mz);  // selected ion m/z (monoisotopic)
      prec.setCharge(first_entry->charge);
      if (!std::isnan(first_entry->intensity))
        prec.setIntensity(static_cast<float>(first_entry->intensity));
      prec.setDriftTime(scalar_im);
      prec.setDriftTimeUnit(DriftTimeUnit::VSSC);

      // Isolation window offsets relative to the selected ion m/z
      double iso_offset_lower = first_entry->isolation_width / 2.0 +
                                (first_entry->monoisotopic_mz - first_entry->isolation_mz);
      double iso_offset_upper = first_entry->isolation_width / 2.0 -
                                (first_entry->monoisotopic_mz - first_entry->isolation_mz);
      prec.setIsolationWindowLowerOffset(std::max(0.0, iso_offset_lower));
      prec.setIsolationWindowUpperOffset(std::max(0.0, iso_offset_upper));

      std::vector<Precursor> precursors;
      precursors.push_back(std::move(prec));
      spec.setPrecursors(std::move(precursors));

      exp.addSpectrum(std::move(spec));
      ++precursor_idx;
      setProgress(ms1_to_emit + precursor_idx);
    }
    endProgress();
    centroider.reportCapSummary(static_cast<size_t>(config.ms1_centroid_max_peaks));
  }

  // =====================================================================
  // loadDIA_: MS1 frames + MS2 frames split by SWATH window
  // =====================================================================
  void BrukerTimsFile::loadDIA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config)
  {
    String tdf_path = handle.get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(std::string(tdf_path), SQLite::OPEN_READONLY);

    // Collect MS1 frame IDs (MS2 frames are mapped via DiaFrameMsMsInfo below)
    std::vector<uint32_t> ms1_frame_ids;

    for (uint32_t fid = handle.min_frame_id(); fid <= handle.max_frame_id(); ++fid)
    {
      if (!frameInRange(fid, config)) continue;
      if (!handle.has_frame(fid)) continue;
      TimsFrame& frame = handle.get_frame(fid);
      if (frame.msms_type == 0)
        ms1_frame_ids.push_back(fid);
    }
    warnPartialMS1AggregationConfig(config, ms1_frame_ids.size());

    // --- MS1 frames ---
    FrameCentroider centroider;
    FrameAggregator ms1_aggregator(0.01);  // finer bin for MS1 (see design spec)
    const size_t ms1_to_emit = config.load_ms1 ? ms1_frame_ids.size() : 0;
    if (config.load_ms1 && config.ms1_n_neighbors > 0) ms1_aggregator.reserve(300'000);
    startProgress(0, ms1_to_emit, "Loading DIA-PASEF MS1 frames");
    for (size_t i = 0; i < ms1_to_emit; ++i)
    {
      MSSpectrum spec;
      if (config.ms1_n_neighbors > 0)
      {
        loadAggregatedMS1Spectrum(handle, ms1_frame_ids, i, config,
                                  ms1_aggregator, centroider, spec);
      }
      else
      {
        TimsFrame& frame = handle.get_frame(ms1_frame_ids[i]);
        loadMS1Spectrum(frame, spec, config, centroider);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
    endProgress();
    if (config.load_ms1)
      centroider.reportCapSummary(static_cast<size_t>(config.ms1_centroid_max_peaks));

    // --- SWATH windows ---
    std::vector<DIAWindow> windows = readDIAWindows(db, *handle.scan2inv_ion_mobility_converter);

    if (windows.empty())
    {
      OPENMS_LOG_WARN << "Warning: DIA dataset detected but no SWATH windows found" << std::endl;
      return;
    }

    // --- MS2 frames (split by SWATH window) ---
    // Use per-WindowGroup iteration for both paths: each diaPASEF frame belongs
    // to exactly one WindowGroup, so we only produce spectra for valid
    // frame/window combinations (not the brute-force frames × all_windows).
    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
    for (auto& [group, frame_ids] : group_to_frames)
    {
      frame_ids.erase(
        std::remove_if(frame_ids.begin(), frame_ids.end(),
                       [&config](uint32_t fid) { return !frameInRange(fid, config); }),
        frame_ids.end());
    }

    // Group DIAWindows by window_group
    std::map<int, std::vector<const DIAWindow*>> group_to_windows;
    for (const auto& w : windows)
    {
      group_to_windows[w.window_group].push_back(&w);
    }

    // Count total work units for progress (frames * windows-per-group)
    Size total_work = 0;
    for (const auto& [group, frames] : group_to_frames)
    {
      auto wit = group_to_windows.find(group);
      if (wit != group_to_windows.end())
        total_work += frames.size() * wit->second.size();
    }

    const bool do_aggregate = config.dia_ms2_n_neighbors > 0;

    if (do_aggregate)
    {
      // === Aggregated path: per-WindowGroup iteration with RT-neighbor summing ===
      startProgress(0, total_work, "Loading DIA-PASEF MS2 frames (aggregated)");
      Size progress_count = 0;

      const int N = config.dia_ms2_n_neighbors;
      FrameAggregator aggregator;

      for (const auto& [group, frame_ids] : group_to_frames)
      {
        auto wit = group_to_windows.find(group);
        if (wit == group_to_windows.end()) continue;
        const auto& dia_windows = wit->second;

        for (const DIAWindow* win : dia_windows)
        {
          for (size_t i = 0; i < frame_ids.size(); ++i)
          {
            setProgress(progress_count++);
            aggregator.clear();
            size_t contributing_frames = 0;

            // Determine neighbor range
            size_t lo = (i >= static_cast<size_t>(N)) ? i - N : 0;
            size_t hi = std::min(i + N, frame_ids.size() - 1);

            // Aggregate peaks from all neighbor frames
            for (size_t ni = lo; ni <= hi; ++ni)
            {
              TimsFrame& nframe = handle.get_frame(frame_ids[ni]);
              if (nframe.num_peaks == 0) continue;

              std::vector<uint32_t> scan_ids(nframe.num_peaks);
              std::vector<uint32_t> intensities(nframe.num_peaks);
              std::vector<double> mzs(nframe.num_peaks);

              nframe.save_to_buffs(nullptr, scan_ids.data(), nullptr, intensities.data(),
                                   mzs.data(), nullptr, nullptr);

              bool frame_contributed = false;
              for (uint32_t p = 0; p < nframe.num_peaks; ++p)
              {
                // Filter by scan bounds (integer comparison, no IM conversion needed)
                if (scan_ids[p] >= win->scan_begin && scan_ids[p] <= win->scan_end)
                {
                  aggregator.addPeak(mzs[p], intensities[p], scan_ids[p]);
                  frame_contributed = true;
                }
              }
              if (frame_contributed) ++contributing_frames;
            }

            // Denoise skipped when the grid was populated from a single frame (or none),
            // because spatial denoising requires neighbor context that only cross-frame
            // aggregation provides. Also skipped when min_support <= 0 (user explicitly
            // disabled the filter).
            bool skip_denoise = (contributing_frames <= 1) || config.dia_ms2_min_support <= 0;
            using CA = BrukerTimsFile::Config::CentroidAlgo;
            const auto ms2_algo = effectiveMS2Algo(config, /*emit_warnings=*/false);
            // For HillBased we keep the unprocessed aggregator output (denoised
            // only) and run the hill helper after; for Greedy2D we delegate to
            // finalizeCentroided which combines smoothing+local-maxima inline.
            auto peaks = (ms2_algo == CA::GREEDY2D)
              ? aggregator.finalizeCentroided(config.dia_ms2_min_support, skip_denoise)
              : aggregator.finalize(config.dia_ms2_min_support, skip_denoise);

            if (peaks.empty()) continue;

            // Build spectrum from center frame metadata
            TimsFrame& center_frame = handle.get_frame(frame_ids[i]);
            MSSpectrum spec;
            spec.setRT(center_frame.time);
            spec.setMSLevel(2);
            spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
            spec.setType(SpectrumSettings::SpectrumType::CENTROID);
            spec.setNativeID("frame=" + String(center_frame.id) + " windowGroup=" + String(win->window_group) + " scan=" + String(win->scan_begin));

            Precursor prec;
            prec.setMZ(win->mz_center);
            prec.setIsolationWindowLowerOffset(win->mz_width / 2.0);
            prec.setIsolationWindowUpperOffset(win->mz_width / 2.0);
            spec.setPrecursors({prec});

            spec.setMetaValue("ion mobility lower limit", win->im_lower);
            spec.setMetaValue("ion mobility upper limit", win->im_upper);

            DataArrays::FloatDataArray im_array;
            IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

            // Convert aggregator output to FrameBuffers once for all algos
            // (HILL/PPIM consume the buffers; GREEDY2D/OFF emit from them).
            // This also enables the isotopic prefilter to run before dispatch.
            FrameBuffers buf;
            aggregatorPeaksToBuffers(peaks, handle, center_frame.id, buf);

            if (config.isotopic_prefilter)
            {
              const size_t before = buf.size();
              const size_t dropped = isotopicPrefilter(buf, config.isotopic_prefilter_tol_ppm);
              OPENMS_LOG_INFO << "Isotopic prefilter (DIA-MS2 aggregated, frame="
                              << center_frame.id << " window=" << win->window_group
                              << "): dropped " << dropped << " / " << before
                              << " peaks ("
                              << (100.0 * dropped / std::max<size_t>(1, before))
                              << "%)" << std::endl;
              if (buf.empty()) continue;
            }

            std::vector<float> b_im_lo, b_im_hi, b_mz_lo, b_mz_hi;
            bool wrote_hill_bounds = false;
            if (ms2_algo == CA::HILL_BASED)
            {
              std::vector<double> out_mz, out_int;
              std::vector<float>  out_im;
              const bool want_bounds = config.expose_hill_bounds;
              HillBoundsOut bounds = want_bounds
                  ? HillBoundsOut{&b_im_lo, &b_im_hi, &b_mz_lo, &b_mz_hi}
                  : HillBoundsOut{nullptr, nullptr, nullptr, nullptr};
              const bool produced = runMS2HillCentroider(
                  config, buf.mz.data(), buf.intensity.data(), buf.im.data(),
                  buf.size(), out_mz, out_int, out_im,
                  buf.scan_ids.data(), bounds);
              if (!produced) continue;
              for (size_t k = 0; k < out_mz.size(); ++k)
              {
                spec.emplace_back(out_mz[k], static_cast<float>(out_int[k]));
                im_array.push_back(out_im[k]);
              }
              spec.setIMPeakType(IMPeakType::IM_CENTROIDED);
              wrote_hill_bounds = want_bounds;
            }
            else
            {
              // GREEDY2D (finalizeCentroided already produced centroids) and OFF
              // (raw passthrough). Emit straight from buf.
              for (size_t k = 0; k < buf.size(); ++k)
              {
                spec.emplace_back(buf.mz[k], static_cast<float>(buf.intensity[k]));
                im_array.push_back(static_cast<float>(buf.im[k]));
              }
              spec.setIMPeakType(ms2_algo == CA::GREEDY2D ? IMPeakType::IM_CENTROIDED
                                                          : IMPeakType::IM_PROFILE);
            }

            spec.getFloatDataArrays().push_back(std::move(im_array));
            if (wrote_hill_bounds)
              attachHillBoundsArrays(spec, b_im_lo, b_im_hi, b_mz_lo, b_mz_hi);
            exp.addSpectrum(std::move(spec));
          }
        }
      }
      endProgress();
    }
    else
    {
      // === Raw path: per-WindowGroup iteration (no RT aggregation) ===
      // HillBased centroiding can still run here — it operates on the single
      // frame's peaks across IM scans, independent of cross-RT summing.
      using CA = BrukerTimsFile::Config::CentroidAlgo;
      const auto ms2_algo_raw = effectiveMS2Algo(config, /*emit_warnings=*/false);
      const bool do_hill_raw = (ms2_algo_raw == CA::HILL_BASED);
      startProgress(0, total_work, do_hill_raw
                    ? "Loading DIA-PASEF MS2 frames (hill, no RT aggregation)"
                    : "Loading DIA-PASEF MS2 frames");
      Size progress_count = 0;

      for (const auto& [group, frame_ids] : group_to_frames)
      {
        auto wit = group_to_windows.find(group);
        if (wit == group_to_windows.end()) continue;
        const auto& dia_windows = wit->second;

        for (const DIAWindow* win : dia_windows)
        {
          for (size_t i = 0; i < frame_ids.size(); ++i)
          {
            setProgress(progress_count++);
            TimsFrame& frame = handle.get_frame(frame_ids[i]);
            if (frame.num_peaks == 0) continue;

            std::vector<uint32_t> raw_scan_ids(frame.num_peaks);
            std::vector<uint32_t> raw_intensities(frame.num_peaks);
            std::vector<double> raw_mzs(frame.num_peaks);
            std::vector<double> raw_ims(frame.num_peaks);

            frame.save_to_buffs(nullptr, raw_scan_ids.data(), nullptr,
                                raw_intensities.data(),
                                raw_mzs.data(), raw_ims.data(), nullptr);

            MSSpectrum spec;
            spec.setRT(frame.time);
            spec.setMSLevel(2);
            spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
            spec.setType(SpectrumSettings::SpectrumType::CENTROID);
            spec.setNativeID("frame=" + String(frame.id) + " windowGroup=" + String(win->window_group) + " scan=" + String(win->scan_begin));

            Precursor prec;
            prec.setMZ(win->mz_center);
            prec.setIsolationWindowLowerOffset(win->mz_width / 2.0);
            prec.setIsolationWindowUpperOffset(win->mz_width / 2.0);
            spec.setPrecursors({prec});

            spec.setMetaValue("ion mobility lower limit", win->im_lower);
            spec.setMetaValue("ion mobility upper limit", win->im_upper);

            DataArrays::FloatDataArray im_array;
            IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

            // Window-filter into FrameBuffers once for all algos. IM is
            // already in 1/K0 from save_to_buffs — no scan_id→IM conversion
            // needed here (unlike the aggregated path). scan_ids carried
            // through for HillBased max_scan_gap and the isotopic prefilter.
            FrameBuffers buf;
            buf.reserve(frame.num_peaks);
            for (uint32_t p = 0; p < frame.num_peaks; ++p)
            {
              if (raw_ims[p] >= win->im_lower && raw_ims[p] <= win->im_upper)
              {
                buf.mz.push_back(raw_mzs[p]);
                buf.intensity.push_back(static_cast<double>(raw_intensities[p]));
                buf.im.push_back(raw_ims[p]);
                buf.scan_ids.push_back(raw_scan_ids[p]);
              }
            }
            if (buf.empty()) continue;

            if (config.isotopic_prefilter)
            {
              const size_t before = buf.size();
              const size_t dropped = isotopicPrefilter(buf, config.isotopic_prefilter_tol_ppm);
              OPENMS_LOG_INFO << "Isotopic prefilter (DIA-MS2 raw, frame="
                              << frame.id << " window=" << win->window_group
                              << "): dropped " << dropped << " / " << before
                              << " peaks ("
                              << (100.0 * dropped / std::max<size_t>(1, before))
                              << "%)" << std::endl;
              if (buf.empty()) continue;
            }

            std::vector<float> b_im_lo, b_im_hi, b_mz_lo, b_mz_hi;
            bool wrote_hill_bounds = false;
            if (do_hill_raw)
            {
              std::vector<double> out_mz, out_int;
              std::vector<float>  out_im;
              const bool want_bounds = config.expose_hill_bounds;
              HillBoundsOut bounds = want_bounds
                  ? HillBoundsOut{&b_im_lo, &b_im_hi, &b_mz_lo, &b_mz_hi}
                  : HillBoundsOut{nullptr, nullptr, nullptr, nullptr};
              const bool produced = runMS2HillCentroider(
                  config, buf.mz.data(), buf.intensity.data(), buf.im.data(),
                  buf.size(), out_mz, out_int, out_im,
                  buf.scan_ids.data(), bounds);
              if (!produced) continue;
              for (size_t k = 0; k < out_mz.size(); ++k)
              {
                spec.emplace_back(out_mz[k], static_cast<float>(out_int[k]));
                im_array.push_back(out_im[k]);
              }
              spec.setIMPeakType(IMPeakType::IM_CENTROIDED);
              wrote_hill_bounds = want_bounds;
            }
            else
            {
              // OFF (or Greedy2D in raw mode — falls through here per
              // effectiveDIAMS2Algo: Greedy2D requires aggregation).
              for (size_t k = 0; k < buf.size(); ++k)
              {
                spec.emplace_back(buf.mz[k], static_cast<float>(buf.intensity[k]));
                im_array.push_back(static_cast<float>(buf.im[k]));
              }
              spec.setIMPeakType(IMPeakType::IM_PROFILE);
            }

            if (!spec.empty())
            {
              spec.getFloatDataArrays().push_back(std::move(im_array));
              if (wrote_hill_bounds)
                attachHillBoundsArrays(spec, b_im_lo, b_im_hi, b_mz_lo, b_mz_hi);
              exp.addSpectrum(std::move(spec));
            }
          }
        }
      }
      endProgress();
    }
  }

  // =====================================================================
  // loadFrames_: raw frames without precursor grouping or SWATH splitting
  // =====================================================================
  void BrukerTimsFile::loadFrames_(TimsDataHandle& handle, MSExperiment& exp, const Config& config)
  {
    // FRAME mode returns raw per-frame spectra; MS1 aggregation would break that
    // contract, so warn and ignore the knob.
    if (config.ms1_n_neighbors > 0)
    {
      OPENMS_LOG_WARN << "Warning: ms1_n_neighbors (=" << config.ms1_n_neighbors
                      << ") is ignored when export_mode=FRAME because FRAME mode returns "
                      << "raw frames. Set export_mode=AUTO or SPECTRUM to enable MS1 "
                      << "frame aggregation." << std::endl;
    }

    // Iterate all frames at each MS level
    for (int level = 1; level <= 2; ++level)
    {
      // Honor load_ms1: skip level==1 entirely when the user opts out of MS1.
      if (level == 1 && !config.load_ms1) continue;

      // Collect frame IDs at this level
      std::vector<uint32_t> frame_ids;
      for (uint32_t fid = handle.min_frame_id(); fid <= handle.max_frame_id(); ++fid)
      {
        if (!frameInRange(fid, config)) continue;
        if (!handle.has_frame(fid)) continue;
        TimsFrame& frame = handle.get_frame(fid);
        if ((level == 1 && frame.msms_type == 0) ||
            (level == 2 && frame.msms_type != 0))
        {
          frame_ids.push_back(fid);
        }
      }

      if (frame_ids.empty())
      {
        if (level == 2)
          OPENMS_LOG_WARN << "Warning: no MS2 frames found (skipping)" << std::endl;
        continue;
      }

      FrameCentroider centroider;
      startProgress(0, frame_ids.size(), String("Loading MS") + String(level) + " frames");
      for (size_t i = 0; i < frame_ids.size(); ++i)
      {
        TimsFrame& frame = handle.get_frame(frame_ids[i]);
        MSSpectrum spec;
        if (level == 1)
          loadMS1Spectrum(frame, spec, config, centroider);
        else
          frameToSpectrum(frame, spec, level);
        exp.addSpectrum(std::move(spec));
        setProgress(i);
      }
      endProgress();
      if (level == 1)
        centroider.reportCapSummary(static_cast<size_t>(config.ms1_centroid_max_peaks));
    }
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
