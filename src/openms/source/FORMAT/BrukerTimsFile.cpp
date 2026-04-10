// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/DATAACCESS/SwathFileConsumer.h>
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

  static constexpr size_t MAX_CENTROID_PEAKS = 10000;

  // Reusable buffer for IM-dimension centroiding of a single frame.
  // Translated from Sage's PeakBuffer (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
  // Buffers persist across frames within a single load() call for memory reuse.
  struct FrameCentroider
  {
    std::vector<ImsPeak> peaks;
    std::vector<size_t> order;       // indices sorted by descending intensity
    std::vector<ImsPeak> agg_buff;   // centroided output

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

    // Centroid the loaded frame by collapsing the IM dimension.
    // Iterates peaks in descending intensity order. For each apex peak, aggregates
    // all unconsumed neighbors within m/z (ppm) and IM (percent) tolerances.
    // Output is capped at MAX_CENTROID_PEAKS entries.
    void centroid(float mz_ppm, float im_pct,
                  std::vector<double>& out_mz,
                  std::vector<double>& out_intensity,
                  std::vector<float>& out_im)
    {
      agg_buff.clear();
      agg_buff.reserve(std::min(peaks.size(), MAX_CENTROID_PEAKS + 1));

      const float utol = mz_ppm / 1e6f;
      const float im_tol_frac = im_pct / 100.0f;
      size_t global_consumed = 0;

      for (size_t idx : order)
      {
        if (peaks[idx].intensity <= 0.0f) continue;  // already consumed

        if (agg_buff.size() > MAX_CENTROID_PEAKS)
        {
          if (peaks[idx].intensity > 200.0f)
          {
            OPENMS_LOG_DEBUG << "FrameCentroider: reached MAX_CENTROID_PEAKS at index "
                             << idx << "/" << peaks.size()
                             << " intensity=" << peaks[idx].intensity << std::endl;
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
  };

  // Helper: check if MS1 centroiding is enabled, warn on partial config
  static bool isCentroidingEnabled(const BrukerTimsFile::Config& config)
  {
    bool has_mz = config.ms1_centroid_mz_ppm > 0.0f;
    bool has_im = config.ms1_centroid_im_pct > 0.0f;
    if (has_mz != has_im)
    {
      OPENMS_LOG_WARN << "Warning: ms1_centroid_mz_ppm and ms1_centroid_im_pct must both be > 0 "
                      << "to enable MS1 centroiding. Centroiding disabled." << std::endl;
      return false;
    }
    return has_mz && has_im;
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

    /// Helper for DIA MS2 frame aggregation and denoising.
    /// Bins peaks from multiple frames onto a sparse (mz_bin, scan_id) grid,
    /// applies spatial denoising, and outputs the surviving peaks.
    class DIAFrameAggregator
    {
    public:
      static constexpr double MZ_BIN_WIDTH = 0.02; // Da — absorbs frame-to-frame m/z jitter

      struct OutputPeak
      {
        double mz;        // intensity-weighted mean m/z
        double intensity;  // summed intensity
        uint32_t scan_id;  // native scan index (for IM conversion)
      };

      /// Add a peak to the grid. Call for every peak from every neighbor frame.
      void addPeak(double mz, uint32_t intensity, uint32_t scan_id)
      {
        int64_t mz_bin = static_cast<int64_t>(std::round(mz / MZ_BIN_WIDTH));
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
                uint64_t nkey = (static_cast<uint64_t>(static_cast<uint32_t>(mz_bin + dm)) << 32)
                              | (scan_id + ds);
                if (grid_.count(nkey)) ++neighbors;
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
        std::unordered_map<uint64_t, Cell> denoised;
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
                uint64_t nkey = (static_cast<uint64_t>(static_cast<uint32_t>(mz_bin + dm)) << 32)
                              | (scan_id + ds);
                if (grid_.count(nkey)) ++neighbors;
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

        std::unordered_map<uint64_t, double> smoothed;
        for (const auto& [key, cell] : denoised)
        {
          uint32_t scan_id = static_cast<uint32_t>(key & 0xFFFFFFFF);
          uint32_t mz_bin = static_cast<uint32_t>(key >> 32);
          double weighted_sum = 0.0;
          for (int dm = -2; dm <= 2; ++dm)
          {
            for (int ds = -2; ds <= 2; ++ds)
            {
              uint64_t nkey = (static_cast<uint64_t>(static_cast<uint32_t>(mz_bin + dm)) << 32)
                            | (scan_id + ds);
              auto it = denoised.find(nkey);
              if (it != denoised.end())
              {
                weighted_sum += it->second.intensity_sum * kernel[dm + 2][ds + 2];
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
              uint64_t nkey = (static_cast<uint64_t>(static_cast<uint32_t>(mz_bin + dm)) << 32)
                            | (scan_id + ds);
              auto it = smoothed.find(nkey);
              if (it != smoothed.end() && it->second > val) is_max = false;
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
              uint64_t nkey = (static_cast<uint64_t>(static_cast<uint32_t>(center_mz_bin + dm)) << 32)
                            | (center_scan + ds);
              auto it = denoised.find(nkey);
              if (it != denoised.end())
              {
                double int_val = it->second.intensity_sum;
                double mz_val = it->second.mz_weighted_sum / it->second.intensity_sum;
                total_intensity += int_val;
                mz_weighted += mz_val * int_val;
                scan_weighted += static_cast<double>(static_cast<uint32_t>(nkey & 0xFFFFFFFF)) * int_val;
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

    private:
      struct Cell
      {
        double intensity_sum = 0.0;
        double mz_weighted_sum = 0.0;
      };
      std::unordered_map<uint64_t, Cell> grid_;
    };

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

    FrameCounts readFrameCounts(SQLite::Database& db)
    {
      FrameCounts counts;

      SQLite::Statement q_total(db, "SELECT COUNT(*) FROM Frames");
      if (q_total.executeStep())
        counts.total = static_cast<uint32_t>(q_total.getColumn(0).getInt());

      SQLite::Statement q_ms1(db, "SELECT COUNT(*) FROM Frames WHERE MsMsType = 0");
      if (q_ms1.executeStep())
        counts.ms1 = static_cast<uint32_t>(q_ms1.getColumn(0).getInt());

      // MS2 = everything that is not MS1 (MsMsType != 0)
      SQLite::Statement q_ms2(db, "SELECT COUNT(*) FROM Frames WHERE MsMsType != 0");
      if (q_ms2.executeStep())
        counts.ms2 = static_cast<uint32_t>(q_ms2.getColumn(0).getInt());

      return counts;
    }

    // Count distinct precursors (for DDA MS2 spectrum count)
    uint32_t countDDAPrecursors(SQLite::Database& db)
    {
      try
      {
        SQLite::Statement q(db, "SELECT COUNT(*) FROM Precursors");
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

    // Extract data from frame using save_to_buffs
    std::vector<uint32_t> intensities(frame.num_peaks);
    std::vector<double> mzs(frame.num_peaks);
    std::vector<double> inv_ion_mobilities(frame.num_peaks);

    frame.save_to_buffs(nullptr, nullptr, nullptr, intensities.data(),
                        mzs.data(), inv_ion_mobilities.data(), nullptr);

    // Run centroiding
    centroider.loadFrame(mzs.data(), intensities.data(), inv_ion_mobilities.data(), frame.num_peaks);

    std::vector<double> cent_mz, cent_intensity;
    std::vector<float> cent_im;
    centroider.centroid(config.ms1_centroid_mz_ppm, config.ms1_centroid_im_pct,
                        cent_mz, cent_intensity, cent_im);

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
    spec.setIMPeakType(IMPeakType::IM_CENTROIDED);
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

    // Count MS1 frames
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ++meta.nr_ms1_spectra;
    }

    // Count MS2 spectra per window (= frames per WindowGroup)
    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);
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

    // --- MS1 frames ---
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;

    std::vector<uint32_t> ms1_frame_ids;
    for (uint32_t fid = handle->min_frame_id(); fid <= handle->max_frame_id(); ++fid)
    {
      if (handle->has_frame(fid) && handle->get_frame(fid).msms_type == 0)
        ms1_frame_ids.push_back(fid);
    }

    startProgress(0, ms1_frame_ids.size(), "Streaming DIA-PASEF MS1 frames");
    for (size_t i = 0; i < ms1_frame_ids.size(); ++i)
    {
      TimsFrame& frame = handle->get_frame(ms1_frame_ids[i]);
      MSSpectrum spec;
      if (do_centroid)
        centroidMS1Frame(frame, spec, config, centroider);
      else
        frameToSpectrum(frame, spec, 1);
      consumer.consumeSpectrum(spec);
      setProgress(i);
    }
    endProgress();

    // --- MS2 frames: raw per-WindowGroup iteration (no aggregation) ---
    auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);
    if (windows.empty())
    {
      OPENMS_LOG_WARN << "Warning: DIA dataset detected but no SWATH windows found" << std::endl;
      return;
    }

    auto group_to_frames = readFrameToWindowGroupMapping(db, windows);

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
    bool is_dia = isDIA_(tdf_path);
    Config::ExportMode mode = config.export_mode;

    if (mode == Config::FRAME)
    {
      loadFrames_(*handle, exp, config);
    }
    else if (mode == Config::SPECTRUM || (mode == Config::AUTO && !is_dia))
    {
      loadDDA_(*handle, exp, config);
    }
    else // AUTO + DIA
    {
      loadDIA_(*handle, exp, config);
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

    bool is_dia = isDIA(db);
    Config::ExportMode mode = config.export_mode;

    // Compute expected size for consumer
    FrameCounts counts = readFrameCounts(db);
    size_t expected = 0;
    if (mode == Config::FRAME)
    {
      expected = counts.total;
    }
    else if (is_dia && mode != Config::SPECTRUM)
    {
      // DIA: MS1 frames + MS2 frames * windows
      auto windows = readDIAWindows(db, *handle->scan2inv_ion_mobility_converter);
      expected = counts.ms1 + static_cast<size_t>(counts.ms2) * windows.size();
    }
    else
    {
      // DDA: MS1 frames + per-precursor MS2 spectra
      uint32_t num_precursors = countDDAPrecursors(db);
      expected = counts.ms1 + num_precursors;
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
      loadFrames_(*handle, exp, config);
    else if (is_dia && mode != Config::SPECTRUM)
      loadDIA_(*handle, exp, config);
    else
      loadDDA_(*handle, exp, config);

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
      if (!handle.has_frame(fid)) continue;
      TimsFrame& frame = handle.get_frame(fid);
      if (frame.msms_type == 0)
        ms1_frame_ids.push_back(fid);
      else
        ms2_frame_ids.push_back(fid);
    }

    // Read DDA precursors from SQL
    std::vector<DDAPrecursorInfo> precursor_entries = readDDAPrecursors(db);

    // Group by precursor ID
    std::map<int, std::vector<const DDAPrecursorInfo*>> precursor_groups;
    for (const auto& entry : precursor_entries)
    {
      precursor_groups[entry.precursor_id].push_back(&entry);
    }

    uint32_t num_ms2 = static_cast<uint32_t>(precursor_groups.size());
    exp.reserveSpaceSpectra(static_cast<unsigned int>(ms1_frame_ids.size()) + num_ms2);

    startProgress(0, ms1_frame_ids.size() + num_ms2, "Loading DDA-PASEF data");

    // --- MS1 frames ---
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    for (size_t i = 0; i < ms1_frame_ids.size(); ++i)
    {
      TimsFrame& frame = handle.get_frame(ms1_frame_ids[i]);
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame(frame, spec, config, centroider);
      }
      else
      {
        frameToSpectrum(frame, spec, 1);
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
        setProgress(ms1_frame_ids.size() + precursor_idx);
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
        setProgress(ms1_frame_ids.size() + precursor_idx);
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
      spec.setNativeID("scan=" + String(prec_id));

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
      setProgress(ms1_frame_ids.size() + precursor_idx);
    }
    endProgress();
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
      if (!handle.has_frame(fid)) continue;
      TimsFrame& frame = handle.get_frame(fid);
      if (frame.msms_type == 0)
        ms1_frame_ids.push_back(fid);
    }

    // --- MS1 frames ---
    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    startProgress(0, ms1_frame_ids.size(), "Loading DIA-PASEF MS1 frames");
    for (size_t i = 0; i < ms1_frame_ids.size(); ++i)
    {
      TimsFrame& frame = handle.get_frame(ms1_frame_ids[i]);
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame(frame, spec, config, centroider);
      }
      else
      {
        frameToSpectrum(frame, spec, 1);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
    endProgress();

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
      DIAFrameAggregator aggregator;

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

              for (uint32_t p = 0; p < nframe.num_peaks; ++p)
              {
                // Filter by scan bounds (integer comparison, no IM conversion needed)
                if (scan_ids[p] >= win->scan_begin && scan_ids[p] <= win->scan_end)
                {
                  aggregator.addPeak(mzs[p], intensities[p], scan_ids[p]);
                }
              }
            }

            // Denoise (skip if only 1 frame in range)
            // TODO: This checks the index range, not the actual number of neighbor
            // frames that contributed peaks to the grid. If neighbor frames exist
            // but have zero peaks passing the IM filter for this window, the grid
            // contains only single-frame data yet denoising still runs — which may
            // remove valid isolated peaks. Consider counting actual contributing
            // frames if this becomes an issue in practice.
            bool skip_denoise = (hi - lo) < 1;
            auto peaks = config.dia_ms2_centroid
              ? aggregator.finalizeCentroided(config.dia_ms2_min_support, skip_denoise)
              : aggregator.finalize(config.dia_ms2_min_support, skip_denoise);

            if (peaks.empty()) continue;

            // Build spectrum from center frame metadata
            TimsFrame& center_frame = handle.get_frame(frame_ids[i]);
            MSSpectrum spec;
            spec.setRT(center_frame.time);
            spec.setMSLevel(2);
            spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
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

            for (const auto& peak : peaks)
            {
              spec.emplace_back(peak.mz, peak.intensity);
              // Convert scan_id → 1/K0 using center frame's calibration
              double im_val = 0.0;
              handle.scan2inv_ion_mobility_converter->convert(
                center_frame.id, &im_val, &peak.scan_id, 1);
              im_array.push_back(static_cast<float>(im_val));
            }

            spec.getFloatDataArrays().push_back(std::move(im_array));
            spec.setIMPeakType(config.dia_ms2_centroid ? IMPeakType::IM_CENTROIDED : IMPeakType::IM_PROFILE);
            exp.addSpectrum(std::move(spec));
          }
        }
      }
      endProgress();
    }
    else
    {
      // === Raw path: per-WindowGroup iteration (no aggregation) ===
      startProgress(0, total_work, "Loading DIA-PASEF MS2 frames");
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
    // Iterate all frames at each MS level
    for (int level = 1; level <= 2; ++level)
    {
      // Collect frame IDs at this level
      std::vector<uint32_t> frame_ids;
      for (uint32_t fid = handle.min_frame_id(); fid <= handle.max_frame_id(); ++fid)
      {
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

      bool do_centroid = (level == 1) && isCentroidingEnabled(config);
      FrameCentroider centroider;
      startProgress(0, frame_ids.size(), String("Loading MS") + String(level) + " frames");
      for (size_t i = 0; i < frame_ids.size(); ++i)
      {
        TimsFrame& frame = handle.get_frame(frame_ids[i]);
        MSSpectrum spec;
        if (do_centroid)
        {
          centroidMS1Frame(frame, spec, config, centroider);
        }
        else
        {
          frameToSpectrum(frame, spec, level);
        }
        exp.addSpectrum(std::move(spec));
        setProgress(i);
      }
      endProgress();
    }
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
