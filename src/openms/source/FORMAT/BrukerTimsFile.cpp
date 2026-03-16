// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <timsrust_cpp_bridge.h>

#include <memory>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

namespace OpenMS
{

  // RAII wrapper for tims_dataset handle
  struct TimsDatasetHandle
  {
    tims_dataset* ptr = nullptr;

    ~TimsDatasetHandle()
    {
      if (ptr) tims_close(ptr);
    }

    TimsDatasetHandle() = default;
    TimsDatasetHandle(const TimsDatasetHandle&) = delete;
    TimsDatasetHandle& operator=(const TimsDatasetHandle&) = delete;
  };

  // RAII wrapper for tims_config
  struct TimsConfigHandle
  {
    tims_config* ptr = nullptr;

    ~TimsConfigHandle()
    {
      if (ptr) tims_config_free(ptr);
    }

    TimsConfigHandle() = default;
    TimsConfigHandle(const TimsConfigHandle&) = delete;
    TimsConfigHandle& operator=(const TimsConfigHandle&) = delete;
  };

  // Helper: expand run-length-encoded scan offsets to per-peak IM values.
  // Translated from Sage's PeakBuffer::expand_mobility_iter().
  // T = float for frameToSpectrum_() / centroiding path,
  // T = double for loadDIA_() SWATH splitting (preserves precision for window boundary comparisons).
  template <typename T>
  static void expandScanOffsets(const uint64_t* scan_offsets, uint32_t num_scans,
                                const double* scan_im, uint32_t num_peaks,
                                std::vector<T>& out_im)
  {
    out_im.resize(num_peaks);
    for (uint32_t scan_idx = 0; scan_idx < num_scans; ++scan_idx)
    {
      uint64_t start = scan_offsets[scan_idx];
      uint64_t end = scan_offsets[scan_idx + 1];
      for (uint64_t p = start; p < end; ++p)
      {
        out_im[p] = static_cast<T>(scan_im[scan_idx]);
      }
    }
  }

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
    // mz_values: converted m/z (double, from tims_convert_tof_to_mz_array)
    // intensities: raw uint32_t from tims_frame::intensities (cast to float internally)
    // im_values: expanded per-peak IM (float, from expandScanOffsets)
    void loadFrame(const double* mz_values, const uint32_t* intensities,
                   const float* im_values, uint32_t count)
    {
      clear();
      peaks.reserve(count);
      for (uint32_t i = 0; i < count; ++i)
      {
        peaks.push_back({static_cast<float>(mz_values[i]),
                         static_cast<float>(intensities[i]),
                         im_values[i]});
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

  // Helper: get error string from handle
  static String getTimsError(tims_dataset* handle)
  {
    if (!handle) return "unknown timsrust error (null handle)";
    char buf[512];
    tims_get_last_error(handle, buf, sizeof(buf));
    return String(buf);
  }

  // Centroid a single MS1 frame: convert TOF->mz, expand scan offsets, run centroider,
  // build MSSpectrum with centroided peaks and IM FloatDataArray.
  // Note: This duplicates the TOF-to-mz and scan-to-IM conversion from frameToSpectrum_()
  // because centroiding needs the intermediate arrays before they are packed into MSSpectrum.
  static void centroidMS1Frame_(tims_dataset* ds, const tims_frame& frame,
                                MSSpectrum& spec, const BrukerTimsFile::Config& config,
                                FrameCentroider& centroider)
  {
    spec.clear(true);
    spec.setRT(frame.rt_seconds);
    spec.setMSLevel(frame.ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    spec.setType(SpectrumSettings::SpectrumType::CENTROID);
    spec.setNativeID("frame=" + String(frame.index));

    if (frame.num_peaks == 0) return;

    // Batch-convert TOF indices to m/z
    std::vector<double> mz_values(frame.num_peaks);
    timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame.tof_indices, frame.num_peaks, mz_values.data());
    if (mz_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "TOF to m/z conversion failed: " + getTimsError(ds));

    // Batch-convert scan indices to IM
    std::vector<uint32_t> scan_indices(frame.num_scans);
    std::iota(scan_indices.begin(), scan_indices.end(), 0u);
    std::vector<double> scan_im(frame.num_scans);
    timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame.num_scans, scan_im.data());
    if (im_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Scan to IM conversion failed: " + getTimsError(ds));

    // Expand scan offsets to per-peak IM
    std::vector<float> im_values;
    expandScanOffsets(frame.scan_offsets, frame.num_scans, scan_im.data(), frame.num_peaks, im_values);

    // Run centroiding
    centroider.loadFrame(mz_values.data(), frame.intensities, im_values.data(), frame.num_peaks);

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
  }

  // Helper: open dataset with optional config
  static void openDataset(const String& path, const BrukerTimsFile::Config& config,
                          TimsDatasetHandle& ds, TimsConfigHandle& cfg)
  {
    bool has_config = (config.smoothing_window != 0 || config.centroiding_window != 0 ||
                       config.calibration_tolerance != 0.0 || config.calibrate);

    if (has_config)
    {
      cfg.ptr = tims_config_create();
      if (!cfg.ptr)
        throw Exception::BaseException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Internal error", "Failed to create timsrust config");

      if (config.smoothing_window != 0)
        tims_config_set_smoothing_window(cfg.ptr, config.smoothing_window);
      if (config.centroiding_window != 0)
        tims_config_set_centroiding_window(cfg.ptr, config.centroiding_window);
      if (config.calibration_tolerance != 0.0)
        tims_config_set_calibration_tolerance(cfg.ptr, config.calibration_tolerance);
      if (config.calibrate)
        tims_config_set_calibrate(cfg.ptr, true);

      timsffi_status status = tims_open_with_config(path.c_str(), cfg.ptr, &ds.ptr);
      if (status != TIMSFFI_OK)
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          path + " (timsrust: " + getTimsError(ds.ptr) + ")");
    }
    else
    {
      timsffi_status status = tims_open(path.c_str(), &ds.ptr);
      if (status != TIMSFFI_OK)
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          path + " (timsrust: " + getTimsError(ds.ptr) + ")");
    }
  }

  bool BrukerTimsFile::isDIA_(void* handle) const
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);
    unsigned int count = 0;
    tims_swath_window* windows = nullptr;
    timsffi_status status = tims_get_swath_windows(ds, &count, &windows);
    if (status == TIMSFFI_OK && windows)
    {
      tims_free_swath_windows(ds, windows);
    }
    return (status == TIMSFFI_OK && count > 0);
  }

  void BrukerTimsFile::load(const String& path, MSExperiment& exp)
  {
    load(path, exp, Config());
  }

  void BrukerTimsFile::load(const String& path, MSExperiment& exp, const Config& config)
  {
    TimsDatasetHandle ds;
    TimsConfigHandle cfg;
    openDataset(path, config, ds, cfg);

    bool is_dia = isDIA_(ds.ptr);
    Config::ExportMode mode = config.export_mode;

    if (mode == Config::FRAME)
    {
      loadFrames_(ds.ptr, exp, config);
    }
    else if (mode == Config::SPECTRUM || (mode == Config::AUTO && !is_dia))
    {
      loadDDA_(ds.ptr, exp, config);  // DDA path (also used for SPECTRUM mode on DIA data)
    }
    else // AUTO + DIA
    {
      loadDIA_(ds.ptr, exp, config);
    }

    // Populate source file metadata
    SourceFile sf;
    sf.setNameOfFile(File::basename(path));
    sf.setPathToFile(File::path(path));
    sf.setFileType("Bruker TDF");
    sf.setNativeIDType("scan number only nativeID format");
    sf.setNativeIDTypeAccession("MS:1000776");
    exp.getSourceFiles().push_back(sf);

    // Sort by RT, interleaved across MS levels
    exp.sortSpectra(true);
    exp.updateRanges();
  }

  void BrukerTimsFile::transform(const String& path, Interfaces::IMSDataConsumer* consumer)
  {
    transform(path, consumer, Config());
  }

  void BrukerTimsFile::transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config)
  {
    TimsDatasetHandle ds;
    TimsConfigHandle cfg;
    openDataset(path, config, ds, cfg);

    tims_file_info_t info;
    timsffi_status status = tims_file_info(ds.ptr, &info);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read file info: " + getTimsError(ds.ptr));

    bool is_dia = isDIA_(ds.ptr);
    Config::ExportMode mode = config.export_mode;

    // Compute expected size
    size_t expected = 0;
    if (mode == Config::FRAME)
    {
      expected = info.num_frames;
    }
    else if (is_dia && mode != Config::SPECTRUM)
    {
      // DIA: MS1 frames + MS2 frames * windows
      unsigned int win_count = 0;
      tims_swath_window* windows = nullptr;
      timsffi_status win_status = tims_get_swath_windows(ds.ptr, &win_count, &windows);
      if (win_status == TIMSFFI_OK && windows) tims_free_swath_windows(ds.ptr, windows);
      expected = info.ms1.count + info.ms2.count * win_count;
    }
    else
    {
      expected = info.ms1.count + info.num_spectra_ms2;
    }

    consumer->setExpectedSize(expected, 0);
    consumer->setExperimentalSettings(ExperimentalSettings());

    // NOTE: This loads into a temporary experiment then feeds to consumer.
    // Not truly constant-memory — a future optimization should iterate
    // frame-by-frame and call consumer->consumeSpectrum() inline.
    OPENMS_LOG_INFO << "BrukerTimsFile::transform(): loading full dataset (streaming optimization pending)" << std::endl;
    MSExperiment exp;
    if (mode == Config::FRAME)
      loadFrames_(ds.ptr, exp, config);
    else if (is_dia && mode != Config::SPECTRUM)
      loadDIA_(ds.ptr, exp, config);
    else
      loadDDA_(ds.ptr, exp, config);

    exp.sortSpectra(true);
    for (auto& spec : exp)
    {
      consumer->consumeSpectrum(spec);
    }
  }

  void BrukerTimsFile::loadDDA_(void* handle, MSExperiment& exp, const Config& config)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (CONCATENATED format) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));

    RAIICleanup ms1_guard([&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); });

    unsigned int num_ms2 = tims_num_spectra(ds);
    exp.reserveSpaceSpectra(ms1_count + num_ms2);

    startProgress(0, ms1_count + num_ms2, "Loading DDA-PASEF data");

    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame_(ds, ms1_frames[i], spec, config, centroider);
      }
      else
      {
        frameToSpectrum_(handle, &ms1_frames[i], spec);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }

    // --- MS2 spectra (one per precursor, scalar IM) ---
    for (unsigned int i = 0; i < num_ms2; ++i)
    {
      tims_spectrum ts;
      status = tims_get_spectrum(ds, i, &ts);
      if (status != TIMSFFI_OK)
      {
        OPENMS_LOG_WARN << "Warning: failed to read MS2 spectrum " << i << ": " << getTimsError(ds) << std::endl;
        continue;
      }

      MSSpectrum spec;
      spec.setRT(ts.rt_seconds);
      spec.setMSLevel(2);
      spec.setDriftTime(ts.im);
      spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      spec.setNativeID("scan=" + String(ts.index));

      // Copy peak data (float -> double widening)
      spec.reserve(ts.num_peaks);
      for (uint32_t p = 0; p < ts.num_peaks; ++p)
      {
        Peak1D peak;
        peak.setMZ(static_cast<double>(ts.mz[p]));
        peak.setIntensity(static_cast<double>(ts.intensity[p]));
        spec.push_back(peak);
      }

      // Precursor metadata
      // Use selected-ion m/z (monoisotopic precursor) for Precursor::getMZ() so
      // search engines use the correct mass for candidate lookup. The quadrupole
      // isolation-window center goes into the isolation window offsets (relative
      // to the selected ion m/z).
      Precursor prec;
      prec.setMZ(ts.precursor_mz);  // selected ion m/z (monoisotopic)
      prec.setCharge(static_cast<int>(ts.charge));
      if (!std::isnan(ts.precursor_intensity))
        prec.setIntensity(static_cast<float>(ts.precursor_intensity));
      prec.setDriftTime(ts.im);
      prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      // Isolation window offsets relative to the selected ion m/z
      double iso_offset_lower = ts.isolation_width / 2.0 + (ts.precursor_mz - ts.isolation_mz);
      double iso_offset_upper = ts.isolation_width / 2.0 - (ts.precursor_mz - ts.isolation_mz);
      prec.setIsolationWindowLowerOffset(std::max(0.0, iso_offset_lower));
      prec.setIsolationWindowUpperOffset(std::max(0.0, iso_offset_upper));

      std::vector<Precursor> precursors;
      precursors.push_back(std::move(prec));
      spec.setPrecursors(std::move(precursors));

      exp.addSpectrum(std::move(spec));
      setProgress(ms1_count + i);
    }
    endProgress();
  }

  void BrukerTimsFile::loadDIA_(void* handle, MSExperiment& exp, const Config& config)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (same as DDA) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));
    RAIICleanup ms1_guard([&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); });

    bool do_centroid = isCentroidingEnabled(config);
    FrameCentroider centroider;
    startProgress(0, ms1_count, "Loading DIA-PASEF MS1 frames");
    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      if (do_centroid)
      {
        centroidMS1Frame_(ds, ms1_frames[i], spec, config, centroider);
      }
      else
      {
        frameToSpectrum_(handle, &ms1_frames[i], spec);
      }
      exp.addSpectrum(std::move(spec));
      setProgress(i);
    }
    endProgress();

    // --- SWATH windows ---
    unsigned int win_count = 0;
    tims_swath_window* windows = nullptr;
    status = tims_get_swath_windows(ds, &win_count, &windows);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read SWATH windows: " + getTimsError(ds));
    RAIICleanup win_guard([&]() { if (windows) tims_free_swath_windows(ds, windows); });

    // --- MS2 frames (split by SWATH window) ---
    unsigned int ms2_count = 0;
    tims_frame* ms2_frames = nullptr;
    status = tims_get_frames_by_level(ds, 2, &ms2_count, &ms2_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS2 frames: " + getTimsError(ds));
    RAIICleanup ms2_guard([&]() { if (ms2_frames) tims_free_frame_array(ds, ms2_frames, ms2_count); });

    startProgress(0, ms2_count, "Loading DIA-PASEF MS2 frames");
    for (unsigned int fi = 0; fi < ms2_count; ++fi)
    {
      setProgress(fi);
      const tims_frame& frame = ms2_frames[fi];
      if (frame.num_peaks == 0) continue;

      // Convert TOF -> m/z for entire frame
      std::vector<double> mz_values(frame.num_peaks);
      timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame.tof_indices, frame.num_peaks, mz_values.data());
      if (mz_status != TIMSFFI_OK)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "TOF to m/z conversion failed: " + getTimsError(ds));

      // Batch scan -> IM conversion
      std::vector<uint32_t> scan_indices(frame.num_scans);
      std::iota(scan_indices.begin(), scan_indices.end(), 0u);
      std::vector<double> scan_im(frame.num_scans);
      timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame.num_scans, scan_im.data());
      if (im_status != TIMSFFI_OK)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "", "Scan to IM conversion failed: " + getTimsError(ds));

      // Build per-peak IM from scan offsets
      std::vector<double> im_values;
      expandScanOffsets<double>(frame.scan_offsets, frame.num_scans, scan_im.data(), frame.num_peaks, im_values);

      // Split peaks by SWATH window IM bounds
      for (unsigned int wi = 0; wi < win_count; ++wi)
      {
        if (windows[wi].is_ms1) continue;

        MSSpectrum spec;
        spec.setRT(frame.rt_seconds);
        spec.setMSLevel(2);
        spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
        spec.setNativeID("frame=" + String(frame.index) + " windowGroup=" + String(wi));

        // Set isolation window
        Precursor prec;
        prec.setMZ(windows[wi].mz_center);
        prec.setIsolationWindowLowerOffset(windows[wi].mz_center - windows[wi].mz_lower);
        prec.setIsolationWindowUpperOffset(windows[wi].mz_upper - windows[wi].mz_center);
        spec.setPrecursors({prec});

        DataArrays::FloatDataArray im_array;
        IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);

        for (uint32_t p = 0; p < frame.num_peaks; ++p)
        {
          // Check if peak falls within this window's IM bounds
          if (im_values[p] >= windows[wi].im_lower && im_values[p] <= windows[wi].im_upper)
          {
            Peak1D peak;
            peak.setMZ(mz_values[p]);
            peak.setIntensity(static_cast<double>(frame.intensities[p]));
            spec.push_back(peak);
            im_array.push_back(static_cast<float>(im_values[p]));
          }
        }

        if (!spec.empty())
        {
          spec.getFloatDataArrays().push_back(std::move(im_array));
          exp.addSpectrum(std::move(spec));
        }
      }
    }
    endProgress();
  }

  void BrukerTimsFile::loadFrames_(void* handle, MSExperiment& exp, const Config& config)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // Load all frames for each MS level
    for (uint8_t level = 1; level <= 2; ++level)
    {
      unsigned int count = 0;
      tims_frame* frames = nullptr;
      timsffi_status status = tims_get_frames_by_level(ds, level, &count, &frames);
      if (status != TIMSFFI_OK)
      {
        // MS2 frame-level access may not be available (e.g., DDA datasets)
        OPENMS_LOG_WARN << "Warning: failed to read frames at MS level " << (int)level
                        << ": " << getTimsError(ds) << " (skipping)" << std::endl;
        continue;
      }
      RAIICleanup frame_guard([&]() { if (frames) tims_free_frame_array(ds, frames, count); });

      bool do_centroid = (level == 1) && isCentroidingEnabled(config);
      FrameCentroider centroider;
      startProgress(0, count, String("Loading MS") + String(level) + " frames");
      for (unsigned int i = 0; i < count; ++i)
      {
        MSSpectrum spec;
        if (do_centroid)
        {
          centroidMS1Frame_(ds, frames[i], spec, config, centroider);
        }
        else
        {
          frameToSpectrum_(handle, &frames[i], spec);
        }
        exp.addSpectrum(std::move(spec));
        setProgress(i);
      }
      endProgress();
    }
  }

  void BrukerTimsFile::frameToSpectrum_(void* handle, const void* frame_ptr, MSSpectrum& spec) const
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);
    const tims_frame* frame = static_cast<const tims_frame*>(frame_ptr);

    spec.clear(true);
    spec.setRT(frame->rt_seconds);
    spec.setMSLevel(frame->ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
    spec.setNativeID("frame=" + String(frame->index));

    if (frame->num_peaks == 0) return;

    // Batch-convert TOF indices to m/z (check return value!)
    std::vector<double> mz_values(frame->num_peaks);
    timsffi_status mz_status = tims_convert_tof_to_mz_array(ds, frame->tof_indices, frame->num_peaks, mz_values.data());
    if (mz_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "TOF to m/z conversion failed: " + getTimsError(ds));

    // Batch-convert scan indices to IM (use batch API for efficiency)
    std::vector<uint32_t> scan_indices(frame->num_scans);
    std::iota(scan_indices.begin(), scan_indices.end(), 0u);
    std::vector<double> scan_im(frame->num_scans);
    timsffi_status im_status = tims_convert_scan_to_im_array(ds, scan_indices.data(), frame->num_scans, scan_im.data());
    if (im_status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Scan to IM conversion failed: " + getTimsError(ds));

    // Build per-peak IM values from scan offsets
    std::vector<float> im_values;
    expandScanOffsets<float>(frame->scan_offsets, frame->num_scans, scan_im.data(), frame->num_peaks, im_values);

    // Fill spectrum peaks
    spec.reserve(frame->num_peaks);
    for (uint32_t i = 0; i < frame->num_peaks; ++i)
    {
      Peak1D peak;
      peak.setMZ(mz_values[i]);
      peak.setIntensity(static_cast<double>(frame->intensities[i]));
      spec.push_back(peak);
    }

    // Set IM float data array with correct CV term
    DataArrays::FloatDataArray im_array;
    im_array.resize(frame->num_peaks);
    std::copy(im_values.begin(), im_values.end(), im_array.begin());
    IMDataConverter::setIMUnit(im_array, DriftTimeUnit::VSSC);
    spec.getFloatDataArrays().push_back(std::move(im_array));
  }

} // namespace OpenMS

#endif // WITH_TIMSRUST
