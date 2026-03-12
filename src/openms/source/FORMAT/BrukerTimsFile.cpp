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

  // Helper: get error string from handle
  static String getTimsError(tims_dataset* handle)
  {
    char buf[512];
    tims_get_last_error(handle, buf, sizeof(buf));
    return String(buf);
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
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
          "timsrust error: " + getTimsError(ds.ptr));
    }
    else
    {
      timsffi_status status = tims_open(path.c_str(), &ds.ptr);
      if (status != TIMSFFI_OK)
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
          "timsrust error: " + getTimsError(ds.ptr));
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

  void BrukerTimsFile::load(const String& path, MSExperiment& exp, const Config& config)
  {
    // Placeholder — will be replaced in Task 8 with full dispatch logic
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::transform(const String& /*path*/, Interfaces::IMSDataConsumer* /*consumer*/, const Config& /*config*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::loadDDA_(void* handle, MSExperiment& exp)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (CONCATENATED format) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));

    // RAII guard for ms1_frames
    auto ms1_guard = [&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); };
    struct Ms1Guard { decltype(ms1_guard)& f; ~Ms1Guard() { f(); } } ms1_g{ms1_guard};

    exp.reserveSpaceSpectra(ms1_count + tims_num_spectra(ds));

    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
    }

    // --- MS2 spectra (one per precursor, scalar IM) ---
    unsigned int num_ms2 = tims_num_spectra(ds);
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
      // Note: isolation_mz is the quadrupole center; precursor_mz is the selected ion.
      // mzML Precursor::setMZ is the isolation target, offsets are relative to it.
      Precursor prec;
      prec.setMZ(ts.isolation_mz);  // isolation window center
      prec.setCharge(static_cast<int>(ts.charge));
      if (!std::isnan(ts.precursor_intensity))
        prec.setIntensity(static_cast<float>(ts.precursor_intensity));
      prec.setDriftTime(ts.im);
      prec.setDriftTimeUnit(DriftTimeUnit::VSSC);
      prec.setIsolationWindowLowerOffset(ts.isolation_width / 2.0);
      prec.setIsolationWindowUpperOffset(ts.isolation_width / 2.0);
      // Store selected ion m/z as user param if different from isolation center
      if (std::abs(ts.precursor_mz - ts.isolation_mz) > 1e-6)
        prec.setMetaValue("selected_ion_mz", ts.precursor_mz);

      std::vector<Precursor> precursors;
      precursors.push_back(std::move(prec));
      spec.setPrecursors(std::move(precursors));

      exp.addSpectrum(std::move(spec));
    }
  }

  void BrukerTimsFile::loadDIA_(void* handle, MSExperiment& exp)
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);

    // --- MS1 frames (same as DDA) ---
    unsigned int ms1_count = 0;
    tims_frame* ms1_frames = nullptr;
    timsffi_status status = tims_get_frames_by_level(ds, 1, &ms1_count, &ms1_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS1 frames: " + getTimsError(ds));
    auto ms1_guard = [&]() { if (ms1_frames) tims_free_frame_array(ds, ms1_frames, ms1_count); };
    struct Ms1Guard { decltype(ms1_guard)& f; ~Ms1Guard() { f(); } } ms1_g{ms1_guard};

    for (unsigned int i = 0; i < ms1_count; ++i)
    {
      MSSpectrum spec;
      frameToSpectrum_(handle, &ms1_frames[i], spec);
      exp.addSpectrum(std::move(spec));
    }

    // --- SWATH windows ---
    unsigned int win_count = 0;
    tims_swath_window* windows = nullptr;
    status = tims_get_swath_windows(ds, &win_count, &windows);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read SWATH windows: " + getTimsError(ds));
    auto win_guard = [&]() { if (windows) tims_free_swath_windows(ds, windows); };
    struct WinGuard { decltype(win_guard)& f; ~WinGuard() { f(); } } win_g{win_guard};

    // --- MS2 frames (split by SWATH window) ---
    unsigned int ms2_count = 0;
    tims_frame* ms2_frames = nullptr;
    status = tims_get_frames_by_level(ds, 2, &ms2_count, &ms2_frames);
    if (status != TIMSFFI_OK)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "", "Failed to read MS2 frames: " + getTimsError(ds));
    auto ms2_guard = [&]() { if (ms2_frames) tims_free_frame_array(ds, ms2_frames, ms2_count); };
    struct Ms2Guard { decltype(ms2_guard)& f; ~Ms2Guard() { f(); } } ms2_g{ms2_guard};

    for (unsigned int fi = 0; fi < ms2_count; ++fi)
    {
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
      std::vector<double> im_values(frame.num_peaks);
      for (uint32_t scan_idx = 0; scan_idx < frame.num_scans; ++scan_idx)
      {
        uint64_t start = frame.scan_offsets[scan_idx];
        uint64_t end = frame.scan_offsets[scan_idx + 1];
        for (uint64_t p = start; p < end; ++p)
        {
          im_values[p] = scan_im[scan_idx];
        }
      }

      // Split peaks by SWATH window IM bounds
      for (unsigned int wi = 0; wi < win_count; ++wi)
      {
        if (windows[wi].is_ms1) continue;

        MSSpectrum spec;
        spec.setRT(frame.rt_seconds);
        spec.setMSLevel(2);
        spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

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
  }

  void BrukerTimsFile::loadFrames_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::frameToSpectrum_(void* handle, const void* frame_ptr, MSSpectrum& spec) const
  {
    tims_dataset* ds = static_cast<tims_dataset*>(handle);
    const tims_frame* frame = static_cast<const tims_frame*>(frame_ptr);

    spec.clear(true);
    spec.setRT(frame->rt_seconds);
    spec.setMSLevel(frame->ms_level);
    spec.setDriftTimeUnit(DriftTimeUnit::VSSC);

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
    std::vector<float> im_values(frame->num_peaks);
    for (uint32_t scan_idx = 0; scan_idx < frame->num_scans; ++scan_idx)
    {
      uint64_t start = frame->scan_offsets[scan_idx];
      uint64_t end = frame->scan_offsets[scan_idx + 1];
      for (uint64_t p = start; p < end; ++p)
      {
        im_values[p] = static_cast<float>(scan_im[scan_idx]);
      }
    }

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
