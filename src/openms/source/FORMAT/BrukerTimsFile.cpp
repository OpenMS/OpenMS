// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/CONCEPT/Exception.h>

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

  void BrukerTimsFile::loadDDA_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  void BrukerTimsFile::loadDIA_(void* /*handle*/, MSExperiment& /*exp*/)
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
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
