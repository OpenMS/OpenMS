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

  void BrukerTimsFile::frameToSpectrum_(void* /*handle*/, const void* /*frame*/, MSSpectrum& /*spec*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

} // namespace OpenMS

#endif // WITH_TIMSRUST
