// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_TIMSRUST

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <cstdint>

namespace OpenMS
{

  /**
   * @brief Reader for Bruker TimsTOF .d directories via timsrust_cpp_bridge.
   *
   * Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
   * Ion mobility data is stored in VSSC (1/K0) units using CONCATENATED format
   * for MS1 and DIA MS2, and scalar drift times for DDA MS2.
   */
  class OPENMS_DLLAPI BrukerTimsFile : public ProgressLogger
  {
  public:
    /// Processing and export configuration
    struct Config
    {
      uint32_t smoothing_window = 0;       ///< 0 = timsrust default; bin count
      uint32_t centroiding_window = 0;     ///< 0 = timsrust default; bin count
      double calibration_tolerance = 0.0;  ///< 0 = timsrust default; m/z tolerance
      bool calibrate = false;              ///< off by default (timsrust 0.4.2 caveat)

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA
    };

    /// Load entire .d directory into MSExperiment
    void load(const String& path, MSExperiment& exp);
    /// @overload with explicit configuration
    void load(const String& path, MSExperiment& exp, const Config& config);

    /// Streaming: read .d and feed spectra to consumer
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer);
    /// @overload with explicit configuration
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config);

  private:
    /// Load DDA-PASEF data: MS1 frames (CONCATENATED) + MS2 spectra (scalar IM)
    void loadDDA_(void* handle, MSExperiment& exp);

    /// Load DIA-PASEF data: MS1 frames (CONCATENATED) + MS2 frames split by SWATH window
    void loadDIA_(void* handle, MSExperiment& exp);

    /// Load raw frames without precursor grouping or SWATH splitting
    void loadFrames_(void* handle, MSExperiment& exp);

    /// Convert a single tims_frame to an MSSpectrum in CONCATENATED format
    void frameToSpectrum_(void* handle, const void* frame, MSSpectrum& spec) const;

    /// Detect DDA vs DIA by checking for SWATH windows
    bool isDIA_(void* handle) const;
  };

} // namespace OpenMS

#endif // WITH_TIMSRUST
