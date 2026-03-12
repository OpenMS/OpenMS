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
   *
   * @section signal_processing Signal processing
   *
   * In AUTO and SPECTRUM export modes, spectra are produced by timsrust's
   * SpectrumReader which always applies TOF-domain smoothing and centroiding
   * to the raw scan data. This processing is the same for both DDA and DIA
   * acquisitions. The @p smoothing_window and @p centroiding_window parameters
   * control the aggressiveness (bin count); 0 uses the library defaults.
   *
   * In FRAME export mode, raw TOF indices and intensities are returned without
   * any signal processing. TOF-to-m/z and scan-to-IM conversions are applied,
   * but no smoothing or centroiding is performed.
   *
   * @section export_modes Export modes
   *
   * - **AUTO** (default): detects DDA vs DIA by querying for SWATH windows.
   *   DDA uses per-precursor MS2 spectra with scalar drift times; DIA uses
   *   per-frame MS2 spectra split by SWATH isolation window with per-peak IM.
   * - **SPECTRUM**: forces the DDA (per-precursor) path regardless of
   *   acquisition type.
   * - **FRAME**: returns raw 4D frames as CONCATENATED spectra (per-peak IM
   *   in FloatDataArray) for both MS1 and MS2, without precursor grouping or
   *   SWATH splitting. MS2 frame-level access may not be available for all
   *   acquisition types.
   */
  class OPENMS_DLLAPI BrukerTimsFile : public ProgressLogger
  {
  public:
    /// Processing and export configuration
    struct Config
    {
      uint32_t smoothing_window = 0;       ///< Smoothing window in TOF bins (0 = timsrust default). Applied in AUTO/SPECTRUM modes only.
      uint32_t centroiding_window = 0;     ///< Centroiding window in TOF bins (0 = timsrust default). Applied in AUTO/SPECTRUM modes only.
      double calibration_tolerance = 0.0;  ///< m/z recalibration tolerance (0 = timsrust default)
      bool calibrate = false;              ///< Enable m/z recalibration (off by default; may fail on some datasets)

      float ms1_centroid_mz_ppm = 0.0f; ///< MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0). Adapted from Sage (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
      float ms1_centroid_im_pct = 0.0f;  ///< MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA; SPECTRUM forces per-precursor; FRAME returns raw 4D frames
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
    void loadDDA_(void* handle, MSExperiment& exp, const Config& config);

    /// Load DIA-PASEF data: MS1 frames (CONCATENATED) + MS2 frames split by SWATH window
    void loadDIA_(void* handle, MSExperiment& exp, const Config& config);

    /// Load raw frames without precursor grouping or SWATH splitting
    void loadFrames_(void* handle, MSExperiment& exp, const Config& config);

    /// Convert a single tims_frame to an MSSpectrum in CONCATENATED format
    void frameToSpectrum_(void* handle, const void* frame, MSSpectrum& spec) const;

    /// Detect DDA vs DIA by checking for SWATH windows
    bool isDIA_(void* handle) const;
  };

} // namespace OpenMS

#endif // WITH_TIMSRUST
