// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <cstdint>
#include <string>

class TimsDataHandle;

namespace OpenMS
{
  class FullSwathFileConsumer;

  /**
   * @brief Reader for Bruker TimsTOF .d directories via opentims.
   *
   * Supports DDA-PASEF, DIA-PASEF, and raw frame-level 4D access.
   * Ion mobility data is stored in VSSC (1/K0) units using IM_PEAK format
   * for MS1 and DIA MS2, and scalar drift times for DDA MS2.
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
   * - **FRAME**: returns raw 4D frames as IM_PEAK spectra (per-peak IM
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
      double calibration_tolerance = 0.0;  ///< m/z recalibration tolerance in Da (0 = default 0.1 Da)
      bool calibrate = false;              ///< Enable m/z recalibration (off by default; may fail on some datasets)

      float ms1_centroid_mz_ppm = 0.0f; ///< MS1 IM-centroiding m/z tolerance in ppm (0 = disabled, suggested: 5.0). Adapted from Sage (Lazear 2023, doi:10.1021/acs.jproteome.3c00486).
      float ms1_centroid_im_pct = 0.0f;  ///< MS1 IM-centroiding ion mobility tolerance in percent (0 = disabled, suggested: 3.0)

      int dia_ms2_n_neighbors = 0;  ///< DIA MS2 frame aggregation: number of adjacent frames on each side (0 = disabled, 1 = 3-frame sum, 2 = 5-frame sum)
      int dia_ms2_min_support = 1;  ///< DIA MS2 denoising: minimum occupied neighbors in 3x3 (mz x IM) grid to retain a point (center excluded)
      bool dia_ms2_centroid = false; ///< DIA MS2 2D peak picking: apply Gaussian smoothing + local maxima detection to produce IM_CENTROIDED spectra

      enum ExportMode { AUTO, SPECTRUM, FRAME };
      ExportMode export_mode = AUTO;       ///< AUTO detects DDA vs DIA; SPECTRUM forces per-precursor; FRAME returns raw 4D frames

      /// Strategy for converting TIMS scan indices to 1/K0 values.
      /// AUTO (default): tries Bruker SDK → rational (TimsCalibration table) → linear.
      enum class TimsCalibrationStrategy { AUTO, BRUKER_SDK, RATIONAL, LINEAR };
      TimsCalibrationStrategy tims_calibration_strategy = TimsCalibrationStrategy::AUTO;

      /// Pressure compensation strategy (only effective with BRUKER_SDK calibration).
      /// Corrects for ambient gas pressure drift during acquisition.
      /// Ignored when using RATIONAL or LINEAR strategies; a warning is logged.
      enum class PressureCompensation { NONE, GLOBAL, PER_FRAME };
      PressureCompensation pressure_compensation = PressureCompensation::NONE;

      /// Path to Bruker SDK library (timsdata.dll / libtimsdata.so).
      /// Empty string (default): discover from OPENMS_BRUKER_SDK_PATH env var.
      std::string bruker_sdk_path;
    };

    /// Metadata for constructing a FullSwathFileConsumer for DIA streaming.
    struct DIAStreamingMetadata
    {
      std::vector<OpenSwath::SwathMap> boundaries;  ///< one SwathMap per DIAWindow (MS2 only)
      int nr_ms1_spectra = 0;                       ///< number of MS1 frames
      std::vector<int> nr_ms2_spectra;              ///< per-window spectrum counts (parallel to boundaries)
    };

    /// Read DIA SWATH boundaries and spectrum counts from a .d directory (SQL only, no peak data).
    /// Also populates exp_settings with source file metadata.
    DIAStreamingMetadata readDIAMetadata(const String& path, ExperimentalSettings& exp_settings);
    /// @overload with explicit configuration
    DIAStreamingMetadata readDIAMetadata(const String& path, ExperimentalSettings& exp_settings,
                                         const Config& config);

    /// Stream DIA spectra to a consumer one-at-a-time without accumulating.
    /// MS2 spectra are always raw (no aggregation/denoising/centroiding).
    /// MS1 centroiding respects the Config settings.
    void loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer);
    /// @overload with explicit configuration
    void loadDIAStreaming(const String& path, FullSwathFileConsumer& consumer,
                          const Config& config);

    /// Load entire .d directory into MSExperiment
    void load(const String& path, MSExperiment& exp);
    /// @overload with explicit configuration
    void load(const String& path, MSExperiment& exp, const Config& config);

    /// Feed spectra from a .d directory to a consumer.
    /// @note Currently loads the full dataset into memory before streaming to the consumer.
    ///       A future optimization should iterate frame-by-frame for true constant-memory operation.
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer);
    /// @overload with explicit configuration
    void transform(const String& path, Interfaces::IMSDataConsumer* consumer, const Config& config);

  private:
    /// Load DDA-PASEF data: MS1 frames (IM_PEAK) + MS2 spectra (scalar IM)
    void loadDDA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Load DIA-PASEF data: MS1 frames (IM_PEAK) + MS2 frames split by SWATH window
    void loadDIA_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Load raw frames without precursor grouping or SWATH splitting
    void loadFrames_(TimsDataHandle& handle, MSExperiment& exp, const Config& config);

    /// Detect DDA vs DIA by checking for SWATH windows
    bool isDIA_(const String& tdf_path) const;

    /// Populate SourceFile metadata from the .d path (no peak data read)
    void loadExperimentalSettings_(const String& path, ExperimentalSettings& settings);
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
