// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/OpenMSConfig.h>
#include <opentims++/tof2mz_converter.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Per-frame TOF index -> m/z converter using the Bruker MzCalibration table.
   *
   * opentims' OpenSourceTof2MzConverter interpolates linearly in sqrt(m/z) between
   * MzAcqRangeLower/Upper and ignores the calibration stored in analysis.tdf. This
   * converter applies the stored calibration instead, with the per-frame temperature
   * correction. The formula follows rustims (rustdf/src/data/calibration.rs, MIT
   * licence), reimplemented here:
   *
   *   tc = 1 + (dC1 * (T1 - Frames.T1) + dC2 * (T2 - Frames.T2)) / 1e6
   *   b  = sqrt(1e12 / (C1 * tc)),  c2 = C2 / tc
   *   t  = tof_index * DigitizerTimebase + DigitizerDelay
   *   solve t = C0 + b*s + c2*s^2 + c3*s^3 for s,  m/z = s^2 - c4
   *
   * T1, T2, dC1, dC2 and C0..C4 come from the MzCalibration row of the frame,
   * Frames.T1 and Frames.T2 from the frame itself.
   *
   * - ModelType 1: c3 = C3 and c4 = C4. Matches the Bruker SDK to better than
   *   1e-4 ppm on the data checked (27 runs, 2017-2026).
   * - ModelType 2: the C3/C4 columns repeat C0/C2 and are ignored (c3 = c4 = 0).
   *   ModelType 2 also has a correction (columns C5..C14) between about m/z 222 and
   *   1320-1525 that is not modelled here, so m/z deviates from the Bruker SDK by up
   *   to a few ppm there (max 1.3-7.8 ppm on the four runs checked). Elsewhere the
   *   result matches the SDK.
   *
   * Thread safety: immutable after construction. convert() and inverse_convert()
   * are safe for concurrent calls.
   */
  class OPENMS_DLLAPI MzCalibrationTof2MzConverter : public Tof2MzConverter
  {
  public:
    /// One row of the MzCalibration table.
    struct Calibration
    {
      int model_type = 0;
      double digitizer_timebase = 0.0;
      double digitizer_delay = 0.0;
      double t1 = 0.0; ///< reference temperature 1
      double t2 = 0.0; ///< reference temperature 2
      double dc1 = 0.0; ///< temperature coefficient for T1 (ppm per degree)
      double dc2 = 0.0; ///< temperature coefficient for T2 (ppm per degree)
      double c0 = 0.0, c1 = 0.0, c2 = 0.0, c3 = 0.0, c4 = 0.0;
    };

    /// Calibration of one frame after the temperature correction.
    struct FrameModel
    {
      double timebase = 0.0;
      double delay = 0.0;
      double c0 = 0.0;
      double b = 0.0;  ///< sqrt(1e12 / (C1 * tc))
      double c2 = 0.0; ///< C2 / tc
      double c3 = 0.0;
      double c4 = 0.0;
    };

    /// Applies the temperature correction for a frame with temperatures @p frame_t1 and @p frame_t2.
    static FrameModel makeFrameModel(const Calibration& cal, double frame_t1, double frame_t2);

    /// @param frame_models  Model for each frame, indexed by frame_id (1-based; index 0 unused). Must not be empty.
    explicit MzCalibrationTof2MzConverter(std::vector<FrameModel> frame_models, std::string description = "");

    void convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size) override;
    void convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size) override;

    /// Returns e.g. "MzCalibrationTof2MzConverter (MzCalibration ModelType 1, per-frame temperature correction)"
    std::string description() override;

    /// TOF index -> m/z for one frame model
    static double tofToMz(const FrameModel& m, double tof_index);

    /// m/z -> TOF index (rounded) for one frame model
    static uint32_t mzToTof(const FrameModel& m, double mz);

  private:
    std::vector<FrameModel> frame_models_; ///< indexed by frame_id (1-based)
    std::string description_;

    /// Model for a frame. frame_id 0 or an unknown frame_id uses the first frame (the latter with a warning).
    const FrameModel& getModel(uint32_t frame_id) const;
  };

  /// Try to create a MzCalibrationTof2MzConverter from the MzCalibration and Frames tables of
  /// the given .d directory. Returns nullptr (and logs why) if the tables or columns are missing,
  /// a frame references an unknown or unsupported (not ModelType 1 or 2) calibration, or the
  /// coefficients are unusable. Frames without temperatures use the reference temperatures.
  OPENMS_DLLAPI std::unique_ptr<Tof2MzConverter> tryCreateMzCalibrationConverter(
    const std::string& tims_dir_path);

} // namespace OpenMS

#endif // WITH_OPENTIMS
