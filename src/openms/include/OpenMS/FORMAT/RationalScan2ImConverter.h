// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/OpenMSConfig.h>
#include <opentims++/scan2inv_ion_mobility_converter.h>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Per-frame rational function converter using Bruker TimsCalibration table.
   *
   * This is the first open-source implementation to read per-frame calibration
   * coefficients from the TimsCalibration table in analysis.tdf. All other
   * open-source tools (OpenTIMS, AlphaTIMS, timsrust) use a simpler linear
   * approximation from GlobalMetadata that is frame-independent.
   *
   * Physical model (ModelType=2):
   *   Inner term maps scan index to TIMS elution voltage:
   *     V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
   *   Outer term maps voltage to inverse reduced ion mobility:
   *     1/K0 = 1 / (c6 + c7 / V)
   *
   * Combined:
   *   1/K0 = 1 / (c6 + c7 / (c2 + ((c3 - c2) / c1) * (scan - c4 - c0)))
   *
   * Limitations vs Bruker SDK:
   *   - No pressure compensation (ambient gas pressure drift correction)
   *   - No recalibrated-state awareness (uses raw calibration from acquisition)
   *   - Only ModelType=2 supported; unknown types should fall back to linear
   *
   * Thread safety: immutable after construction. convert() and inverse_convert()
   * are safe for concurrent calls.
   */
  class OPENMS_DLLAPI RationalScan2ImConverter : public Scan2InvIonMobilityConverter
  {
  public:
    /// Calibration coefficients from one row of the TimsCalibration table.
    struct Coefficients
    {
      double c0, c1, c2, c3, c4, c5, c6, c7, c8, c9;
      // c5, c8, c9: stored but unused in ModelType=2. No public documentation
      // assigns them meaning. Retained for forward compatibility.
    };

    /// @param calibrations  Map of calibration ID -> coefficients (from TimsCalibration table)
    /// @param frame_to_cal  Calibration ID for each frame, indexed by frame_id (1-based; index 0 unused)
    RationalScan2ImConverter(
      std::unordered_map<uint32_t, Coefficients> calibrations,
      std::vector<uint32_t> frame_to_cal);

    void convert(uint32_t frame_id, double* inv_ion_mobilities,
                 const double* scans, uint32_t size) override;
    void convert(uint32_t frame_id, double* inv_ion_mobilities,
                 const uint32_t* scans, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, uint32_t* scans,
                         const double* inv_ion_mobilities, uint32_t size) override;

    /// Returns e.g. "RationalScan2ImConverter (2 calibration segments)"
    std::string description() const override;

  private:
    std::unordered_map<uint32_t, Coefficients> calibrations_;
    std::vector<uint32_t> frame_to_cal_;  ///< indexed by frame_id (1-based)

    /// Look up calibration for a frame. Falls back to first calibration for
    /// out-of-range frame_id with a warning.
    const Coefficients& getCalibration(uint32_t frame_id) const;

    /// V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
    /// 1/K0 = 1.0 / (c6 + c7 / V)
    static double applyFormula(const Coefficients& c, double scan);

    /// scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1.0/inv_k0 - c6) - c2)
    static double invertFormula(const Coefficients& c, double inv_k0);
  };

  /// Try to create a RationalScan2ImConverter by reading the TimsCalibration
  /// table from the given .d directory. Returns nullptr if the table is missing,
  /// has unknown ModelType, or any frame has a NULL/invalid TimsCalibration FK.
  OPENMS_DLLAPI std::unique_ptr<Scan2InvIonMobilityConverter> tryCreateRationalConverter(
    const std::string& tims_dir_path);

} // namespace OpenMS

#endif // WITH_OPENTIMS
