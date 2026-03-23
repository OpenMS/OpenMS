// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <tof2mz_converter.h>
#include <scan2inv_ion_mobility_converter.h>
#include <string>
#include <cstdint>

namespace OpenMS
{

  /**
   * @brief Open-source TOF-to-m/z converter using linear-in-sqrt calibration.
   *
   * Implements the same calibration model as timsrust: the relationship between
   * TOF index and m/z is linear in sqrt(m/z) space.
   *
   * @code
   *   sqrt(mz) = intercept + slope * tof_index
   *   mz       = (intercept + slope * tof_index)^2
   * @endcode
   *
   * Calibration parameters are derived from the GlobalMetadata table in
   * analysis.tdf: MzAcqRangeLower, MzAcqRangeUpper, DigitizerNumSamples.
   */
  class OPENMS_DLLAPI OpenSourceTof2MzConverter : public Tof2MzConverter
  {
  public:
    /**
     * @brief Construct from acquisition range metadata.
     *
     * @param mz_min         MzAcqRangeLower from GlobalMetadata
     * @param mz_max         MzAcqRangeUpper from GlobalMetadata
     * @param tof_max_index  DigitizerNumSamples from GlobalMetadata
     * @param is_otof_control  True if AcquisitionSoftware == "Bruker otofControl"
     *                         (extends range by +/-5 Da)
     */
    OpenSourceTof2MzConverter(double mz_min, double mz_max, uint32_t tof_max_index,
                              bool is_otof_control);

    void convert(uint32_t frame_id, double* mzs, const double* tofs, uint32_t size) override;
    void convert(uint32_t frame_id, double* mzs, const uint32_t* tofs, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, uint32_t* tofs, const double* mzs, uint32_t size) override;
    std::string description() override;

    /// Update calibration parameters (used by OLS recalibration)
    void updateCalibration(double new_intercept, double new_slope);
    double intercept() const { return intercept_; }
    double slope() const { return slope_; }

  private:
    double intercept_;
    double slope_;
  };

  /**
   * @brief Open-source scan-to-inverse-ion-mobility converter (linear model).
   *
   * The scan number and 1/K0 are related by a simple linear model:
   *
   * @code
   *   1/K0 = intercept + slope * scan_index
   * @endcode
   *
   * where intercept = OneOverK0AcqRangeUpper and slope is derived from the
   * acquisition range and maximum number of scans per frame.
   */
  class OPENMS_DLLAPI OpenSourceScan2ImConverter : public Scan2InvIonMobilityConverter
  {
  public:
    /**
     * @brief Construct from ion mobility acquisition range metadata.
     *
     * @param im_min          OneOverK0AcqRangeLower from GlobalMetadata
     * @param im_max          OneOverK0AcqRangeUpper from GlobalMetadata
     * @param scan_max_index  MAX(NumScans) from the Frames table
     */
    OpenSourceScan2ImConverter(double im_min, double im_max, uint32_t scan_max_index);

    void convert(uint32_t frame_id, double* inv_ion_mobilities, const double* scans, uint32_t size) override;
    void convert(uint32_t frame_id, double* inv_ion_mobilities, const uint32_t* scans, uint32_t size) override;
    void inverse_convert(uint32_t frame_id, uint32_t* scans, const double* inv_ion_mobilities, uint32_t size) override;
    std::string description() const override;

  private:
    double intercept_;
    double slope_;
  };

  /**
   * @brief Factory for OpenSourceTof2MzConverter.
   *
   * Reads calibration metadata from the analysis.tdf SQLite database inside
   * the .d directory and constructs an OpenSourceTof2MzConverter.
   */
  class OPENMS_DLLAPI OpenSourceTof2MzConverterFactory : public Tof2MzConverterFactory
  {
  public:
    std::unique_ptr<Tof2MzConverter> produce(TimsDataHandle& TDH,
      pressure_compensation_strategy pcs = NoPressureCompensation) override;
  };

  /**
   * @brief Factory for OpenSourceScan2ImConverter.
   *
   * Reads ion mobility calibration metadata from the analysis.tdf SQLite
   * database and the Frames table to construct an OpenSourceScan2ImConverter.
   */
  class OPENMS_DLLAPI OpenSourceScan2ImConverterFactory : public Scan2InvIonMobilityConverterFactory
  {
  public:
    std::unique_ptr<Scan2InvIonMobilityConverter> produce(TimsDataHandle& TDH,
      pressure_compensation_strategy pcs = NoPressureCompensation) override;
  };

} // namespace OpenMS

#endif // WITH_OPENTIMS
