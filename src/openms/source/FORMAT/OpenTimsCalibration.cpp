// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/OpenTimsCalibration.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <opentims.h>
#include <SQLiteCpp/SQLiteCpp.h>
#include <cmath>

namespace OpenMS
{

  // --- OpenSourceTof2MzConverter ---

  OpenSourceTof2MzConverter::OpenSourceTof2MzConverter(
    double mz_min, double mz_max, uint32_t tof_max_index, bool is_otof_control)
  {
    if (is_otof_control)
    {
      mz_min -= 5.0;
      mz_max += 5.0;
    }
    intercept_ = std::sqrt(mz_min);
    slope_ = (std::sqrt(mz_max) - std::sqrt(mz_min)) / static_cast<double>(tof_max_index);
  }

  void OpenSourceTof2MzConverter::convert(uint32_t /*frame_id*/, double* mzs,
    const double* tofs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = intercept_ + slope_ * tofs[i];
      mzs[i] = val * val;
    }
  }

  void OpenSourceTof2MzConverter::convert(uint32_t /*frame_id*/, double* mzs,
    const uint32_t* tofs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = intercept_ + slope_ * static_cast<double>(tofs[i]);
      mzs[i] = val * val;
    }
  }

  void OpenSourceTof2MzConverter::inverse_convert(uint32_t /*frame_id*/, uint32_t* tofs,
    const double* mzs, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      tofs[i] = static_cast<uint32_t>((std::sqrt(mzs[i]) - intercept_) / slope_);
    }
  }

  std::string OpenSourceTof2MzConverter::description()
  {
    return "OpenSourceTof2MzConverter (linear-in-sqrt, OpenMS)";
  }

  void OpenSourceTof2MzConverter::updateCalibration(double new_intercept, double new_slope)
  {
    intercept_ = new_intercept;
    slope_ = new_slope;
  }

  // --- OpenSourceScan2ImConverter ---

  OpenSourceScan2ImConverter::OpenSourceScan2ImConverter(
    double im_min, double im_max, uint32_t scan_max_index)
  {
    intercept_ = im_max;
    slope_ = (im_min - im_max) / static_cast<double>(scan_max_index);
  }

  void OpenSourceScan2ImConverter::convert(uint32_t /*frame_id*/, double* inv_ion_mobilities,
    const double* scans, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = intercept_ + slope_ * scans[i];
    }
  }

  void OpenSourceScan2ImConverter::convert(uint32_t /*frame_id*/, double* inv_ion_mobilities,
    const uint32_t* scans, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = intercept_ + slope_ * static_cast<double>(scans[i]);
    }
  }

  void OpenSourceScan2ImConverter::inverse_convert(uint32_t /*frame_id*/, uint32_t* scans,
    const double* inv_ion_mobilities, uint32_t size)
  {
    for (uint32_t i = 0; i < size; ++i)
    {
      scans[i] = static_cast<uint32_t>((inv_ion_mobilities[i] - intercept_) / slope_);
    }
  }

  std::string OpenSourceScan2ImConverter::description() const
  {
    return "OpenSourceScan2ImConverter (linear, OpenMS)";
  }

  // --- Factories ---

  std::unique_ptr<Tof2MzConverter> OpenSourceTof2MzConverterFactory::produce(
    TimsDataHandle& TDH, pressure_compensation_strategy /*pcs*/)
  {
    std::string tdf_path = TDH.get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

    // Read calibration parameters from GlobalMetadata
    double mz_min = 0, mz_max = 0;
    uint32_t tof_max = 0;
    bool is_otof = false;

    SQLite::Statement query(db, "SELECT Key, Value FROM GlobalMetadata "
      "WHERE Key IN ('MzAcqRangeLower','MzAcqRangeUpper','DigitizerNumSamples','AcquisitionSoftware')");

    while (query.executeStep())
    {
      std::string key = query.getColumn(0).getString();
      std::string val = query.getColumn(1).getString();
      if (key == "MzAcqRangeLower") mz_min = std::stod(val);
      else if (key == "MzAcqRangeUpper") mz_max = std::stod(val);
      else if (key == "DigitizerNumSamples") tof_max = static_cast<uint32_t>(std::stoul(val));
      else if (key == "AcquisitionSoftware") is_otof = (val == "Bruker otofControl");
    }

    if (mz_min <= 0 || mz_max <= 0 || tof_max == 0)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        tdf_path, "OpenSourceTof2MzConverterFactory: missing calibration metadata (MzAcqRangeLower/MzAcqRangeUpper/DigitizerNumSamples)");

    return std::make_unique<OpenSourceTof2MzConverter>(mz_min, mz_max, tof_max, is_otof);
  }

  std::unique_ptr<Scan2InvIonMobilityConverter> OpenSourceScan2ImConverterFactory::produce(
    TimsDataHandle& TDH, pressure_compensation_strategy /*pcs*/)
  {
    std::string tdf_path = TDH.get_tims_dir_path() + "/analysis.tdf";
    SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

    double im_min = 0, im_max = 0;

    SQLite::Statement meta_query(db, "SELECT Key, Value FROM GlobalMetadata "
      "WHERE Key IN ('OneOverK0AcqRangeLower','OneOverK0AcqRangeUpper')");
    while (meta_query.executeStep())
    {
      std::string key = meta_query.getColumn(0).getString();
      std::string val = meta_query.getColumn(1).getString();
      if (key == "OneOverK0AcqRangeLower") im_min = std::stod(val);
      else if (key == "OneOverK0AcqRangeUpper") im_max = std::stod(val);
    }

    SQLite::Statement scan_query(db, "SELECT MAX(NumScans) FROM Frames");
    uint32_t scan_max = 0;
    if (scan_query.executeStep())
    {
      scan_max = static_cast<uint32_t>(scan_query.getColumn(0).getInt());
    }

    if (im_min <= 0 || im_max <= 0 || scan_max == 0)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        tdf_path, "OpenSourceScan2ImConverterFactory: missing calibration metadata (OneOverK0AcqRangeLower/OneOverK0AcqRangeUpper/NumScans)");

    return std::make_unique<OpenSourceScan2ImConverter>(im_min, im_max, scan_max);
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
