// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/RationalScan2ImConverter.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <SQLiteCpp/SQLiteCpp.h>
#include <cmath>
#include <sstream>

namespace OpenMS
{

  RationalScan2ImConverter::RationalScan2ImConverter(
    std::unordered_map<uint32_t, Coefficients> calibrations,
    std::vector<uint32_t> frame_to_cal)
    : calibrations_(std::move(calibrations)),
      frame_to_cal_(std::move(frame_to_cal))
  {
  }

  const RationalScan2ImConverter::Coefficients&
  RationalScan2ImConverter::getCalibration(uint32_t frame_id) const
  {
    // frame_to_cal_ is 1-based (index 0 unused)
    if (frame_id > 0 && frame_id < frame_to_cal_.size())
    {
      auto it = calibrations_.find(frame_to_cal_[frame_id]);
      if (it != calibrations_.end())
      {
        return it->second;
      }
    }
    // Fallback: use first calibration entry (should not happen in valid data)
    if (frame_id != 0)
    {
      OPENMS_LOG_WARN << "RationalScan2ImConverter: frame_id " << frame_id
                      << " out of range, using first calibration" << std::endl;
    }
    return calibrations_.begin()->second;
  }

  double RationalScan2ImConverter::applyFormula(const Coefficients& c, double scan)
  {
    // V = c2 + ((c3 - c2) / c1) * (scan - c4 - c0)
    double V = c.c2 + ((c.c3 - c.c2) / c.c1) * (scan - c.c4 - c.c0);
    if (V == 0.0) V = 1e-10; // guard against division by zero (physically unreachable)
    // 1/K0 = 1 / (c6 + c7 / V)
    double denom = c.c6 + c.c7 / V;
    if (denom == 0.0) denom = 1e-10; // guard
    return 1.0 / denom;
  }

  double RationalScan2ImConverter::invertFormula(const Coefficients& c, double inv_k0)
  {
    // scan = c0 + c4 + (c1 / (c3 - c2)) * (c7 / (1/inv_k0 - c6) - c2)
    double recip = 1.0 / inv_k0;
    double inner_denom = recip - c.c6;
    if (inner_denom == 0.0) inner_denom = 1e-10; // guard
    double slope_denom = c.c3 - c.c2;
    if (slope_denom == 0.0) slope_denom = 1e-10; // guard
    return c.c0 + c.c4 + (c.c1 / slope_denom) * (c.c7 / inner_denom - c.c2);
  }

  void RationalScan2ImConverter::convert(uint32_t frame_id,
    double* inv_ion_mobilities, const double* scans, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = applyFormula(c, scans[i]);
    }
  }

  void RationalScan2ImConverter::convert(uint32_t frame_id,
    double* inv_ion_mobilities, const uint32_t* scans, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      inv_ion_mobilities[i] = applyFormula(c, static_cast<double>(scans[i]));
    }
  }

  void RationalScan2ImConverter::inverse_convert(uint32_t frame_id,
    uint32_t* scans, const double* inv_ion_mobilities, uint32_t size)
  {
    const auto& c = getCalibration(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      double val = invertFormula(c, inv_ion_mobilities[i]);
      scans[i] = val > 0.0 ? static_cast<uint32_t>(val + 0.5) : 0;
    }
  }

  std::string RationalScan2ImConverter::description() const
  {
    return "RationalScan2ImConverter (" + std::to_string(calibrations_.size())
           + " calibration segment" + (calibrations_.size() != 1 ? "s" : "") + ")";
  }

  // =========================================================================
  // Factory function: try to build a RationalScan2ImConverter from SQLite
  // =========================================================================

  std::unique_ptr<Scan2InvIonMobilityConverter> tryCreateRationalConverter(
    const std::string& tims_dir_path)
  {
    std::string tdf_path = tims_dir_path + "/analysis.tdf";

    try
    {
      SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

      // 1. Read TimsCalibration table
      std::unordered_map<uint32_t, RationalScan2ImConverter::Coefficients> calibrations;
      {
        SQLite::Statement q(db,
          "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 "
          "FROM TimsCalibration");

        while (q.executeStep())
        {
          uint32_t id = static_cast<uint32_t>(q.getColumn(0).getInt());
          int model_type = q.getColumn(1).getInt();

          if (model_type != 2)
          {
            OPENMS_LOG_WARN << "TimsCalibration ModelType=" << model_type
                            << " unsupported (only ModelType=2 implemented), "
                            << "falling back to linear converter" << std::endl;
            return nullptr;
          }

          RationalScan2ImConverter::Coefficients c;
          c.c0 = q.getColumn(2).getDouble();
          c.c1 = q.getColumn(3).getDouble();
          c.c2 = q.getColumn(4).getDouble();
          c.c3 = q.getColumn(5).getDouble();
          c.c4 = q.getColumn(6).getDouble();
          c.c5 = q.getColumn(7).getDouble();
          c.c6 = q.getColumn(8).getDouble();
          c.c7 = q.getColumn(9).getDouble();
          c.c8 = q.getColumn(10).getDouble();
          c.c9 = q.getColumn(11).getDouble();

          calibrations[id] = c;
        }
      }

      if (calibrations.empty())
      {
        OPENMS_LOG_DEBUG << "TimsCalibration table empty, "
                         << "falling back to linear converter" << std::endl;
        return nullptr;
      }

      // 2. Read frame-to-calibration mapping from Frames table
      //    (separate try/catch: column may not exist in old TDF versions)
      std::vector<uint32_t> frame_to_cal;
      try
      {
        SQLite::Statement q(db, "SELECT Id, TimsCalibration FROM Frames ORDER BY Id");
        while (q.executeStep())
        {
          uint32_t frame_id = static_cast<uint32_t>(q.getColumn(0).getInt());
          if (q.getColumn(1).isNull())
          {
            OPENMS_LOG_WARN << "Frame " << frame_id
                            << " has NULL TimsCalibration FK, "
                            << "falling back to linear converter" << std::endl;
            return nullptr;
          }
          uint32_t cal_id = static_cast<uint32_t>(q.getColumn(1).getInt());

          // Ensure vector is large enough (frame IDs are 1-based)
          if (frame_to_cal.size() <= frame_id)
          {
            frame_to_cal.resize(frame_id + 1, 0);
          }
          frame_to_cal[frame_id] = cal_id;

          // Verify calibration ID exists
          if (calibrations.find(cal_id) == calibrations.end())
          {
            OPENMS_LOG_WARN << "Frame " << frame_id
                            << " references unknown TimsCalibration Id=" << cal_id
                            << ", falling back to linear converter" << std::endl;
            return nullptr;
          }
        }
      }
      catch (const SQLite::Exception&)
      {
        // TimsCalibration column does not exist in Frames table (old TDF format)
        OPENMS_LOG_DEBUG << "Frames table has no TimsCalibration column, "
                         << "falling back to linear converter" << std::endl;
        return nullptr;
      }

      if (frame_to_cal.size() <= 1) // only index 0 (unused)
      {
        OPENMS_LOG_DEBUG << "No frames found, falling back to linear converter"
                         << std::endl;
        return nullptr;
      }

      OPENMS_LOG_INFO << "TIMS calibration: rational (TimsCalibration table, "
                      << calibrations.size() << " segment"
                      << (calibrations.size() != 1 ? "s" : "")
                      << " for " << (frame_to_cal.size() - 1) << " frames)"
                      << std::endl;

      return std::make_unique<RationalScan2ImConverter>(
        std::move(calibrations), std::move(frame_to_cal));
    }
    catch (const SQLite::Exception&)
    {
      // TimsCalibration table does not exist
      OPENMS_LOG_DEBUG << "TimsCalibration table not found in " << tdf_path
                       << ", falling back to linear converter" << std::endl;
      return nullptr;
    }
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
