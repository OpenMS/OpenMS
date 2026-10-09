// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/MzCalibrationTof2MzConverter.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <SQLiteCpp/SQLiteCpp.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

namespace OpenMS
{

  namespace
  {
    /// ModelType 2 correction at base-curve m/z @p m: polynomial inside [low, high], 0 outside
    double correctionAt(const MzCalibrationTof2MzConverter::Correction& c, double m)
    {
      if (c.n <= 0 || m < c.low || m > c.high) return 0.0;
      double value = 0.0;
      for (int k = c.n - 1; k >= 0; --k)
      {
        value = value * m + c.coefficients[k];
      }
      return value;
    }
  }

  MzCalibrationTof2MzConverter::FrameModel MzCalibrationTof2MzConverter::makeFrameModel(
    const Calibration& cal, double frame_t1, double frame_t2)
  {
    const double tc = 1.0 + (cal.dc1 * (cal.t1 - frame_t1) + cal.dc2 * (cal.t2 - frame_t2)) / 1.0e6;
    FrameModel m;
    m.timebase = cal.digitizer_timebase;
    m.delay = cal.digitizer_delay;
    m.c0 = cal.c0;
    m.b = std::sqrt(1.0e12 / (cal.c1 * tc));
    m.c2 = cal.c2 / tc;
    // In ModelType 2 the C3/C4 columns repeat C0/C2; they are no cubic term or mass offset there
    if (cal.model_type == 1)
    {
      m.c3 = cal.c3;
      m.c4 = cal.c4;
    }
    else if (cal.model_type == 2)
    {
      m.correction = cal.correction;
    }
    return m;
  }

  MzCalibrationTof2MzConverter::MzCalibrationTof2MzConverter(
    std::vector<FrameModel> frame_models, std::string description)
    : frame_models_(std::move(frame_models)),
      description_(std::move(description))
  {
  }

  const MzCalibrationTof2MzConverter::FrameModel&
  MzCalibrationTof2MzConverter::getModel(uint32_t frame_id) const
  {
    // frame_models_ is 1-based (index 0 is a copy of the first frame); b == 0 marks a frame_id gap
    if (frame_id < frame_models_.size() && frame_models_[frame_id].b != 0.0)
    {
      return frame_models_[frame_id];
    }
    OPENMS_LOG_WARN << "MzCalibrationTof2MzConverter: frame_id " << frame_id
                    << " out of range, using first frame" << std::endl;
    return frame_models_[frame_models_.size() > 1 ? 1 : 0];
  }

  double MzCalibrationTof2MzConverter::tofToMz(const FrameModel& m, double tof_index)
  {
    const double t = tof_index * m.timebase + m.delay;
    // Linear-in-sqrt(m/z) estimate; exact when c2 = c3 = 0
    double s = (t - m.c0) / m.b;
    if (m.c3 != 0.0)
    {
      // cubic: Newton iterations from the linear estimate
      for (int i = 0; i < 8; ++i)
      {
        const double f = m.c0 + s * (m.b + s * (m.c2 + s * m.c3)) - t;
        const double df = m.b + s * (2.0 * m.c2 + 3.0 * m.c3 * s);
        if (df == 0.0) break;
        const double step = f / df;
        s -= step;
        if (std::abs(step) < 1e-12) break;
      }
    }
    else if (m.c2 != 0.0)
    {
      // quadratic c2*s^2 + b*s + (c0 - t) = 0; this form of the root avoids cancellation for small c2
      const double disc = m.b * m.b - 4.0 * m.c2 * (m.c0 - t);
      if (disc >= 0.0)
      {
        const double q = -0.5 * (m.b + std::sqrt(disc));
        s = (m.c0 - t) / q;
      }
    }
    const double mz = s * s - m.c4;
    return mz - correctionAt(m.correction, mz);
  }

  uint32_t MzCalibrationTof2MzConverter::mzToTof(const FrameModel& m, double mz)
  {
    // Undo the ModelType 2 correction: solve base - correction(base) = mz. The correction is a few
    // mDa and changes by far less than 1 per m/z unit, so the fixed-point iteration converges fast.
    double base = mz;
    for (int i = 0; i < 4; ++i)
    {
      base = mz + correctionAt(m.correction, base);
    }
    const double s = std::sqrt(std::max(base + m.c4, 0.0));
    const double t = m.c0 + s * (m.b + s * (m.c2 + s * m.c3));
    const double tof = (t - m.delay) / m.timebase;
    return tof > 0.0 ? static_cast<uint32_t>(tof + 0.5) : 0;
  }

  void MzCalibrationTof2MzConverter::convert(uint32_t frame_id,
    double* mzs, const double* tofs, uint32_t size)
  {
    const auto& m = getModel(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      mzs[i] = tofToMz(m, tofs[i]);
    }
  }

  void MzCalibrationTof2MzConverter::convert(uint32_t frame_id,
    double* mzs, const uint32_t* tofs, uint32_t size)
  {
    const auto& m = getModel(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      mzs[i] = tofToMz(m, static_cast<double>(tofs[i]));
    }
  }

  void MzCalibrationTof2MzConverter::inverse_convert(uint32_t frame_id,
    uint32_t* tofs, const double* mzs, uint32_t size)
  {
    const auto& m = getModel(frame_id);
    for (uint32_t i = 0; i < size; ++i)
    {
      tofs[i] = mzToTof(m, mzs[i]);
    }
  }

  std::string MzCalibrationTof2MzConverter::description()
  {
    return description_.empty() ? std::string("MzCalibrationTof2MzConverter") : description_;
  }

  // =========================================================================
  // Factory function: try to build a MzCalibrationTof2MzConverter from SQLite
  // =========================================================================

  std::unique_ptr<Tof2MzConverter> tryCreateMzCalibrationConverter(
    const std::string& tims_dir_path)
  {
    const std::string tdf_path = tims_dir_path + "/analysis.tdf";
    using Calibration = MzCalibrationTof2MzConverter::Calibration;

    try
    {
      SQLite::Database db(tdf_path, SQLite::OPEN_READONLY);

      // 1. Read the MzCalibration table. Rows that no frame uses may be unsupported or
      //    incomplete (older files carry such rows); only rows used by a frame are checked.
      std::unordered_map<uint32_t, Calibration> calibrations;
      std::unordered_map<uint32_t, std::string> unusable; // calibration id -> reason
      {
        SQLite::Statement q(db,
          "SELECT Id, ModelType, DigitizerTimebase, DigitizerDelay, T1, T2, dC1, dC2, "
          "C0, C1, C2, C3, C4 FROM MzCalibration");
        while (q.executeStep())
        {
          const uint32_t id = static_cast<uint32_t>(q.getColumn(0).getInt());
          const int model_type = q.getColumn(1).getInt();
          if (model_type != 1 && model_type != 2)
          {
            unusable[id] = "ModelType " + std::to_string(model_type) + " is not supported (only 1 and 2)";
            continue;
          }
          bool has_null = false;
          for (int col = 2; col <= 10; ++col) // DigitizerTimebase .. C2
          {
            has_null = has_null || q.getColumn(col).isNull();
          }
          if (has_null)
          {
            unusable[id] = "NULL coefficients";
            continue;
          }
          Calibration c;
          c.model_type = model_type;
          c.digitizer_timebase = q.getColumn(2).getDouble();
          c.digitizer_delay = q.getColumn(3).getDouble();
          c.t1 = q.getColumn(4).getDouble();
          c.t2 = q.getColumn(5).getDouble();
          c.dc1 = q.getColumn(6).getDouble();
          c.dc2 = q.getColumn(7).getDouble();
          c.c0 = q.getColumn(8).getDouble();
          c.c1 = q.getColumn(9).getDouble();
          c.c2 = q.getColumn(10).getDouble();
          c.c3 = q.getColumn(11).isNull() ? 0.0 : q.getColumn(11).getDouble();
          c.c4 = q.getColumn(12).isNull() ? 0.0 : q.getColumn(12).getDouble();
          if (!(c.c1 > 0.0) || !(c.digitizer_timebase > 0.0))
          {
            unusable[id] = "C1 and DigitizerTimebase must be positive";
            continue;
          }
          calibrations[id] = c;
        }
      }

      // ModelType 2: correction polynomial in C5..C14 (these columns only exist in files with ModelType 2)
      size_t n_without_correction = 0;
      for (auto& [id, cal] : calibrations)
      {
        if (cal.model_type != 2) continue;
        try
        {
          SQLite::Statement q(db, "SELECT C5, C6, C7, C8, C9, C10, C11, C12, C13, C14 FROM MzCalibration WHERE Id = ?");
          q.bind(1, static_cast<int64_t>(id));
          if (!q.executeStep()) continue;
          bool usable = !q.getColumn(0).isNull() && !q.getColumn(1).isNull() && !q.getColumn(2).isNull();
          MzCalibrationTof2MzConverter::Correction corr;
          if (usable)
          {
            corr.low = q.getColumn(0).getDouble();
            corr.high = q.getColumn(1).getDouble();
            corr.n = q.getColumn(2).getInt();
            usable = corr.low < corr.high && corr.n >= 1 && corr.n <= static_cast<int>(corr.coefficients.size());
          }
          for (int k = 0; usable && k < corr.n; ++k)
          {
            usable = !q.getColumn(3 + k).isNull();
            if (usable) corr.coefficients[k] = q.getColumn(3 + k).getDouble();
          }
          if (usable)
          {
            cal.correction = corr;
          }
          else
          {
            ++n_without_correction;
          }
        }
        catch (const SQLite::Exception&)
        {
          ++n_without_correction; // no C5..C14 columns
        }
      }

      // 2. One model per frame, with the frame's temperatures
      std::vector<MzCalibrationTof2MzConverter::FrameModel> frame_models;
      std::set<int> model_types;
      size_t n_frames = 0;
      size_t n_without_temperature = 0;
      {
        SQLite::Statement q(db, "SELECT Id, MzCalibration, T1, T2 FROM Frames ORDER BY Id");
        while (q.executeStep())
        {
          const uint32_t frame_id = static_cast<uint32_t>(q.getColumn(0).getInt());
          if (q.getColumn(1).isNull())
          {
            OPENMS_LOG_WARN << "Frame " << frame_id << " has no MzCalibration, "
                            << "falling back to the linear m/z converter" << std::endl;
            return nullptr;
          }
          const uint32_t cal_id = static_cast<uint32_t>(q.getColumn(1).getInt());
          auto it = calibrations.find(cal_id);
          if (it == calibrations.end())
          {
            auto bad = unusable.find(cal_id);
            OPENMS_LOG_WARN << "Frame " << frame_id << " uses MzCalibration Id=" << cal_id << " ("
                            << (bad != unusable.end() ? bad->second : std::string("unknown Id"))
                            << "), falling back to the linear m/z converter" << std::endl;
            return nullptr;
          }
          const Calibration& cal = it->second;
          double frame_t1 = cal.t1;
          double frame_t2 = cal.t2;
          if (q.getColumn(2).isNull() || q.getColumn(3).isNull())
          {
            ++n_without_temperature;
          }
          else
          {
            frame_t1 = q.getColumn(2).getDouble();
            frame_t2 = q.getColumn(3).getDouble();
          }
          if (frame_models.size() <= frame_id)
          {
            frame_models.resize(frame_id + 1);
          }
          frame_models[frame_id] = MzCalibrationTof2MzConverter::makeFrameModel(cal, frame_t1, frame_t2);
          model_types.insert(cal.model_type);
          ++n_frames;
        }
      }

      if (n_frames == 0)
      {
        OPENMS_LOG_DEBUG << "No frames found, falling back to the linear m/z converter" << std::endl;
        return nullptr;
      }
      // index 0 is used for frame_id 0: same model as the first frame
      const auto first = std::find_if(frame_models.begin() + 1, frame_models.end(),
                                      [](const auto& m) { return m.b != 0.0; });
      if (first != frame_models.end() && frame_models[0].b == 0.0)
      {
        frame_models[0] = *first;
      }

      if (n_without_temperature > 0)
      {
        OPENMS_LOG_WARN << n_without_temperature << " frame(s) have no T1/T2; their m/z uses the "
                        << "reference temperatures of the MzCalibration table" << std::endl;
      }

      std::string types;
      for (int t : model_types)
      {
        types += (types.empty() ? "" : "/") + std::to_string(t);
      }
      std::string description = "MzCalibrationTof2MzConverter (MzCalibration ModelType " + types
                                + ", per-frame temperature correction)";
      OPENMS_LOG_INFO << "TIMS m/z calibration: MzCalibration table (ModelType " << types << ", "
                      << n_frames << " frames, per-frame temperature correction)" << std::endl;
      if (n_without_correction > 0 && model_types.count(2))
      {
        OPENMS_LOG_WARN << "TIMS m/z calibration: " << n_without_correction << " ModelType 2 calibration(s) "
                        << "without usable correction terms (C5..C14); m/z uses the base curve only, which "
                        << "can deviate from the Bruker SDK by a few ppm" << std::endl;
      }

      return std::make_unique<MzCalibrationTof2MzConverter>(std::move(frame_models), std::move(description));
    }
    catch (const SQLite::Exception& e)
    {
      // MzCalibration table or the Frames T1/T2/MzCalibration columns do not exist
      OPENMS_LOG_DEBUG << "MzCalibration not usable in " << tdf_path << " (" << e.what()
                       << "), falling back to the linear m/z converter" << std::endl;
      return nullptr;
    }
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
