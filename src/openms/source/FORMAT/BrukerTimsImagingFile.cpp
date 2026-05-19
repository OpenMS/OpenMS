// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/FORMAT/BrukerTimsImagingFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/SYSTEM/File.h>

#include <sqlite3.h>

#include <algorithm>
#include <limits>
#include <unordered_map>

namespace OpenMS
{
  namespace
  {
    // Parse the integer frame id from a Bruker TDF native ID of the form
    // "frame=<id>" or "frame=<id> ...". Returns -1 if the prefix is missing
    // or the integer cannot be parsed.
    Int parseFrameId(const String& native_id)
    {
      const std::string prefix = "frame=";
      const std::size_t pos = native_id.find(prefix);
      if (pos != 0) { return -1; }
      const std::size_t start = prefix.size();
      std::size_t end = native_id.find(' ', start);
      if (end == std::string::npos) { end = native_id.size(); }
      try
      {
        return String(native_id.substr(start, end - start)).toInt();
      }
      catch (...)
      {
        return -1;
      }
    }
  }

  String BrukerTimsImagingFile::resolveTdfPath_(const String& d_folder)
  {
    const String tdf = d_folder + "/analysis.tdf";
    if (!File::exists(tdf))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tdf);
    }
    return tdf;
  }

  bool BrukerTimsImagingFile::isImagingDataset(const String& path)
  {
    try
    {
      return readGlobalMetadataValue(path, "MaldiApplicationType") == "Imaging";
    }
    catch (...)
    {
      return false;
    }
  }

  String BrukerTimsImagingFile::readGlobalMetadataValue(const String& d_folder, const String& key)
  {
    const String tdf = resolveTdfPath_(d_folder);
    SqliteConnector conn(tdf, SqliteConnector::SqlOpenMode::READ_ONLY);
    if (!conn.tableExists("GlobalMetadata")) { return ""; }

    sqlite3_stmt* stmt = nullptr;
    conn.prepareStatement(&stmt, "SELECT Value FROM GlobalMetadata WHERE Key = ?1 LIMIT 1;");
    sqlite3_bind_text(stmt, 1, key.c_str(), -1, SQLITE_TRANSIENT);

    String result;
    if (sqlite3_step(stmt) == SQLITE_ROW)
    {
      const unsigned char* txt = sqlite3_column_text(stmt, 0);
      if (txt != nullptr) { result = reinterpret_cast<const char*>(txt); }
    }
    sqlite3_finalize(stmt);
    return result;
  }

  std::vector<BrukerTimsImagingFile::MaldiPixel>
  BrukerTimsImagingFile::readMaldiFrameInfo(const String& d_folder)
  {
    const String tdf = resolveTdfPath_(d_folder);
    SqliteConnector conn(tdf, SqliteConnector::SqlOpenMode::READ_ONLY);
    if (!conn.tableExists("MaldiFrameInfo"))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  tdf, "MaldiFrameInfo table missing");
    }

    sqlite3_stmt* stmt = nullptr;
    conn.prepareStatement(&stmt, "SELECT Frame, XIndexPos, YIndexPos FROM MaldiFrameInfo ORDER BY Frame;");

    std::vector<MaldiPixel> rows;
    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
      MaldiPixel p;
      p.frame_id = sqlite3_column_int(stmt, 0);
      p.x        = sqlite3_column_int(stmt, 1);
      p.y        = sqlite3_column_int(stmt, 2);
      rows.push_back(p);
    }
    sqlite3_finalize(stmt);

    if (rows.empty())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  tdf, "MaldiFrameInfo is empty");
    }
    return rows;
  }

  void BrukerTimsImagingFile::computeBoundingBox_(const std::vector<MaldiPixel>& rows,
                                                  Int& x_min, Int& y_min,
                                                  UInt& width, UInt& height)
  {
    Int x_max = std::numeric_limits<Int>::min();
    Int y_max = std::numeric_limits<Int>::min();
    x_min = y_min = std::numeric_limits<Int>::max();

    for (const MaldiPixel& p : rows)
    {
      x_min = std::min(x_min, p.x); x_max = std::max(x_max, p.x);
      y_min = std::min(y_min, p.y); y_max = std::max(y_max, p.y);
    }
    width  = static_cast<UInt>(x_max - x_min + 1);
    height = static_cast<UInt>(y_max - y_min + 1);
  }

  void BrukerTimsImagingFile::load(const String& path, MSImagingExperiment& exp)
  {
    load(path, exp, Config{});
  }

  void BrukerTimsImagingFile::load(const String& path, MSImagingExperiment& exp, const Config& config)
  {
    if (!File::exists(path))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path);
    }

    // Resolve analysis.tdf eagerly so an I/O error surfaces as FileNotReadable
    // here, rather than being masked by isImagingDataset() (which swallows
    // exceptions to behave as a non-throwing probe).
    const String tdf = resolveTdfPath_(path);

    if (config.strict_imaging_only && !isImagingDataset(path))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Not a MALDI imaging dataset (MaldiApplicationType != 'Imaging')", path);
    }

    // Reject multi-section datasets — geometry is strictly 2D, callers should
    // split sections into separate MSImagingExperiment objects.
    {
      SqliteConnector conn(tdf, SqliteConnector::SqlOpenMode::READ_ONLY);
      if (conn.columnExists("MaldiFrameInfo", "ZIndexPos"))
      {
        sqlite3_stmt* stmt = nullptr;
        conn.prepareStatement(&stmt, "SELECT COUNT(DISTINCT ZIndexPos) FROM MaldiFrameInfo;");
        int distinct_z = 0;
        if (sqlite3_step(stmt) == SQLITE_ROW) { distinct_z = sqlite3_column_int(stmt, 0); }
        sqlite3_finalize(stmt);
        if (distinct_z > 1)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Multi-section MALDI imaging (>1 distinct ZIndexPos) is out of scope "
                                        "for the 2D MSImagingGeometry; split into per-section .d folders or "
                                        "filter MaldiFrameInfo before loading", path);
        }
      }
    }

    const std::vector<MaldiPixel> rows = readMaldiFrameInfo(path);

    BrukerTimsFile::Config inner_config = config.inner_config;
    inner_config.export_mode = BrukerTimsFile::Config::ExportMode::FRAME;
    inner_config.load_ms1 = true;

    MSExperiment ms_exp;
    BrukerTimsFile inner;
    inner.setLogType(this->getLogType());
    inner.load(path, ms_exp, inner_config);

    std::unordered_map<Int, Size> frame_to_index;
    frame_to_index.reserve(ms_exp.size());
    for (Size i = 0; i < ms_exp.size(); ++i)
    {
      const Int frame_id = parseFrameId(ms_exp[i].getNativeID());
      if (frame_id < 0) { continue; }
      frame_to_index[frame_id] = i;
    }

    Int x_min = 0, y_min = 0;
    UInt width = 0, height = 0;
    computeBoundingBox_(rows, x_min, y_min, width, height);

    MSImagingGeometry geometry;
    geometry.setDimensions(width, height);

    String spot_size_str = readGlobalMetadataValue(path, "MaldiSpotSize_um");
    if (spot_size_str.empty()) { spot_size_str = readGlobalMetadataValue(path, "ImagingAreaSpotSize_um"); }
    double spot_size = 1.0;
    if (!spot_size_str.empty())
    {
      try { spot_size = spot_size_str.toDouble(); }
      catch (...) { spot_size = 1.0; }
    }
    geometry.setPixelSize(spot_size, spot_size, "micrometer");

    Size dropped = 0;
    for (const MaldiPixel& row : rows)
    {
      const auto it = frame_to_index.find(row.frame_id);
      if (it == frame_to_index.end()) { ++dropped; continue; }
      geometry.addPixel(static_cast<UInt>(row.x - x_min),
                        static_cast<UInt>(row.y - y_min),
                        it->second);
    }
    if (dropped > 0)
    {
      OPENMS_LOG_WARN << "BrukerTimsImagingFile: dropped " << dropped
                      << " pixel(s) whose Frame ID was not present in the loaded spectra." << std::endl;
    }

    exp.setMSExperiment(std::move(ms_exp));
    exp.setGeometry(std::move(geometry));
    exp.validate();
  }

} // namespace OpenMS

#endif // WITH_OPENTIMS
