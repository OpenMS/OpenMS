// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/BrukerTimsImagingFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/SYSTEM/File.h>

#include <filesystem>

using namespace OpenMS;
using namespace std;

namespace
{
  // Build a fake .d directory with only an analysis.tdf SQLite file inside.
  // Optional MaldiApplicationType + MaldiFrameInfo rows.
  // Returns the path to the .d folder (a subdir of @p parent).
  // Build a fake .d directory with only an analysis.tdf SQLite file inside.
  // Each row in @p xyz_rows is (frame, x, y, z); pass z<0 to omit the
  // ZIndexPos column from the schema entirely.
  String createFakeD(const String& parent,
                     const String& name,
                     const String& maldi_application_type,
                     const std::vector<std::tuple<int,int,int,int>>& xyz_rows,
                     bool with_z_column = false,
                     const String& spot_size_um = "")
  {
    const String d_path = parent + "/" + name + ".d";
    std::filesystem::create_directories(d_path.c_str());
    const String tdf_path = d_path + "/analysis.tdf";

    SqliteConnector conn(tdf_path, SqliteConnector::SqlOpenMode::READWRITE_OR_CREATE);

    conn.executeStatement("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT);");
    if (!maldi_application_type.empty())
    {
      conn.executeStatement("INSERT INTO GlobalMetadata (Key, Value) VALUES ('MaldiApplicationType', '"
                             + maldi_application_type + "');");
    }
    if (!spot_size_um.empty())
    {
      conn.executeStatement("INSERT INTO GlobalMetadata (Key, Value) VALUES ('MaldiSpotSize_um', '"
                             + spot_size_um + "');");
    }

    if (with_z_column)
    {
      conn.executeStatement("CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER, ZIndexPos INTEGER);");
    }
    else
    {
      conn.executeStatement("CREATE TABLE MaldiFrameInfo (Frame INTEGER, XIndexPos INTEGER, YIndexPos INTEGER);");
    }
    for (const auto& r : xyz_rows)
    {
      const int frame = std::get<0>(r);
      const int x     = std::get<1>(r);
      const int y     = std::get<2>(r);
      const int z     = std::get<3>(r);
      String sql = "INSERT INTO MaldiFrameInfo VALUES (" + String(frame) + ", "
                   + String(x) + ", " + String(y);
      if (with_z_column) { sql += ", " + String(z); }
      sql += ");";
      conn.executeStatement(sql);
    }

    return d_path;
  }
}

START_TEST(BrukerTimsImagingFile, "$Id$")

File::TempDir tmp(false);

START_SECTION(static bool isImagingDataset(const String& path))
{
  // Imaging dataset
  const String d_img = createFakeD(tmp.getPath(), "imaging", "Imaging",
                                    {{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 1, 1, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_img), true)

  // SingleSpectra dataset (dried droplet)
  const String d_dd  = createFakeD(tmp.getPath(), "droplet", "SingleSpectra",
                                    {{1, 0, 0, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_dd), false)

  // Missing folder
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(tmp.getPath() + "/does_not_exist.d"), false)
}
END_SECTION

START_SECTION(static String readGlobalMetadataValue(const String& d_folder, const String& key))
{
  const String d = createFakeD(tmp.getPath(), "gm", "Imaging",
                                {{1, 0, 0, 0}}, /*with_z_column=*/false, /*spot=*/"50");
  TEST_EQUAL(BrukerTimsImagingFile::readGlobalMetadataValue(d, "MaldiApplicationType"), "Imaging")
  TEST_EQUAL(BrukerTimsImagingFile::readGlobalMetadataValue(d, "MaldiSpotSize_um"), "50")
  TEST_EQUAL(BrukerTimsImagingFile::readGlobalMetadataValue(d, "UnknownKey"), "")

  TEST_EXCEPTION(Exception::FileNotReadable,
                 BrukerTimsImagingFile::readGlobalMetadataValue(tmp.getPath() + "/nope.d", "MaldiApplicationType"))
}
END_SECTION

START_SECTION(static std::vector<MaldiPixel> readMaldiFrameInfo(const String& d_folder))
{
  // 2x2 grid with deliberate offset (XIndexPos in {10, 11}, YIndexPos in {5, 6}).
  // readMaldiFrameInfo does not normalize; it returns raw coords.
  const String d = createFakeD(tmp.getPath(), "grid2x2", "Imaging",
                                {{10, 10, 5, 0}, {11, 11, 5, 0}, {12, 10, 6, 0}, {13, 11, 6, 0}});
  auto rows = BrukerTimsImagingFile::readMaldiFrameInfo(d);
  TEST_EQUAL(rows.size(), 4)
  TEST_EQUAL(rows[0].frame_id, 10)
  TEST_EQUAL(rows[0].x, 10)
  TEST_EQUAL(rows[0].y, 5)
  TEST_EQUAL(rows[3].frame_id, 13)
  TEST_EQUAL(rows[3].x, 11)
  TEST_EQUAL(rows[3].y, 6)

  // ZIndexPos column is present but ignored by the 2D reader.
  const String d_with_z = createFakeD(tmp.getPath(), "withz", "Imaging",
                                       {{1, 3, 4, 0}, {2, 3, 5, 0}}, /*with_z_column=*/true);
  auto rows_z = BrukerTimsImagingFile::readMaldiFrameInfo(d_with_z);
  TEST_EQUAL(rows_z.size(), 2)
  TEST_EQUAL(rows_z[0].x, 3)
  TEST_EQUAL(rows_z[0].y, 4)
  TEST_EQUAL(rows_z[1].y, 5)

  // Empty MaldiFrameInfo -> ParseError
  const String d_empty = createFakeD(tmp.getPath(), "empty", "Imaging", {});
  TEST_EXCEPTION(Exception::ParseError, BrukerTimsImagingFile::readMaldiFrameInfo(d_empty))
}
END_SECTION

START_SECTION(void load(const String& path, MSImagingExperiment& exp))
{
  // Non-existent path -> FileNotReadable.
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  TEST_EXCEPTION(Exception::FileNotReadable, f.load(tmp.getPath() + "/nope.d", exp))

  // .d folder exists but analysis.tdf is missing -> FileNotReadable
  // (must not be masked as InvalidValue by the strict-imaging check).
  const String d_no_tdf = tmp.getPath() + "/no_tdf.d";
  std::filesystem::create_directories(d_no_tdf.c_str());
  TEST_EXCEPTION(Exception::FileNotReadable, f.load(d_no_tdf, exp))
}
END_SECTION

START_SECTION([SingleSpectra] strict_imaging_only rejects dried-droplet datasets)
{
  const String d_dd = createFakeD(tmp.getPath(), "rej", "SingleSpectra",
                                   {{1, 0, 0, 0}, {2, 1, 0, 0}});
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  TEST_EXCEPTION(Exception::InvalidValue, f.load(d_dd, exp))
}
END_SECTION

START_SECTION([Multi-section] rejects datasets with multiple distinct ZIndexPos)
{
  // Two slices: z=0 and z=1 -> 2D geometry cannot represent this.
  const String d_multi = createFakeD(tmp.getPath(), "multiz", "Imaging",
                                      {{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 0, 0, 1}, {4, 1, 0, 1}},
                                      /*with_z_column=*/true);
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  TEST_EXCEPTION(Exception::InvalidValue, f.load(d_multi, exp))
}
END_SECTION

END_TEST

#else

#include <OpenMS/CONCEPT/ClassTest.h>
START_TEST(BrukerTimsImagingFile, "$Id$")
END_TEST

#endif // WITH_OPENTIMS
