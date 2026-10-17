// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/config.h>

#ifdef WITH_OPENTIMS

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/BrukerTimsFile.h>
#include <OpenMS/FORMAT/BrukerTimsImagingFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <filesystem>
#include <fstream>

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
  std::string createFakeD(const std::string& parent,
                     const std::string& name,
                     const std::string& maldi_application_type,
                     const std::vector<std::tuple<int,int,int,int>>& xyz_rows,
                     bool with_z_column = false,
                     const std::string& spot_size_um = "")
  {
    const std::string d_path = parent + "/" + name + ".d";
    std::filesystem::create_directories(d_path.c_str());
    const std::string tdf_path = d_path + "/analysis.tdf";

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
      std::string sql = "INSERT INTO MaldiFrameInfo VALUES (" + StringUtils::toStr(frame) + ", "
                   + StringUtils::toStr(x) + ", " + StringUtils::toStr(y);
      if (with_z_column) { sql += ", " + StringUtils::toStr(z); }
      sql += ");";
      conn.executeStatement(sql);
    }

    return d_path;
  }
}

START_TEST(BrukerTimsImagingFile, "$Id$")

TempDir tmp(false);

START_SECTION(static bool isImagingDataset(const String& path))
{
  // Imaging dataset
  const std::string d_img = createFakeD(tmp.getPath(), "imaging", "Imaging",
                                    {{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 1, 1, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_img), true)

  // SingleSpectra dataset (dried droplet)
  const std::string d_dd  = createFakeD(tmp.getPath(), "droplet", "SingleSpectra",
                                    {{1, 0, 0, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_dd), false)

  // Missing folder
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(tmp.getPath() + "/does_not_exist.d"), false)
}
END_SECTION

START_SECTION(static String readGlobalMetadataValue(const String& d_folder, const String& key))
{
  const std::string d = createFakeD(tmp.getPath(), "gm", "Imaging",
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
  const std::string d = createFakeD(tmp.getPath(), "grid2x2", "Imaging",
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
  const std::string d_with_z = createFakeD(tmp.getPath(), "withz", "Imaging",
                                       {{1, 3, 4, 0}, {2, 3, 5, 0}}, /*with_z_column=*/true);
  auto rows_z = BrukerTimsImagingFile::readMaldiFrameInfo(d_with_z);
  TEST_EQUAL(rows_z.size(), 2)
  TEST_EQUAL(rows_z[0].x, 3)
  TEST_EQUAL(rows_z[0].y, 4)
  TEST_EQUAL(rows_z[1].y, 5)

  // Empty MaldiFrameInfo -> ParseError
  const std::string d_empty = createFakeD(tmp.getPath(), "empty", "Imaging", {});
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
  const std::string d_no_tdf = tmp.getPath() + "/no_tdf.d";
  std::filesystem::create_directories(d_no_tdf.c_str());
  TEST_EXCEPTION(Exception::FileNotReadable, f.load(d_no_tdf, exp))

  // TSF data (MALDI without ion mobility, written for most MALDI imaging runs) cannot be read: said
  // as such (BrukerTimsFile::checkDataFormat()), not as a missing analysis.tdf.
  const std::string d_tsf = tmp.getPath() + "/tsf.d";
  std::filesystem::create_directories(d_tsf.c_str());
  { std::ofstream os((d_tsf + "/analysis.tsf").c_str()); os << "no data"; }
  TEST_EXCEPTION(Exception::ParseError, f.load(d_tsf, exp))
}
END_SECTION

START_SECTION([SingleSpectra] strict_imaging_only rejects dried-droplet datasets)
{
  const std::string d_dd = createFakeD(tmp.getPath(), "rej", "SingleSpectra",
                                   {{1, 0, 0, 0}, {2, 1, 0, 0}});
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  TEST_EXCEPTION(Exception::InvalidValue, f.load(d_dd, exp))
}
END_SECTION

START_SECTION([EXTRA] '.d.zip' archives)
{
  // The metadata is read from the analysis.tdf of the .d directory that BrukerTimsFile reads from the archive
  // (#10596): it is extracted alone, so the fake .d folders, which hold no frame data, suffice.
  // Zips the directory @p dir (holding .d folders) into an archive next to it
  auto zip = [](const std::string& dir)
  {
    const std::string archive = dir + ".d.zip";
    ZipArchiveFile::zipDirectory(dir, archive);
    return archive;
  };
  const std::string img_dir = tmp.getPath() + "/zip_imaging";
  createFakeD(img_dir, "run", "Imaging", {{1, 0, 0, 0}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 1, 1, 0}}, false, "50");
  const std::string img_zip = zip(img_dir);
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(img_zip), true)
  TEST_EQUAL(BrukerTimsImagingFile::readGlobalMetadataValue(img_zip, "MaldiSpotSize_um"), "50")
  TEST_EQUAL(BrukerTimsImagingFile::readMaldiFrameInfo(img_zip).size(), 4)

  const std::string dd_dir = tmp.getPath() + "/zip_droplet";
  createFakeD(dd_dir, "run", "SingleSpectra", {{1, 0, 0, 0}});
  const std::string dd_zip = zip(dd_dir);
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(dd_zip), false)
  // load() checks the metadata in the archive before any frame is read
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  TEST_EXCEPTION(Exception::InvalidValue, f.load(dd_zip, exp))
  const std::string multi_dir = tmp.getPath() + "/zip_multiz";
  createFakeD(multi_dir, "run", "Imaging", {{1, 0, 0, 0}, {2, 0, 0, 1}}, /*with_z_column=*/true);
  TEST_EXCEPTION(Exception::InvalidValue, f.load(zip(multi_dir), exp))

  // As BrukerTimsFile: the shallowest .d (the first by name among equally deep ones), none below __MACOSX
  const std::string nested_dir = tmp.getPath() + "/zip_nested";
  createFakeD(nested_dir + "/__MACOSX", "a", "SingleSpectra", {{1, 0, 0, 0}});
  createFakeD(nested_dir + "/outer", "b", "Imaging", {{1, 0, 0, 0}});
  createFakeD(nested_dir + "/outer", "c", "SingleSpectra", {{1, 0, 0, 0}});
  createFakeD(nested_dir + "/outer/deeper", "a", "SingleSpectra", {{1, 0, 0, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(zip(nested_dir)), true) // outer/b.d
  const std::string top_dir = tmp.getPath() + "/zip_top";
  createFakeD(top_dir, "z", "SingleSpectra", {{1, 0, 0, 0}});
  createFakeD(top_dir + "/outer", "a", "Imaging", {{1, 0, 0, 0}});
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(zip(top_dir)), false) // z.d
  // An entry './a.d/analysis.tdf' is extracted to a.d, as deep as b.d, so a.d (the first by name) is read
  const std::string dot_dir = tmp.getPath() + "/zip_dot";
  const std::string dot_a = createFakeD(dot_dir, "a", "Imaging", {{1, 0, 0, 0}, {2, 1, 0, 0}});
  const std::string dot_b = createFakeD(dot_dir, "b", "SingleSpectra", {{1, 0, 0, 0}});
  const std::string dot_zip = dot_dir + ".d.zip";
  ZipArchiveFile::addOrReplaceFromFile(dot_zip, "./a.d/analysis.tdf", dot_a + "/analysis.tdf");
  ZipArchiveFile::addOrReplaceFromFile(dot_zip, "b.d/analysis.tdf", dot_b + "/analysis.tdf");
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(dot_zip), true)
  TEST_EQUAL(BrukerTimsImagingFile::readMaldiFrameInfo(dot_zip).size(), 2)
  // BrukerTimsFile, which reads the frames, takes a.d as well: it refuses the imaging run outside export mode FRAME
  {
    BrukerTimsFile frames;
    MSExperiment frames_exp;
    TEST_EXCEPTION(Exception::InvalidFileType, frames.load(dot_zip, frames_exp))
  }

  // No .d directory with analysis.tdf in the archive
  const std::string empty_dir = tmp.getPath() + "/zip_empty";
  File::makeDir(empty_dir);
  { std::ofstream os((empty_dir + "/readme.txt").c_str()); os << "no Bruker data"; }
  const std::string empty_zip = zip(empty_dir);
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(empty_zip), false)
  TEST_EXCEPTION(Exception::ParseError, BrukerTimsImagingFile::readGlobalMetadataValue(empty_zip, "MaldiApplicationType"))
  TEST_EXCEPTION(Exception::ParseError, f.load(empty_zip, exp))
  TEST_EXCEPTION(Exception::FileNotReadable, f.load(tmp.getPath() + "/missing.d.zip", exp))
}
END_SECTION

START_SECTION([Multi-section] rejects datasets with multiple distinct ZIndexPos)
{
  // Two slices: z=0 and z=1 -> 2D geometry cannot represent this.
  const std::string d_multi = createFakeD(tmp.getPath(), "multiz", "Imaging",
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
