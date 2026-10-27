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
#include <OpenMS/FORMAT/ZstdCompression.h>
#include <OpenMS/IMAGING/MSImagingExperiment.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>

using namespace OpenMS;
using namespace std;

namespace
{
  // Build a fake .d directory with only an analysis.tdf SQLite file inside.
  // Optional MaldiApplicationType + MaldiFrameInfo rows.
  // Returns the path to the .d folder (a subdir of @p parent).
  // Build a fake .d directory with only an analysis.tdf SQLite file inside
  // (@p db_file "analysis.tsf": TSF data, timsTOF without ion mobility).
  // Each row in @p xyz_rows is (frame, x, y, z); pass z<0 to omit the
  // ZIndexPos column from the schema entirely.
  std::string createFakeD(const std::string& parent,
                     const std::string& name,
                     const std::string& maldi_application_type,
                     const std::vector<std::tuple<int,int,int,int>>& xyz_rows,
                     bool with_z_column = false,
                     const std::string& spot_size_um = "",
                     const std::string& db_file = "analysis.tdf")
  {
    const std::string d_path = parent + "/" + name + ".d";
    std::filesystem::create_directories(d_path.c_str());
    const std::string tdf_path = d_path + "/" + db_file;

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

  // Adds frames 1..@p n_frames to the TSF .d directory @p d_path (made by createFakeD() with "analysis.tsf"): the
  // MzCalibration and Frames tables, and analysis.tsf_bin with one peak per frame, of intensity 100 * frame id, at a TOF
  // index whose m/z is 379.09864007179004 by the Bruker SDK (MassIVE MSV000088438, rs17e_1mgml_no_tims_1_0_C2_MS.d)
  void addTsfFrames(const std::string& d_path, int n_frames)
  {
    SqliteConnector conn(d_path + "/analysis.tsf", SqliteConnector::SqlOpenMode::READWRITE);
    conn.executeStatement("INSERT INTO GlobalMetadata (Key, Value) VALUES ('TimsCompressionType', '3');");
    conn.executeStatement("CREATE TABLE MzCalibration (Id INTEGER PRIMARY KEY, ModelType INTEGER, DigitizerTimebase REAL, "
                          "DigitizerDelay REAL, T1 REAL, T2 REAL, dC1 REAL, dC2 REAL, C0 REAL, C1 REAL, C2 REAL, C3 REAL, C4 REAL);");
    conn.executeStatement("INSERT INTO MzCalibration VALUES (1, 1, 0.19999999999999998, 25627.6, 25.661657387366994, "
                          "27.37440317917224, 20.0, 0.0, 314.2024623155012, 156059.18973537948, -0.0020587980315838962, 0.0, 0.0);");
    conn.executeStatement("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, Polarity CHAR(1), MsMsType INTEGER, TimsId INTEGER, "
                          "NumPeaks INTEGER, MzCalibration INTEGER, T1 REAL, T2 REAL);");
    std::ofstream bin((d_path + "/analysis.tsf_bin").c_str(), std::ios::binary);
    uint64_t offset = 0;
    for (int id = 1; id <= n_frames; ++id)
    {
      const double index = 119863.65367102693;
      const float intensity = 100.0f * static_cast<float>(id);
      std::string line(sizeof(double) + sizeof(float), '\0');
      std::memcpy(&line[0], &index, sizeof(double));
      std::memcpy(&line[sizeof(double)], &intensity, sizeof(float));
      std::string line_zstd;
      ZstdCompression::compressData(line.data(), line.size(), line_zstd);
      const uint32_t header[2] = {static_cast<uint32_t>(8 + line_zstd.size()), static_cast<uint32_t>(line_zstd.size())};
      bin.write(reinterpret_cast<const char*>(header), sizeof(header));
      bin.write(line_zstd.data(), line_zstd.size());
      conn.executeStatement("INSERT INTO Frames VALUES (" + std::to_string(id) + ", " + std::to_string(id) + ".5, '+', 0, " +
                            std::to_string(offset) + ", 1, 1, 25.69411772046709, 27.491792202523506);");
      offset += header[0];
    }
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

  // TSF data (MALDI without ion mobility, written for most MALDI imaging runs): the pixels of MaldiFrameInfo in
  // analysis.tsf, a spectrum per frame (one peak each, see addTsfFrames()), without ion mobility
  const std::string d_tsf = createFakeD(tmp.getPath(), "tsf_load", "Imaging", {{1, 12, 11, 0}, {2, 13, 11, 0}, {3, 12, 12, 0}},
                                        false, "20", "analysis.tsf");
  addTsfFrames(d_tsf, 3);
  f.load(d_tsf, exp);
  TEST_EQUAL(exp.getMSExperiment().size(), 3)
  TEST_EQUAL(exp.getGeometry().getWidth(), 2)
  TEST_EQUAL(exp.getGeometry().getHeight(), 2)
  TEST_REAL_SIMILAR(exp.getGeometry().getPixelSizeX(), 20.0)
  const auto& pixels = exp.getGeometry().getPixels();
  TEST_EQUAL(pixels.size(), 3)
  for (const auto& pixel : pixels)
  {
    const MSSpectrum& spec = exp.getMSExperiment()[pixel.spectrum_index];
    TEST_EQUAL(spec.size(), 1)
    TEST_EQUAL(spec.containsIMData(), false)
    if (spec.size() != 1) continue;
    TEST_TRUE(std::abs(spec[0].getMZ() - 379.09864007179004) / 379.09864007179004 * 1e6 < 1e-4) // ppm
    // frame 3 is pixel (12, 12), at (0, 1) after the shift to zero; frame 2 is (13, 11) -> (1, 0)
    const std::string frame = spec.getNativeID();
    if (frame == "frame=3") TEST_TRUE(pixel.x == 0 && pixel.y == 1 && spec[0].getIntensity() == 300.0)
    if (frame == "frame=2") TEST_TRUE(pixel.x == 1 && pixel.y == 0 && spec[0].getIntensity() == 200.0)
  }
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

START_SECTION([EXTRA] TSF data: the metadata is read from analysis.tsf)
{
  // TSF data (timsTOF without ion mobility, as most MALDI imaging runs) has the GlobalMetadata and MaldiFrameInfo
  // tables of TDF data in analysis.tsf
  const std::string d = createFakeD(tmp.getPath(), "tsf_imaging", "Imaging", {{1, 12, 11, 0}, {2, 13, 11, 0}, {3, 12, 12, 0}},
                                    false, "20", "analysis.tsf");
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d), true)
  TEST_EQUAL(BrukerTimsImagingFile::readGlobalMetadataValue(d, "MaldiSpotSize_um"), "20")
  auto rows = BrukerTimsImagingFile::readMaldiFrameInfo(d);
  TEST_EQUAL(rows.size(), 3)
  TEST_EQUAL(rows[2].frame_id, 3)
  TEST_EQUAL(rows[2].x, 12)
  TEST_EQUAL(rows[2].y, 12)
  const std::string d_dd = createFakeD(tmp.getPath(), "tsf_droplet", "SingleSpectra", {{1, 0, 0, 0}}, false, "", "analysis.tsf");
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_dd), false)

  // A .d with both: analysis.tdf is read, as BrukerTimsFile::detectFormat() reports TDF data for it
  const std::string d_both = createFakeD(tmp.getPath(), "tdf_and_tsf", "SingleSpectra", {{1, 0, 0, 0}});
  std::filesystem::copy_file(d + "/analysis.tsf", d_both + "/analysis.tsf");
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(d_both), false)

  // '.d.zip': the shallowest .d with TSF data, as for TDF data
  const std::string zip_dir = tmp.getPath() + "/zip_tsf";
  createFakeD(zip_dir + "/outer", "a", "SingleSpectra", {{1, 0, 0, 0}}, false, "", "analysis.tsf");
  createFakeD(zip_dir, "b", "Imaging", {{1, 0, 0, 0}, {2, 1, 0, 0}}, false, "", "analysis.tsf");
  const std::string archive = zip_dir + ".d.zip";
  ZipArchiveFile::zipDirectory(zip_dir, archive);
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(archive), true) // b.d
  TEST_EQUAL(BrukerTimsImagingFile::readMaldiFrameInfo(archive).size(), 2)
  // An archive with TDF data in one .d and TSF data in another: the TDF data is read, even from a deeper .d
  const std::string mixed_dir = tmp.getPath() + "/zip_tdf_tsf";
  createFakeD(mixed_dir, "a", "SingleSpectra", {{1, 0, 0, 0}}, false, "", "analysis.tsf");
  createFakeD(mixed_dir + "/outer", "b", "Imaging", {{1, 0, 0, 0}});
  const std::string mixed = mixed_dir + ".d.zip";
  ZipArchiveFile::zipDirectory(mixed_dir, mixed);
  TEST_EQUAL(BrukerTimsImagingFile::isImagingDataset(mixed), true) // outer/b.d
}
END_SECTION

#ifdef OPENTIMS_TSF_IMAGING_TEST_DATA
START_SECTION([EXTRA] TSF data: a MALDI imaging run, as the Bruker SDK reads it)
{
  // MassIVE MSV000088438, 20210920_vc_rugose_1mMTCA_gordon.d (timsTOF fleX). TIMSCONVERT 2.0.0 (Bruker SDK) gives 276
  // spectra with 74188 peaks, at the XIndexPos/YIndexPos of MaldiFrameInfo (x 12..39, y 11..34); frame 1 has 478 peaks,
  // the most intense at m/z 379.0929137613323 (335138 counts), and is spot R00X019Y012.
  BrukerTimsImagingFile f;
  MSImagingExperiment exp;
  f.load(OPENTIMS_TSF_IMAGING_TEST_DATA, exp);
  TEST_EQUAL(exp.getMSExperiment().size(), 276)
  TEST_EQUAL(exp.getGeometry().getPixels().size(), 276)
  TEST_EQUAL(exp.getGeometry().getWidth(), 28)
  TEST_EQUAL(exp.getGeometry().getHeight(), 24)
  Size peaks = 0;
  for (const MSSpectrum& spec : exp.getMSExperiment()) peaks += spec.size();
  TEST_EQUAL(peaks, 74188)
  for (const auto& pixel : exp.getGeometry().getPixels())
  {
    const MSSpectrum& spec = exp.getMSExperiment()[pixel.spectrum_index];
    if (spec.getNativeID() != "frame=1") continue;
    TEST_TRUE(pixel.x == 19 - 12 && pixel.y == 12 - 11)
    TEST_STRING_EQUAL(spec.getMetaValue("maldi_spot_id").toString(), "R00X019Y012")
    TEST_EQUAL(spec.size(), 478)
    if (spec.empty()) continue;
    const auto top = std::max_element(spec.begin(), spec.end(),
      [](const Peak1D& a, const Peak1D& b) { return a.getIntensity() < b.getIntensity(); });
    TEST_REAL_SIMILAR(top->getIntensity(), 335138.0)
    TEST_TRUE(std::abs(top->getMZ() - 379.0929137613323) / 379.0929137613323 * 1e6 < 1e-4) // ppm
  }
}
END_SECTION
#endif

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
