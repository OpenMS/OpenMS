// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathOSWParquetRoundTrip, "$Id$")

START_SECTION(void round-trip write/read .oswpq archive using RAF path)
{
  // Build a reference LightTargetedExperiment from a TraML file
  const String input_file = OPENMS_GET_TEST_DATA_PATH("MRMAssay_detectingTransistionCompound_input.TraML");
  TraMLFile traml;
  TargetedExperiment targeted_exp;
  traml.load(input_file, targeted_exp);

  OpenSwath::LightTargetedExperiment light_exp;
  OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, light_exp);
  TEST_EQUAL(light_exp.compounds.size() > 0, true)

  // Write to a single .oswpq archive (do NOT create a directory) to exercise archive writer path
  File::TempDir tmp_dir;
  const String out_archive = tmp_dir.getPath() + "/roundtrip.oswpq";

  OpenSwathOSWParquetWriter writer;
  FeatureMap empty_map;
  writer.write(out_archive, light_exp, empty_map, 1, String("test_input"), false);

  // Archive and embedded sidecar should exist (sidecar is written inside the zip)
  TEST_EQUAL(File::exists(out_archive), true)
  {
    auto entries = ZipArchiveFile::listEntries(out_archive);
    bool found = false;
    for (const auto& e : entries) if (e == "roundtrip.oswpq.idx.json") { found = true; break; }
    TEST_EQUAL(found, true)
  }

  // Verify the RAF path works: ZipRandomAccessFile::Open should succeed directly on the archive
  {
    std::unique_ptr<File::TempDir> raf_tmp;
    auto ra_res = ZipRandomAccessFile::Open(out_archive, "library/precursors.parquet", raf_tmp);
#if __has_include(<zip.h>)
    TEST_EQUAL(ra_res.ok(), true)
#else
    TEST_EQUAL(ra_res.status().IsNotImplemented(), true)
#endif
  }

  // Read back using TransitionParquetFile and verify the round-trip data
  TransitionParquetFile reader;
  OpenSwath::LightTargetedExperiment roundtrip_exp;
  reader.convertParquetToTargetedExperiment(out_archive, roundtrip_exp);

  TEST_EQUAL(roundtrip_exp.compounds.size(), light_exp.compounds.size())
  TEST_EQUAL(roundtrip_exp.transitions.size(), light_exp.transitions.size())
}
END_SECTION

END_TEST
