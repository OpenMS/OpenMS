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

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathOSWParquetRoundTrip, "$Id$")

START_SECTION(void round-trip write/read .oswpq archive using RAF path)
{
#ifdef WITH_PARQUET
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

  // Archive and sidecar should exist
  TEST_EQUAL(File::exists(out_archive), true)
  TEST_EQUAL(File::exists(out_archive + ".idx.json"), true)

  // Read back using TransitionParquetFile which should prefer RAF-based reads when available
  // Reset test extraction counter and assert no extraction occurs (i.e., RAF path used)
  ZipArchiveFile::testResetExtractionCount();
  TransitionParquetFile reader;
  OpenSwath::LightTargetedExperiment roundtrip_exp;
  reader.convertParquetToTargetedExperiment(out_archive, roundtrip_exp);

  TEST_EQUAL(ZipArchiveFile::testGetExtractionCount(), 0)

  TEST_EQUAL(roundtrip_exp.compounds.size(), light_exp.compounds.size())
  TEST_EQUAL(roundtrip_exp.transitions.size(), light_exp.transitions.size())
#else
  NOT_TESTABLE
#endif
}
END_SECTION

END_TEST
