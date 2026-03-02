// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace OpenMS;

START_TEST(OpenSwathOSWParquetWriter, "$Id$")

OpenSwathOSWParquetWriter writer;

START_SECTION(void write(const String&, const OpenSwath::LightTargetedExperiment&, const FeatureMap&, UInt64, const String&, bool))
  OpenSwath::LightTargetedExperiment lib;
  OpenSwath::LightCompound c1; c1.id = "P1"; c1.charge = 2; c1.rt = 10.0; lib.compounds.push_back(c1);
  OpenSwath::LightCompound c2; c2.id = "7"; c2.charge = 1; c2.rt = 20.0; lib.compounds.push_back(c2);
  OpenSwath::LightTransition t1; t1.transition_name = "tr1"; t1.peptide_ref = "P1"; t1.product_mz = 100.0; lib.transitions.push_back(t1);
  OpenSwath::LightTransition t2; t2.transition_name = "tr2"; t2.peptide_ref = "7"; t2.product_mz = 200.0; lib.transitions.push_back(t2);

  // Build a feature map with two features referencing P1 and 7
  FeatureMap fmap;
  Feature f1; f1.setUniqueId(1); f1.setRT(10.0); f1.setIntensity(100.0); f1.setMetaValue("PeptideRef", "P1");
  Feature sub1; sub1.setRT(10.0); sub1.setIntensity(10.0); f1.getSubordinates().push_back(sub1);
  Feature f2; f2.setUniqueId(2); f2.setRT(20.0); f2.setIntensity(200.0); f2.setMetaValue("PeptideRef", "7");
  Feature sub2; sub2.setRT(20.0); sub2.setIntensity(20.0); f2.getSubordinates().push_back(sub2);
  fmap.push_back(f1);
  fmap.push_back(f2);

  // Use File::TempDir so the temporary directory is removed automatically
  // even if the test aborts early. Create a subdirectory for the .oswpq
  // content so it is contained inside the TempDir.
  File::TempDir tmp_dir;
  String base = tmp_dir.getPath() + "/oswpq";
  if (File::exists(base)) File::removeDirRecursively(base);
  File::makeDir(base);

  // write library and run
  writer.write(base, lib, fmap, 90, "input.mzML", false);

  const String lib_prec = base + "/library/precursors.parquet";
  const String run_features = base + "/runs/run_id=90/features.parquet";
  TEST_EQUAL(File::exists(lib_prec), true)
  TEST_EQUAL(File::exists(run_features), true)

  // read precursor ids from library
  auto prec_table = ParquetFile::readTable(lib_prec);
  auto prec_id_col = ParquetFile::getColumn(prec_table, "precursor_id");
  std::unordered_set<int64_t> lib_ids;
  for (int64_t r = 0; r < prec_table->num_rows(); ++r)
  {
    lib_ids.insert(ParquetFile::getInt64(prec_id_col, r, 0, false));
  }

  // read run feature precursor ids and ensure they exist in library
  auto features_table = ParquetFile::readTable(run_features);
  auto feature_prec_col = ParquetFile::getColumn(features_table, "precursor_id");
  for (int64_t r = 0; r < features_table->num_rows(); ++r)
  {
    int64_t pid = ParquetFile::getInt64(feature_prec_col, r, 0, false);
    TEST_EQUAL(lib_ids.find(pid) != lib_ids.end(), true)
  }

  // TempDir destructor will remove the temporary directory and its contents.
END_SECTION

END_TEST
