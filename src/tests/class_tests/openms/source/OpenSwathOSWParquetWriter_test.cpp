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
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/DATAACCESS/MSChromatogramParquetConsumer.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <unordered_map>
#include <unordered_set>

using namespace OpenMS;

START_TEST(OpenSwathOSWParquetWriter, "$Id$")

OpenSwathOSWParquetWriter writer;

START_SECTION(void write(const std::string&, const OpenSwath::LightTargetedExperiment&, const FeatureMap&, UInt64, const std::string&, bool))
  OpenSwath::LightTargetedExperiment lib;
  OpenSwath::LightCompound c1; c1.id = "0"; c1.charge = 2; c1.rt = 10.0; c1.sequence = "PEPTIDEA"; lib.compounds.push_back(c1);
  OpenSwath::LightCompound c2; c2.id = "7"; c2.charge = 1; c2.rt = 20.0; c2.sequence = "PEPTIDEB"; lib.compounds.push_back(c2);
  OpenSwath::LightTransition t1; t1.transition_name = "3"; t1.peptide_ref = "0"; t1.precursor_mz = 500.0; t1.product_mz = 100.0; lib.transitions.push_back(t1);
  OpenSwath::LightTransition t2; t2.transition_name = "100"; t2.peptide_ref = "7"; t2.precursor_mz = 600.0; t2.product_mz = 200.0; lib.transitions.push_back(t2);

  // Build a feature map with two features referencing canonical precursor IDs 0 and 7.
  FeatureMap fmap;
  Feature f1; f1.setUniqueId(1); f1.setRT(10.0); f1.setIntensity(100.0); f1.setMetaValue("PeptideRef", "0");
  Feature sub1; sub1.setRT(10.0); sub1.setIntensity(10.0); sub1.setMetaValue("FeatureLevel", "MS2"); sub1.setMetaValue("native_id", "3"); f1.getSubordinates().push_back(sub1);
  Feature f2; f2.setUniqueId(2); f2.setRT(20.0); f2.setIntensity(200.0); f2.setMetaValue("PeptideRef", "7");
  Feature sub2; sub2.setRT(20.0); sub2.setIntensity(20.0); sub2.setMetaValue("FeatureLevel", "MS2"); sub2.setMetaValue("native_id", "100"); f2.getSubordinates().push_back(sub2);
  fmap.push_back(f1);
  fmap.push_back(f2);

  // Use File::TempDir so the temporary directory is removed automatically
  // even if the test aborts early. Create a subdirectory for the .oswpq
  // content so it is contained inside the TempDir.
  File::TempDir tmp_dir;
  std::string base = tmp_dir.getPath() + "/oswpq";
  if (File::exists(base)) File::removeDirRecursively(base);
  File::makeDir(base);

  // write library and run
  writer.write(base, lib, fmap, 90, "input.mzML", false);

  const std::string lib_prec = base + "/library/precursors.parquet";
  const std::string run_features = base + "/runs/run_id=90/features.parquet";
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

  TEST_EQUAL(lib_ids.size(), 2)
  TEST_EQUAL(lib_ids.count(0), 1)
  TEST_EQUAL(lib_ids.count(7), 1)

  // The library transition table must preserve zero/sparse canonical IDs too.
  const std::string lib_trans = base + "/library/transitions.parquet";
  auto transition_table = ParquetFile::readTable(lib_trans);
  auto transition_id_col = ParquetFile::getColumn(transition_table, "transition_id");
  auto transition_precursor_col = ParquetFile::getColumn(transition_table, "precursor_id");
  std::unordered_map<int64_t, int64_t> library_transition_to_precursor;
  for (int64_t r = 0; r < transition_table->num_rows(); ++r)
  {
    library_transition_to_precursor.emplace(
      ParquetFile::getInt64(transition_id_col, r, 0, false),
      ParquetFile::getInt64(transition_precursor_col, r, 0, false));
  }
  TEST_EQUAL(library_transition_to_precursor.size(), 2)
  TEST_EQUAL(library_transition_to_precursor.at(3), 0)
  TEST_EQUAL(library_transition_to_precursor.at(100), 7)

  // Run-level transition foreign keys must remain in that same canonical domain.
  const std::string feature_transition_path = base + "/runs/run_id=90/feature_transition.parquet";
  TEST_EQUAL(File::exists(feature_transition_path), true)
  auto feature_transition_table = ParquetFile::readTable(feature_transition_path);
  auto feature_transition_id_col = ParquetFile::getColumn(feature_transition_table, "transition_id");
  std::unordered_set<int64_t> feature_transition_ids;
  for (int64_t r = 0; r < feature_transition_table->num_rows(); ++r)
  {
    feature_transition_ids.insert(ParquetFile::getInt64(feature_transition_id_col, r, -1, false));
  }
  TEST_EQUAL(feature_transition_ids.size(), 2)
  TEST_EQUAL(feature_transition_ids.count(3), 1)
  TEST_EQUAL(feature_transition_ids.count(100), 1)

  // Cross-output regression: XIC and OSWPQ written from the same canonical
  // LightTargetedExperiment must use exactly the same precursor/transition
  // foreign keys. Previously the OSWPQ writer shifted these to a separate
  // one-based domain while XIC preserved zero.
  const std::string xic_path = tmp_dir.getPath() + "/canonical_ids.xic";
  {
    MSChromatogramParquetConsumer xic_writer(xic_path, 90, "input.mzML", lib);
    MSChromatogram tr3;
    tr3.setNativeID("3");
    xic_writer.consumeChromatogram(tr3);
    MSChromatogram tr100;
    tr100.setNativeID("100");
    xic_writer.consumeChromatogram(tr100);
    xic_writer.finalize();
  }

  auto xic_table = ParquetFile::readTable(xic_path);
  auto xic_precursor_col = ParquetFile::getColumn(xic_table, XICSchema::PRECURSOR_ID);
  auto xic_transition_col = ParquetFile::getColumn(xic_table, XICSchema::TRANSITION_ID);
  TEST_EQUAL(xic_table->num_rows(), 2)
  for (int64_t r = 0; r < xic_table->num_rows(); ++r)
  {
    const int64_t precursor_id = ParquetFile::getInt64(xic_precursor_col, r, -1, false);
    const int64_t transition_id = ParquetFile::getInt64(xic_transition_col, r, -1, false);
    TEST_EQUAL(library_transition_to_precursor.count(transition_id), 1)
    TEST_EQUAL(library_transition_to_precursor.at(transition_id), precursor_id)
  }

  // Appending must also reject an existing archive with the same row counts but
  // a different canonical ID domain (for example an archive written by the former
  // one-based allocator). A count-only compatibility check would miss this.
  OpenSwath::LightTargetedExperiment incompatible_lib = lib;
  incompatible_lib.compounds[0].id = "1";
  incompatible_lib.compounds[1].id = "8";
  incompatible_lib.transitions[0].transition_name = "4";
  incompatible_lib.transitions[0].peptide_ref = "1";
  incompatible_lib.transitions[1].transition_name = "101";
  incompatible_lib.transitions[1].peptide_ref = "8";
  OpenSwathLibraryIDNormalizer::validateCanonicalIDs(incompatible_lib);

  const std::string incompatible_base = tmp_dir.getPath() + "/incompatible_oswpq";
  File::makeDir(incompatible_base);
  TransitionParquetFile().convertLightTargetedExperimentToParquet(incompatible_base, incompatible_lib);
  TEST_EXCEPTION(Exception::InvalidValue,
                 writer.write(incompatible_base, lib, fmap, 91, "input.mzML", false))

  // TempDir destructor will remove the temporary directory and its contents.
END_SECTION

END_TEST
