// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h>
///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/FORMAT/ParquetFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathParquetExporter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void writeFeatureScores(const std::string& filename, const OpenSwathFeatureScoreTable& table)))
{
  OpenSwathFeatureScoreRow r;
  r.protein_id = 1; r.peptide_id = 2; r.precursor_id = 3;
  r.protein_accession = "PROT1";
  r.unmodified_sequence = "PEPTIDEK";
  r.modified_sequence = "PEPTIDEK(UniMod:1)";
  r.precursor_mz = 450.7; r.precursor_charge = 2; r.run_id = 99; r.feature_id = 7777;
  r.exp_rt = 120.5; r.norm_rt = 118.0; r.delta_rt = 2.5; r.left_width = 118.0; r.right_width = 123.0;

  OpenSwathFeatureScoreRow r2 = r;
  r2.feature_id = 8888;
  r2.modified_sequence = "ACDEFGHIK";

  OpenSwathFeatureScoreTable table; // no dynamic FEATURE_MS1/MS2 score columns
  table.rows = {r, r2};

  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  OpenSwathParquetExporter::writeFeatureScores(fn, table);

  // read the parquet back and check the round-tripped key columns
  auto tbl = ParquetFile::readTable(fn);
  TEST_NOT_EQUAL(tbl, nullptr)
  TEST_EQUAL(tbl->num_rows(), 2)

  auto fid = ParquetFile::getColumn(tbl, "FEATURE_ID");
  TEST_EQUAL(ParquetFile::getInt64(fid, 0, -1, false), 7777)
  TEST_EQUAL(ParquetFile::getInt64(fid, 1, -1, false), 8888)
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "PROTEIN_ACCESSION"), 0), "PROT1")
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "MODIFIED_SEQUENCE"), 0), "PEPTIDEK(UniMod:1)")
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "MODIFIED_SEQUENCE"), 1), "ACDEFGHIK")
  TEST_REAL_SIMILAR(ParquetFile::getDouble(ParquetFile::getColumn(tbl, "PRECURSOR_MZ"), 0, -1, false), 450.7)
  TEST_REAL_SIMILAR(ParquetFile::getDouble(ParquetFile::getColumn(tbl, "EXP_RT"), 0, -1, false), 120.5)
}
END_SECTION

START_SECTION((static void writeTransitionScores(const std::string& filename, const OpenSwathTransitionScoreTable& table)))
{
  OpenSwathTransitionScoreRow r;
  r.precursor_id = 5; r.transition_id = 42; r.transition_traml_id = "tr42";
  r.product_mz = 500.25; r.transition_charge = 1; r.transition_type = "y";
  r.transition_ordinal = 3; r.annotation = "y3";
  r.transition_detecting = true; r.transition_decoy = false;

  OpenSwathTransitionScoreRow r2 = r;
  r2.transition_id = 43; r2.annotation = "y4"; r2.transition_decoy = true;

  OpenSwathTransitionScoreTable table;
  table.rows = {r, r2};

  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  OpenSwathParquetExporter::writeTransitionScores(fn, table);

  auto tbl = ParquetFile::readTable(fn);
  TEST_NOT_EQUAL(tbl, nullptr)
  TEST_EQUAL(tbl->num_rows(), 2)

  auto tid = ParquetFile::getColumn(tbl, "TRANSITION_ID");
  TEST_EQUAL(ParquetFile::getInt64(tid, 0, -1, false), 42)
  TEST_EQUAL(ParquetFile::getInt64(tid, 1, -1, false), 43)
  TEST_REAL_SIMILAR(ParquetFile::getDouble(ParquetFile::getColumn(tbl, "PRODUCT_MZ"), 0, -1, false), 500.25)
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "ANNOTATION"), 1), "y4")
  // decoy flag round-trips per row
  auto dec = ParquetFile::getColumn(tbl, "TRANSITION_DECOY");
  TEST_EQUAL(ParquetFile::getBool(dec, 0, true, false), false)
  TEST_EQUAL(ParquetFile::getBool(dec, 1, false, false), true)

  // an empty table is allowed (schema-only output) and writes a 0-row file
  OpenSwathTransitionScoreTable empty;
  std::string fn2;
  NEW_TMP_FILE(fn2)
  fn2 += ".parquet";
  OpenSwathParquetExporter::writeTransitionScores(fn2, empty);
  TEST_EQUAL(ParquetFile::rowCount(fn2), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
