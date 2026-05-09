// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <sstream>

using namespace OpenMS;

namespace
{
  void copySharedInferenceFixture_(const String& filename)
  {
    File::remove(filename);
    TEST_EQUAL(File::copy(OPENMS_GET_TEST_DATA_PATH("PyProphet_inference_test.osw"), filename), true)
  }

  Size countRows_(const String& filename, const String& table_name)
  {
    SqliteConnector conn(filename, SqliteConnector::SqlOpenMode::READ_ONLY);
    return conn.countTableRows(table_name);
  }

  Size countFilteredRows_(const String& filename, const String& table_name, const String& where_clause)
  {
    SqliteConnector conn(filename, SqliteConnector::SqlOpenMode::READWRITE);
    const String temp_view = "__oswfile_inference_test_count_view";
    conn.executeStatement("DROP VIEW IF EXISTS " + temp_view + ";");
    conn.executeStatement("CREATE TEMP VIEW " + temp_view + " AS SELECT * FROM " + table_name + " WHERE " + where_clause + ";");
    return conn.countTableRows(temp_view);
  }

  void createLargeTransitionOSW_(const String& filename)
  {
    SqliteConnector conn(filename);
    conn.executeStatement(
      "CREATE TABLE TRANSITION(ID INTEGER PRIMARY KEY NOT NULL, TYPE TEXT NOT NULL, DECOY INTEGER NOT NULL);"
      "CREATE TABLE SCORE_TRANSITION(FEATURE_ID INTEGER NOT NULL, TRANSITION_ID INTEGER NOT NULL, PEP DOUBLE NOT NULL);"
      "CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INTEGER NOT NULL, PEPTIDE_ID INTEGER NOT NULL);"
    );

    std::ostringstream transition_sql;
    transition_sql << "INSERT INTO TRANSITION (ID, TYPE, DECOY) VALUES ";
    for (Int id = 1; id <= 6; ++id)
    {
      if (id > 1)
      {
        transition_sql << ", ";
      }
      transition_sql << "(" << (1000 + id) << ", 'y', 0)";
    }
    transition_sql << ";";
    conn.executeStatement(transition_sql.str());

    conn.executeStatement(
      "INSERT INTO TRANSITION_PEPTIDE_MAPPING (TRANSITION_ID, PEPTIDE_ID) VALUES "
      "(1001, 101), (1002, 101), (1003, 101), (1004, 102), (1005, 102), (1006, 102);"
    );

    std::ostringstream score_sql;
    score_sql << "INSERT INTO SCORE_TRANSITION (FEATURE_ID, TRANSITION_ID, PEP) VALUES ";
    bool first = true;
    for (Int feature_id = 1; feature_id <= 40; ++feature_id)
    {
      const double pep = feature_id <= 20 ? 0.2 : 0.8;
      for (Int transition_id = 1001; transition_id <= 1006; ++transition_id)
      {
        if (!first)
        {
          score_sql << ", ";
        }
        first = false;
        score_sql << "(" << feature_id << ", " << transition_id << ", " << pep << ")";
      }
    }
    score_sql << ";";
    conn.executeStatement(score_sql.str());
  }
} // namespace

START_TEST(OSWFileInference, "$Id$")

START_SECTION(OSWFile IPF readers)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedInferenceFixture_(tmp_osw);

  OSWFile osw(tmp_osw);

  PeptidoformInferenceConfig config;
  config.ipf_ms1_scoring = false;
  config.ipf_ms2_scoring = false;
  config.ipf_max_transition_pep = 0.7;
  const auto precursor_rows = osw.readIPFPrecursorData(config);
  TEST_EQUAL(precursor_rows.size(), 389)
  auto feature1 = std::find_if(precursor_rows.begin(), precursor_rows.end(), [](const auto& row) { return row.feature_id == -9078977811506172301LL; });
  TEST_EQUAL(feature1 != precursor_rows.end(), true)
  TEST_REAL_SIMILAR(feature1->ms2_peakgroup_pep, 0.00314231908146708)
  TEST_EQUAL(feature1->ms1_precursor_pep.has_value(), false)
  TEST_EQUAL(feature1->ms2_precursor_pep.has_value(), false)

  const auto transition_rows = osw.readIPFTransitionData(config);
  TEST_EQUAL(transition_rows.size(), 9431)
  auto support_row = std::find_if(transition_rows.begin(), transition_rows.end(),
                                  [](const auto& row) { return row.feature_id == -9078977811506172301LL && row.transition_id == 3275 && row.peptide_id == 305; });
  auto nonsupport_row = std::find_if(transition_rows.begin(), transition_rows.end(),
                                     [](const auto& row) { return row.feature_id == -9078977811506172301LL && row.transition_id == 3275 && row.peptide_id == -1; });
  TEST_EQUAL(support_row != transition_rows.end(), true)
  TEST_EQUAL(nonsupport_row != transition_rows.end(), true)
  TEST_EQUAL(support_row->bmask, 1)
  TEST_EQUAL(nonsupport_row->bmask, 0)
  TEST_EQUAL(support_row->num_peptidoforms, 1)
  // The shared PyProphet test fixture is single-run and intentionally does not
  // contain alignment tables, so across-run propagation is covered elsewhere.
}
END_SECTION

START_SECTION(OSWFile level-context readers)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedInferenceFixture_(tmp_osw);

  OSWFile osw(tmp_osw);

  const auto peptide_global = osw.readLevelContextData(InferenceLevel::Peptide, InferenceContext::Global);
  TEST_EQUAL(peptide_global.size(), 682)
  for (const auto& row : peptide_global)
  {
    TEST_EQUAL(row.run_id.has_value(), false)
  }

  const auto peptide_run = osw.readLevelContextData(InferenceLevel::Peptide, InferenceContext::RunSpecific);
  TEST_EQUAL(peptide_run.size(), 682)
  for (const auto& row : peptide_run)
  {
    TEST_EQUAL(row.run_id.has_value(), true)
  }

  const auto protein_global = osw.readLevelContextData(InferenceLevel::Protein, InferenceContext::Global);
  TEST_EQUAL(protein_global.size(), 32)

  const auto gene_global = osw.readLevelContextData(InferenceLevel::Gene, InferenceContext::Global);
  TEST_EQUAL(gene_global.size(), 32)
}
END_SECTION

START_SECTION(OSWFile writers)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedInferenceFixture_(tmp_osw);

  String out_osw;
  NEW_TMP_FILE(out_osw);
  File::remove(out_osw);

  OSWFile osw(tmp_osw);
  osw.writeIPFResults(out_osw, {{1, 101, 0.2, 0.1, 0.1}, {5, 102, 0.3, 0.2, 0.2}});
  TEST_EQUAL(countRows_(out_osw, "SCORE_IPF"), 2)

  osw.writeIPFResults(out_osw, {{1, 101, 0.25, 0.15, 0.12}});
  TEST_EQUAL(countRows_(out_osw, "SCORE_IPF"), 1)

  OSWFile out_writer(out_osw);
  out_writer.writeLevelContextResults("", InferenceLevel::Peptide, InferenceContext::RunSpecific,
                                      {{1, 101, InferenceContext::RunSpecific, 5.0, 0.1, 0.1, 0.1}});
  out_writer.writeLevelContextResults("", InferenceLevel::Peptide, InferenceContext::Global,
                                      {{std::nullopt, 101, InferenceContext::Global, 6.0, 0.2, 0.2, 0.2}});
  TEST_EQUAL(countRows_(out_osw, "SCORE_PEPTIDE"), 2)

  out_writer.writeLevelContextResults("", InferenceLevel::Peptide, InferenceContext::Global,
                                      {
                                        {std::nullopt, 101, InferenceContext::Global, 6.5, 0.25, 0.25, 0.25},
                                        {std::nullopt, 102, InferenceContext::Global, 5.5, 0.35, 0.35, 0.35}
                                      });
  TEST_EQUAL(countRows_(out_osw, "SCORE_PEPTIDE"), 3)
  TEST_EQUAL(countFilteredRows_(out_osw, "SCORE_PEPTIDE", "CONTEXT = 'run-specific'"), 1)
  TEST_EQUAL(countFilteredRows_(out_osw, "SCORE_PEPTIDE", "CONTEXT = 'global'"), 2)
}
END_SECTION

START_SECTION(OSWFile transition reader scales with thresholded evidence)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  createLargeTransitionOSW_(tmp_osw);

  OSWFile osw(tmp_osw);
  PeptidoformInferenceConfig config;
  config.ipf_max_transition_pep = 0.5;
  config.ipf_h0 = true;

  const auto transition_rows = osw.readIPFTransitionData(config);
  TEST_EQUAL(transition_rows.size(), 360)
  TEST_EQUAL(std::count_if(transition_rows.begin(), transition_rows.end(),
                           [](const auto& row) { return row.feature_id > 20; }),
             0)
  TEST_EQUAL(std::count_if(transition_rows.begin(), transition_rows.end(),
                           [](const auto& row) { return row.peptide_id == -1; }),
             120)
}
END_SECTION

END_TEST
