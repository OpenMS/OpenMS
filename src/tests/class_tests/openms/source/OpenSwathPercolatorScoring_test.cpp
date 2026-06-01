// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <sqlite3.h>

#include <memory>
#include <vector>

using namespace OpenMS;

namespace
{
  void copyFixture_(const String& fixture_name, const String& target_path)
  {
    File::remove(target_path);
    TEST_EQUAL(File::copy(OPENMS_GET_TEST_DATA_PATH(fixture_name), target_path), true)
  }

  Int64 querySqliteInt64_(const String& filename, const String& query)
  {
    sqlite3* db = nullptr;
    TEST_EQUAL(sqlite3_open_v2(filename.c_str(), &db, SQLITE_OPEN_READONLY, nullptr), SQLITE_OK)
    if (db == nullptr)
    {
      return 0;
    }

    sqlite3_stmt* stmt = nullptr;
    TEST_EQUAL(sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, nullptr), SQLITE_OK)
    Int64 out = 0;
    if (stmt != nullptr)
    {
      TEST_EQUAL(sqlite3_step(stmt), SQLITE_ROW)
      out = sqlite3_column_int64(stmt, 0);
      sqlite3_finalize(stmt);
    }
    sqlite3_close(db);
    return out;
  }

  struct ExtractedArchive
  {
    String base_dir;
    std::unique_ptr<File::TempDir> temp_dir;
  };

  ExtractedArchive unzipArchive_(const String& archive_path)
  {
    ExtractedArchive archive;
    archive.base_dir = ZipArchiveFile::unzipDirectory(archive_path, archive.temp_dir);
    return archive;
  }

  std::vector<Int64> readRunIds_(const String& base_dir)
  {
    const auto runs_table = ParquetFile::readTable(base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    std::vector<Int64> run_ids;
    run_ids.reserve(runs_table->num_rows());
    for (int64_t row = 0; row < runs_table->num_rows(); ++row)
    {
      run_ids.push_back(ParquetFile::getInt64(run_id_col, row, 0, false));
    }
    return run_ids;
  }

  Size countNonNull_(const std::shared_ptr<arrow::Array>& array)
  {
    return array == nullptr ? 0 : static_cast<Size>(array->length() - array->null_count());
  }

  void checkParquetFeatureColumns_(const String& base_dir,
                                   const std::vector<String>& expected_columns,
                                   const String& non_null_score_column)
  {
    Size total_non_null = 0;
    for (const Int64 run_id : readRunIds_(base_dir))
    {
      const auto table = ParquetFile::readTable(base_dir + "/runs/run_id=" + String(run_id) + "/features.parquet");
      TEST_EQUAL(table->num_rows() > 0, true)
      for (const auto& column : expected_columns)
      {
        TEST_NOT_EQUAL(ParquetFile::getOptionalColumn(table, column), nullptr)
      }
      total_non_null += countNonNull_(ParquetFile::getOptionalColumn(table, non_null_score_column));
    }
    TEST_EQUAL(total_non_null > 0, true)
  }

  void checkParquetTransitionColumns_(const String& base_dir,
                                      const std::vector<String>& expected_columns,
                                      const String& non_null_score_column)
  {
    Size total_non_null = 0;
    for (const Int64 run_id : readRunIds_(base_dir))
    {
      const auto table = ParquetFile::readTable(base_dir + "/runs/run_id=" + String(run_id) + "/feature_transition.parquet");
      TEST_EQUAL(table->num_rows() > 0, true)
      for (const auto& column : expected_columns)
      {
        TEST_NOT_EQUAL(ParquetFile::getOptionalColumn(table, column), nullptr)
      }
      total_non_null += countNonNull_(ParquetFile::getOptionalColumn(table, non_null_score_column));
    }
    TEST_EQUAL(total_non_null > 0, true)
  }
} // namespace

START_TEST(OpenSwathPercolatorScoring, "$Id$")

START_SECTION(OpenSwathPercolatorScoring())
{
  OpenSwathPercolatorScoring* ptr = nullptr;
  ptr = new OpenSwathPercolatorScoring();
  TEST_NOT_EQUAL(ptr, nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(ScoreSummary score(const String& input_path, Level level, const String& output_path))
{
  OpenSwathPercolatorScoring scorer;

  String ms1_osw;
  NEW_TMP_FILE(ms1_osw);
  copyFixture_("PyProphet_inference_test.osw", ms1_osw);
  const auto ms1_summary = scorer.score(ms1_osw, OpenSwathPercolatorScoring::Level::MS1);
  TEST_EQUAL(ms1_summary.total_rows > 0, true)
  TEST_EQUAL(ms1_summary.target_rows > 0, true)
  TEST_EQUAL(ms1_summary.decoy_rows > 0, true)
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1;"),
             querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM FEATURE_MS1;"))
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1 WHERE RANK < 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  String ms2_osw;
  NEW_TMP_FILE(ms2_osw);
  copyFixture_("PyProphet_inference_test.osw", ms2_osw);
  const auto ms2_summary = scorer.score(ms2_osw, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(ms2_summary.total_rows > 0, true)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2;"),
             querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM FEATURE_MS2;"))
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE RANK < 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  const auto transition_summary = scorer.score(ms2_osw, OpenSwathPercolatorScoring::Level::TRANSITION);
  TEST_EQUAL(transition_summary.total_rows > 0, true)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION;") > 0, true)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION WHERE RANK != 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  String ms1ms2_osw;
  NEW_TMP_FILE(ms1ms2_osw);
  copyFixture_("PyProphet_inference_test.osw", ms1ms2_osw);
  const auto ms1ms2_summary = scorer.score(ms1ms2_osw, OpenSwathPercolatorScoring::Level::MS1MS2);
  TEST_EQUAL(ms1ms2_summary.total_rows > 0, true)
  TEST_EQUAL(querySqliteInt64_(ms1ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2;"),
             querySqliteInt64_(ms1ms2_osw, "SELECT COUNT(*) FROM FEATURE_MS2;"))

  String ms1_oswpq;
  NEW_TMP_FILE(ms1_oswpq);
  ms1_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms1_oswpq);
  const auto ms1_oswpq_summary = scorer.score(ms1_oswpq, OpenSwathPercolatorScoring::Level::MS1);
  TEST_EQUAL(ms1_oswpq_summary.total_rows > 0, true)
  const auto ms1_archive = unzipArchive_(ms1_oswpq);
  checkParquetFeatureColumns_(ms1_archive.base_dir,
                              {"score_ms1_score", "score_ms1_peak_group_rank", "score_ms1_qvalue", "score_ms1_pep"},
                              "score_ms1_score");

  String ms2_oswpq;
  NEW_TMP_FILE(ms2_oswpq);
  ms2_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms2_oswpq);
  const auto ms2_oswpq_summary = scorer.score(ms2_oswpq, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(ms2_oswpq_summary.total_rows > 0, true)
  const auto transition_oswpq_summary = scorer.score(ms2_oswpq, OpenSwathPercolatorScoring::Level::TRANSITION);
  TEST_EQUAL(transition_oswpq_summary.total_rows > 0, true)
  const auto transition_archive = unzipArchive_(ms2_oswpq);
  checkParquetFeatureColumns_(transition_archive.base_dir,
                              {"score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_qvalue", "score_ms2_pep"},
                              "score_ms2_score");
  checkParquetTransitionColumns_(transition_archive.base_dir,
                                 {"score_transition_score", "score_transition_rank", "score_transition_qvalue", "score_transition_pep"},
                                 "score_transition_score");

  String ms1ms2_oswpq;
  NEW_TMP_FILE(ms1ms2_oswpq);
  ms1ms2_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms1ms2_oswpq);
  const auto ms1ms2_oswpq_summary = scorer.score(ms1ms2_oswpq, OpenSwathPercolatorScoring::Level::MS1MS2);
  TEST_EQUAL(ms1ms2_oswpq_summary.total_rows > 0, true)
  const auto ms1ms2_archive = unzipArchive_(ms1ms2_oswpq);
  checkParquetFeatureColumns_(ms1ms2_archive.base_dir,
                              {"score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_qvalue", "score_ms2_pep"},
                              "score_ms2_score");
}
END_SECTION

END_TEST
