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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <sqlite3.h>

#include <memory>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  void copyFixture_(const std::string& fixture_name, const std::string& target_path)
  {
    File::remove(target_path);
    TEST_TRUE(File::copy(OPENMS_GET_TEST_DATA_PATH(fixture_name), target_path))
  }

  Int64 querySqliteInt64_(const std::string& filename, const std::string& query)
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

  void executeSqlite_(const std::string& filename, const std::string& query)
  {
    sqlite3* db = nullptr;
    TEST_EQUAL(sqlite3_open_v2(filename.c_str(), &db, SQLITE_OPEN_READWRITE, nullptr), SQLITE_OK)
    if (db == nullptr)
    {
      return;
    }

    char* error_message = nullptr;
    const int rc = sqlite3_exec(db, query.c_str(), nullptr, nullptr, &error_message);
    TEST_EQUAL(rc, SQLITE_OK)
    if (error_message != nullptr)
    {
      sqlite3_free(error_message);
    }
    sqlite3_close(db);
  }

  struct ExtractedArchive
  {
    std::string base_dir;
    std::unique_ptr<File::TempDir> temp_dir;
  };

  ExtractedArchive unzipArchive_(const std::string& archive_path)
  {
    ExtractedArchive archive;
    archive.base_dir = ZipArchiveFile::unzipDirectory(archive_path, archive.temp_dir);
    return archive;
  }

  std::vector<Int64> readRunIds_(const std::string& base_dir)
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

  void checkParquetFeatureColumns_(const std::string& base_dir,
                                   const std::vector<std::string>& expected_columns,
                                   const std::string& non_null_score_column)
  {
    Size total_rows = 0;
    Size total_non_null = 0;
    for (const Int64 run_id : readRunIds_(base_dir))
    {
      const auto table = ParquetFile::readTable(base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet");
      total_rows += static_cast<Size>(table->num_rows());
      for (const auto& column : expected_columns)
      {
        TEST_NOT_EQUAL(ParquetFile::getOptionalColumn(table, column), nullptr)
      }
      total_non_null += countNonNull_(ParquetFile::getOptionalColumn(table, non_null_score_column));
    }
    TEST_TRUE(total_rows > 0)
    TEST_TRUE(total_non_null > 0)
  }

  void checkParquetTransitionColumns_(const std::string& base_dir,
                                      const std::vector<std::string>& expected_columns,
                                      const std::string& non_null_score_column)
  {
    Size total_rows = 0;
    Size total_non_null = 0;
    for (const Int64 run_id : readRunIds_(base_dir))
    {
      const auto table = ParquetFile::readTable(base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet");
      total_rows += static_cast<Size>(table->num_rows());
      for (const auto& column : expected_columns)
      {
        TEST_NOT_EQUAL(ParquetFile::getOptionalColumn(table, column), nullptr)
      }
      total_non_null += countNonNull_(ParquetFile::getOptionalColumn(table, non_null_score_column));
    }
    TEST_TRUE(total_rows > 0)
    TEST_TRUE(total_non_null > 0)
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

START_SECTION(ScoreSummary score(const std::string& input_path, Level level, const std::string& output_path))
{
  OpenSwathPercolatorScoring scorer;

  std::string ms1_osw;
  NEW_TMP_FILE(ms1_osw);
  ms1_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", ms1_osw);
  const auto ms1_summary = scorer.score(ms1_osw, OpenSwathPercolatorScoring::Level::MS1);
  TEST_TRUE(ms1_summary.total_rows > 0)
  TEST_TRUE(ms1_summary.target_rows > 0)
  TEST_TRUE(ms1_summary.decoy_rows > 0)
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1;"),
             querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM FEATURE_MS1;"))
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1 WHERE RANK < 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms1_osw, "SELECT COUNT(*) FROM SCORE_MS1 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  std::string ms2_osw;
  NEW_TMP_FILE(ms2_osw);
  ms2_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", ms2_osw);
  const auto ms2_summary = scorer.score(ms2_osw, OpenSwathPercolatorScoring::Level::MS2);
  TEST_TRUE(ms2_summary.total_rows > 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2;"),
             querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM FEATURE_MS2;"))
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE RANK < 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  std::string ms2_nullfiltered_osw;
  NEW_TMP_FILE(ms2_nullfiltered_osw);
  ms2_nullfiltered_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", ms2_nullfiltered_osw);
  const Int64 feature_ms2_count = querySqliteInt64_(ms2_nullfiltered_osw, "SELECT COUNT(*) FROM FEATURE_MS2;");
  executeSqlite_(ms2_nullfiltered_osw,
                 "ALTER TABLE FEATURE_MS2 ADD COLUMN VAR_NULLFILTER_TEST REAL; "
                 "UPDATE FEATURE_MS2 SET VAR_NULLFILTER_TEST = CAST(FEATURE_ID AS REAL); "
                 "UPDATE FEATURE_MS2 SET VAR_NULLFILTER_TEST = NULL "
                 "WHERE FEATURE_ID = (SELECT FEATURE_ID FROM FEATURE_MS2 ORDER BY FEATURE_ID LIMIT 1);");
  const auto ms2_nullfiltered_summary = scorer.score(ms2_nullfiltered_osw, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(ms2_nullfiltered_summary.total_rows, static_cast<Size>(feature_ms2_count - 1))
  TEST_EQUAL(querySqliteInt64_(ms2_nullfiltered_osw, "SELECT COUNT(*) FROM SCORE_MS2;"), feature_ms2_count - 1)
  TEST_EQUAL(querySqliteInt64_(ms2_nullfiltered_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE RANK < 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms2_nullfiltered_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  const auto transition_summary = scorer.score(ms2_osw, OpenSwathPercolatorScoring::Level::TRANSITION);
  TEST_TRUE(transition_summary.total_rows > 0)
  TEST_TRUE(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION;") > 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION WHERE RANK != 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(ms2_osw, "SELECT COUNT(*) FROM SCORE_TRANSITION WHERE QVALUE < 0 OR QVALUE > 1;"), 0)

  std::string ms1ms2_osw;
  NEW_TMP_FILE(ms1ms2_osw);
  ms1ms2_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", ms1ms2_osw);
  const auto ms1ms2_summary = scorer.score(ms1ms2_osw, OpenSwathPercolatorScoring::Level::MS1MS2);
  TEST_TRUE(ms1ms2_summary.total_rows > 0)
  TEST_EQUAL(querySqliteInt64_(ms1ms2_osw, "SELECT COUNT(*) FROM SCORE_MS2;"),
             querySqliteInt64_(ms1ms2_osw, "SELECT COUNT(*) FROM FEATURE_MS2;"))

  std::string ms1_oswpq;
  NEW_TMP_FILE(ms1_oswpq);
  ms1_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms1_oswpq);
  TEST_EXCEPTION(Exception::InvalidValue, scorer.score(ms1_oswpq, OpenSwathPercolatorScoring::Level::MS1))

  std::string ms2_oswpq;
  NEW_TMP_FILE(ms2_oswpq);
  ms2_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms2_oswpq);
  const auto ms2_oswpq_summary = scorer.score(ms2_oswpq, OpenSwathPercolatorScoring::Level::MS2);
  TEST_TRUE(ms2_oswpq_summary.total_rows > 0)
  const auto transition_archive = unzipArchive_(ms2_oswpq);
  checkParquetFeatureColumns_(transition_archive.base_dir,
                              {"score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_qvalue", "score_ms2_pep"},
                              "score_ms2_score");
  std::string invalid_dir_output;
  NEW_TMP_FILE(invalid_dir_output);
  invalid_dir_output += ".oswpq.zip";
  TEST_EXCEPTION(Exception::InvalidValue,
                 scorer.score(transition_archive.base_dir, OpenSwathPercolatorScoring::Level::MS2, invalid_dir_output))
  TEST_EXCEPTION(Exception::InvalidValue, scorer.score(ms2_oswpq, OpenSwathPercolatorScoring::Level::TRANSITION))

  std::string ms1ms2_oswpq;
  NEW_TMP_FILE(ms1ms2_oswpq);
  ms1ms2_oswpq += ".oswpq";
  copyFixture_("OpenSwathWorkflow_tworuns_1_17.output.oswpq", ms1ms2_oswpq);
  const auto ms1ms2_oswpq_summary = scorer.score(ms1ms2_oswpq, OpenSwathPercolatorScoring::Level::MS1MS2);
  TEST_TRUE(ms1ms2_oswpq_summary.total_rows > 0)
  const auto ms1ms2_archive = unzipArchive_(ms1ms2_oswpq);
  checkParquetFeatureColumns_(ms1ms2_archive.base_dir,
                              {"score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_qvalue", "score_ms2_pep"},
                              "score_ms2_score");
}
END_SECTION

END_TEST
