// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <arrow/api.h>
#include <arrow/array/concatenate.h>
#include <sqlite3.h>

#include <algorithm>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <tuple>
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

  // Reverse every row and remap all generated identifiers without changing scientific values.
  void permuteParquet_(const std::string& path)
  {
    const auto table = ParquetFile::readTable(path);
    if (table->num_rows() == 0)
    {
      return;
    }
    std::vector<std::shared_ptr<arrow::Array>> columns;
    for (int col = 0; col < table->num_columns(); ++col)
    {
      const auto array = ParquetFile::getColumn(table, table->field(col)->name());
      const auto& name = table->field(col)->name();
      if (name == "run_id" || name == "precursor_id" || name == "feature_id")
      {
        arrow::Int64Builder builder;
        for (int64_t row = table->num_rows(); row-- > 0;)
        {
          ParquetFile::appendOrThrow(builder.Append(-1 - ParquetFile::getInt64(array, row, 0, false)), name);
        }
        columns.push_back(ParquetFile::finishArray(builder, name));
      }
      else
      {
        arrow::ArrayVector slices;
        for (int64_t row = table->num_rows(); row-- > 0;)
        {
          slices.push_back(array->Slice(row, 1));
        }
        const auto result = arrow::Concatenate(slices);
        TEST_TRUE(result.ok())
        columns.push_back(result.ValueOrDie());
      }
    }
    ParquetFile::writeTable(arrow::Table::Make(table->schema(), columns), path);
  }

  using ScientificScores = std::map<std::tuple<std::string, std::string, Int64, bool, double, double>, std::multiset<std::vector<double>>>;

  ScientificScores sqliteScientificScores_(const std::string& path)
  {
    ScientificScores scores;
    sqlite3* db = nullptr;
    TEST_EQUAL(sqlite3_open_v2(path.c_str(), &db, SQLITE_OPEN_READONLY, nullptr), SQLITE_OK)
    sqlite3_stmt* stmt = nullptr;
    const char* query = "SELECT RUN.FILENAME, PEPTIDE.MODIFIED_SEQUENCE, PRECURSOR.CHARGE, PRECURSOR.DECOY, "
      "FEATURE.EXP_RT, SCORE_MS2.SCORE, SCORE_MS2.RANK, SCORE_MS2.QVALUE, SCORE_MS2.PEP, FEATURE_MS2.VAR_XCORR_SHAPE "
      "FROM FEATURE JOIN RUN ON RUN.ID = FEATURE.RUN_ID JOIN PRECURSOR ON PRECURSOR.ID = FEATURE.PRECURSOR_ID "
      "JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "
      "JOIN PEPTIDE ON PEPTIDE.ID = PEPTIDE_ID JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
      "JOIN FEATURE_MS2 ON FEATURE.ID = FEATURE_MS2.FEATURE_ID;";
    TEST_EQUAL(sqlite3_prepare_v2(db, query, -1, &stmt, nullptr), SQLITE_OK)
    while (sqlite3_step(stmt) == SQLITE_ROW)
    {
      const auto text = [&](int col) { return std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt, col))); };
      const auto key = std::make_tuple(text(0), text(1), sqlite3_column_int64(stmt, 2), sqlite3_column_int(stmt, 3) != 0,
                                       sqlite3_column_double(stmt, 4), sqlite3_column_double(stmt, 9));
      scores[key].insert({sqlite3_column_double(stmt, 5), sqlite3_column_double(stmt, 6),
                          sqlite3_column_double(stmt, 7), sqlite3_column_double(stmt, 8)});
    }
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return scores;
  }

  struct ExtractedArchive
  {
    std::string base_dir;
    std::unique_ptr<TempDir> temp_dir;
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

  /// Peak groups of one run in which, in a fraction of the target groups, an interference peak beats the true peak
  /// on var_library_corr but not on var_xcorr_shape, var_xcorr_coelution and var_norm_rt_score (review of #10610).
  /// Decoy groups get interference peaks at the same rate. Portable: std::mt19937 is specified by the standard and
  /// the features are sums of uniforms (no <random> distributions, no libm), so every platform draws the same data.
  struct ControlledCase
  {
    RescoreInput input;
    std::vector<int> truth; // 1: the true peak of its group
  };

  double approxGauss_(std::mt19937& rng, const double mean, const double sd)
  {
    double s = 0; // sum of 4 uniforms, variance 4/12
    for (int k = 0; k < 4; ++k) s += (static_cast<double>(rng()) + 0.5) / 4294967296.0;
    return mean + sd * (s - 2.0) * 1.7320508075688772;
  }

  /// n_groups target and n_groups decoy groups of 5 peak groups; a fraction frac_true of the groups has a true
  /// peak (targets only) and, of those, a fraction frac_interference an interference peak
  ControlledCase makeControlledCase_(const int n_groups, const double frac_true, const double frac_interference,
                                     const unsigned seed)
  {
    ControlledCase c;
    c.input.feature_names = {"var_library_corr", "var_xcorr_shape", "var_xcorr_coelution", "var_norm_rt_score", "var_noise"};
    std::mt19937 rng(seed);
    auto uniform = [&rng]() { return (static_cast<double>(rng()) + 0.5) / 4294967296.0; };
    int key = 0;
    for (int decoy = 0; decoy < 2; ++decoy)
    {
      for (int g = 0; g < n_groups; ++g, ++key)
      {
        const bool has_true = uniform() < frac_true;
        const bool has_interference = has_true && uniform() < frac_interference;
        std::vector<std::pair<std::vector<double>, int>> peaks;
        if (has_true && !decoy) // true peak: good on everything
        {
          peaks.push_back({{approxGauss_(rng, 0.65, 0.2), approxGauss_(rng, 0.8, 0.1), approxGauss_(rng, 1.8, 1.0),
                            approxGauss_(rng, 1.8, 1.0), approxGauss_(rng, 0, 1)}, 1});
        }
        if (has_interference) // best library correlation, poor shape, coelution and RT
        {
          peaks.push_back({{approxGauss_(rng, 0.8, 0.1), approxGauss_(rng, 0.45, 0.15), approxGauss_(rng, 0.3, 1.0),
                            approxGauss_(rng, 0.0, 1.0), approxGauss_(rng, 0, 1)}, 0});
        }
        while (peaks.size() < 5) // noise
        {
          peaks.push_back({{approxGauss_(rng, 0.2, 0.25), approxGauss_(rng, 0.4, 0.15), approxGauss_(rng, 0.0, 1.0),
                            approxGauss_(rng, 0.0, 1.0), approxGauss_(rng, 0, 1)}, 0});
        }
        for (Size i = peaks.size() - 1; i > 0; --i) // portable shuffle
        {
          std::swap(peaks[i], peaks[rng() % (i + 1)]);
        }
        for (auto& peak : peaks)
        {
          c.input.features.push_back(peak.first);
          c.input.is_decoy.push_back(decoy == 1);
          c.input.cv_group_keys.push_back(key);
          c.truth.push_back(peak.second);
        }
      }
    }
    return c;
  }

  int trueInRows_(const ControlledCase& c, const std::vector<Size>& rows)
  {
    int n = 0;
    for (const Size i : rows) n += c.truth[i];
    return n;
  }

  /// Trains Percolator on @p rows (as OpenSwathPercolatorScoring::score does with the training subset), scores all
  /// rows with the raw model and counts the target groups at 1% FDR (target-decoy on the group tops) whose top
  /// pick is the true peak. @p trained is set to false if the model weights only one feature, i.e. if
  /// Percolator's sanity check replaced the trained weights by the initial direction.
  int trueTopPicksAtOnePercent_(const ControlledCase& c, const std::vector<Size>& rows, bool& trained)
  {
    RescoreInput train;
    train.feature_names = c.input.feature_names;
    for (const Size i : rows)
    {
      train.features.push_back(c.input.features[i]);
      train.is_decoy.push_back(c.input.is_decoy[i]);
      train.cv_group_keys.push_back(c.input.cv_group_keys[i]);
    }
    Percolator percolator;
    Param p = percolator.getParameters();
    p.setValue("subset_max_train", 0);
    percolator.setParameters(p);
    const PercolatorModel model = percolator.train(train);
    trained = std::count_if(model.weights.begin(), model.weights.end() - 1, [](double w) { return w != 0.0; }) > 1;

    const int n_groups = c.input.cv_group_keys.back() + 1;
    std::vector<double> score(c.input.features.size());
    std::vector<long> top(n_groups, -1);
    for (Size i = 0; i < score.size(); ++i)
    {
      double v = model.weights.back();
      for (Size j = 0; j < c.input.features[i].size(); ++j) v += model.weights[j] * c.input.features[i][j];
      score[i] = v;
      long& t = top[c.input.cv_group_keys[i]];
      if (t < 0 || score[i] > score[t]) t = static_cast<long>(i);
    }
    std::vector<int> order(n_groups);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return score[top[a]] > score[top[b]]; });
    int n_target = 0, n_decoy = 0, cut = 0; // the longest prefix with decoys <= 1% of targets
    for (int k = 0; k < n_groups; ++k)
    {
      (c.input.is_decoy[top[order[k]]] ? n_decoy : n_target)++;
      if (n_decoy <= 0.01 * n_target) cut = k + 1;
    }
    int n_true = 0;
    for (int k = 0; k < cut; ++k)
    {
      if (!c.input.is_decoy[top[order[k]]] && c.truth[top[order[k]]] == 1) ++n_true;
    }
    return n_true;
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

START_SECTION([EXTRA] training on the best peak group per precursor plus a sample of the others)
{
  // a cap far below the number of peak groups: training uses a subset, all peak groups are scored
  std::string capped_osw;
  NEW_TMP_FILE(capped_osw);
  capped_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", capped_osw);
  const Int64 n_features = querySqliteInt64_(capped_osw, "SELECT COUNT(*) FROM FEATURE_MS2;");
  const Int64 n_precursors = querySqliteInt64_(capped_osw, "SELECT COUNT(DISTINCT PRECURSOR_ID) FROM FEATURE;");
  TEST_TRUE(n_features > n_precursors + 5) // otherwise nothing would be left out

  OpenSwathPercolatorScoring capped;
  Param p = capped.getParameters();
  p.setValue("percolator:subset_max_train", 5);
  capped.setParameters(p);
  const auto capped_summary = capped.score(capped_osw, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(static_cast<Int64>(capped_summary.total_rows), n_features)
  TEST_EQUAL(querySqliteInt64_(capped_osw, "SELECT COUNT(*) FROM SCORE_MS2;"), n_features)
  TEST_EQUAL(querySqliteInt64_(capped_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE QVALUE < 0 OR QVALUE > 1;"), 0)
  TEST_EQUAL(querySqliteInt64_(capped_osw, "SELECT COUNT(*) FROM SCORE_MS2 WHERE RANK < 1;"), 0)

  // the same seed draws the same sample: the scores are reproducible
  std::string capped_again_osw;
  NEW_TMP_FILE(capped_again_osw);
  capped_again_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", capped_again_osw);
  capped.score(capped_again_osw, OpenSwathPercolatorScoring::Level::MS2);
  executeSqlite_(capped_osw, "ATTACH DATABASE '" + capped_again_osw + "' AS again; "
                 "CREATE TABLE score_differences AS SELECT COUNT(*) AS N FROM SCORE_MS2 a "
                 "JOIN again.SCORE_MS2 b ON a.FEATURE_ID = b.FEATURE_ID WHERE a.SCORE != b.SCORE;");
  TEST_EQUAL(querySqliteInt64_(capped_osw, "SELECT N FROM score_differences;"), 0)

  // a main score that is not among the features: all peak groups are used, no error
  std::string unknown_score_osw;
  NEW_TMP_FILE(unknown_score_osw);
  unknown_score_osw += ".osw";
  copyFixture_("PyProphet_inference_test.osw", unknown_score_osw);
  p.setValue("training_main_score", "var_no_such_score");
  capped.setParameters(p);
  const auto unknown_summary = capped.score(unknown_score_osw, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(static_cast<Int64>(unknown_summary.total_rows), n_features)
}
END_SECTION

START_SECTION([EXTRA] scientific input is invariant to generated IDs and row order)
{
  std::string original;
  NEW_TMP_FILE(original);
  original += ".osw";
  std::string permuted;
  NEW_TMP_FILE(permuted);
  permuted += ".osw";
  copyFixture_("PyProphet_inference_test.osw", original);
  copyFixture_("PyProphet_inference_test.osw", permuted);
  // Force equal retention times. Distinct score vectors, not feature IDs, must break these ties.
  executeSqlite_(original, "UPDATE FEATURE SET EXP_RT = 100;");
  executeSqlite_(permuted, "UPDATE FEATURE SET EXP_RT = 100;");
  // Reverse identifier order and physical score-table order. Preserve all scientific input.
  executeSqlite_(permuted,
    "PRAGMA foreign_keys=OFF; UPDATE RUN SET ID=-1-ID; UPDATE FEATURE SET RUN_ID=-1-RUN_ID; "
    "UPDATE PRECURSOR SET ID=-1-ID; UPDATE FEATURE SET PRECURSOR_ID=-1-PRECURSOR_ID; "
    "UPDATE PRECURSOR_PEPTIDE_MAPPING SET PRECURSOR_ID=-1-PRECURSOR_ID; "
    "UPDATE FEATURE SET ID=-1-ID; UPDATE FEATURE_MS1 SET FEATURE_ID=-1-FEATURE_ID; "
    "UPDATE FEATURE_MS2 SET FEATURE_ID=-1-FEATURE_ID; "
    "CREATE TABLE reversed AS SELECT * FROM FEATURE_MS2 ORDER BY FEATURE_ID; "
    "DELETE FROM FEATURE_MS2; INSERT INTO FEATURE_MS2 SELECT * FROM reversed; DROP TABLE reversed;");
  for (const int cap : {0, 5})
  {
    OpenSwathPercolatorScoring scorer;
    Param p = scorer.getParameters();
    p.setValue("percolator:subset_max_train", cap);
    scorer.setParameters(p);
    const auto first = scorer.score(original, OpenSwathPercolatorScoring::Level::MS2);
    const auto second = scorer.score(permuted, OpenSwathPercolatorScoring::Level::MS2);
    TEST_EQUAL(first.total_rows, second.total_rows)
    const auto expected = sqliteScientificScores_(original);
    const auto actual = sqliteScientificScores_(permuted);
    TEST_TRUE(!expected.empty())
    TEST_EQUAL(expected.size(), actual.size())
    TEST_TRUE(expected == actual)
    TEST_EQUAL(querySqliteInt64_(original, "SELECT COUNT(*) FROM SCORE_MS2 WHERE RANK=1;"),
               querySqliteInt64_(original, "SELECT COUNT(*) FROM (SELECT DISTINCT RUN_ID, PRECURSOR_ID FROM FEATURE);"))
  }
  executeSqlite_(permuted, "DELETE FROM PRECURSOR_PEPTIDE_MAPPING;");
  OpenSwathPercolatorScoring scorer;
  TEST_EXCEPTION(Exception::MissingInformation, scorer.score(permuted, OpenSwathPercolatorScoring::Level::MS2))
  executeSqlite_(original, "UPDATE FEATURE SET EXP_RT = 1e999;");
  TEST_EXCEPTION(Exception::InvalidValue, scorer.score(original, OpenSwathPercolatorScoring::Level::MS2))
}
END_SECTION

START_SECTION([EXTRA] Parquet scientific input is invariant to generated IDs and row order)
{
  const auto original = unzipArchive_(OPENMS_GET_TEST_DATA_PATH("OpenSwathWorkflow_tworuns_1_17.output.oswpq"));
  const auto permuted = unzipArchive_(OPENMS_GET_TEST_DATA_PATH("OpenSwathWorkflow_tworuns_1_17.output.oswpq"));
  const auto run_ids = readRunIds_(original.base_dir);
  permuteParquet_(permuted.base_dir + "/runs/runs.parquet");
  permuteParquet_(permuted.base_dir + "/library/precursors.parquet");
  for (const Int64 id : run_ids)
  {
    const auto old_dir = permuted.base_dir + "/runs/run_id=" + StringUtils::toStr(id);
    const auto new_dir = permuted.base_dir + "/runs/run_id=" + StringUtils::toStr(-1 - id);
    permuteParquet_(old_dir + "/features.parquet");
    TEST_TRUE(File::rename(old_dir, new_dir))
  }
  OpenSwathPercolatorScoring scorer;
  const auto first = scorer.score(original.base_dir, OpenSwathPercolatorScoring::Level::MS2);
  const auto second = scorer.score(permuted.base_dir, OpenSwathPercolatorScoring::Level::MS2);
  TEST_EQUAL(first.total_rows, second.total_rows)
  for (const Int64 id : run_ids)
  {
    const auto expected = ParquetFile::readTable(original.base_dir + "/runs/run_id=" + StringUtils::toStr(id) + "/features.parquet");
    const auto actual = ParquetFile::readTable(permuted.base_dir + "/runs/run_id=" + StringUtils::toStr(-1 - id) + "/features.parquet");
    TEST_EQUAL(expected->num_rows(), actual->num_rows())
    for (const std::string column : {"score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_qvalue", "score_ms2_pep"})
    {
      const auto left = ParquetFile::getColumn(expected, column);
      const auto right = ParquetFile::getColumn(actual, column);
      for (int64_t row = 0; row < expected->num_rows(); ++row)
      {
        // The permutation reversed every row. Null score rows must also remain null.
        TEST_TRUE(left->Slice(row, 1)->Equals(right->Slice(actual->num_rows() - 1 - row, 1)))
      }
    }
  }
}
END_SECTION

START_SECTION([EXTRA] native SQL molecules and retained metadata validation)
{
  std::string path;
  NEW_TMP_FILE(path);
  path += ".osw";
  copyFixture_("PyProphet_inference_test.osw", path);
  // Actual native PQP/OSW molecule tables, with no peptide mapping at all.
  executeSqlite_(path,
    "DROP TABLE IF EXISTS PRECURSOR_COMPOUND_MAPPING; DROP TABLE IF EXISTS COMPOUND; "
    "CREATE TABLE COMPOUND(ID INT PRIMARY KEY NOT NULL, COMPOUND_NAME TEXT NOT NULL, "
    "SUM_FORMULA TEXT NOT NULL, SMILES TEXT NOT NULL, ADDUCTS TEXT NOT NULL, DECOY INT NOT NULL); "
    "CREATE TABLE PRECURSOR_COMPOUND_MAPPING(PRECURSOR_ID INT NOT NULL, COMPOUND_ID INT NOT NULL); "
    "INSERT INTO COMPOUND SELECT ID, MODIFIED_SEQUENCE, 'C6H12O6', 'C(O)C', '[M+H]+', 0 FROM PEPTIDE; "
    "INSERT INTO PRECURSOR_COMPOUND_MAPPING SELECT PRECURSOR_ID, PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING; "
    "DROP TABLE PRECURSOR_PEPTIDE_MAPPING; DROP TABLE PEPTIDE;");
  OpenSwathPercolatorScoring scorer;
  const auto summary = scorer.score(path, OpenSwathPercolatorScoring::Level::MS2);
  TEST_TRUE(summary.total_rows > 0)
  TEST_EQUAL(querySqliteInt64_(path, "SELECT COUNT(*) FROM SCORE_MS2;"), summary.total_rows)
  // The bundled legacy native OSW schema predates the optional ADDUCTS column.
  executeSqlite_(path,
    "CREATE TABLE legacy_compound AS SELECT ID,COMPOUND_NAME,SUM_FORMULA,SMILES,DECOY FROM COMPOUND; "
    "DROP TABLE COMPOUND; ALTER TABLE legacy_compound RENAME TO COMPOUND;");
  TEST_EQUAL(scorer.score(path, OpenSwathPercolatorScoring::Level::MS2).total_rows, summary.total_rows)
  executeSqlite_(path, "DELETE FROM PRECURSOR_COMPOUND_MAPPING;");
  TEST_EXCEPTION(Exception::MissingInformation, scorer.score(path, OpenSwathPercolatorScoring::Level::MS2))
}
END_SECTION

START_SECTION([EXTRA] incomplete SQL rows and unused NULL metadata do not abort)
{
  std::string path;
  NEW_TMP_FILE(path);
  path += ".osw";
  copyFixture_("PyProphet_inference_test.osw", path);
  // Recreate only metadata tables to permit explicit NULLs in this regression fixture.
  executeSqlite_(path,
    "CREATE TABLE nullable_run AS SELECT * FROM RUN; DROP TABLE RUN; ALTER TABLE nullable_run RENAME TO RUN; "
    "CREATE TABLE nullable_precursor AS SELECT * FROM PRECURSOR; DROP TABLE PRECURSOR; ALTER TABLE nullable_precursor RENAME TO PRECURSOR; "
    "CREATE TABLE nullable_peptide AS SELECT * FROM PEPTIDE; DROP TABLE PEPTIDE; ALTER TABLE nullable_peptide RENAME TO PEPTIDE; "
    "CREATE TABLE nullable_ms2 AS SELECT * FROM FEATURE_MS2; DROP TABLE FEATURE_MS2; ALTER TABLE nullable_ms2 RENAME TO FEATURE_MS2; "
    "INSERT INTO RUN(ID,FILENAME) VALUES(-999,NULL); "
    "INSERT INTO PRECURSOR(ID,CHARGE,DECOY) VALUES(-999,NULL,NULL); "
    "INSERT INTO PEPTIDE(ID,MODIFIED_SEQUENCE) VALUES(-999,NULL); "
    "INSERT INTO PRECURSOR_PEPTIDE_MAPPING VALUES(-999,-999); "
    "INSERT INTO FEATURE(ID,RUN_ID,PRECURSOR_ID,EXP_RT,NORM_RT,DELTA_RT,LEFT_WIDTH,RIGHT_WIDTH) "
    "VALUES(-999,-999,-999,1e999,0,0,0,0); "
    "INSERT INTO FEATURE_MS2(FEATURE_ID) VALUES(-999);");
  OpenSwathPercolatorScoring scorer;
  const auto summary = scorer.score(path, OpenSwathPercolatorScoring::Level::MS2);
  TEST_TRUE(summary.total_rows > 0)
  TEST_EQUAL(querySqliteInt64_(path, "SELECT COUNT(*) FROM SCORE_MS2 WHERE FEATURE_ID=-999;"), 0)
  // Retained missing metadata must report MissingInformation, not a raw SQL extraction error.
  executeSqlite_(path, "UPDATE PRECURSOR SET CHARGE=NULL WHERE ID=(SELECT PRECURSOR_ID FROM FEATURE WHERE ID<>-999 LIMIT 1);");
  TEST_EXCEPTION(Exception::MissingInformation, scorer.score(path, OpenSwathPercolatorScoring::Level::MS2))
  executeSqlite_(path, "UPDATE PRECURSOR SET CHARGE=2, DECOY=NULL WHERE ID=(SELECT PRECURSOR_ID FROM FEATURE WHERE ID<>-999 LIMIT 1);");
  TEST_EXCEPTION(Exception::MissingInformation, scorer.score(path, OpenSwathPercolatorScoring::Level::MS2))
}
END_SECTION

START_SECTION([EXTRA] missing SQL feature tables precede metadata reads)
{
  for (const auto level : {OpenSwathPercolatorScoring::Level::MS1, OpenSwathPercolatorScoring::Level::MS2,
                           OpenSwathPercolatorScoring::Level::MS1MS2, OpenSwathPercolatorScoring::Level::TRANSITION})
  {
    std::string path;
    NEW_TMP_FILE(path);
    path += ".osw";
    copyFixture_("PyProphet_inference_test.osw", path);
    executeSqlite_(path, "DROP TABLE FEATURE; DROP TABLE RUN; DROP TABLE PEPTIDE;");
    OpenSwathPercolatorScoring scorer;
    TEST_EXCEPTION(Exception::Precondition, scorer.score(path, level))
  }
}
END_SECTION

START_SECTION([EXTRA] native Parquet molecules with empty modified sequence)
{
  const auto archive = unzipArchive_(OPENMS_GET_TEST_DATA_PATH("OpenSwathWorkflow_tworuns_1_17.output.oswpq"));
  TransitionParquetFile library;
  OpenSwath::LightTargetedExperiment experiment;
  OpenSwathLibraryIDNormalizer::SourceIDMapping source_ids;
  library.convertParquetToTargetedExperiment(archive.base_dir, experiment, &source_ids);
  // The workflow archive includes unused library precursors with no transitions. Native library
  // writing requires a transition carrying precursor_mz for every compound being written.
  std::set<std::string> referenced_precursors;
  for (const auto& transition : experiment.transitions)
  {
    referenced_precursors.insert(transition.peptide_ref);
  }
  std::erase_if(experiment.compounds, [&](const auto& compound)
  {
    return !referenced_precursors.contains(compound.id);
  });
  TEST_TRUE(!experiment.compounds.empty())
  for (auto& compound : experiment.compounds)
  {
    compound.sequence.clear();
  }
  // Exercise the native writer rather than inventing compound columns in Parquet.
  library.convertLightTargetedExperimentToParquet(archive.base_dir, experiment, &source_ids);
  const auto precursors = ParquetFile::readTable(archive.base_dir + "/library/precursors.parquet");
  const auto sequences = ParquetFile::getColumn(precursors, "modified_sequence");
  for (int64_t row = 0; row < precursors->num_rows(); ++row)
  {
    TEST_TRUE(ParquetFile::getString(sequences, row).empty())
  }
  // An unused native library row with NULL metadata must not abort shared feature readers.
  std::vector<std::shared_ptr<arrow::Array>> with_unused;
  for (const auto& field : precursors->schema()->fields())
  {
    const auto values = ParquetFile::getColumn(precursors, field->name());
    std::shared_ptr<arrow::Array> extra;
    if (field->name() == "precursor_id")
    {
      arrow::Int64Builder builder;
      ParquetFile::appendOrThrow(builder.Append(-999), "precursor_id");
      extra = ParquetFile::finishArray(builder, "precursor_id");
    }
    else
    {
      const auto null_array = arrow::MakeArrayOfNull(field->type(), 1);
      TEST_TRUE(null_array.ok())
      extra = null_array.ValueOrDie();
    }
    const auto combined = arrow::Concatenate({values, extra});
    TEST_TRUE(combined.ok())
    with_unused.push_back(combined.ValueOrDie());
  }
  ParquetFile::writeTable(arrow::Table::Make(precursors->schema(), with_unused), archive.base_dir + "/library/precursors.parquet");
  OpenSwathPercolatorScoring scorer;
  TEST_TRUE(scorer.score(archive.base_dir, OpenSwathPercolatorScoring::Level::MS2).total_rows > 0)
  // A retained NULL charge is missing scientific identity, not an eager reader error.
  const int charge_index = precursors->schema()->GetFieldIndex("charge");
  const auto null_charges = arrow::MakeArrayOfNull(precursors->field(charge_index)->type(), precursors->num_rows() + 1);
  TEST_TRUE(null_charges.ok())
  with_unused[charge_index] = null_charges.ValueOrDie();
  ParquetFile::writeTable(arrow::Table::Make(precursors->schema(), with_unused), archive.base_dir + "/library/precursors.parquet");
  TEST_EXCEPTION(Exception::MissingInformation, scorer.score(archive.base_dir, OpenSwathPercolatorScoring::Level::MS2))
}
END_SECTION

START_SECTION((static std::vector<Size> trainingRows(const RescoreInput& input, const std::string& main_score, Size max_others, UInt64 seed)))
{
  // 6 groups of 4 rows; the best row by var_a is row g % 4 of group g, the last group is a tie (all equal)
  RescoreInput input;
  input.feature_names = {"var_a", "var_b"};
  const std::vector<int> keys = {10, 11, 12, 13, 14, 15}; // keys need not be consecutive
  for (Size g = 0; g < keys.size(); ++g)
  {
    for (Size r = 0; r < 4; ++r)
    {
      const double a = (g == 5) ? 0.5 : (r == g % 4 ? 1.0 : 0.1 * static_cast<double>(r));
      input.features.push_back({a, static_cast<double>(g)});
      input.is_decoy.push_back(g % 2 == 1);
      input.cv_group_keys.push_back(keys[g]);
    }
  }
  const Size n_rows = input.features.size();
  std::vector<Size> all_rows(n_rows);
  std::iota(all_rows.begin(), all_rows.end(), Size(0));

  // everything: no cap, a cap the other rows fit into, an unknown main score, no group keys
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "var_a", 0, 1) == all_rows, true)
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "var_a", 18, 1) == all_rows, true)
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "var_no_such_score", 7, 1) == all_rows, true)
  RescoreInput no_keys = input;
  no_keys.cv_group_keys.clear();
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(no_keys, "var_a", 7, 1) == all_rows, true)

  // which rows enter training: the best of every group, and complete groups while their other rows fit
  auto rowsPerGroup = [&](const std::vector<Size>& rows)
  {
    std::map<int, std::vector<Size>> per_group;
    for (const Size i : rows) per_group[input.cv_group_keys[i]].push_back(i);
    return per_group;
  };
  const std::vector<Size> rows = OpenSwathPercolatorScoring::trainingRows(input, "var_a", 7, 1);
  TEST_EQUAL(std::is_sorted(rows.begin(), rows.end()), true)
  TEST_EQUAL(rows.size(), 12) // 6 best + 2 complete groups (3 others each); a third would exceed 7
  Size n_complete = 0;
  for (const auto& [key, group_rows] : rowsPerGroup(rows))
  {
    const Size g = static_cast<Size>(key - 10);
    TEST_EQUAL(group_rows.size() == 1 || group_rows.size() == 4, true)
    if (group_rows.size() == 1)
    {
      TEST_EQUAL(group_rows[0], 4 * g + (g == 5 ? 0 : g % 4)) // the best; a tie keeps the first row
    }
    n_complete += group_rows.size() == 4;
  }
  TEST_EQUAL(rowsPerGroup(rows).size(), 6)
  TEST_EQUAL(n_complete, 2)
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "var_a", 7, 1) == rows, true) // same seed, same rows

  // a cap below the other rows of one group: only the best rows
  std::vector<Size> best_rows;
  for (Size g = 0; g < 6; ++g) best_rows.push_back(4 * g + (g == 5 ? 0 : g % 4));
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "var_a", 2, 1) == best_rows, true)

  // '-' selects the lowest value: row 1 of groups 0 and 4 (0.1), row 0 of groups 1-3 (0.0), row 0 of the tie group
  const std::vector<Size> lowest = {1, 4, 8, 12, 17, 20};
  TEST_EQUAL(OpenSwathPercolatorScoring::trainingRows(input, "-var_a", 2, 1) == lowest, true)
}
END_SECTION

START_SECTION([EXTRA] training subset when an interference peak beats the true peak on the main score)
{
  // 5,000 target and 5,000 decoy groups of 5 peak groups; 30% of the target groups have a true peak and, of those,
  // 40% an interference peak with a better var_library_corr. Training subset: the best peak group of each group
  // plus at most 500 others (125 complete groups).
  const ControlledCase c = makeControlledCase_(5000, 0.3, 0.4, 1);
  std::vector<Size> all_rows(c.input.features.size());
  std::iota(all_rows.begin(), all_rows.end(), Size(0));
  const int n_true = trueInRows_(c, all_rows);
  TEST_EQUAL(n_true, 1525)

  const std::vector<Size> by_shape = OpenSwathPercolatorScoring::trainingRows(c.input, "var_xcorr_shape", 500, 1);
  const std::vector<Size> by_corr = OpenSwathPercolatorScoring::trainingRows(c.input, "var_library_corr", 500, 1);
  TEST_EQUAL(by_shape.size(), 10500)
  TEST_EQUAL(by_corr.size(), 10500)
  // which candidates enter training: selecting by var_xcorr_shape (the default) keeps 1,447 of the 1,525 true
  // peaks, selecting by var_library_corr 825 (the others lose to the interference or a noise peak)
  TEST_EQUAL(trueInRows_(c, by_shape) >= 0.9 * n_true, true)
  TEST_EQUAL(trueInRows_(c, by_corr) <= 0.6 * n_true, true)

  // subset training against training on all rows (subset_max_train = 0); all three train (no fallback to one
  // feature). True top picks at 1% FDR here: all rows 1,234, by var_xcorr_shape 1,258, by var_library_corr 1,236
  // (over data seeds 1-5, the subsets reach 0.967-1.019 of all rows); allow 5% for platform differences.
  bool trained_all = false, trained_shape = false, trained_corr = false;
  const int with_all = trueTopPicksAtOnePercent_(c, all_rows, trained_all);
  const int with_shape = trueTopPicksAtOnePercent_(c, by_shape, trained_shape);
  const int with_corr = trueTopPicksAtOnePercent_(c, by_corr, trained_corr);
  TEST_EQUAL(trained_all, true)
  TEST_EQUAL(trained_shape, true)
  TEST_EQUAL(trained_corr, true)
  TEST_EQUAL(with_all >= 1100, true)
  TEST_EQUAL(with_shape >= 0.95 * with_all, true)
  TEST_EQUAL(with_corr >= 0.95 * with_all, true)
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
