// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <sqlite3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{
  const std::array<std::string, static_cast<Size>(OpenSwathPercolatorScoring::Level::SIZE_OF_LEVEL)>
    OpenSwathPercolatorScoring::names_of_level = {"ms1", "ms2", "ms1ms2", "transition"};

  struct OpenSwathPercolatorScoring::Impl
  {
    Percolator percolator;
    Int32 transition_peakgroup_max_rank = 1;
    double transition_peakgroup_max_pep = 0.7;
    double transition_max_isotope_overlap = 0.5;
    double transition_min_log_sn = 0.0;
  };

  namespace
  {
    namespace Sql = Internal::SqliteHelper;
    using Level = OpenSwathPercolatorScoring::Level;
    using ScoreSummary = OpenSwathPercolatorScoring::ScoreSummary;

    struct ScoreRow
    {
      Int64 feature_id = -1;
      Int64 run_id = -1;
      Int64 precursor_id = -1;
      std::optional<Int64> transition_id;
      bool decoy = false;
      int cv_group_key = 0;
      std::vector<double> feature_values;
    };

    struct PreparedInput
    {
      RescoreInput input;
      std::vector<ScoreRow> rows;
      bool transition_level = false;
    };

    struct ScoreResult
    {
      double score = -100.0;
      Int32 rank = 1;
      std::optional<double> pvalue;
      double qvalue = 1.0;
      double pep = 1.0;
    };

    struct ParquetWorkspace
    {
      String base_dir;
      String output_path;
      bool output_is_archive = false;
      std::unique_ptr<File::TempDir> temp_dir;
    };

    template <typename K>
    struct PairHash
    {
      Size operator()(const K& key) const noexcept
      {
        return std::hash<typename K::first_type>{}(key.first) ^
               (std::hash<typename K::second_type>{}(key.second) << 1);
      }
    };

    String levelName_(const Level level)
    {
      return OpenSwathPercolatorScoring::names_of_level.at(static_cast<Size>(level));
    }

    String canonicalFeatureNameSqlite_(const String& column_name,
                                       const Level level,
                                       const bool from_ms1_table)
    {
      String out = column_name;
      out.toLower();
      if (level == Level::MS1MS2 && from_ms1_table)
      {
        if (out.hasPrefix("var_"))
        {
          out = "var_ms1_" + out.substr(4);
        }
      }
      return out;
    }

    String canonicalFeatureNameParquet_(const String& column_name,
                                        const Level level,
                                        const bool from_ms1_table)
    {
      (void)from_ms1_table;
      String out = column_name;
      if (out.hasPrefix("var_ms2_"))
      {
        out = "var_" + out.substr(8);
      }
      else if (out.hasPrefix("var_ms1_"))
      {
        if (level == Level::MS1MS2)
        {
          out = out;
        }
        else
        {
          out = "var_" + out.substr(8);
        }
      }
      return out.toLower();
    }

    std::vector<String> sqliteColumnsWithPrefix_(SqliteConnector& conn,
                                                 const String& table_name,
                                                 const String& prefix)
    {
      sqlite3_stmt* stmt = nullptr;
      conn.prepareStatement(&stmt, "PRAGMA table_info('" + table_name + "');");
      std::vector<String> columns;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        const String col_name = Sql::extractString(stmt, 1);
        if (col_name.hasPrefix(prefix))
        {
          columns.push_back(col_name);
        }
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);
      std::sort(columns.begin(), columns.end());
      return columns;
    }

    void checkSqliteRc_(sqlite3* db, const int rc, const String& action)
    {
      if (rc != SQLITE_OK && rc != SQLITE_DONE)
      {
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            action + ": " + sqlite3_errmsg(db));
      }
    }

    void bindNullableDouble_(sqlite3* db, sqlite3_stmt* stmt, const int pos, const std::optional<double>& value)
    {
      const int rc = value.has_value() ? sqlite3_bind_double(stmt, pos, *value) :
                                         sqlite3_bind_null(stmt, pos);
      checkSqliteRc_(db, rc, "Failed to bind nullable REAL");
    }

    PreparedInput finalizeInput_(std::vector<ScoreRow>&& rows,
                                 std::vector<String>&& feature_names,
                                 const bool transition_level)
    {
      if (rows.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No OpenSWATH rows available for Percolator scoring.");
      }
      if (feature_names.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No feature columns available for OpenSWATH Percolator scoring.", "");
      }

      const Size n_rows = rows.size();
      const Size n_cols = feature_names.size();
      std::vector<bool> keep_column(n_cols, false);
      constexpr double eps = 1e-12;

      for (Size col = 0; col < n_cols; ++col)
      {
        double first = std::numeric_limits<double>::quiet_NaN();
        bool first_seen = false;
        for (Size row = 0; row < n_rows; ++row)
        {
          double value = rows[row].feature_values[col];
          if (!std::isfinite(value))
          {
            value = 0.0;
            rows[row].feature_values[col] = 0.0;
          }
          if (!first_seen)
          {
            first = value;
            first_seen = true;
          }
          else if (std::fabs(value - first) > eps)
          {
            keep_column[col] = true;
            break;
          }
        }
      }

      StringList kept_names;
      for (Size col = 0; col < n_cols; ++col)
      {
        if (keep_column[col])
        {
          kept_names.push_back(feature_names[col]);
        }
      }

      if (kept_names.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "All OpenSWATH feature columns are constant after null-to-zero normalization.", "");
      }

      PreparedInput prepared;
      prepared.rows = std::move(rows);
      prepared.transition_level = transition_level;
      prepared.input.feature_names = kept_names;
      prepared.input.features.reserve(n_rows);
      prepared.input.is_decoy.reserve(n_rows);
      prepared.input.cv_group_keys.reserve(n_rows);

      for (const auto& row : prepared.rows)
      {
        std::vector<double> filtered;
        filtered.reserve(kept_names.size());
        for (Size col = 0; col < n_cols; ++col)
        {
          if (keep_column[col])
          {
            filtered.push_back(std::isfinite(row.feature_values[col]) ? row.feature_values[col] : 0.0);
          }
        }
        prepared.input.features.push_back(std::move(filtered));
        prepared.input.is_decoy.push_back(row.decoy);
        prepared.input.cv_group_keys.push_back(row.cv_group_key);
      }

      return prepared;
    }

    PreparedInput extractSqliteRows_(const String& osw_path,
                                     const Level level,
                                     const Int32 transition_peakgroup_max_rank,
                                     const double transition_peakgroup_max_pep,
                                     const double transition_max_isotope_overlap,
                                     const double transition_min_log_sn)
    {
      SqliteConnector conn(osw_path, SqliteConnector::SqlOpenMode::READ_ONLY);
      sqlite3_stmt* stmt = nullptr;
      std::vector<ScoreRow> rows;
      std::vector<String> feature_names;

      if (level == Level::MS1 || level == Level::MS2 || level == Level::MS1MS2)
      {
        const String ms2_table = "FEATURE_MS2";
        const String ms1_table = "FEATURE_MS1";
        if (level != Level::MS1 && !conn.tableExists(ms2_table))
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "OpenSWATH Percolator scoring requires FEATURE_MS2 for ms2 / ms1ms2 scoring.");
        }
        if ((level == Level::MS1 || level == Level::MS1MS2) && !conn.tableExists(ms1_table))
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "OpenSWATH Percolator scoring requires FEATURE_MS1 for ms1 / ms1ms2 scoring.");
        }

        std::vector<String> ms2_cols;
        std::vector<String> ms1_cols;
        if (level == Level::MS2 || level == Level::MS1MS2)
        {
          ms2_cols = sqliteColumnsWithPrefix_(conn, ms2_table, "VAR_");
        }
        if (level == Level::MS1 || level == Level::MS1MS2)
        {
          ms1_cols = sqliteColumnsWithPrefix_(conn, ms1_table, "VAR_");
        }

        for (const auto& c : ms2_cols)
        {
          feature_names.push_back(canonicalFeatureNameSqlite_(c, level, false));
        }
        for (const auto& c : ms1_cols)
        {
          feature_names.push_back(canonicalFeatureNameSqlite_(c, level, true));
        }

        String select_sql =
          "SELECT FEATURE.ID AS FEATURE_ID, FEATURE.RUN_ID, FEATURE.PRECURSOR_ID, PRECURSOR.DECOY";
        for (const auto& c : ms2_cols)
        {
          select_sql += ", FEATURE_MS2." + c;
        }
        for (const auto& c : ms1_cols)
        {
          select_sql += ", FEATURE_MS1." + c;
        }
        select_sql += " FROM FEATURE ";
        select_sql += "INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID ";
        if (level == Level::MS2)
        {
          select_sql += "INNER JOIN FEATURE_MS2 ON FEATURE.ID = FEATURE_MS2.FEATURE_ID ";
        }
        else if (level == Level::MS1)
        {
          select_sql += "INNER JOIN FEATURE_MS1 ON FEATURE.ID = FEATURE_MS1.FEATURE_ID ";
        }
        else
        {
          select_sql += "INNER JOIN FEATURE_MS2 ON FEATURE.ID = FEATURE_MS2.FEATURE_ID ";
          select_sql += "LEFT JOIN FEATURE_MS1 ON FEATURE.ID = FEATURE_MS1.FEATURE_ID ";
        }
        select_sql += "ORDER BY FEATURE.RUN_ID, FEATURE.PRECURSOR_ID, FEATURE.EXP_RT, FEATURE.ID;";

        conn.prepareStatement(&stmt, select_sql);
        Sql::SqlState state = Sql::nextRow(stmt);
        std::unordered_map<std::pair<Int64, Int64>, int, PairHash<std::pair<Int64, Int64>>> group_to_key;
        int next_group_key = 0;
        while (state == Sql::SqlState::SQL_ROW)
        {
          ScoreRow row;
          int col = 0;
          row.feature_id = Sql::extractInt64(stmt, col++);
          row.run_id = Sql::extractInt64(stmt, col++);
          row.precursor_id = Sql::extractInt64(stmt, col++);
          row.decoy = Sql::extractInt(stmt, col++) != 0;
          row.feature_values.reserve(feature_names.size());
          for (Size i = 0; i < ms2_cols.size(); ++i, ++col)
          {
            row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ? 0.0 : sqlite3_column_double(stmt, col));
          }
          for (Size i = 0; i < ms1_cols.size(); ++i, ++col)
          {
            row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ? 0.0 : sqlite3_column_double(stmt, col));
          }
          const std::pair<Int64, Int64> group{row.run_id, row.precursor_id};
          auto it = group_to_key.find(group);
          if (it == group_to_key.end())
          {
            it = group_to_key.emplace(group, next_group_key++).first;
          }
          row.cv_group_key = it->second;
          rows.push_back(std::move(row));
          state = Sql::nextRow(stmt, state);
        }
        sqlite3_finalize(stmt);
        return finalizeInput_(std::move(rows), std::move(feature_names), false);
      }

      if (!conn.tableExists("SCORE_MS2"))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH Percolator scoring requires prior ms2/ms1ms2 scoring (missing SCORE_MS2).");
      }
      if (!conn.tableExists("FEATURE_TRANSITION"))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH Percolator scoring requires FEATURE_TRANSITION.");
      }

      const auto transition_cols = sqliteColumnsWithPrefix_(conn, "FEATURE_TRANSITION", "VAR_");
      feature_names = transition_cols;

      String select_sql =
        "SELECT FEATURE.RUN_ID, FEATURE.ID AS FEATURE_ID, FEATURE.PRECURSOR_ID, "
        "       FEATURE_TRANSITION.TRANSITION_ID, TRANSITION.DECOY";
      for (const auto& c : transition_cols)
      {
        select_sql += ", FEATURE_TRANSITION." + c;
      }
      select_sql +=
        " FROM FEATURE_TRANSITION "
        "INNER JOIN FEATURE ON FEATURE_TRANSITION.FEATURE_ID = FEATURE.ID "
        "INNER JOIN PRECURSOR ON FEATURE.PRECURSOR_ID = PRECURSOR.ID "
        "INNER JOIN SCORE_MS2 ON FEATURE.ID = SCORE_MS2.FEATURE_ID "
        "INNER JOIN TRANSITION ON FEATURE_TRANSITION.TRANSITION_ID = TRANSITION.ID "
        "WHERE SCORE_MS2.RANK <= " + String(transition_peakgroup_max_rank) +
        "  AND SCORE_MS2.PEP <= " + String(transition_peakgroup_max_pep) +
        "  AND FEATURE_TRANSITION.VAR_ISOTOPE_OVERLAP_SCORE <= " + String(transition_max_isotope_overlap) +
        "  AND FEATURE_TRANSITION.VAR_LOG_SN_SCORE > " + String(transition_min_log_sn) +
        "  AND PRECURSOR.DECOY = 0 "
        "ORDER BY FEATURE.RUN_ID, FEATURE.PRECURSOR_ID, FEATURE.EXP_RT, FEATURE_TRANSITION.TRANSITION_ID;";

      conn.prepareStatement(&stmt, select_sql);
      Sql::SqlState state = Sql::nextRow(stmt);
      int next_group_key = 0;
      while (state == Sql::SqlState::SQL_ROW)
      {
        ScoreRow row;
        int col = 0;
        row.run_id = Sql::extractInt64(stmt, col++);
        row.feature_id = Sql::extractInt64(stmt, col++);
        row.precursor_id = Sql::extractInt64(stmt, col++);
        row.transition_id = Sql::extractInt64(stmt, col++);
        row.decoy = Sql::extractInt(stmt, col++) != 0;
        row.feature_values.reserve(feature_names.size());
        for (Size i = 0; i < transition_cols.size(); ++i, ++col)
        {
          row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ? 0.0 : sqlite3_column_double(stmt, col));
        }
        row.cv_group_key = next_group_key++;
        rows.push_back(std::move(row));
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      if (rows.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH Percolator scoring found no candidate transitions. "
                                      "Run ms2/ms1ms2 scoring first and check the transition candidate filters.");
      }

      return finalizeInput_(std::move(rows), std::move(feature_names), true);
    }

    PreparedInput extractParquetRows_(const String& oswpq_path,
                                      const Level level,
                                      const Int32 transition_peakgroup_max_rank,
                                      const double transition_peakgroup_max_pep,
                                      const double transition_max_isotope_overlap,
                                      const double transition_min_log_sn)
    {
      OpenSwathOSWParquetReader reader;
      std::vector<ScoreRow> rows;
      std::vector<String> feature_names;

      if (level == Level::MS2)
      {
        const auto result = reader.fetchPeakGroupFeatures(oswpq_path, "ms2");
        std::vector<Size> selected_cols;
        for (Size i = 0; i < result.ms2_columns.size(); ++i)
        {
          if (result.ms2_columns[i].hasPrefix("var_ms2_"))
          {
            selected_cols.push_back(i);
            feature_names.push_back(canonicalFeatureNameParquet_(result.ms2_columns[i], level, false));
          }
        }
        std::unordered_map<std::pair<Int64, Int64>, int, PairHash<std::pair<Int64, Int64>>> group_to_key;
        int next_group_key = 0;
        for (Size i = 0; i < result.feature_id.size(); ++i)
        {
          ScoreRow row;
          row.feature_id = result.feature_id[i];
          row.run_id = result.run_id[i];
          row.precursor_id = result.precursor_id[i];
          row.decoy = result.decoy[i];
          row.feature_values.reserve(selected_cols.size());
          for (const Size col_idx : selected_cols)
          {
            const auto& values = result.ms2_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ? values[i] : 0.0);
          }
          const std::pair<Int64, Int64> group{row.run_id, row.precursor_id};
          auto it = group_to_key.find(group);
          if (it == group_to_key.end())
          {
            it = group_to_key.emplace(group, next_group_key++).first;
          }
          row.cv_group_key = it->second;
          rows.push_back(std::move(row));
        }
        return finalizeInput_(std::move(rows), std::move(feature_names), false);
      }

      if (level == Level::MS1 || level == Level::MS1MS2)
      {
        const auto result = reader.fetchPeakGroupFeatures(oswpq_path, "ms1ms2");
        std::vector<Size> ms2_selected;
        std::vector<Size> ms1_selected;
        if (level == Level::MS1MS2)
        {
          for (Size i = 0; i < result.ms2_columns.size(); ++i)
          {
            if (result.ms2_columns[i].hasPrefix("var_ms2_"))
            {
              ms2_selected.push_back(i);
              feature_names.push_back(canonicalFeatureNameParquet_(result.ms2_columns[i], level, false));
            }
          }
        }
        for (Size i = 0; i < result.ms1_columns.size(); ++i)
        {
          if (result.ms1_columns[i].hasPrefix("var_ms1_"))
          {
            ms1_selected.push_back(i);
            feature_names.push_back(canonicalFeatureNameParquet_(result.ms1_columns[i], level, true));
          }
        }

        if (level == Level::MS1 && ms1_selected.empty())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "OpenSWATH OSWPQ ms1 scoring requires var_ms1_* feature columns.");
        }
        if (level == Level::MS1MS2 && (ms1_selected.empty() || feature_names.empty()))
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "OpenSWATH OSWPQ ms1ms2 scoring requires both var_ms2_* and var_ms1_* feature columns.");
        }

        std::unordered_map<std::pair<Int64, Int64>, int, PairHash<std::pair<Int64, Int64>>> group_to_key;
        int next_group_key = 0;
        for (Size i = 0; i < result.feature_id.size(); ++i)
        {
          ScoreRow row;
          row.feature_id = result.feature_id[i];
          row.run_id = result.run_id[i];
          row.precursor_id = result.precursor_id[i];
          row.decoy = result.decoy[i];
          row.feature_values.reserve(feature_names.size());
          for (const Size col_idx : ms2_selected)
          {
            const auto& values = result.ms2_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ? values[i] : 0.0);
          }
          for (const Size col_idx : ms1_selected)
          {
            const auto& values = result.ms1_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ? values[i] : 0.0);
          }
          const std::pair<Int64, Int64> group{row.run_id, row.precursor_id};
          auto it = group_to_key.find(group);
          if (it == group_to_key.end())
          {
            it = group_to_key.emplace(group, next_group_key++).first;
          }
          row.cv_group_key = it->second;
          rows.push_back(std::move(row));
        }
        return finalizeInput_(std::move(rows), std::move(feature_names), false);
      }

      const auto transition_result = reader.fetchTransitionFeatures(oswpq_path);
      feature_names = transition_result.transition_var_columns;
      for (auto& feature_name : feature_names)
      {
        feature_name = feature_name.toLower();
      }

      const auto precursors_table = ParquetFile::readTable(oswpq_path + "/library/precursors.parquet");
      const auto precursor_id_array = ParquetFile::getColumn(precursors_table, "precursor_id");
      const auto precursor_decoy_array = ParquetFile::getOptionalColumn(precursors_table, "decoy");
      std::unordered_map<Int64, bool> precursor_decoy;
      precursor_decoy.reserve(precursors_table->num_rows());
      for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
      {
        precursor_decoy.emplace(
          ParquetFile::getInt64(precursor_id_array, row, 0, false),
          ParquetFile::getBool(precursor_decoy_array, row, false, true));
      }

      using FeatureKey = std::pair<Int64, Int64>;
      std::unordered_map<FeatureKey, std::tuple<Int64, Int32, double>, PairHash<FeatureKey>> feature_to_score;
      feature_to_score.reserve(transition_result.feature_id.size());

      const auto runs_table = ParquetFile::readTable(oswpq_path + "/runs/runs.parquet");
      const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");
      for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
      {
        const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
        const String run_dir = oswpq_path + "/runs/run_id=" + String(run_id);
        const auto features_table = ParquetFile::readTable(run_dir + "/features.parquet");
        const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
        const auto feature_precursor_col = ParquetFile::getColumn(features_table, "precursor_id");
        auto score_ms2_rank_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_peak_group_rank");
        if (!score_ms2_rank_col)
        {
          score_ms2_rank_col = ParquetFile::getOptionalColumn(features_table, "SCORE_MS2_PEAK_GROUP_RANK");
        }
        auto score_ms2_pep_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_pep");
        if (!score_ms2_pep_col)
        {
          score_ms2_pep_col = ParquetFile::getOptionalColumn(features_table, "SCORE_MS2_PEP");
        }
        if (!score_ms2_rank_col || !score_ms2_pep_col)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Transition-level OpenSWATH OSWPQ Percolator scoring requires prior ms2/ms1ms2 scoring "
                                        "(missing score_ms2_peak_group_rank / score_ms2_pep columns in features.parquet).");
        }

        for (int64_t feat_row = 0; feat_row < features_table->num_rows(); ++feat_row)
        {
          const Int64 feature_id = ParquetFile::getInt64(feature_id_col, feat_row, 0, false);
          const Int64 precursor_id = ParquetFile::getInt64(feature_precursor_col, feat_row, 0, false);
          const Int32 rank = static_cast<Int32>(ParquetFile::getInt64(score_ms2_rank_col, feat_row, 0, true));
          const double pep = ParquetFile::getDouble(score_ms2_pep_col, feat_row, std::numeric_limits<double>::quiet_NaN(), true);
          feature_to_score.emplace(FeatureKey{run_id, feature_id}, std::make_tuple(precursor_id, rank, pep));
        }
      }

      const auto isotope_it = std::find(feature_names.begin(), feature_names.end(), "var_isotope_overlap_score");
      const Size isotope_col_idx = isotope_it != feature_names.end() ?
        static_cast<Size>(std::distance(feature_names.begin(), isotope_it)) :
        std::numeric_limits<Size>::max();
      const auto log_sn_it = std::find(feature_names.begin(), feature_names.end(), "var_log_sn_score");
      const Size log_sn_col_idx = log_sn_it != feature_names.end() ?
        static_cast<Size>(std::distance(feature_names.begin(), log_sn_it)) :
        std::numeric_limits<Size>::max();

      int next_group_key = 0;
      for (Size row_idx = 0; row_idx < transition_result.feature_id.size(); ++row_idx)
      {
        const FeatureKey feature_key{transition_result.run_id[row_idx], transition_result.feature_id[row_idx]};
        const auto score_it = feature_to_score.find(feature_key);
        if (score_it == feature_to_score.end())
        {
          continue;
        }

        const auto& [precursor_id, peakgroup_rank, peakgroup_pep] = score_it->second;
        if (peakgroup_rank <= 0 ||
            peakgroup_rank > transition_peakgroup_max_rank ||
            !std::isfinite(peakgroup_pep) ||
            peakgroup_pep > transition_peakgroup_max_pep)
        {
          continue;
        }
        const auto precursor_it = precursor_decoy.find(precursor_id);
        if (precursor_it != precursor_decoy.end() && precursor_it->second)
        {
          continue;
        }

        ScoreRow out;
        out.feature_id = transition_result.feature_id[row_idx];
        out.run_id = transition_result.run_id[row_idx];
        out.precursor_id = precursor_id;
        out.transition_id = transition_result.transition_id[row_idx];
        out.decoy = transition_result.decoy[row_idx];
        out.feature_values.reserve(feature_names.size());
        for (Size col_idx = 0; col_idx < transition_result.transition_var_values.size(); ++col_idx)
        {
          out.feature_values.push_back(transition_result.transition_var_values[col_idx][row_idx]);
        }
        if (isotope_col_idx < out.feature_values.size())
        {
          const double isotope_overlap = out.feature_values[isotope_col_idx];
          if (!std::isfinite(isotope_overlap) || isotope_overlap > transition_max_isotope_overlap)
          {
            continue;
          }
        }
        if (log_sn_col_idx < out.feature_values.size())
        {
          const double log_sn = out.feature_values[log_sn_col_idx];
          if (!std::isfinite(log_sn) || !(log_sn > transition_min_log_sn))
          {
            continue;
          }
        }
        out.cv_group_key = next_group_key++;
        rows.push_back(std::move(out));
      }

      if (rows.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH OSWPQ Percolator scoring found no candidate transitions. "
                                      "Run ms2/ms1ms2 scoring first and check the transition candidate filters.");
      }

      return finalizeInput_(std::move(rows), std::move(feature_names), true);
    }

    std::vector<ScoreResult> computeRanks_(const PreparedInput& prepared,
                                           const RescoreOutput& output)
    {
      if (output.scores.size() != prepared.rows.size() ||
          output.q_values.size() != prepared.rows.size() ||
          output.peps.size() != prepared.rows.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Percolator output size does not match OpenSWATH input size.", "");
      }

      std::vector<ScoreResult> results(prepared.rows.size());
      for (Size i = 0; i < prepared.rows.size(); ++i)
      {
        results[i].score = output.scores[i];
        results[i].qvalue = output.q_values[i];
        results[i].pep = output.peps[i];
      }

      if (prepared.transition_level)
      {
        for (auto& result : results)
        {
          result.rank = 1;
        }
        return results;
      }

      std::unordered_map<int, std::vector<Size>> groups;
      for (Size i = 0; i < prepared.rows.size(); ++i)
      {
        groups[prepared.rows[i].cv_group_key].push_back(i);
      }

      for (auto& [_, indices] : groups)
      {
        std::stable_sort(indices.begin(), indices.end(),
                         [&](const Size lhs, const Size rhs)
                         {
                           return results[lhs].score > results[rhs].score;
                         });
        for (Size rank = 0; rank < indices.size(); ++rank)
        {
          results[indices[rank]].rank = static_cast<Int32>(rank + 1);
        }
      }
      return results;
    }

    String prepareSqliteOutput_(const String& input_path, const String& output_path)
    {
      const String target_path = output_path.empty() ? input_path : output_path;
      if (target_path == input_path)
      {
        return target_path;
      }
      if (File::exists(target_path) && !File::remove(target_path))
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, target_path);
      }
      if (!File::copy(input_path, target_path))
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, target_path);
      }
      return target_path;
    }

    void validateOutputPath_(const String& input_path, const String& output_path)
    {
      if (output_path.empty())
      {
        return;
      }

      const bool input_is_osw = input_path.hasSuffix(".osw");
      if (input_is_osw)
      {
        if (!output_path.hasSuffix(".osw"))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "SQLite OSW input can only be written to another .osw output path.",
                                        output_path);
        }
        return;
      }

      if (output_path.hasSuffix(".osw"))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "OSWPQ input can only be written to another OSWPQ path.",
                                      output_path);
      }

      if (File::isDirectory(input_path))
      {
        if (File::exists(output_path) && !File::isDirectory(output_path))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Directory OSWPQ input can only be written to another directory path.",
                                        output_path);
        }
      }
      else
      {
        if (File::exists(output_path) && File::isDirectory(output_path))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Archive OSWPQ input can only be written to another archive path.",
                                        output_path);
        }
        if (!(output_path.hasSuffix(".oswpq") || output_path.hasSuffix(".oswpq.zip")))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Archive OSWPQ input can only be written to an '.oswpq' or '.oswpq.zip' output path.",
                                        output_path);
        }
      }
    }

    void writeSqliteScores_(const String& osw_path,
                            const Level level,
                            const PreparedInput& prepared,
                            const std::vector<ScoreResult>& results)
    {
      SqliteConnector conn(osw_path, SqliteConnector::SqlOpenMode::READWRITE);
      const bool transition_level = level == Level::TRANSITION;
      const String table_name =
        level == Level::MS1 ? "SCORE_MS1" :
        level == Level::TRANSITION ? "SCORE_TRANSITION" :
        "SCORE_MS2";

      String create_sql = "DROP TABLE IF EXISTS " + table_name + "; ";
      if (transition_level)
      {
        create_sql +=
          "CREATE TABLE " + table_name + "("
          "FEATURE_ID INTEGER NOT NULL, "
          "TRANSITION_ID INTEGER NOT NULL, "
          "SCORE REAL NOT NULL, "
          "RANK INTEGER NOT NULL, "
          "PVALUE REAL, "
          "QVALUE REAL NOT NULL, "
          "PEP REAL NOT NULL);";
      }
      else
      {
        create_sql +=
          "CREATE TABLE " + table_name + "("
          "FEATURE_ID INTEGER NOT NULL, "
          "SCORE REAL NOT NULL, "
          "RANK INTEGER NOT NULL, "
          "PVALUE REAL, "
          "QVALUE REAL NOT NULL, "
          "PEP REAL NOT NULL);";
      }
      conn.executeStatement(create_sql);
      if (table_name == "SCORE_MS2")
      {
        conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);");
      }
      else if (table_name == "SCORE_TRANSITION")
      {
        conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_transition_feature_id ON SCORE_TRANSITION (FEATURE_ID);");
        conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_transition_transition_id ON SCORE_TRANSITION (TRANSITION_ID);");
      }

      sqlite3* db = conn.getDB();
      sqlite3_stmt* insert_stmt = nullptr;
      if (transition_level)
      {
        conn.prepareStatement(&insert_stmt,
          "INSERT INTO SCORE_TRANSITION (FEATURE_ID, TRANSITION_ID, SCORE, RANK, PVALUE, QVALUE, PEP) "
          "VALUES (?, ?, ?, ?, ?, ?, ?);");
      }
      else
      {
        conn.prepareStatement(&insert_stmt,
          "INSERT INTO " + table_name + " (FEATURE_ID, SCORE, RANK, PVALUE, QVALUE, PEP) "
          "VALUES (?, ?, ?, ?, ?, ?);");
      }

      conn.executeStatement("BEGIN TRANSACTION;");
      for (Size i = 0; i < prepared.rows.size(); ++i)
      {
        const auto& row = prepared.rows[i];
        const auto& result = results[i];
        int rc = sqlite3_bind_int64(insert_stmt, 1, row.feature_id);
        checkSqliteRc_(db, rc, "Failed to bind FEATURE_ID");
        int bind_col = 2;
        if (transition_level)
        {
          rc = sqlite3_bind_int64(insert_stmt, bind_col++, *row.transition_id);
          checkSqliteRc_(db, rc, "Failed to bind TRANSITION_ID");
        }
        rc = sqlite3_bind_double(insert_stmt, bind_col++, result.score);
        checkSqliteRc_(db, rc, "Failed to bind SCORE");
        rc = sqlite3_bind_int(insert_stmt, bind_col++, result.rank);
        checkSqliteRc_(db, rc, "Failed to bind RANK");
        bindNullableDouble_(db, insert_stmt, bind_col++, result.pvalue);
        rc = sqlite3_bind_double(insert_stmt, bind_col++, result.qvalue);
        checkSqliteRc_(db, rc, "Failed to bind QVALUE");
        rc = sqlite3_bind_double(insert_stmt, bind_col++, result.pep);
        checkSqliteRc_(db, rc, "Failed to bind PEP");
        rc = sqlite3_step(insert_stmt);
        if (rc != SQLITE_DONE)
        {
          sqlite3_finalize(insert_stmt);
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Failed to insert OpenSWATH Percolator score row: " + String(sqlite3_errmsg(db)));
        }
        sqlite3_reset(insert_stmt);
        sqlite3_clear_bindings(insert_stmt);
      }
      sqlite3_finalize(insert_stmt);
      conn.executeStatement("END TRANSACTION;");
    }

    ParquetWorkspace prepareParquetWorkspace_(const String& input_path, const String& output_path)
    {
      const bool input_is_dir = File::isDirectory(input_path);
      ParquetWorkspace workspace;
      workspace.output_is_archive = !input_is_dir;
      workspace.output_path = output_path.empty() ? input_path : output_path;

      if (input_is_dir)
      {
        if (!output_path.empty() && File::exists(workspace.output_path) && !File::isDirectory(workspace.output_path))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Directory OSWPQ input can only be written to another directory path.", workspace.output_path);
        }
        if (!output_path.empty() && workspace.output_path != input_path)
        {
          if (File::exists(workspace.output_path))
          {
            File::removeDirRecursively(workspace.output_path);
          }
          if (!File::copyDirRecursively(input_path, workspace.output_path))
          {
            throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, workspace.output_path);
          }
          workspace.base_dir = workspace.output_path;
        }
        else
        {
          workspace.base_dir = input_path;
        }
      }
      else
      {
        if (!output_path.empty() && File::exists(workspace.output_path) && File::isDirectory(workspace.output_path))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Archive OSWPQ input can only be written to another archive path.", workspace.output_path);
        }
        workspace.base_dir = ZipArchiveFile::unzipDirectory(input_path, workspace.temp_dir);
      }
      return workspace;
    }

    std::vector<String> removalPrefixesForLevel_(const Level level)
    {
      if (level == Level::MS1)
      {
        return {"score_ms1_", "SCORE_MS1_"};
      }
      if (level == Level::TRANSITION)
      {
        return {"score_transition_", "SCORE_TRANSITION_"};
      }
      return {"score_ms2_", "SCORE_MS2_"};
    }

    bool hasAnyPrefix_(const String& value, const std::vector<String>& prefixes)
    {
      for (const auto& prefix : prefixes)
      {
        if (value.hasPrefix(prefix))
        {
          return true;
        }
      }
      return false;
    }

    void appendFeatureScoreColumns_(const String& file_path,
                                    const std::vector<String>& prefixes_to_remove,
                                    const std::vector<std::shared_ptr<arrow::Field>>& extra_fields,
                                    const std::vector<std::shared_ptr<arrow::Array>>& extra_arrays)
    {
      auto table = ParquetFile::readTable(file_path);
      std::vector<std::shared_ptr<arrow::Field>> fields;
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      fields.reserve(table->num_columns() + extra_fields.size());
      arrays.reserve(table->num_columns() + extra_arrays.size());

      for (int i = 0; i < table->num_columns(); ++i)
      {
        const auto field = table->field(i);
        const String field_name = field->name();
        if (hasAnyPrefix_(field_name, prefixes_to_remove))
        {
          continue;
        }
        fields.push_back(field);
        arrays.push_back(table->column(i)->chunk(0));
      }
      for (Size i = 0; i < extra_fields.size(); ++i)
      {
        fields.push_back(extra_fields[i]);
        arrays.push_back(extra_arrays[i]);
      }
      ParquetFile::writeTable(arrow::Table::Make(arrow::schema(fields), arrays), file_path);
    }

    void writeParquetScores_(const ParquetWorkspace& workspace,
                             const Level level,
                             const PreparedInput& prepared,
                             const std::vector<ScoreResult>& results)
    {
      const auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
      const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");
      const auto prefixes_to_remove = removalPrefixesForLevel_(level);

      if (level == Level::TRANSITION)
      {
        using FeatureTransitionKey = std::pair<Int64, Int64>;
        std::unordered_map<Int64, std::unordered_map<FeatureTransitionKey, Size, PairHash<FeatureTransitionKey>>> run_feature_transition_to_index;
        for (Size i = 0; i < prepared.rows.size(); ++i)
        {
          run_feature_transition_to_index[prepared.rows[i].run_id][FeatureTransitionKey{prepared.rows[i].feature_id, *prepared.rows[i].transition_id}] = i;
        }

        for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
        {
          const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
          const String ft_path = workspace.base_dir + "/runs/run_id=" + String(run_id) + "/feature_transition.parquet";
          auto table = ParquetFile::readTable(ft_path);
          const auto feature_id_col = ParquetFile::getColumn(table, "feature_id");
          const auto transition_id_col = ParquetFile::getColumn(table, "transition_id");

          arrow::DoubleBuilder score_builder;
          arrow::Int32Builder rank_builder;
          arrow::DoubleBuilder pvalue_builder;
          arrow::DoubleBuilder qvalue_builder;
          arrow::DoubleBuilder pep_builder;
          for (int64_t row = 0; row < table->num_rows(); ++row)
          {
            const Int64 feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
            const Int64 transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
            const auto run_it = run_feature_transition_to_index.find(run_id);
            if (run_it != run_feature_transition_to_index.end())
            {
              const auto row_it = run_it->second.find(FeatureTransitionKey{feature_id, transition_id});
              if (row_it != run_it->second.end())
              {
                const auto& result = results[row_it->second];
                ParquetFile::appendOrThrow(score_builder.Append(result.score), "score_transition_score");
                ParquetFile::appendOrThrow(rank_builder.Append(result.rank), "score_transition_rank");
                ParquetFile::appendOrThrow(pvalue_builder.AppendNull(), "score_transition_pvalue");
                ParquetFile::appendOrThrow(qvalue_builder.Append(result.qvalue), "score_transition_qvalue");
                ParquetFile::appendOrThrow(pep_builder.Append(result.pep), "score_transition_pep");
                continue;
              }
            }
            ParquetFile::appendOrThrow(score_builder.AppendNull(), "score_transition_score");
            ParquetFile::appendOrThrow(rank_builder.AppendNull(), "score_transition_rank");
            ParquetFile::appendOrThrow(pvalue_builder.AppendNull(), "score_transition_pvalue");
            ParquetFile::appendOrThrow(qvalue_builder.AppendNull(), "score_transition_qvalue");
            ParquetFile::appendOrThrow(pep_builder.AppendNull(), "score_transition_pep");
          }

          appendFeatureScoreColumns_(ft_path,
                                     prefixes_to_remove,
                                     {
                                       arrow::field("score_transition_score", arrow::float64(), true),
                                       arrow::field("score_transition_rank", arrow::int32(), true),
                                       arrow::field("score_transition_pvalue", arrow::float64(), true),
                                       arrow::field("score_transition_qvalue", arrow::float64(), true),
                                       arrow::field("score_transition_pep", arrow::float64(), true)
                                     },
                                     {
                                       ParquetFile::finishArray(score_builder, "score_transition_score"),
                                       ParquetFile::finishArray(rank_builder, "score_transition_rank"),
                                       ParquetFile::finishArray(pvalue_builder, "score_transition_pvalue"),
                                       ParquetFile::finishArray(qvalue_builder, "score_transition_qvalue"),
                                       ParquetFile::finishArray(pep_builder, "score_transition_pep")
                                     });
        }
      }
      else
      {
        std::unordered_map<Int64, std::unordered_map<Int64, Size>> run_feature_to_index;
        for (Size i = 0; i < prepared.rows.size(); ++i)
        {
          run_feature_to_index[prepared.rows[i].run_id][prepared.rows[i].feature_id] = i;
        }

        const String score_name = level == Level::MS1 ? "score_ms1_score" : "score_ms2_score";
        const String rank_name = level == Level::MS1 ? "score_ms1_peak_group_rank" : "score_ms2_peak_group_rank";
        const String pvalue_name = level == Level::MS1 ? "score_ms1_pvalue" : "score_ms2_pvalue";
        const String qvalue_name = level == Level::MS1 ? "score_ms1_qvalue" : "score_ms2_qvalue";
        const String pep_name = level == Level::MS1 ? "score_ms1_pep" : "score_ms2_pep";

        for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
        {
          const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
          const String features_path = workspace.base_dir + "/runs/run_id=" + String(run_id) + "/features.parquet";
          auto table = ParquetFile::readTable(features_path);
          const auto feature_id_col = ParquetFile::getColumn(table, "feature_id");

          arrow::DoubleBuilder score_builder;
          arrow::Int32Builder rank_builder;
          arrow::DoubleBuilder pvalue_builder;
          arrow::DoubleBuilder qvalue_builder;
          arrow::DoubleBuilder pep_builder;
          for (int64_t row = 0; row < table->num_rows(); ++row)
          {
            const Int64 feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
            const auto run_it = run_feature_to_index.find(run_id);
            if (run_it == run_feature_to_index.end())
            {
              ParquetFile::appendOrThrow(score_builder.AppendNull(), score_name);
              ParquetFile::appendOrThrow(rank_builder.AppendNull(), rank_name);
              ParquetFile::appendOrThrow(pvalue_builder.AppendNull(), pvalue_name);
              ParquetFile::appendOrThrow(qvalue_builder.AppendNull(), qvalue_name);
              ParquetFile::appendOrThrow(pep_builder.AppendNull(), pep_name);
              continue;
            }
            const auto result_it = run_it->second.find(feature_id);
            if (result_it == run_it->second.end())
            {
              ParquetFile::appendOrThrow(score_builder.AppendNull(), score_name);
              ParquetFile::appendOrThrow(rank_builder.AppendNull(), rank_name);
              ParquetFile::appendOrThrow(pvalue_builder.AppendNull(), pvalue_name);
              ParquetFile::appendOrThrow(qvalue_builder.AppendNull(), qvalue_name);
              ParquetFile::appendOrThrow(pep_builder.AppendNull(), pep_name);
              continue;
            }
            const auto& result = results[result_it->second];
            ParquetFile::appendOrThrow(score_builder.Append(result.score), score_name);
            ParquetFile::appendOrThrow(rank_builder.Append(result.rank), rank_name);
            ParquetFile::appendOrThrow(pvalue_builder.AppendNull(), pvalue_name);
            ParquetFile::appendOrThrow(qvalue_builder.Append(result.qvalue), qvalue_name);
            ParquetFile::appendOrThrow(pep_builder.Append(result.pep), pep_name);
          }

          appendFeatureScoreColumns_(features_path,
                                     prefixes_to_remove,
                                     {
                                       arrow::field(std::string(score_name), arrow::float64(), true),
                                       arrow::field(std::string(rank_name), arrow::int32(), true),
                                       arrow::field(std::string(pvalue_name), arrow::float64(), true),
                                       arrow::field(std::string(qvalue_name), arrow::float64(), true),
                                       arrow::field(std::string(pep_name), arrow::float64(), true)
                                     },
                                     {
                                       ParquetFile::finishArray(score_builder, score_name),
                                       ParquetFile::finishArray(rank_builder, rank_name),
                                       ParquetFile::finishArray(pvalue_builder, pvalue_name),
                                       ParquetFile::finishArray(qvalue_builder, qvalue_name),
                                       ParquetFile::finishArray(pep_builder, pep_name)
                                     });
        }
      }

      if (workspace.output_is_archive)
      {
        ZipArchiveFile::zipDirectory(workspace.base_dir, workspace.output_path);
      }
      ZipArchiveFile::writeSidecarIndex(workspace.output_path);
    }

    ScoreSummary makeSummary_(const PreparedInput& prepared)
    {
      ScoreSummary summary;
      summary.total_rows = prepared.rows.size();
      std::set<Int64> features;
      for (const auto& row : prepared.rows)
      {
        if (row.decoy) ++summary.decoy_rows;
        else ++summary.target_rows;
        features.insert(row.feature_id);
      }
      summary.feature_count = features.size();
      return summary;
    }
  } // namespace

  OpenSwathPercolatorScoring::OpenSwathPercolatorScoring() :
    DefaultParamHandler("OpenSwathPercolatorScoring"),
    impl_(std::make_unique<Impl>())
  {
    defaults_.setValue("transition:peakgroup_max_rank", impl_->transition_peakgroup_max_rank,
                       "Transition scoring: maximum MS2 peak-group rank retained for transition rescoring.");
    defaults_.setMinInt("transition:peakgroup_max_rank", 1);

    defaults_.setValue("transition:peakgroup_max_pep", impl_->transition_peakgroup_max_pep,
                       "Transition scoring: maximum MS2 peak-group PEP retained for transition rescoring.");
    defaults_.setMinFloat("transition:peakgroup_max_pep", 0.0);
    defaults_.setMaxFloat("transition:peakgroup_max_pep", 1.0);

    defaults_.setValue("transition:max_isotope_overlap", impl_->transition_max_isotope_overlap,
                       "Transition scoring: maximum transition isotope-overlap score retained for rescoring.");
    defaults_.setMinFloat("transition:max_isotope_overlap", 0.0);

    defaults_.setValue("transition:min_log_sn", impl_->transition_min_log_sn,
                       "Transition scoring: minimum log signal-to-noise score retained for rescoring.");

    defaults_.insert("percolator:", impl_->percolator.getDefaults());
    defaultsToParam_();
    updateMembers_();
  }

  OpenSwathPercolatorScoring::~OpenSwathPercolatorScoring() = default;

  void OpenSwathPercolatorScoring::updateMembers_()
  {
    impl_->transition_peakgroup_max_rank = static_cast<Int32>(param_.getValue("transition:peakgroup_max_rank"));
    impl_->transition_peakgroup_max_pep = static_cast<double>(param_.getValue("transition:peakgroup_max_pep"));
    impl_->transition_max_isotope_overlap = static_cast<double>(param_.getValue("transition:max_isotope_overlap"));
    impl_->transition_min_log_sn = static_cast<double>(param_.getValue("transition:min_log_sn"));
    impl_->percolator.setParameters(param_.copy("percolator:", true));
  }

  OpenSwathPercolatorScoring::ScoreSummary OpenSwathPercolatorScoring::score(const String& input_path,
                                                                              const Level level,
                                                                              const String& output_path)
  {
    if (!File::exists(input_path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_path);
    }
    validateOutputPath_(input_path, output_path);

    OPENMS_LOG_INFO << "OpenSWATH Percolator scoring (" << levelName_(level)
                    << ") on '" << input_path << "'\n";

    PreparedInput prepared;
    if (input_path.hasSuffix(".osw"))
    {
      const String target_path = prepareSqliteOutput_(input_path, output_path);
      prepared = extractSqliteRows_(target_path, level,
                                    impl_->transition_peakgroup_max_rank,
                                    impl_->transition_peakgroup_max_pep,
                                    impl_->transition_max_isotope_overlap,
                                    impl_->transition_min_log_sn);
      const RescoreOutput output = impl_->percolator.rescore(prepared.input);
      const auto results = computeRanks_(prepared, output);
      writeSqliteScores_(target_path, level, prepared, results);
      return makeSummary_(prepared);
    }

    ParquetWorkspace workspace = prepareParquetWorkspace_(input_path, output_path);
    prepared = extractParquetRows_(workspace.base_dir, level,
                                   impl_->transition_peakgroup_max_rank,
                                   impl_->transition_peakgroup_max_pep,
                                   impl_->transition_max_isotope_overlap,
                                   impl_->transition_min_log_sn);
    const RescoreOutput output = impl_->percolator.rescore(prepared.input);
    const auto results = computeRanks_(prepared, output);
    writeParquetScores_(workspace, level, prepared, results);
    return makeSummary_(prepared);
  }
} // namespace OpenMS
