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
#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/SqliteConnector_impl.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <arrow/api.h>
#include <sqlite3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{
  using String = std::string;

  const std::array<std::string, static_cast<Size>(OpenSwathPercolatorScoring::Level::SIZE_OF_LEVEL)>
    OpenSwathPercolatorScoring::names_of_level = {"ms1", "ms2", "ms1ms2", "transition"};

  struct OpenSwathPercolatorScoring::Impl
  {
    Percolator percolator;
    Int32 transition_peakgroup_max_rank = 1;
    double transition_peakgroup_max_pep = 0.7;
    double transition_max_isotope_overlap = 0.5;
    double transition_min_log_sn = 0.0;
    std::string training_main_score = "var_xcorr_shape";
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
      double exp_rt = 0.0;
      std::optional<Int64> transition_id;
      bool decoy = false;
      Size run_rank = 0;
      Size precursor_rank = 0;
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
      std::unique_ptr<TempDir> temp_dir;
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
      StringUtils::toLower(out);
      if (level == Level::MS1MS2 && from_ms1_table)
      {
        if (StringUtils::hasPrefix(out, "var_"))
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
      if (StringUtils::hasPrefix(out, "var_ms2_"))
      {
        out = "var_" + out.substr(8);
      }
      else if (StringUtils::hasPrefix(out, "var_ms1_"))
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
      return StringUtils::toLower(out);
    }

    std::vector<String> sqliteColumnsWithPrefix_(SqliteConnector& conn,
                                                 const String& table_name,
                                                 const String& prefix)
    {
      sqlite3_stmt* stmt = nullptr;
      Internal::SqliteHelper::prepareStatement(conn, &stmt, "PRAGMA table_info('" + table_name + "');");
      std::vector<String> columns;
      Sql::SqlState state = Sql::nextRow(stmt);
      while (state == Sql::SqlState::SQL_ROW)
      {
        const String col_name = Sql::extractString(stmt, 1);
        if (StringUtils::hasPrefix(col_name, prefix))
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

    std::pair<double, double> meanAndSampleStdDev_(const std::vector<double>& values)
    {
      if (values.empty())
      {
        return {
          std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN()
        };
      }

      double mean = 0.0;
      for (const double value : values)
      {
        mean += value;
      }
      mean /= static_cast<double>(values.size());

      if (values.size() == 1)
      {
        return {mean, 0.0};
      }

      double squared_sum = 0.0;
      for (const double value : values)
      {
        const double diff = value - mean;
        squared_sum += diff * diff;
      }
      const double variance = squared_sum / static_cast<double>(values.size() - 1);
      return {mean, std::sqrt(std::max(0.0, variance))};
    }

    using ScientificPrecursor = std::tuple<std::vector<std::vector<String>>, Int64, bool>;

    struct ScientificKeys
    {
      std::map<Int64, String> runs;
      std::map<Int64, ScientificPrecursor> precursors;
      std::unordered_map<Int64, bool> precursor_decoys;
    };

    ScientificKeys scientificKeysSqlite_(SqliteConnector& conn)
    {
      ScientificKeys keys;
      sqlite3_stmt* stmt = nullptr;
      if (!conn.tableExists("RUN") || !conn.tableExists("PRECURSOR")) return keys;
      Sql::prepareStatement(conn, &stmt, "SELECT ID, FILENAME FROM RUN;");
      for (auto state = Sql::nextRow(stmt); state == Sql::SqlState::SQL_ROW; state = Sql::nextRow(stmt, state))
      {
        if (sqlite3_column_type(stmt, 1) != SQLITE_NULL)
        {
          keys.runs.emplace(Sql::extractInt64(stmt, 0), Sql::extractString(stmt, 1));
        }
      }
      sqlite3_finalize(stmt);
      std::set<Int64> incomplete_precursors;
      auto readIdentities = [&](const String& query, const int identity_columns)
      {
        Sql::prepareStatement(conn, &stmt, query);
        for (auto state = Sql::nextRow(stmt); state == Sql::SqlState::SQL_ROW; state = Sql::nextRow(stmt, state))
        {
          // NULL metadata is irrelevant until a retained feature needs this identity.
          if ((identity_columns == 1 && sqlite3_column_type(stmt, 1) == SQLITE_NULL) ||
              sqlite3_column_type(stmt, identity_columns + 1) == SQLITE_NULL ||
              sqlite3_column_type(stmt, identity_columns + 2) == SQLITE_NULL)
          {
            incomplete_precursors.insert(Sql::extractInt64(stmt, 0));
            continue;
          }
          std::vector<String> identity;
          for (int col = 1; col <= identity_columns; ++col)
          {
            identity.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ? "" : Sql::extractString(stmt, col));
          }
          if (std::all_of(identity.begin(), identity.end(), [](const String& value) { return value.empty(); }))
          {
            incomplete_precursors.insert(Sql::extractInt64(stmt, 0));
            continue;
          }
          identity.insert(identity.begin(), identity_columns == 1 ? "peptide" : "compound");
          auto& key = keys.precursors[Sql::extractInt64(stmt, 0)];
          std::get<0>(key).push_back(std::move(identity));
          std::get<1>(key) = Sql::extractInt64(stmt, identity_columns + 1);
          std::get<2>(key) = sqlite3_column_int(stmt, identity_columns + 2) != 0;
        }
        sqlite3_finalize(stmt);
      };
      if (conn.tableExists("PEPTIDE") && conn.tableExists("PRECURSOR_PEPTIDE_MAPPING"))
      {
        readIdentities("SELECT PRECURSOR.ID, PEPTIDE.MODIFIED_SEQUENCE, PRECURSOR.CHARGE, PRECURSOR.DECOY "
          "FROM PRECURSOR JOIN PRECURSOR_PEPTIDE_MAPPING ON PRECURSOR.ID = PRECURSOR_ID "
          "JOIN PEPTIDE ON PEPTIDE.ID = PEPTIDE_ID;", 1);
      }
      if (conn.tableExists("COMPOUND") && conn.tableExists("PRECURSOR_COMPOUND_MAPPING"))
      {
        // Older native OSW libraries omit ADDUCTS, as supported by TransitionPQPFile.
        const String adducts = Sql::columnExists(Sql::getNativeHandle(conn), "COMPOUND", "ADDUCTS") ?
          "COMPOUND.ADDUCTS" : "''";
        readIdentities("SELECT PRECURSOR.ID, COMPOUND.COMPOUND_NAME, COMPOUND.SUM_FORMULA, COMPOUND.SMILES, " +
          adducts + ", PRECURSOR.CHARGE, PRECURSOR.DECOY "
          "FROM PRECURSOR JOIN PRECURSOR_COMPOUND_MAPPING ON PRECURSOR.ID = PRECURSOR_ID "
          "JOIN COMPOUND ON COMPOUND.ID = COMPOUND_ID;", 4);
      }
      for (const auto id : incomplete_precursors) keys.precursors.erase(id);
      // Mapping row order and peptide IDs must not affect precursor identity.
      for (auto& [id, key] : keys.precursors)
      {
        auto& sequences = std::get<0>(key);
        std::sort(sequences.begin(), sequences.end());
        sequences.erase(std::unique(sequences.begin(), sequences.end()), sequences.end());
      }
      return keys;
    }

    ScientificKeys scientificKeysParquet_(const String& path)
    {
      ScientificKeys keys;
      const auto runs = ParquetFile::readTable(path + "/runs/runs.parquet");
      const auto run_ids = ParquetFile::getColumn(runs, "run_id");
      const auto filenames = ParquetFile::getOptionalColumn(runs, "filename");
      for (int64_t i = 0; i < runs->num_rows(); ++i)
      {
        keys.runs.emplace(ParquetFile::getInt64(run_ids, i, 0, false), ParquetFile::getString(filenames, i));
      }
      const auto precursors = ParquetFile::readTable(path + "/library/precursors.parquet");
      const auto ids = ParquetFile::getColumn(precursors, "precursor_id");
      const auto sequences = ParquetFile::getOptionalColumn(precursors, "modified_sequence");
      const auto charges = ParquetFile::getOptionalColumn(precursors, "charge");
      const auto decoys = ParquetFile::getOptionalColumn(precursors, "decoy");
      const auto names = ParquetFile::getOptionalColumn(precursors, "traml_id");
      const auto masses = ParquetFile::getOptionalColumn(precursors, "precursor_mz");
      for (int64_t i = 0; i < precursors->num_rows(); ++i)
      {
        const Int64 id = ParquetFile::getInt64(ids, i, 0, false);
        const bool decoy = ParquetFile::getBool(decoys, i, false, true);
        // Keep labels available for transition filtering even when identity is missing.
        keys.precursor_decoys.emplace(id, decoy);
        if (!charges || charges->IsNull(i) || (decoys && decoys->IsNull(i))) continue;
        const String sequence = ParquetFile::getString(sequences, i);
        std::vector<String> identity;
        if (!sequence.empty())
        {
          identity = {"peptide", sequence};
        }
        else
        {
          const String name = ParquetFile::getString(names, i);
          const double mz = ParquetFile::getDouble(masses, i, std::numeric_limits<double>::quiet_NaN(), true);
          if (name.empty() || !std::isfinite(mz)) continue;
          // Native molecule libraries carry the source compound name in traml_id, not a sequence.
          identity = {"compound", name, StringUtils::toStr(mz)};
        }
        keys.precursors.emplace(id, ScientificPrecursor{{std::move(identity)},
          ParquetFile::getInt64(charges, i, 0, false), decoy});
      }
      return keys;
    }

    PreparedInput finalizeInput_(std::vector<ScoreRow>&& rows,
                                 std::vector<String>&& feature_names,
                                 const bool transition_level,
                                 const ScientificKeys& keys)
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
      std::vector<bool> active_column(n_cols, false);
      constexpr double eps = 1e-12;

      for (Size col = 0; col < n_cols; ++col)
      {
        for (Size row = 0; row < n_rows; ++row)
        {
          if (std::isfinite(rows[row].feature_values[col]))
          {
            active_column[col] = true;
            break;
          }
        }
      }

      const Size all_null_columns = static_cast<Size>(std::count(active_column.begin(), active_column.end(), false));
      if (all_null_columns == n_cols)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No feature columns available for OpenSWATH Percolator scoring after dropping all-null columns.", "");
      }

      std::vector<ScoreRow> filtered_rows;
      filtered_rows.reserve(n_rows);
      for (auto& row : rows)
      {
        bool valid_row = true;
        for (Size col = 0; col < n_cols; ++col)
        {
          if (active_column[col] && !std::isfinite(row.feature_values[col]))
          {
            valid_row = false;
            break;
          }
        }
        if (valid_row)
        {
          filtered_rows.push_back(std::move(row));
        }
      }

      if (filtered_rows.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "No OpenSWATH rows remained after dropping rows with missing active feature values.", "");
      }

      // Assign equal integer ranks to equal scientific identities once, outside the row comparator.
      std::map<String, Size> run_ranks;
      std::map<ScientificPrecursor, Size> precursor_ranks;
      for (const auto& row : filtered_rows)
      {
        if (!keys.runs.contains(row.run_id) || keys.runs.at(row.run_id).empty() ||
            !keys.precursors.contains(row.precursor_id) || std::get<0>(keys.precursors.at(row.precursor_id)).empty())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Deterministic OpenSWATH Percolator scoring requires a run filename and peptide or compound identity "
            "for run " + StringUtils::toStr(row.run_id) + ", precursor " + StringUtils::toStr(row.precursor_id) + ".");
        }
        if (!std::isfinite(row.exp_rt))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Deterministic OpenSWATH Percolator scoring requires finite retention times.", StringUtils::toStr(row.exp_rt));
        }
        run_ranks.emplace(keys.runs.at(row.run_id), 0);
        precursor_ranks.emplace(keys.precursors.at(row.precursor_id), 0);
      }
      Size rank = 0;
      for (auto& [key, value] : run_ranks) value = rank++;
      rank = 0;
      for (auto& [key, value] : precursor_ranks) value = rank++;
      for (auto& row : filtered_rows)
      {
        row.run_rank = run_ranks.at(keys.runs.at(row.run_id));
        row.precursor_rank = precursor_ranks.at(keys.precursors.at(row.precursor_id));
      }

      std::vector<bool> keep_column(n_cols, false);
      for (Size col = 0; col < n_cols; ++col)
      {
        if (!active_column[col])
        {
          continue;
        }

        const double first = filtered_rows[0].feature_values[col];
        for (Size row = 1; row < filtered_rows.size(); ++row)
        {
          if (std::fabs(filtered_rows[row].feature_values[col] - first) > eps)
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
                                      "All OpenSWATH feature columns are constant after dropping incomplete rows.", "");
      }

      PreparedInput prepared;
      prepared.rows = std::move(filtered_rows);
      prepared.transition_level = transition_level;
      prepared.input.feature_names = kept_names;

      std::stable_sort(prepared.rows.begin(), prepared.rows.end(),
                       [&keep_column](const ScoreRow& lhs, const ScoreRow& rhs)
                       {
                         const auto left = std::tie(lhs.run_rank, lhs.precursor_rank);
                         const auto right = std::tie(rhs.run_rank, rhs.precursor_rank);
                         if (left != right) return left < right;
                         if (lhs.exp_rt != rhs.exp_rt) return lhs.exp_rt < rhs.exp_rt;
                         if (lhs.decoy != rhs.decoy) return lhs.decoy < rhs.decoy;
                         // Only finite active scores participate. Generated feature/transition IDs do not.
                         for (Size col = 0; col < lhs.feature_values.size(); ++col)
                         {
                           if (keep_column[col] && lhs.feature_values[col] != rhs.feature_values[col])
                           {
                             return lhs.feature_values[col] < rhs.feature_values[col];
                           }
                         }
                         return false;
                       });

      if (transition_level)
      {
        for (Size i = 0; i < prepared.rows.size(); ++i)
        {
          prepared.rows[i].cv_group_key = static_cast<int>(i);
        }
      }
      else
      {
        // Keep the original one-run/one-precursor partition. Scientific aliases must not merge.
        // ponytail: distinct assay aliases with identical scientific sort keys keep their partitions,
        // but alias-specific CV assignment is not guaranteed invariant without a stable assay identity.
        std::map<std::pair<Int64, Int64>, int> group_to_key;
        for (auto& row : prepared.rows)
        {
          const auto group = std::make_pair(row.run_id, row.precursor_id);
          const auto [it, inserted] = group_to_key.emplace(group, static_cast<int>(group_to_key.size()));
          row.cv_group_key = it->second;
        }
      }

      prepared.input.features.reserve(prepared.rows.size());
      prepared.input.is_decoy.reserve(prepared.rows.size());
      prepared.input.cv_group_keys.reserve(prepared.rows.size());

      for (const auto& row : prepared.rows)
      {
        std::vector<double> filtered;
        filtered.reserve(kept_names.size());
        for (Size col = 0; col < n_cols; ++col)
        {
          if (keep_column[col])
          {
            filtered.push_back(row.feature_values[col]);
          }
        }
        prepared.input.features.push_back(std::move(filtered));
        prepared.input.is_decoy.push_back(row.decoy);
        prepared.input.cv_group_keys.push_back(row.cv_group_key);
      }

      OPENMS_LOG_INFO << "OpenSWATH Percolator input cleanup: retained "
                      << prepared.rows.size() << " / " << n_rows
                      << " rows after dropping rows with missing active scores; removed "
                      << all_null_columns << " all-null and "
                      << static_cast<Size>(std::count(keep_column.begin(), keep_column.end(), false)) - all_null_columns
                      << " constant feature columns.\n";

      return prepared;
    }

    std::vector<double> computeRawModelScores_(const PreparedInput& prepared,
                                               const PercolatorModel& model)
    {
      if (model.weights.size() != prepared.input.feature_names.size() + 1)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Percolator model weights do not match the filtered OpenSWATH feature matrix.", "");
      }

      std::vector<double> raw_scores(prepared.input.features.size(), model.weights.back());
      for (Size row_idx = 0; row_idx < prepared.input.features.size(); ++row_idx)
      {
        for (Size col_idx = 0; col_idx < prepared.input.features[row_idx].size(); ++col_idx)
        {
          raw_scores[row_idx] += prepared.input.features[row_idx][col_idx] * model.weights[col_idx];
        }
      }
      return raw_scores;
    }

    std::vector<double> computeOpenSwathDScores_(const PreparedInput& prepared,
                                                 const std::vector<double>& raw_scores)
    {
      if (raw_scores.size() != prepared.rows.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "OpenSWATH raw-score vector does not match the prepared row count.", "");
      }

      std::unordered_map<int, std::vector<Size>> groups;
      for (Size i = 0; i < prepared.rows.size(); ++i)
      {
        groups[prepared.rows[i].cv_group_key].push_back(i);
      }

      std::vector<double> top_decoy_scores;
      top_decoy_scores.reserve(groups.size());
      for (auto& [_, indices] : groups)
      {
        const auto best_it = std::max_element(indices.begin(), indices.end(),
                                              [&](const Size lhs, const Size rhs)
                                              {
                                                return raw_scores[lhs] < raw_scores[rhs];
                                              });
        if (best_it != indices.end())
        {
          const Size top_index = *best_it;
          if (prepared.rows[top_index].decoy && std::isfinite(raw_scores[top_index]))
          {
            top_decoy_scores.push_back(raw_scores[top_index]);
          }
        }
      }

      const auto [decoy_mean, decoy_stddev] = meanAndSampleStdDev_(top_decoy_scores);
      if (!(std::isfinite(decoy_mean) && std::isfinite(decoy_stddev) && decoy_stddev > 0.0))
      {
        OPENMS_LOG_WARN << "OpenSWATH Percolator scoring could not standardize the classifier score with the "
                        << "top-decoy distribution; writing raw Percolator classifier scores instead.\n";
        return raw_scores;
      }

      std::vector<double> d_scores(raw_scores.size());
      for (Size i = 0; i < raw_scores.size(); ++i)
      {
        d_scores[i] = (raw_scores[i] - decoy_mean) / decoy_stddev;
      }
      return d_scores;
    }

    /**
      @brief Reads the peak-group rows for ms1/ms2/ms1ms2 scoring with one scan per table.

      Equivalent to joining FEATURE, PRECURSOR, FEATURE_MS2 and FEATURE_MS1 in SQL and ordering by
      (RUN_ID, PRECURSOR_ID, EXP_RT, FEATURE.ID), but joined and sorted here: for millions of
      peak groups SQLite answers that query with an index lookup per row, an index built on the
      fly for FEATURE_MS1 and a temporary B-tree for the order, which takes several times as
      long as reading the tables front to back.

      @return false if a feature has more than one row in a score table (the SQL join would
              return it once per combination); @p rows is then unusable and the caller runs the query.
    */
    bool readFeatureRowsByScans_(SqliteConnector& conn,
                                 const Level level,
                                 const std::vector<String>& ms2_cols,
                                 const std::vector<String>& ms1_cols,
                                 std::vector<ScoreRow>& rows)
    {
      const double nan = std::numeric_limits<double>::quiet_NaN();
      const Size n_values = ms2_cols.size() + ms1_cols.size();
      sqlite3_stmt* stmt = nullptr;

      std::unordered_map<Int64, bool> precursor_decoy;
      Sql::prepareStatement(conn, &stmt, "SELECT ID, DECOY FROM PRECURSOR;");
      for (Sql::SqlState state = Sql::nextRow(stmt); state == Sql::SqlState::SQL_ROW; state = Sql::nextRow(stmt, state))
      {
        precursor_decoy[Sql::extractInt64(stmt, 0)] = sqlite3_column_int(stmt, 1) != 0;
      }
      sqlite3_finalize(stmt);

      rows.clear();
      std::unordered_map<Int64, Size> row_of_feature;
      Sql::prepareStatement(conn, &stmt, "SELECT ID, RUN_ID, PRECURSOR_ID, EXP_RT FROM FEATURE;");
      for (Sql::SqlState state = Sql::nextRow(stmt); state == Sql::SqlState::SQL_ROW; state = Sql::nextRow(stmt, state))
      {
        ScoreRow row;
        row.feature_id = Sql::extractInt64(stmt, 0);
        row.run_id = Sql::extractInt64(stmt, 1);
        row.precursor_id = Sql::extractInt64(stmt, 2);
        row.exp_rt = sqlite3_column_type(stmt, 3) == SQLITE_NULL ? nan : sqlite3_column_double(stmt, 3);
        const auto precursor = precursor_decoy.find(row.precursor_id);
        if (precursor == precursor_decoy.end())
        {
          continue; // inner join with PRECURSOR
        }
        row.decoy = precursor->second;
        if (!row_of_feature.emplace(row.feature_id, rows.size()).second)
        {
          sqlite3_finalize(stmt);
          return false;
        }
        rows.push_back(std::move(row));
      }
      sqlite3_finalize(stmt);

      // Fills the columns [offset, offset + cols.size()) of every feature that has a row in `table`;
      // found[i] tells whether row i has one.
      auto readScores = [&](const String& table, const std::vector<String>& cols, const Size offset, std::vector<char>& found) -> bool
      {
        String sql = "SELECT FEATURE_ID";
        for (const auto& c : cols)
        {
          sql += ", " + c;
        }
        sql += " FROM " + table + ";";
        found.assign(rows.size(), 0);
        sqlite3_stmt* score_stmt = nullptr;
        Sql::prepareStatement(conn, &score_stmt, sql);
        for (Sql::SqlState state = Sql::nextRow(score_stmt); state == Sql::SqlState::SQL_ROW; state = Sql::nextRow(score_stmt, state))
        {
          const auto feature = row_of_feature.find(Sql::extractInt64(score_stmt, 0));
          if (feature == row_of_feature.end())
          {
            continue;
          }
          if (found[feature->second])
          {
            sqlite3_finalize(score_stmt);
            return false;
          }
          found[feature->second] = 1;
          std::vector<double>& values = rows[feature->second].feature_values;
          if (values.empty())
          {
            values.assign(n_values, nan);
          }
          for (Size i = 0; i < cols.size(); ++i)
          {
            const int col = static_cast<int>(i) + 1;
            values[offset + i] = sqlite3_column_type(score_stmt, col) == SQLITE_NULL ? nan : sqlite3_column_double(score_stmt, col);
          }
        }
        sqlite3_finalize(score_stmt);
        return true;
      };

      // The score table that is joined with INNER JOIN decides which features remain:
      // FEATURE_MS2 for ms2 and ms1ms2, FEATURE_MS1 for ms1. In ms1ms2, FEATURE_MS1 is a LEFT JOIN.
      std::vector<char> has_required, has_optional;
      if (level == Level::MS1)
      {
        if (!readScores("FEATURE_MS1", ms1_cols, 0, has_required)) return false;
      }
      else
      {
        if (!readScores("FEATURE_MS2", ms2_cols, 0, has_required)) return false;
        if (level == Level::MS1MS2 && !readScores("FEATURE_MS1", ms1_cols, ms2_cols.size(), has_optional)) return false;
      }
      row_of_feature.clear();

      std::vector<ScoreRow> joined;
      joined.reserve(rows.size());
      for (Size i = 0; i < rows.size(); ++i)
      {
        if (has_required[i])
        {
          joined.push_back(std::move(rows[i]));
        }
      }
      // finalizeInput_ establishes scientific order after discarding incomplete score rows.

      std::unordered_map<std::pair<Int64, Int64>, int, PairHash<std::pair<Int64, Int64>>> group_to_key;
      int next_group_key = 0;
      for (ScoreRow& row : joined)
      {
        const std::pair<Int64, Int64> group{row.run_id, row.precursor_id};
        auto it = group_to_key.find(group);
        if (it == group_to_key.end())
        {
          it = group_to_key.emplace(group, next_group_key++).first;
        }
        row.cv_group_key = it->second;
      }
      rows = std::move(joined);
      return true;
    }

    PreparedInput extractSqliteRows_(const String& osw_path,
                                     const Level level,
                                     const Int32 transition_peakgroup_max_rank,
                                     const double transition_peakgroup_max_pep,
                                     const double transition_max_isotope_overlap,
                                     const double transition_min_log_sn)
    {
      SqliteConnector conn(osw_path, SqliteConnector::SqlOpenMode::READ_ONLY);
      for (const String& table : {String("FEATURE"), level == Level::MS1 ? String("FEATURE_MS1") :
           level == Level::TRANSITION ? String("FEATURE_TRANSITION") : String("FEATURE_MS2")})
      {
        if (!conn.tableExists(table))
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "OpenSWATH Percolator scoring requires " + table + ".");
        }
      }
      if (level == Level::MS1MS2 && !conn.tableExists("FEATURE_MS1"))
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OpenSWATH Percolator scoring requires FEATURE_MS1.");
      }
      const auto keys = scientificKeysSqlite_(conn);
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

        if (readFeatureRowsByScans_(conn, level, ms2_cols, ms1_cols, rows))
        {
          return finalizeInput_(std::move(rows), std::move(feature_names), false, keys);
        }
        rows.clear(); // fall back to the query below

        String select_sql =
          "SELECT FEATURE.ID AS FEATURE_ID, FEATURE.RUN_ID, FEATURE.PRECURSOR_ID, FEATURE.EXP_RT, PRECURSOR.DECOY";
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

        Internal::SqliteHelper::prepareStatement(conn, &stmt, select_sql);
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
          row.exp_rt = sqlite3_column_type(stmt, col) == SQLITE_NULL ? std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col);
          ++col;
          row.decoy = sqlite3_column_int(stmt, col++) != 0;
          row.feature_values.reserve(feature_names.size());
          for (Size i = 0; i < ms2_cols.size(); ++i, ++col)
          {
            row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ?
              std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col));
          }
          for (Size i = 0; i < ms1_cols.size(); ++i, ++col)
          {
            row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ?
              std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col));
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
        return finalizeInput_(std::move(rows), std::move(feature_names), false, keys);
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
      const auto isotope_it = std::find(feature_names.begin(), feature_names.end(), "VAR_ISOTOPE_OVERLAP_SCORE");
      const Size isotope_col_idx = isotope_it != feature_names.end() ?
        static_cast<Size>(std::distance(feature_names.begin(), isotope_it)) :
        std::numeric_limits<Size>::max();
      const auto log_sn_it = std::find(feature_names.begin(), feature_names.end(), "VAR_LOG_SN_SCORE");
      const Size log_sn_col_idx = log_sn_it != feature_names.end() ?
        static_cast<Size>(std::distance(feature_names.begin(), log_sn_it)) :
        std::numeric_limits<Size>::max();

      String select_sql =
        "SELECT FEATURE.RUN_ID, FEATURE.ID AS FEATURE_ID, FEATURE.PRECURSOR_ID, FEATURE.EXP_RT, "
        "       FEATURE_TRANSITION.TRANSITION_ID, TRANSITION.DECOY, "
        "       SCORE_MS2.RANK, SCORE_MS2.PEP";
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
        "WHERE SCORE_MS2.RANK <= " + StringUtils::toStr(transition_peakgroup_max_rank) +
        "  AND PRECURSOR.DECOY = 0 "
        "ORDER BY FEATURE.RUN_ID, FEATURE.PRECURSOR_ID, FEATURE.EXP_RT, FEATURE_TRANSITION.TRANSITION_ID;";

      Internal::SqliteHelper::prepareStatement(conn, &stmt, select_sql);
      Sql::SqlState state = Sql::nextRow(stmt);
      int next_group_key = 0;
      std::vector<ScoreRow> base_rows;
      std::vector<ScoreRow> feature_filtered_rows;
      std::vector<ScoreRow> relaxed_rows;
      std::vector<ScoreRow> strict_rows;
      while (state == Sql::SqlState::SQL_ROW)
      {
        ScoreRow row;
        int col = 0;
        row.run_id = Sql::extractInt64(stmt, col++);
        row.feature_id = Sql::extractInt64(stmt, col++);
        row.precursor_id = Sql::extractInt64(stmt, col++);
        row.exp_rt = sqlite3_column_type(stmt, col) == SQLITE_NULL ? std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col);
        ++col;
        row.transition_id = Sql::extractInt64(stmt, col++);
        row.decoy = sqlite3_column_int(stmt, col++) != 0;
        (void)Sql::extractInt(stmt, col++);
        const double peakgroup_pep = sqlite3_column_type(stmt, col) == SQLITE_NULL ?
          std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col);
        ++col;
        row.feature_values.reserve(feature_names.size());
        for (Size i = 0; i < transition_cols.size(); ++i, ++col)
        {
          row.feature_values.push_back(sqlite3_column_type(stmt, col) == SQLITE_NULL ?
            std::numeric_limits<double>::quiet_NaN() : sqlite3_column_double(stmt, col));
        }
        row.cv_group_key = next_group_key++;
        base_rows.push_back(row);

        bool passes_transition_feature_filters = true;
        if (isotope_col_idx < row.feature_values.size())
        {
          const double isotope_overlap = row.feature_values[isotope_col_idx];
          if (!std::isfinite(isotope_overlap) || isotope_overlap > transition_max_isotope_overlap)
          {
            passes_transition_feature_filters = false;
          }
        }
        if (passes_transition_feature_filters && log_sn_col_idx < row.feature_values.size())
        {
          const double log_sn = row.feature_values[log_sn_col_idx];
          if (!std::isfinite(log_sn) || !(log_sn > transition_min_log_sn))
          {
            passes_transition_feature_filters = false;
          }
        }

        if (passes_transition_feature_filters)
        {
          feature_filtered_rows.push_back(row);
          relaxed_rows.push_back(row);
          if (std::isfinite(peakgroup_pep) && peakgroup_pep <= transition_peakgroup_max_pep)
          {
            strict_rows.push_back(std::move(row));
          }
        }
        state = Sql::nextRow(stmt, state);
      }
      sqlite3_finalize(stmt);

      if (!strict_rows.empty())
      {
        rows = std::move(strict_rows);
      }
      else if (!feature_filtered_rows.empty())
      {
        rows = std::move(relaxed_rows);
        OPENMS_LOG_WARN << "Transition-level OpenSWATH Percolator scoring found no candidates after applying the "
                        << "peak-group PEP filter (" << transition_peakgroup_max_pep << "); "
                        << "falling back to rank- and transition-feature-filtered candidates only.\n";
      }
      else
      {
        rows = std::move(base_rows);
        if (!rows.empty())
        {
          OPENMS_LOG_WARN << "Transition-level OpenSWATH Percolator scoring found no candidates after applying the "
                          << "transition feature filters; falling back to rank-filtered precursor transitions only.\n";
        }
      }

      if (rows.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH Percolator scoring found no candidate transitions. "
                                      "Run ms2/ms1ms2 scoring first and check the transition candidate filters.");
      }

      return finalizeInput_(std::move(rows), std::move(feature_names), true, keys);
    }

    PreparedInput extractParquetRows_(const String& oswpq_path,
                                      const Level level,
                                      const Int32 transition_peakgroup_max_rank,
                                      const double transition_peakgroup_max_pep,
                                      const double transition_max_isotope_overlap,
                                      const double transition_min_log_sn)
    {
      const auto keys = scientificKeysParquet_(oswpq_path);
      OpenSwathOSWParquetReader reader;
      std::vector<ScoreRow> rows;
      std::vector<String> feature_names;

      if (level == Level::MS2)
      {
        const auto result = reader.fetchPeakGroupFeatures(oswpq_path, "ms2");
        std::vector<Size> selected_cols;
        for (Size i = 0; i < result.ms2_columns.size(); ++i)
        {
          if (StringUtils::hasPrefix(result.ms2_columns[i], "var_ms2_"))
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
          row.exp_rt = result.exp_rt[i];
          row.decoy = result.decoy[i];
          row.feature_values.reserve(selected_cols.size());
          for (const Size col_idx : selected_cols)
          {
            const auto& values = result.ms2_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ?
              values[i] : std::numeric_limits<double>::quiet_NaN());
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
        return finalizeInput_(std::move(rows), std::move(feature_names), false, keys);
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
            if (StringUtils::hasPrefix(result.ms2_columns[i], "var_ms2_"))
            {
              ms2_selected.push_back(i);
              feature_names.push_back(canonicalFeatureNameParquet_(result.ms2_columns[i], level, false));
            }
          }
        }
        for (Size i = 0; i < result.ms1_columns.size(); ++i)
        {
          if (StringUtils::hasPrefix(result.ms1_columns[i], "var_ms1_"))
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
          row.exp_rt = result.exp_rt[i];
          row.decoy = result.decoy[i];
          row.feature_values.reserve(feature_names.size());
          for (const Size col_idx : ms2_selected)
          {
            const auto& values = result.ms2_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ?
              values[i] : std::numeric_limits<double>::quiet_NaN());
          }
          for (const Size col_idx : ms1_selected)
          {
            const auto& values = result.ms1_values[col_idx];
            row.feature_values.push_back(i < values.size() && std::isfinite(values[i]) ?
              values[i] : std::numeric_limits<double>::quiet_NaN());
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
        return finalizeInput_(std::move(rows), std::move(feature_names), false, keys);
      }

      const auto transition_result = reader.fetchTransitionFeatures(oswpq_path);
      feature_names.reserve(transition_result.transition_var_columns.size());
      for (const auto& feature_name : transition_result.transition_var_columns)
      {
        feature_names.push_back(StringUtils::toLowered(feature_name));
      }

      using FeatureKey = std::pair<Int64, Int64>;
      std::unordered_map<FeatureKey, std::tuple<Int64, Int32, double>, PairHash<FeatureKey>> feature_to_score;
      feature_to_score.reserve(transition_result.feature_id.size());

      for (const auto& [run_id, filename] : keys.runs)
      {
        const String run_dir = oswpq_path + "/runs/run_id=" + StringUtils::toStr(run_id);
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
          feature_to_score[std::make_pair(run_id, feature_id)] = std::make_tuple(precursor_id, rank, pep);
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
      std::vector<ScoreRow> base_rows;
      std::vector<ScoreRow> feature_filtered_rows;
      std::vector<ScoreRow> relaxed_rows;
      std::vector<ScoreRow> strict_rows;
      for (Size row_idx = 0; row_idx < transition_result.feature_id.size(); ++row_idx)
      {
        const FeatureKey feature_key = std::make_pair(transition_result.run_id[row_idx],
                                                      transition_result.feature_id[row_idx]);
        const auto score_it = feature_to_score.find(feature_key);
        if (score_it == feature_to_score.end())
        {
          continue;
        }

        const auto& [precursor_id, peakgroup_rank, peakgroup_pep] = score_it->second;
        if (peakgroup_rank <= 0 ||
            peakgroup_rank > transition_peakgroup_max_rank)
        {
          continue;
        }
        const auto precursor_it = keys.precursor_decoys.find(precursor_id);
        if (precursor_it != keys.precursor_decoys.end() && precursor_it->second)
        {
          continue;
        }

        ScoreRow out;
        out.feature_id = transition_result.feature_id[row_idx];
        out.run_id = transition_result.run_id[row_idx];
        out.precursor_id = precursor_id;
        out.exp_rt = transition_result.exp_rt[row_idx];
        out.transition_id = transition_result.transition_id[row_idx];
        out.decoy = transition_result.decoy[row_idx];
        out.cv_group_key = next_group_key++;
        out.feature_values.reserve(feature_names.size());
        for (Size col_idx = 0; col_idx < transition_result.transition_var_values.size(); ++col_idx)
        {
          out.feature_values.push_back(transition_result.transition_var_values[col_idx][row_idx]);
        }
        base_rows.push_back(out);
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
        feature_filtered_rows.push_back(out);
        relaxed_rows.push_back(out);
        if (std::isfinite(peakgroup_pep) && peakgroup_pep <= transition_peakgroup_max_pep)
        {
          strict_rows.push_back(std::move(out));
        }
      }

      if (!strict_rows.empty())
      {
        rows = std::move(strict_rows);
      }
      else if (!feature_filtered_rows.empty())
      {
        rows = std::move(relaxed_rows);
        OPENMS_LOG_WARN << "Transition-level OpenSWATH OSWPQ Percolator scoring found no candidates after applying the "
                        << "peak-group PEP filter (" << transition_peakgroup_max_pep << "); "
                        << "falling back to rank- and transition-feature-filtered candidates only.\n";
      }
      else
      {
        rows = std::move(base_rows);
        if (!rows.empty())
        {
          OPENMS_LOG_WARN << "Transition-level OpenSWATH OSWPQ Percolator scoring found no candidates after applying the "
                          << "transition feature filters; falling back to rank-filtered precursor transitions only.\n";
        }
      }

      if (rows.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Transition-level OpenSWATH OSWPQ Percolator scoring found no candidate transitions. "
                                      "Run ms2/ms1ms2 scoring first and check the transition candidate filters.");
      }

      return finalizeInput_(std::move(rows), std::move(feature_names), true, keys);
    }

    std::vector<ScoreResult> computePeakGroupResults_(const PreparedInput& prepared,
                                                      const RescoreOutput& output,
                                                      const std::vector<double>& exported_scores)
    {
      if (output.scores.size() != prepared.rows.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Percolator output size does not match OpenSWATH input size.", "");
      }
      if (output.q_values.size() != prepared.rows.size() || output.peps.size() != prepared.rows.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Percolator q-value / PEP output size does not match OpenSWATH input size.", "");
      }
      if (exported_scores.size() != prepared.rows.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Exported OpenSWATH score vector size does not match the prepared row count.", "");
      }

      std::vector<ScoreResult> results(prepared.rows.size());
      for (Size i = 0; i < prepared.rows.size(); ++i)
      {
        results[i].score = exported_scores[i];
        results[i].qvalue = output.q_values[i];
        results[i].pep = output.peps[i];
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

      const bool input_is_osw = StringUtils::hasSuffix(input_path, ".osw");
      if (input_is_osw)
      {
        if (!StringUtils::hasSuffix(output_path, ".osw"))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "SQLite OSW input can only be written to another .osw output path.",
                                        output_path);
        }
        return;
      }

      if (StringUtils::hasSuffix(output_path, ".osw"))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "OSWPQ input can only be written to another OSWPQ path.",
                                      output_path);
      }

      if (File::isDirectory(input_path))
      {
        if (StringUtils::hasSuffix(output_path, ".oswpq.zip"))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Directory OSWPQ input can only be written to another directory path.",
                                        output_path);
        }
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
        if (!(StringUtils::hasSuffix(output_path, ".oswpq") || StringUtils::hasSuffix(output_path, ".oswpq.zip")))
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
      sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
      sqlite3_stmt* insert_stmt = nullptr;
      bool transaction_started = false;
      conn.executeStatement("BEGIN IMMEDIATE TRANSACTION;");
      transaction_started = true;
      try
      {
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

        if (transition_level)
        {
          Internal::SqliteHelper::prepareStatement(conn, &insert_stmt,
            "INSERT INTO SCORE_TRANSITION (FEATURE_ID, TRANSITION_ID, SCORE, RANK, PVALUE, QVALUE, PEP) "
            "VALUES (?, ?, ?, ?, ?, ?, ?);");
        }
        else
        {
          Internal::SqliteHelper::prepareStatement(conn, &insert_stmt,
            "INSERT INTO " + table_name + " (FEATURE_ID, SCORE, RANK, PVALUE, QVALUE, PEP) "
            "VALUES (?, ?, ?, ?, ?, ?);");
        }

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
            throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "Failed to insert OpenSWATH Percolator score row: " + std::string(sqlite3_errmsg(db)));
          }
          rc = sqlite3_reset(insert_stmt);
          checkSqliteRc_(db, rc, "Failed to reset insert statement");
          rc = sqlite3_clear_bindings(insert_stmt);
          checkSqliteRc_(db, rc, "Failed to clear insert statement bindings");
        }
        sqlite3_finalize(insert_stmt);
        insert_stmt = nullptr;
        conn.executeStatement("COMMIT;");
        transaction_started = false;
      }
      catch (...)
      {
        if (insert_stmt != nullptr)
        {
          sqlite3_finalize(insert_stmt);
          insert_stmt = nullptr;
        }
        if (transaction_started)
        {
          try
          {
            conn.executeStatement("ROLLBACK;");
          }
          catch (...)
          {
          }
        }
        throw;
      }
    }

    ParquetWorkspace prepareParquetWorkspace_(const String& input_path, const String& output_path)
    {
      const bool input_is_dir = File::isDirectory(input_path);
      ParquetWorkspace workspace;
      workspace.output_is_archive = !input_is_dir;
      workspace.output_path = output_path.empty() ? input_path : output_path;

      if (input_is_dir)
      {
        if (!output_path.empty() && StringUtils::hasSuffix(workspace.output_path, ".oswpq.zip"))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Directory OSWPQ input can only be written to another directory path.", workspace.output_path);
        }
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
        if (StringUtils::hasPrefix(value, prefix))
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
          run_feature_transition_to_index[prepared.rows[i].run_id][std::make_pair(prepared.rows[i].feature_id, *prepared.rows[i].transition_id)] = i;
        }

        for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
        {
          const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
          const String ft_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet";
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
              const auto row_it = run_it->second.find(std::make_pair(feature_id, transition_id));
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
          const String features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
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

    defaults_.setValue("training_main_score", impl_->training_main_score,
                       "Peak-group score that selects the best peak group of each precursor for the Percolator training "
                       "set of ms1, ms2 and ms1ms2 rescoring (prefix '-' if lower values are better). The training set holds "
                       "these peak groups plus all peak groups of precursors drawn at random, up to "
                       "percolator:subset_max_train further peak groups (0: train on all peak groups); the model then "
                       "scores all peak groups. If the score is not among the features, all peak groups are used.");
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
    impl_->training_main_score = param_.getValue("training_main_score").toString();
    impl_->percolator.setParameters(param_.copy("percolator:", true));
  }

  std::vector<Size> OpenSwathPercolatorScoring::trainingRows(const RescoreInput& input,
                                                             const std::string& main_score,
                                                             const Size max_others,
                                                             const UInt64 seed)
  {
    const Size n_rows = input.features.size();
    std::vector<Size> all_rows(n_rows);
    std::iota(all_rows.begin(), all_rows.end(), Size(0));
    if (max_others == 0 || input.cv_group_keys.size() != n_rows)
    {
      return all_rows;
    }
    String main_name = main_score;
    double direction = 1.0;
    if (StringUtils::hasPrefix(main_name, "-"))
    {
      direction = -1.0;
      main_name = main_name.substr(1);
    }
    const auto name_it = std::find(input.feature_names.begin(), input.feature_names.end(), main_name);
    if (name_it == input.feature_names.end())
    {
      return all_rows;
    }
    const Size main_col = static_cast<Size>(name_it - input.feature_names.begin());

    // the rows of every group (a group is one precursor in one run), groups in the order of their first row
    std::unordered_map<int, Size> group_of_key;
    std::vector<std::vector<Size>> group_rows;
    for (Size i = 0; i < n_rows; ++i)
    {
      const auto inserted = group_of_key.emplace(input.cv_group_keys[i], group_rows.size());
      if (inserted.second)
      {
        group_rows.emplace_back();
      }
      group_rows[inserted.first->second].push_back(i);
    }

    // the best row of every group; ties keep the first
    std::vector<char> take(n_rows, 0);
    Size n_others = 0;
    for (const auto& rows : group_rows)
    {
      Size best = rows.front();
      for (const Size row : rows)
      {
        if (direction * input.features[row][main_col] > direction * input.features[best][main_col])
        {
          best = row;
        }
      }
      take[best] = 1;
      n_others += rows.size() - 1;
    }
    if (n_others <= max_others)
    {
      return all_rows; // everything would be taken
    }

    // complete groups in random order while their other rows fit; RandomShuffler draws the same
    // shuffle on every platform (std::shuffle does not)
    std::vector<Size> group_order(group_rows.size());
    std::iota(group_order.begin(), group_order.end(), Size(0));
    Math::RandomShuffler shuffler;
    shuffler.seed(seed);
    shuffler.portable_random_shuffle(group_order.begin(), group_order.end());
    Size n_added = 0;
    for (const Size group : group_order)
    {
      const Size n_group_others = group_rows[group].size() - 1;
      if (n_added + n_group_others > max_others)
      {
        break;
      }
      for (const Size row : group_rows[group])
      {
        take[row] = 1;
      }
      n_added += n_group_others;
    }

    std::vector<Size> rows;
    rows.reserve(group_rows.size() + n_added);
    for (Size i = 0; i < n_rows; ++i)
    {
      if (take[i])
      {
        rows.push_back(i);
      }
    }
    return rows;
  }

  namespace
  {
    /**
      @brief The training input of peak-group (ms1/ms2/ms1ms2) rescoring: the rows of
      OpenSwathPercolatorScoring::trainingRows().

      @return @p input itself if all rows are used, otherwise @p subset
    */
    const RescoreInput& trainingInput_(const RescoreInput& input, const String& main_score, const Size max_others,
                                       const UInt64 seed, RescoreInput& subset)
    {
      const Size n_rows = input.features.size();
      if (max_others > 0 && input.cv_group_keys.size() == n_rows)
      {
        String main_name = main_score;
        if (StringUtils::hasPrefix(main_name, "-"))
        {
          main_name = main_name.substr(1);
        }
        if (std::find(input.feature_names.begin(), input.feature_names.end(), main_name) == input.feature_names.end())
        {
          OPENMS_LOG_INFO << "OpenSWATH Percolator training: '" << main_name
                          << "' is not among the features, training on all " << n_rows << " rows.\n";
        }
      }
      const std::vector<Size> rows = OpenSwathPercolatorScoring::trainingRows(input, main_score, max_others, seed);
      if (rows.size() == n_rows)
      {
        return input;
      }

      subset = RescoreInput();
      subset.feature_names = input.feature_names;
      subset.features.reserve(rows.size());
      std::unordered_map<int, Size> rows_of_group;
      for (const Size i : rows) // keeps the order of the input
      {
        subset.features.push_back(input.features[i]);
        subset.is_decoy.push_back(input.is_decoy[i]);
        subset.cv_group_keys.push_back(input.cv_group_keys[i]);
        if (!input.scan_numbers.empty()) subset.scan_numbers.push_back(input.scan_numbers[i]);
        if (!input.spec_file_numbers.empty()) subset.spec_file_numbers.push_back(input.spec_file_numbers[i]);
        if (!input.exp_masses.empty()) subset.exp_masses.push_back(input.exp_masses[i]);
        if (!input.calc_masses.empty()) subset.calc_masses.push_back(input.calc_masses[i]);
        ++rows_of_group[input.cv_group_keys[i]];
      }
      const Size n_drawn = static_cast<Size>(std::count_if(rows_of_group.begin(), rows_of_group.end(),
                                                           [](const auto& group) { return group.second > 1; }));
      OPENMS_LOG_INFO << "OpenSWATH Percolator training on " << rows.size() << " of " << n_rows
                      << " peak groups: the best of each of " << rows_of_group.size() << " precursors by " << main_score
                      << ", and all peak groups of " << n_drawn << " of these precursors drawn at random.\n";
      return subset;
    }
  }

  OpenSwathPercolatorScoring::ScoreSummary OpenSwathPercolatorScoring::score(const std::string& input_path,
                                                                              const Level level,
                                                                              const std::string& output_path)
  {
    if (!File::exists(input_path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_path);
    }
    validateOutputPath_(input_path, output_path);

    OPENMS_LOG_INFO << "OpenSWATH Percolator scoring (" << levelName_(level)
                    << ") on '" << input_path << "'\n";

    // Peak groups are trained on a subset (see trainingInput_); transitions on all rows.
    const Param percolator_param = impl_->percolator.getParameters();
    const Size max_others = static_cast<Size>(static_cast<int>(percolator_param.getValue("subset_max_train")));
    const UInt64 seed = static_cast<UInt64>(static_cast<int>(percolator_param.getValue("seed")));
    auto trainingInputForLevel_ = [&](const PreparedInput& input, RescoreInput& subset) -> const RescoreInput&
    {
      if (level == Level::TRANSITION)
      {
        return input.input;
      }
      return trainingInput_(input.input, impl_->training_main_score, max_others, seed, subset);
    };

    PreparedInput prepared;
    if (StringUtils::hasSuffix(input_path, ".osw"))
    {
      const String target_path = prepareSqliteOutput_(input_path, output_path);
      prepared = extractSqliteRows_(target_path, level,
                                    impl_->transition_peakgroup_max_rank,
                                    impl_->transition_peakgroup_max_pep,
                                    impl_->transition_max_isotope_overlap,
                                    impl_->transition_min_log_sn);
      RescoreInput training_subset;
      const PercolatorModel model = impl_->percolator.train(trainingInputForLevel_(prepared, training_subset));
      const RescoreOutput output = impl_->percolator.score(prepared.input, model);
      const auto exported_scores = computeOpenSwathDScores_(prepared, computeRawModelScores_(prepared, model));
      const auto results = computePeakGroupResults_(prepared, output, exported_scores);
      writeSqliteScores_(target_path, level, prepared, results);
      return makeSummary_(prepared);
    }

    ParquetWorkspace workspace = prepareParquetWorkspace_(input_path, output_path);
    prepared = extractParquetRows_(workspace.base_dir, level,
                                   impl_->transition_peakgroup_max_rank,
                                   impl_->transition_peakgroup_max_pep,
                                   impl_->transition_max_isotope_overlap,
                                   impl_->transition_min_log_sn);
    RescoreInput training_subset;
    const PercolatorModel model = impl_->percolator.train(trainingInputForLevel_(prepared, training_subset));
    const RescoreOutput output = impl_->percolator.score(prepared.input, model);
    const auto exported_scores = computeOpenSwathDScores_(prepared, computeRawModelScores_(prepared, model));
    const auto results = computePeakGroupResults_(prepared, output, exported_scores);
    writeParquetScores_(workspace, level, prepared, results);
    return makeSummary_(prepared);
  }
} // namespace OpenMS
