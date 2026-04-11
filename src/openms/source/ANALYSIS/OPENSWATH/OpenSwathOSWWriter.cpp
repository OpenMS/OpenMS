// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>

#include <OpenMS/FORMAT/SqliteConnector.h>

#include <sqlite3.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <sstream>

namespace OpenMS
{
  namespace
  {
    const std::vector<const char*>& featureColumns()
    {
      static const std::vector<const char*> columns = {
        "ID", "RUN_ID", "PRECURSOR_ID", "EXP_RT", "EXP_IM", "NORM_RT", "DELTA_RT",
        "LEFT_WIDTH", "RIGHT_WIDTH", "EXP_IM_LEFTWIDTH", "EXP_IM_RIGHTWIDTH"
      };
      return columns;
    }

    const std::vector<const char*>& featureMS1Columns()
    {
      static const std::vector<const char*> columns = {
        "FEATURE_ID", "AREA_INTENSITY", "APEX_INTENSITY", "EXP_IM", "DELTA_IM",
        "VAR_MASSDEV_SCORE", "VAR_IM_MS1_DELTA_SCORE", "VAR_MI_SCORE",
        "VAR_MI_CONTRAST_SCORE", "VAR_MI_COMBINED_SCORE", "VAR_ISOTOPE_CORRELATION_SCORE",
        "VAR_ISOTOPE_OVERLAP_SCORE", "VAR_XCORR_COELUTION", "VAR_XCORR_COELUTION_CONTRAST",
        "VAR_XCORR_COELUTION_COMBINED", "VAR_XCORR_SHAPE", "VAR_XCORR_SHAPE_CONTRAST",
        "VAR_XCORR_SHAPE_COMBINED"
      };
      return columns;
    }

    const std::vector<const char*>& featureMS2Columns()
    {
      static const std::vector<const char*> columns = {
        "FEATURE_ID", "AREA_INTENSITY", "TOTAL_AREA_INTENSITY", "APEX_INTENSITY",
        "EXP_IM", "EXP_IM_LEFTWIDTH", "EXP_IM_RIGHTWIDTH", "DELTA_IM", "TOTAL_MI",
        "VAR_BSERIES_SCORE", "VAR_DOTPROD_SCORE", "VAR_INTENSITY_SCORE",
        "VAR_ISOTOPE_CORRELATION_SCORE", "VAR_ISOTOPE_OVERLAP_SCORE", "VAR_LIBRARY_CORR",
        "VAR_LIBRARY_DOTPROD", "VAR_LIBRARY_MANHATTAN", "VAR_LIBRARY_RMSD",
        "VAR_LIBRARY_ROOTMEANSQUARE", "VAR_LIBRARY_SANGLE", "VAR_LOG_SN_SCORE",
        "VAR_MANHATTAN_SCORE", "VAR_MASSDEV_SCORE", "VAR_MASSDEV_SCORE_WEIGHTED",
        "VAR_MI_SCORE", "VAR_MI_WEIGHTED_SCORE", "VAR_MI_RATIO_SCORE", "VAR_NORM_RT_SCORE",
        "VAR_XCORR_COELUTION", "VAR_XCORR_COELUTION_WEIGHTED", "VAR_XCORR_SHAPE",
        "VAR_XCORR_SHAPE_WEIGHTED", "VAR_YSERIES_SCORE", "VAR_ELUTION_MODEL_FIT_SCORE",
        "VAR_IM_XCORR_SHAPE", "VAR_IM_XCORR_COELUTION", "VAR_IM_DELTA_SCORE",
        "VAR_IM_LOG_INTENSITY"
      };
      return columns;
    }

    const std::vector<const char*>& featurePrecursorColumns()
    {
      static const std::vector<const char*> columns = {
        "FEATURE_ID", "ISOTOPE", "AREA_INTENSITY", "APEX_INTENSITY"
      };
      return columns;
    }

    const std::vector<const char*>& featureTransitionColumns()
    {
      static const std::vector<const char*> columns = {
        "FEATURE_ID", "TRANSITION_ID", "AREA_INTENSITY", "TOTAL_AREA_INTENSITY",
        "APEX_RT", "APEX_INTENSITY", "RT_FWHM", "MASSERROR_PPM", "TOTAL_MI",
        "VAR_INTENSITY_SCORE", "VAR_INTENSITY_RATIO_SCORE", "VAR_LOG_INTENSITY",
        "VAR_XCORR_COELUTION", "VAR_XCORR_SHAPE", "VAR_LOG_SN_SCORE", "VAR_MASSDEV_SCORE",
        "VAR_MI_SCORE", "VAR_MI_RATIO_SCORE", "VAR_ISOTOPE_CORRELATION_SCORE",
        "VAR_ISOTOPE_OVERLAP_SCORE", "EXP_IM", "EXP_IM_LEFTWIDTH", "EXP_IM_RIGHTWIDTH",
        "DELTA_IM", "VAR_IM_DELTA_SCORE", "VAR_IM_LOG_INTENSITY",
        "VAR_IM_XCORR_COELUTION_CONTRAST", "VAR_IM_XCORR_SHAPE_CONTRAST",
        "VAR_IM_XCORR_COELUTION_COMBINED", "VAR_IM_XCORR_SHAPE_COMBINED",
        "START_POSITION_AT_5", "END_POSITION_AT_5", "START_POSITION_AT_10",
        "END_POSITION_AT_10", "START_POSITION_AT_50", "END_POSITION_AT_50",
        "TOTAL_WIDTH", "TAILING_FACTOR", "ASYMMETRY_FACTOR", "SLOPE_OF_BASELINE",
        "BASELINE_DELTA_2_HEIGHT", "POINTS_ACROSS_BASELINE", "POINTS_ACROSS_HALF_HEIGHT"
      };
      return columns;
    }

    bool isNullLiteral(const String& value)
    {
      return value.empty() || value == "NULL" || value == "nan" || value == "-nan" ||
             value == "NaN" || value == "-NaN";
    }

    String oswValue(const String& value)
    {
      return isNullLiteral(value) ? String("NULL") : value;
    }

    String oswValue(const DataValue& value)
    {
      if (value.isEmpty())
      {
        return "NULL";
      }
      return oswValue(value.toString());
    }

    String oswValue(const double value)
    {
      if (std::isnan(value))
      {
        return "NULL";
      }
      return String(value);
    }

    String oswValue(const UInt64 value)
    {
      return String(static_cast<Int64>(value));
    }

    String oswValueAt(const std::vector<String>& values, Size index)
    {
      if (index >= values.size())
      {
        return "NULL";
      }
      return oswValue(values[index]);
    }

    std::vector<String> makeTransitionRow(const String& feature_id)
    {
      std::vector<String> row(featureTransitionColumns().size(), "NULL");
      row[0] = feature_id;
      return row;
    }

    String makeInsertStatement(const char* table, const std::vector<const char*>& columns)
    {
      std::stringstream sql;
      sql << "INSERT INTO " << table << " (";
      for (Size i = 0; i < columns.size(); ++i)
      {
        if (i != 0)
        {
          sql << ", ";
        }
        sql << columns[i];
      }
      sql << ") VALUES (";
      for (Size i = 0; i < columns.size(); ++i)
      {
        if (i != 0)
        {
          sql << ", ";
        }
        sql << "?";
      }
      sql << ");";
      return String(sql.str());
    }

    void throwSQLiteError(sqlite3* db, const String& operation)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        operation + " failed: " + sqlite3_errmsg(db));
    }

    void bindRow(sqlite3* db, sqlite3_stmt* stmt, const std::vector<String>& row, Size expected_columns)
    {
      if (row.size() != expected_columns)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OSW row has " + String(row.size()) + " values, expected " + String(expected_columns));
      }

      for (Size i = 0; i < row.size(); ++i)
      {
        int rc = SQLITE_OK;
        if (isNullLiteral(row[i]))
        {
          rc = sqlite3_bind_null(stmt, static_cast<int>(i + 1));
        }
        else
        {
          rc = sqlite3_bind_text(stmt, static_cast<int>(i + 1), row[i].c_str(),
                                 static_cast<int>(row[i].size()), SQLITE_TRANSIENT);
        }
        if (rc != SQLITE_OK)
        {
          throwSQLiteError(db, "sqlite3_bind");
        }
      }

      int rc = sqlite3_step(stmt);
      if (rc != SQLITE_DONE)
      {
        throwSQLiteError(db, "sqlite3_step");
      }
      rc = sqlite3_reset(stmt);
      if (rc != SQLITE_OK)
      {
        throwSQLiteError(db, "sqlite3_reset");
      }
      rc = sqlite3_clear_bindings(stmt);
      if (rc != SQLITE_OK)
      {
        throwSQLiteError(db, "sqlite3_clear_bindings");
      }
    }

    void writeTableRows(sqlite3* db, const char* table, const std::vector<const char*>& columns,
                        const std::vector<std::vector<String>>& rows)
    {
      if (rows.empty())
      {
        return;
      }

      sqlite3_stmt* stmt = nullptr;
      SqliteConnector::prepareStatement(db, &stmt, makeInsertStatement(table, columns));
      try
      {
        for (const auto& row : rows)
        {
          bindRow(db, stmt, row, columns.size());
        }
      }
      catch (...)
      {
        sqlite3_finalize(stmt);
        throw;
      }
      sqlite3_finalize(stmt);
    }

    void appendRows(std::vector<std::vector<String>>& target, std::vector<std::vector<String>>& source)
    {
      target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
      source.clear();
    }

    void appendInsertSQL(std::stringstream& sql, const char* table, const std::vector<const char*>& columns,
                         const std::vector<std::vector<String>>& rows)
    {
      for (const auto& row : rows)
      {
        sql << "INSERT INTO " << table << " (";
        for (Size i = 0; i < columns.size(); ++i)
        {
          if (i != 0)
          {
            sql << ", ";
          }
          sql << columns[i];
        }
        sql << ") VALUES (";
        for (Size i = 0; i < row.size(); ++i)
        {
          if (i != 0)
          {
            sql << ", ";
          }
          sql << oswValue(row[i]);
        }
        sql << "); ";
      }
    }
  }

  void OpenSwathOSWWriter::OSWData::reserve(Size row_count)
  {
    feature_rows.reserve(row_count);
    feature_ms1_rows.reserve(row_count);
    feature_ms2_rows.reserve(row_count);
    feature_precursor_rows.reserve(row_count);
    feature_transition_rows.reserve(row_count);
  }

  void OpenSwathOSWWriter::OSWData::append(OSWData&& rhs)
  {
    appendRows(feature_rows, rhs.feature_rows);
    appendRows(feature_ms1_rows, rhs.feature_ms1_rows);
    appendRows(feature_ms2_rows, rhs.feature_ms2_rows);
    appendRows(feature_precursor_rows, rhs.feature_precursor_rows);
    appendRows(feature_transition_rows, rhs.feature_transition_rows);
  }

  bool OpenSwathOSWWriter::OSWData::empty() const
  {
    return feature_rows.empty() && feature_ms1_rows.empty() && feature_ms2_rows.empty() &&
           feature_precursor_rows.empty() && feature_transition_rows.empty();
  }

  OpenSwathOSWWriter::OpenSwathOSWWriter(const String& output_filename, bool uis_scores) :
    output_filename_(output_filename),
    doWrite_(!output_filename.empty()),
    enable_uis_scoring_(uis_scores)
  {}

  bool OpenSwathOSWWriter::isActive() const
  {
    return doWrite_;
  }

  void OpenSwathOSWWriter::writeHeader()
  {
    // Open database
    SqliteConnector conn(output_filename_);

    // Create SQL structure
    const char * create_sql =
      "CREATE TABLE RUN(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "FILENAME TEXT NOT NULL); " \

      "CREATE TABLE FEATURE(" \
      "ID INT PRIMARY KEY NOT NULL," \
      "RUN_ID INT NOT NULL," \
      "PRECURSOR_ID INT NOT NULL," \
      "EXP_RT REAL NOT NULL," \
      "EXP_IM REAL, " \
      "NORM_RT REAL NOT NULL," \
      "DELTA_RT REAL NOT NULL," \
      "LEFT_WIDTH REAL NOT NULL," \
      "RIGHT_WIDTH REAL NOT NULL," \
      "EXP_IM_LEFTWIDTH REAL," \
      "EXP_IM_RIGHTWIDTH REAL); " \

      // MS1-level scores
      "CREATE TABLE FEATURE_MS1(" \
      "FEATURE_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "EXP_IM REAL," \
      "DELTA_IM REAL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_CONTRAST_SCORE REAL NULL," \
      "VAR_MI_COMBINED_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL," \
      "VAR_IM_MS1_DELTA_SCORE REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_COELUTION_CONTRAST REAL NULL," \
      "VAR_XCORR_COELUTION_COMBINED REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_XCORR_SHAPE_CONTRAST REAL NULL," \
      "VAR_XCORR_SHAPE_COMBINED REAL NULL); " \

      // MS2-level scores
      "CREATE TABLE FEATURE_MS2(" \
      "FEATURE_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "TOTAL_AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "EXP_IM REAL," \
      "EXP_IM_LEFTWIDTH REAL," \
      "EXP_IM_RIGHTWIDTH REAL," \
      "DELTA_IM REAL," \
      "TOTAL_MI REAL NULL," \
      "VAR_BSERIES_SCORE REAL NULL," \
      "VAR_DOTPROD_SCORE REAL NULL," \
      "VAR_INTENSITY_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL," \
      "VAR_LIBRARY_CORR REAL NULL," \
      "VAR_LIBRARY_DOTPROD REAL NULL," \
      "VAR_LIBRARY_MANHATTAN REAL NULL," \
      "VAR_LIBRARY_RMSD REAL NULL," \
      "VAR_LIBRARY_ROOTMEANSQUARE REAL NULL," \
      "VAR_LIBRARY_SANGLE REAL NULL," \
      "VAR_LOG_SN_SCORE REAL NULL," \
      "VAR_MANHATTAN_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE_WEIGHTED REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_WEIGHTED_SCORE REAL NULL," \
      "VAR_MI_RATIO_SCORE REAL NULL," \
      "VAR_NORM_RT_SCORE REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_COELUTION_WEIGHTED REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_XCORR_SHAPE_WEIGHTED REAL NULL," \
      "VAR_YSERIES_SCORE REAL NULL," \
      "VAR_ELUTION_MODEL_FIT_SCORE REAL NULL," \

      "VAR_IM_XCORR_SHAPE REAL NULL," \
      "VAR_IM_XCORR_COELUTION REAL NULL," \
      "VAR_IM_DELTA_SCORE REAL NULL," \
      "VAR_IM_LOG_INTENSITY REAL NULL);" \

      "CREATE TABLE FEATURE_PRECURSOR(" \
      "FEATURE_ID INT NOT NULL," \
      "ISOTOPE INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "APEX_INTENSITY REAL NOT NULL);" \

      // Transition-level scores
      "CREATE TABLE FEATURE_TRANSITION(" \
      "FEATURE_ID INT NOT NULL," \
      "TRANSITION_ID INT NOT NULL," \
      "AREA_INTENSITY REAL NOT NULL," \
      "TOTAL_AREA_INTENSITY REAL NOT NULL," \
      "APEX_RT REAL NULL," \
      "APEX_INTENSITY REAL NOT NULL," \
      "RT_FWHM REAL NOT NULL," \
      "MASSERROR_PPM REAL NULL,"
      "TOTAL_MI REAL NULL," \
      "VAR_INTENSITY_SCORE REAL NULL," \
      "VAR_INTENSITY_RATIO_SCORE REAL NULL," \
      "VAR_LOG_INTENSITY REAL NULL," \
      "VAR_XCORR_COELUTION REAL NULL," \
      "VAR_XCORR_SHAPE REAL NULL," \
      "VAR_LOG_SN_SCORE REAL NULL," \
      "VAR_MASSDEV_SCORE REAL NULL," \
      "VAR_MI_SCORE REAL NULL," \
      "VAR_MI_RATIO_SCORE REAL NULL," \
      "VAR_ISOTOPE_CORRELATION_SCORE REAL NULL," \
      "VAR_ISOTOPE_OVERLAP_SCORE REAL NULL, " \
      "EXP_IM REAL NULL," \
      "EXP_IM_LEFTWIDTH REAL," \
      "EXP_IM_RIGHTWIDTH REAL," \
      "DELTA_IM REAL NULL," \
      "VAR_IM_DELTA_SCORE REAL NULL," \
      "VAR_IM_LOG_INTENSITY REAL NULL," \
      "VAR_IM_XCORR_COELUTION_CONTRAST, " \
      "VAR_IM_XCORR_SHAPE_CONTRAST, " \
      "VAR_IM_XCORR_COELUTION_COMBINED, " \
      "VAR_IM_XCORR_SHAPE_COMBINED, " \
      "START_POSITION_AT_5 REAL NULL, " \
      "END_POSITION_AT_5 REAL NULL, " \
      "START_POSITION_AT_10 REAL NULL, " \
      "END_POSITION_AT_10 REAL NULL, " \
      "START_POSITION_AT_50 REAL NULL, " \
      "END_POSITION_AT_50 REAL NULL, " \
      "TOTAL_WIDTH REAL NULL, " \
      "TAILING_FACTOR REAL NULL, " \
      "ASYMMETRY_FACTOR REAL NULL, " \
      "SLOPE_OF_BASELINE REAL NULL, " \
      "BASELINE_DELTA_2_HEIGHT REAL NULL, " \
      "POINTS_ACROSS_BASELINE REAL NULL, " \
      "POINTS_ACROSS_HALF_HEIGHT REAL NULL); ";


    // Execute SQL create statement
    conn.executeStatement(create_sql);
  }


  void OpenSwathOSWWriter::addRun(const UInt64 run_id, const String& input_filename)
  {
    if (!doWrite_) return;
    SqliteConnector conn(output_filename_);
    const UInt64 rid = Internal::SqliteHelper::clearSignBit(run_id);
    // Bind RUN.ID explicitly as int64 so SQLite stores it with INTEGER storage
    // class. executeBindStatement() unconditionally uses sqlite3_bind_blob for
    // every parameter, which stores integers as BLOB and breaks
    // "FEATURE.RUN_ID = RUN.ID" joins.
    sqlite3_stmt* stmt = nullptr;
    SqliteConnector::prepareStatement(conn.getDB(), &stmt,
        "INSERT INTO RUN (ID, FILENAME) VALUES (?, ?);");
    int rc = sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(rid));
    if (rc != SQLITE_OK)
    {
      sqlite3_finalize(stmt);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("sqlite3_bind_int64 failed: ") + sqlite3_errmsg(conn.getDB()));
    }

    rc = sqlite3_bind_text(stmt, 2, input_filename.c_str(),
      static_cast<int>(input_filename.size()), SQLITE_TRANSIENT);
    if (rc != SQLITE_OK)
    {
      sqlite3_finalize(stmt);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("sqlite3_bind_text failed: ") + sqlite3_errmsg(conn.getDB()));
    }

    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    if (rc != SQLITE_DONE)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        String("sqlite3_step failed: ") + sqlite3_errmsg(conn.getDB()));
    }
    run_id_ = rid;
  }

  void OpenSwathOSWWriter::setRunId(const UInt64 run_id)
  {
    run_id_ = Internal::SqliteHelper::clearSignBit(run_id);
  }

  String OpenSwathOSWWriter::getScore(const Feature& feature, const std::string& score_name) const
  {
    const DataValue& dv = feature.getMetaValue(score_name);
    if (dv.isEmpty()) return "NULL";

    if (dv.valueType() == DataValue::DOUBLE_VALUE)
    {
      if (std::isnan(static_cast<double>(dv))) return "NULL";
    }

    String score = dv.toString();
    if (score == "nan" || score == "-nan" || score == "NaN" || score == "-NaN") return "NULL";
    return score;
  }

  std::vector<String> OpenSwathOSWWriter::getSeparateScore(const Feature& feature, const std::string& score_name) const
  {
    std::vector<String> separated_scores;

    const DataValue& dv = feature.getMetaValue(score_name);
    if (!dv.isEmpty())
    {
      if (dv.valueType() == DataValue::STRING_LIST)
      {
        separated_scores = dv.toStringList();
        std::transform(separated_scores.begin(), separated_scores.end(), separated_scores.begin(),
                       [](const String& value) { return oswValue(value); });
      }
      else if (dv.valueType() == DataValue::INT_LIST)
      {
        std::vector<int> int_separated_scores = dv.toIntList();
        std::transform(int_separated_scores.begin(), int_separated_scores.end(), std::back_inserter(separated_scores), [](const int& num) { return String(num); });
      }
      else if (dv.valueType() == DataValue::DOUBLE_LIST)
      {
        std::vector<double> double_separated_scores = dv.toDoubleList();
        std::transform(double_separated_scores.begin(), double_separated_scores.end(), std::back_inserter(separated_scores), [](const double& num) { return oswValue(num); });
      }
      else
      {
        separated_scores.push_back(oswValue(dv.toString()));
      }
    }

    return separated_scores;
  }

  OpenSwathOSWWriter::OSWData OpenSwathOSWWriter::prepareRows(const OpenSwath::LightCompound& /* pep */,
                                                              const OpenSwath::LightTransition* /* transition */,
                                                              const FeatureMap& output,
                                                              const String& id) const
  {
    OSWData rows;
    rows.reserve(output.size());
    std::vector<std::vector<String>> ms2_transition_rows;
    std::vector<std::vector<String>> uis_transition_rows;
    ms2_transition_rows.reserve(output.size());
    uis_transition_rows.reserve(output.size());

    for (const auto& feature_it : output)
    {
      const String feature_id = oswValue(Internal::SqliteHelper::clearSignBit(feature_it.getUniqueId()));

      const auto& masserror_ppm = feature_it.metaValueExists("masserror_ppm") ? getSeparateScore(feature_it, "masserror_ppm") : std::vector<String>();

      const auto& subordinates = feature_it.getSubordinates();
      for (Size i=0; i < subordinates.size(); i++)
      {
        const auto& sub_it = subordinates[i];
        if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2")
        {
          String total_mi = oswValue(sub_it.getMetaValue("total_mi")); // total_mi is not guaranteed to be set
          String masserror_ppm_query = oswValueAt(masserror_ppm, i); // masserror_ppm is not guaranteed to be set

          std::vector<String> transition_row = makeTransitionRow(feature_id);
          transition_row[1] = oswValue(sub_it.getMetaValue("native_id"));
          transition_row[2] = oswValue(sub_it.getIntensity());
          transition_row[3] = oswValue(sub_it.getMetaValue("total_xic"));
          transition_row[4] = oswValue(sub_it.getMetaValue("peak_apex_position"));
          transition_row[5] = oswValue(sub_it.getMetaValue("peak_apex_int"));
          transition_row[6] = oswValue(sub_it.getMetaValue("width_at_50"));
          transition_row[7] = oswValue(masserror_ppm_query);
          transition_row[8] = oswValue(total_mi);

          if (sub_it.metaValueExists("start_position_at_5"))
          {
            transition_row[30] = oswValue(sub_it.getMetaValue("start_position_at_5"));
            transition_row[31] = oswValue(sub_it.getMetaValue("end_position_at_5"));
            transition_row[32] = oswValue(sub_it.getMetaValue("start_position_at_10"));
            transition_row[33] = oswValue(sub_it.getMetaValue("end_position_at_10"));
            transition_row[34] = oswValue(sub_it.getMetaValue("start_position_at_50"));
            transition_row[35] = oswValue(sub_it.getMetaValue("end_position_at_50"));
            transition_row[36] = oswValue(sub_it.getMetaValue("total_width"));
            transition_row[37] = oswValue(sub_it.getMetaValue("tailing_factor"));
            transition_row[38] = oswValue(sub_it.getMetaValue("asymmetry_factor"));
            transition_row[39] = oswValue(sub_it.getMetaValue("slope_of_baseline"));
            transition_row[40] = oswValue(sub_it.getMetaValue("baseline_delta_2_height"));
            transition_row[41] = oswValue(sub_it.getMetaValue("points_across_baseline"));
            transition_row[42] = oswValue(sub_it.getMetaValue("points_across_half_height"));
          }
          ms2_transition_rows.push_back(std::move(transition_row));
        }
        else if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS1" && sub_it.getIntensity() > 0.0)
        {
          std::vector<String> precursor_id;
          oswValue(sub_it.getMetaValue("native_id")).split(OpenMS::String("Precursor_i"), precursor_id);
          rows.feature_precursor_rows.push_back({
            feature_id,
            precursor_id.size() > 1 ? oswValue(precursor_id[1]) : String("NULL"),
            oswValue(sub_it.getIntensity()),
            oswValue(sub_it.getMetaValue("peak_apex_int"))
          });
        }
      }

      // these will be missing if RT scoring is disabled
      double norm_rt = -1, delta_rt = -1;
      if (feature_it.metaValueExists("norm_RT") ) norm_rt = feature_it.getMetaValue("norm_RT");
      if (feature_it.metaValueExists("delta_rt") ) delta_rt = feature_it.getMetaValue("delta_rt");

      rows.feature_rows.push_back({
        feature_id,
        oswValue(run_id_),
        oswValue(id),
        oswValue(feature_it.getRT()),
        getScore(feature_it, "im_drift"),
        oswValue(norm_rt),
        oswValue(delta_rt),
        oswValue(feature_it.getMetaValue("leftWidth")),
        oswValue(feature_it.getMetaValue("rightWidth")),
        getScore(feature_it, "im_drift_left"),
        getScore(feature_it, "im_drift_right")
      });

      rows.feature_ms2_rows.push_back({
        feature_id,
        oswValue(feature_it.getIntensity()),
        getScore(feature_it, "total_xic"),
        getScore(feature_it, "peak_apices_sum"),
        getScore(feature_it, "im_drift"),
        getScore(feature_it, "im_drift_left"),
        getScore(feature_it, "im_drift_right"),
        getScore(feature_it, "im_delta"),
        getScore(feature_it, "total_mi"),
        getScore(feature_it, "var_bseries_score"),
        getScore(feature_it, "var_dotprod_score"),
        getScore(feature_it, "var_intensity_score"),
        getScore(feature_it, "var_isotope_correlation_score"),
        getScore(feature_it, "var_isotope_overlap_score"),
        getScore(feature_it, "var_library_corr"),
        getScore(feature_it, "var_library_dotprod"),
        getScore(feature_it, "var_library_manhattan"),
        getScore(feature_it, "var_library_rmsd"),
        getScore(feature_it, "var_library_rootmeansquare"),
        getScore(feature_it, "var_library_sangle"),
        getScore(feature_it, "var_log_sn_score"),
        getScore(feature_it, "var_manhatt_score"),
        getScore(feature_it, "var_massdev_score"),
        getScore(feature_it, "var_massdev_score_weighted"),
        getScore(feature_it, "var_mi_score"),
        getScore(feature_it, "var_mi_weighted_score"),
        getScore(feature_it, "var_mi_ratio_score"),
        getScore(feature_it, "var_norm_rt_score"),
        getScore(feature_it, "var_xcorr_coelution"),
        getScore(feature_it, "var_xcorr_coelution_weighted"),
        getScore(feature_it, "var_xcorr_shape"),
        getScore(feature_it, "var_xcorr_shape_weighted"),
        getScore(feature_it, "var_yseries_score"),
        getScore(feature_it, "var_elution_model_fit_score"),
        getScore(feature_it, "var_im_xcorr_shape"),
        getScore(feature_it, "var_im_xcorr_coelution"),
        getScore(feature_it, "var_im_delta_score"),
        getScore(feature_it, "im_log_intensity")
      });

      bool enable_ms1 = feature_it.metaValueExists("var_ms1_ppm_diff");
      if (enable_ms1) // only write MS1 scores if they are present
      {
        rows.feature_ms1_rows.push_back({
          feature_id,
          getScore(feature_it, "ms1_area_intensity"),
          getScore(feature_it, "ms1_apex_intensity"),
          getScore(feature_it, "im_ms1_drift"),
          getScore(feature_it, "im_ms1_delta"),
          getScore(feature_it, "var_ms1_ppm_diff"),
          getScore(feature_it, "var_im_ms1_delta_score"),
          getScore(feature_it, "var_ms1_mi_score"),
          getScore(feature_it, "var_ms1_mi_contrast_score"),
          getScore(feature_it, "var_ms1_mi_combined_score"),
          getScore(feature_it, "var_ms1_isotope_correlation"),
          getScore(feature_it, "var_ms1_isotope_overlap"),
          getScore(feature_it, "var_ms1_xcorr_coelution"),
          getScore(feature_it, "var_ms1_xcorr_coelution_contrast"),
          getScore(feature_it, "var_ms1_xcorr_coelution_combined"),
          getScore(feature_it, "var_ms1_xcorr_shape"),
          getScore(feature_it, "var_ms1_xcorr_shape_contrast"),
          getScore(feature_it, "var_ms1_xcorr_shape_combined")
        });
      }

      if (enable_uis_scoring_)
      {
        auto id_target_transition_names = getSeparateScore(feature_it, "id_target_transition_names");
        auto id_target_area_intensity = getSeparateScore(feature_it, "id_target_area_intensity");
        auto id_target_total_area_intensity = getSeparateScore(feature_it, "id_target_total_area_intensity");
        auto id_target_apex_intensity = getSeparateScore(feature_it, "id_target_apex_intensity");
        auto id_target_peak_apex_position = getSeparateScore(feature_it, "id_target_peak_apex_position");
        auto id_target_peak_fwhm = getSeparateScore(feature_it, "id_target_width_at_50");
        auto id_target_total_mi = getSeparateScore(feature_it, "id_target_total_mi");
        auto id_target_intensity_score = getSeparateScore(feature_it, "id_target_intensity_score");
        auto id_target_intensity_ratio_score = getSeparateScore(feature_it, "id_target_intensity_ratio_score");
        auto id_target_log_intensity = getSeparateScore(feature_it, "id_target_ind_log_intensity");
        auto id_target_ind_xcorr_coelution = getSeparateScore(feature_it, "id_target_ind_xcorr_coelution");
        auto id_target_ind_xcorr_shape = getSeparateScore(feature_it, "id_target_ind_xcorr_shape");
        auto id_target_ind_log_sn_score = getSeparateScore(feature_it, "id_target_ind_log_sn_score");
        auto id_target_ind_massdev_score = getSeparateScore(feature_it, "id_target_ind_massdev_score");
        auto id_target_ind_mi_score = getSeparateScore(feature_it, "id_target_ind_mi_score");
        auto id_target_ind_mi_ratio_score = getSeparateScore(feature_it, "id_target_ind_mi_ratio_score");
        auto id_target_ind_isotope_correlation = getSeparateScore(feature_it, "id_target_ind_isotope_correlation");
        auto id_target_ind_isotope_overlap = getSeparateScore(feature_it, "id_target_ind_isotope_overlap");
        // Ion Mobility scores
        auto id_target_ind_im_drift = getSeparateScore(feature_it, "id_target_ind_im_drift");
        auto id_target_ind_im_drift_left = getSeparateScore(feature_it, "id_target_ind_im_drift_left");
        auto id_target_ind_im_drift_right = getSeparateScore(feature_it, "id_target_ind_im_drift_right");
        auto id_target_ind_im_delta = getSeparateScore(feature_it, "id_target_ind_im_delta");
        auto id_target_ind_im_delta_score = getSeparateScore(feature_it, "id_target_ind_im_delta_score");
        auto id_target_ind_im_log_intensity = getSeparateScore(feature_it, "id_target_ind_im_log_intensity");
        auto id_target_ind_im_contrast_coelution = getSeparateScore(feature_it, "id_target_ind_im_contrast_coelution");
        auto id_target_ind_im_contrast_shape = getSeparateScore(feature_it, "id_target_ind_im_contrast_shape");
        auto id_target_ind_im_sum_contrast_coelution = getSeparateScore(feature_it, "id_target_ind_im_sum_contrast_coelution");
        auto id_target_ind_im_sum_contrast_shape = getSeparateScore(feature_it, "id_target_ind_im_sum_contrast_shape");


        // check if there are compute_peak_shape_metrics scores
        auto id_target_ind_start_position_at_5 = getSeparateScore(feature_it, "id_target_ind_start_position_at_5");
        bool enable_compute_peak_shape_metrics = id_target_ind_start_position_at_5.size() > 0 && id_target_ind_start_position_at_5[0] != "0";
        
        // get scores for peak shape metrics will just be empty vector if not present
        auto start_position_at_5 = getSeparateScore(feature_it, "id_target_ind_start_position_at_5");
        auto end_position_at_5 = getSeparateScore(feature_it, "id_target_ind_end_position_at_5");
        auto start_position_at_10 = getSeparateScore(feature_it, "id_target_ind_start_position_at_10");
        auto end_position_at_10 = getSeparateScore(feature_it, "id_target_ind_end_position_at_10");
        auto start_position_at_50 = getSeparateScore(feature_it, "id_target_ind_start_position_at_50");
        auto end_position_at_50 = getSeparateScore(feature_it, "id_target_ind_end_position_at_50");
        auto total_width = getSeparateScore(feature_it, "id_target_ind_total_width");
        auto tailing_factor = getSeparateScore(feature_it, "id_target_ind_tailing_factor");
        auto asymmetry_factor = getSeparateScore(feature_it, "id_target_ind_asymmetry_factor");
        auto slope_of_baseline = getSeparateScore(feature_it, "id_target_ind_slope_of_baseline");
        auto baseline_delta_2_height = getSeparateScore(feature_it, "id_target_ind_baseline_delta_2_height");
        auto points_across_baseline = getSeparateScore(feature_it, "id_target_ind_points_across_baseline");
        auto points_across_half_height = getSeparateScore(feature_it, "id_target_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_target_num_transitions"))
        {
          int id_target_num_transitions = feature_it.getMetaValue("id_target_num_transitions");

          for (int i = 0; i < id_target_num_transitions; ++i)
          {
            std::vector<String> transition_row = makeTransitionRow(feature_id);
            transition_row[1] = oswValueAt(id_target_transition_names, i);
            transition_row[2] = oswValueAt(id_target_area_intensity, i);
            transition_row[3] = oswValueAt(id_target_total_area_intensity, i);
            transition_row[4] = oswValueAt(id_target_peak_apex_position, i);
            transition_row[5] = oswValueAt(id_target_apex_intensity, i);
            transition_row[6] = oswValueAt(id_target_peak_fwhm, i);
            transition_row[7] = oswValueAt(id_target_ind_massdev_score, i);
            transition_row[8] = oswValueAt(id_target_total_mi, i);
            transition_row[9] = oswValueAt(id_target_intensity_score, i);
            transition_row[10] = oswValueAt(id_target_intensity_ratio_score, i);
            transition_row[11] = oswValueAt(id_target_log_intensity, i);
            transition_row[12] = oswValueAt(id_target_ind_xcorr_coelution, i);
            transition_row[13] = oswValueAt(id_target_ind_xcorr_shape, i);
            transition_row[14] = oswValueAt(id_target_ind_log_sn_score, i);
            transition_row[15] = oswValueAt(id_target_ind_massdev_score, i);
            transition_row[16] = oswValueAt(id_target_ind_mi_score, i);
            transition_row[17] = oswValueAt(id_target_ind_mi_ratio_score, i);
            transition_row[18] = oswValueAt(id_target_ind_isotope_correlation, i);
            transition_row[19] = oswValueAt(id_target_ind_isotope_overlap, i);
            transition_row[20] = oswValueAt(id_target_ind_im_drift, i);
            transition_row[21] = oswValueAt(id_target_ind_im_drift_left, i);
            transition_row[22] = oswValueAt(id_target_ind_im_drift_right, i);
            transition_row[23] = oswValueAt(id_target_ind_im_delta, i);
            transition_row[24] = oswValueAt(id_target_ind_im_delta_score, i);
            transition_row[25] = oswValueAt(id_target_ind_im_log_intensity, i);
            transition_row[26] = oswValueAt(id_target_ind_im_contrast_coelution, i);
            transition_row[27] = oswValueAt(id_target_ind_im_contrast_shape, i);
            transition_row[28] = oswValueAt(id_target_ind_im_sum_contrast_coelution, i);
            transition_row[29] = oswValueAt(id_target_ind_im_sum_contrast_shape, i);
            if (enable_compute_peak_shape_metrics)
            {
              transition_row[30] = oswValueAt(start_position_at_5, i);
              transition_row[31] = oswValueAt(end_position_at_5, i);
              transition_row[32] = oswValueAt(start_position_at_10, i);
              transition_row[33] = oswValueAt(end_position_at_10, i);
              transition_row[34] = oswValueAt(start_position_at_50, i);
              transition_row[35] = oswValueAt(end_position_at_50, i);
              transition_row[36] = oswValueAt(total_width, i);
              transition_row[37] = oswValueAt(tailing_factor, i);
              transition_row[38] = oswValueAt(asymmetry_factor, i);
              transition_row[39] = oswValueAt(slope_of_baseline, i);
              transition_row[40] = oswValueAt(baseline_delta_2_height, i);
              transition_row[41] = oswValueAt(points_across_baseline, i);
              transition_row[42] = oswValueAt(points_across_half_height, i);
            }
            uis_transition_rows.push_back(std::move(transition_row));

          }
        }

        auto id_decoy_transition_names = getSeparateScore(feature_it, "id_decoy_transition_names");
        auto id_decoy_area_intensity = getSeparateScore(feature_it, "id_decoy_area_intensity");
        auto id_decoy_total_area_intensity = getSeparateScore(feature_it, "id_decoy_total_area_intensity");
        auto id_decoy_apex_intensity = getSeparateScore(feature_it, "id_decoy_apex_intensity");
        auto id_decoy_peak_apex_position = getSeparateScore(feature_it, "id_decoy_peak_apex_position");
        auto id_decoy_peak_fwhm = getSeparateScore(feature_it, "id_decoy_width_at_50");
        auto id_decoy_total_mi = getSeparateScore(feature_it, "id_decoy_total_mi");
        auto id_decoy_intensity_score = getSeparateScore(feature_it, "id_decoy_intensity_score");
        auto id_decoy_intensity_ratio_score = getSeparateScore(feature_it, "id_decoy_intensity_ratio_score");
        auto id_decoy_log_intensity = getSeparateScore(feature_it, "id_decoy_ind_log_intensity");
        auto id_decoy_ind_xcorr_coelution = getSeparateScore(feature_it, "id_decoy_ind_xcorr_coelution");
        auto id_decoy_ind_xcorr_shape = getSeparateScore(feature_it, "id_decoy_ind_xcorr_shape");
        auto id_decoy_ind_log_sn_score = getSeparateScore(feature_it, "id_decoy_ind_log_sn_score");
        auto id_decoy_ind_massdev_score = getSeparateScore(feature_it, "id_decoy_ind_massdev_score");
        auto id_decoy_ind_mi_score = getSeparateScore(feature_it, "id_decoy_ind_mi_score");
        auto id_decoy_ind_mi_ratio_score = getSeparateScore(feature_it, "id_decoy_ind_mi_ratio_score");
        auto id_decoy_ind_isotope_correlation = getSeparateScore(feature_it, "id_decoy_ind_isotope_correlation");
        auto id_decoy_ind_isotope_overlap = getSeparateScore(feature_it, "id_decoy_ind_isotope_overlap");

        // Ion Mobility scores
        auto id_decoy_ind_im_drift = getSeparateScore(feature_it, "id_decoy_ind_im_drift");
        auto id_decoy_ind_im_drift_left = getSeparateScore(feature_it, "id_decoy_ind_im_drift_left");
        auto id_decoy_ind_im_drift_right = getSeparateScore(feature_it, "id_decoy_ind_im_drift_right");
        auto id_decoy_ind_im_delta = getSeparateScore(feature_it, "id_decoy_ind_im_delta");
        auto id_decoy_ind_ind_im_delta_score = getSeparateScore(feature_it, "id_decoy_ind_im_delta_score");
        auto id_decoy_ind_log_intensity = getSeparateScore(feature_it, "id_decoy_ind_im_log_intensity");
        auto id_decoy_ind_im_contrast_coelution = getSeparateScore(feature_it, "id_decoy_ind_im_contrast_coelution");
        auto id_decoy_ind_im_contrast_shape = getSeparateScore(feature_it, "id_decoy_ind_im_contrast_shape");
        auto id_decoy_ind_im_sum_contrast_coelution = getSeparateScore(feature_it, "id_decoy_ind_im_sum_contrast_coelution");
        auto id_decoy_ind_im_sum_contrast_shape = getSeparateScore(feature_it, "id_decoy_ind_im_sum_contrast_shape");

        // get scores for peak shape metrics will just be empty vector if not present
        auto decoy_start_position_at_5 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_5");
        auto decoy_end_position_at_5 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_5");
        auto decoy_start_position_at_10 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_10");
        auto decoy_end_position_at_10 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_10");
        auto decoy_start_position_at_50 = getSeparateScore(feature_it, "id_decoy_ind_start_position_at_50");
        auto decoy_end_position_at_50 = getSeparateScore(feature_it, "id_decoy_ind_end_position_at_50");
        auto decoy_total_width = getSeparateScore(feature_it, "id_decoy_ind_total_width");
        auto decoy_tailing_factor = getSeparateScore(feature_it, "id_decoy_ind_tailing_factor");
        auto decoy_asymmetry_factor = getSeparateScore(feature_it, "id_decoy_ind_asymmetry_factor");
        auto decoy_slope_of_baseline = getSeparateScore(feature_it, "id_decoy_ind_slope_of_baseline");
        auto decoy_baseline_delta_2_height = getSeparateScore(feature_it, "id_decoy_ind_baseline_delta_2_height");
        auto decoy_points_across_baseline = getSeparateScore(feature_it, "id_decoy_ind_points_across_baseline");
        auto decoy_points_across_half_height = getSeparateScore(feature_it, "id_decoy_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_decoy_num_transitions"))
        {
          int id_decoy_num_transitions = feature_it.getMetaValue("id_decoy_num_transitions");

          for (int i = 0; i < id_decoy_num_transitions; ++i)
          {
            std::vector<String> transition_row = makeTransitionRow(feature_id);
            transition_row[1] = oswValueAt(id_decoy_transition_names, i);
            transition_row[2] = oswValueAt(id_decoy_area_intensity, i);
            transition_row[3] = oswValueAt(id_decoy_total_area_intensity, i);
            transition_row[4] = oswValueAt(id_decoy_peak_apex_position, i);
            transition_row[5] = oswValueAt(id_decoy_apex_intensity, i);
            transition_row[6] = oswValueAt(id_decoy_peak_fwhm, i);
            transition_row[7] = oswValueAt(id_decoy_ind_massdev_score, i);
            transition_row[8] = oswValueAt(id_decoy_total_mi, i);
            transition_row[9] = oswValueAt(id_decoy_intensity_score, i);
            transition_row[10] = oswValueAt(id_decoy_intensity_ratio_score, i);
            transition_row[11] = oswValueAt(id_decoy_log_intensity, i);
            transition_row[12] = oswValueAt(id_decoy_ind_xcorr_coelution, i);
            transition_row[13] = oswValueAt(id_decoy_ind_xcorr_shape, i);
            transition_row[14] = oswValueAt(id_decoy_ind_log_sn_score, i);
            transition_row[15] = oswValueAt(id_decoy_ind_massdev_score, i);
            transition_row[16] = oswValueAt(id_decoy_ind_mi_score, i);
            transition_row[17] = oswValueAt(id_decoy_ind_mi_ratio_score, i);
            transition_row[18] = oswValueAt(id_decoy_ind_isotope_correlation, i);
            transition_row[19] = oswValueAt(id_decoy_ind_isotope_overlap, i);
            transition_row[20] = oswValueAt(id_decoy_ind_im_drift, i);
            transition_row[21] = oswValueAt(id_decoy_ind_im_drift_left, i);
            transition_row[22] = oswValueAt(id_decoy_ind_im_drift_right, i);
            transition_row[23] = oswValueAt(id_decoy_ind_im_delta, i);
            transition_row[24] = oswValueAt(id_decoy_ind_ind_im_delta_score, i);
            transition_row[25] = oswValueAt(id_decoy_ind_log_intensity, i);
            transition_row[26] = oswValueAt(id_decoy_ind_im_contrast_coelution, i);
            transition_row[27] = oswValueAt(id_decoy_ind_im_contrast_shape, i);
            transition_row[28] = oswValueAt(id_decoy_ind_im_sum_contrast_coelution, i);
            transition_row[29] = oswValueAt(id_decoy_ind_im_sum_contrast_shape, i);

            if (enable_compute_peak_shape_metrics)
            {
              transition_row[30] = oswValueAt(decoy_start_position_at_5, i);
              transition_row[31] = oswValueAt(decoy_end_position_at_5, i);
              transition_row[32] = oswValueAt(decoy_start_position_at_10, i);
              transition_row[33] = oswValueAt(decoy_end_position_at_10, i);
              transition_row[34] = oswValueAt(decoy_start_position_at_50, i);
              transition_row[35] = oswValueAt(decoy_end_position_at_50, i);
              transition_row[36] = oswValueAt(decoy_total_width, i);
              transition_row[37] = oswValueAt(decoy_tailing_factor, i);
              transition_row[38] = oswValueAt(decoy_asymmetry_factor, i);
              transition_row[39] = oswValueAt(decoy_slope_of_baseline, i);
              transition_row[40] = oswValueAt(decoy_baseline_delta_2_height, i);
              transition_row[41] = oswValueAt(decoy_points_across_baseline, i);
              transition_row[42] = oswValueAt(decoy_points_across_half_height, i);
            }
            uis_transition_rows.push_back(std::move(transition_row));
          }
        }
      }
    }

    rows.feature_transition_rows = (enable_uis_scoring_ && !uis_transition_rows.empty()) ?
      std::move(uis_transition_rows) : std::move(ms2_transition_rows);

    return rows;
  }

  String OpenSwathOSWWriter::prepareLine(const OpenSwath::LightCompound& pep,
                                         const OpenSwath::LightTransition* transition,
                                         const FeatureMap& output,
                                         const String& id) const
  {
    const OSWData rows = prepareRows(pep, transition, output, id);
    std::stringstream sql;
    appendInsertSQL(sql, "FEATURE", featureColumns(), rows.feature_rows);
    appendInsertSQL(sql, "FEATURE_MS1", featureMS1Columns(), rows.feature_ms1_rows);
    appendInsertSQL(sql, "FEATURE_PRECURSOR", featurePrecursorColumns(), rows.feature_precursor_rows);
    appendInsertSQL(sql, "FEATURE_MS2", featureMS2Columns(), rows.feature_ms2_rows);
    appendInsertSQL(sql, "FEATURE_TRANSITION", featureTransitionColumns(), rows.feature_transition_rows);
    return String(sql.str());
  }

  void OpenSwathOSWWriter::writeRows(const OSWData& osw_output)
  {
    if (osw_output.empty())
    {
      return;
    }

    SqliteConnector conn(output_filename_);
    conn.executeStatement("BEGIN TRANSACTION");
    try
    {
      writeTableRows(conn.getDB(), "FEATURE", featureColumns(), osw_output.feature_rows);
      writeTableRows(conn.getDB(), "FEATURE_MS1", featureMS1Columns(), osw_output.feature_ms1_rows);
      writeTableRows(conn.getDB(), "FEATURE_PRECURSOR", featurePrecursorColumns(), osw_output.feature_precursor_rows);
      writeTableRows(conn.getDB(), "FEATURE_MS2", featureMS2Columns(), osw_output.feature_ms2_rows);
      writeTableRows(conn.getDB(), "FEATURE_TRANSITION", featureTransitionColumns(), osw_output.feature_transition_rows);
    }
    catch (...)
    {
      try
      {
        conn.executeStatement("ROLLBACK TRANSACTION");
      }
      catch (...)
      {
      }
      throw;
    }
    conn.executeStatement("END TRANSACTION");
  }

  void OpenSwathOSWWriter::writeLines(const std::vector<String>& to_osw_output)
  {
    SqliteConnector conn(output_filename_);
    conn.executeStatement("BEGIN TRANSACTION");
    for (Size i = 0; i < to_osw_output.size(); i++)
    {
      conn.executeStatement(to_osw_output[i]);
    }
    conn.executeStatement("END TRANSACTION");
  }
}
