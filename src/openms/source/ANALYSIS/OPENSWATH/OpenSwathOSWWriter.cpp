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
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <new>
#include <sstream>
#include <utility>

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

    bool hasOnlyTrailingWhitespace(const char* position)
    {
      while (*position != '\0')
      {
        if (std::isspace(static_cast<unsigned char>(*position)) == 0)
        {
          return false;
        }
        ++position;
      }
      return true;
    }

    bool parseInt64(const String& value, Int64& result)
    {
      errno = 0;
      char* end = nullptr;
      const long long parsed = std::strtoll(value.c_str(), &end, 10);
      if (end == value.c_str() || errno == ERANGE || !hasOnlyTrailingWhitespace(end))
      {
        return false;
      }
      result = static_cast<Int64>(parsed);
      return true;
    }

    bool parseDouble(const String& value, double& result)
    {
      errno = 0;
      char* end = nullptr;
      const double parsed = std::strtod(value.c_str(), &end);
      if (end == value.c_str() || errno == ERANGE || !hasOnlyTrailingWhitespace(end))
      {
        return false;
      }
      result = parsed;
      return true;
    }

    OpenSwathOSWWriter::OSWValue makeOSWValueFromString(const String& value)
    {
      if (isNullLiteral(value))
      {
        return OpenSwathOSWWriter::OSWValue::null();
      }

      Int64 int_value = 0;
      if (parseInt64(value, int_value))
      {
        return OpenSwathOSWWriter::OSWValue(int_value);
      }

      double double_value = 0.0;
      if (parseDouble(value, double_value))
      {
        if (std::isnan(double_value))
        {
          return OpenSwathOSWWriter::OSWValue::null();
        }
        return OpenSwathOSWWriter::OSWValue(double_value);
      }

      return OpenSwathOSWWriter::OSWValue::text(value);
    }

    String oswValueToString(const OpenSwathOSWWriter::OSWValue& value)
    {
      switch (value.type)
      {
        case OpenSwathOSWWriter::OSWValue::Type::Null:
          return "NULL";
        case OpenSwathOSWWriter::OSWValue::Type::Int64:
          return String(value.asInt());
        case OpenSwathOSWWriter::OSWValue::Type::Double:
          return std::isnan(value.asDouble()) ? String("NULL") : String(value.asDouble());
        case OpenSwathOSWWriter::OSWValue::Type::Text:
          return value.asText();
      }
      return "NULL";
    }

    String oswValueToSQLLiteral(const OpenSwathOSWWriter::OSWValue& value)
    {
      if (value.type != OpenSwathOSWWriter::OSWValue::Type::Text)
      {
        return oswValueToString(value);
      }

      String quoted = value.asText();
      quoted.substitute("'", "''");
      return String("'") + quoted + "'";
    }

    OpenSwathOSWWriter::OSWValue oswValue(const String& value)
    {
      return makeOSWValueFromString(value);
    }

    OpenSwathOSWWriter::OSWValue oswValue(const DataValue& value)
    {
      if (value.isEmpty())
      {
        return OpenSwathOSWWriter::OSWValue::null();
      }

      if (value.valueType() == DataValue::INT_VALUE)
      {
        return OpenSwathOSWWriter::OSWValue(static_cast<Int64>(static_cast<long long>(value)));
      }
      if (value.valueType() == DataValue::DOUBLE_VALUE)
      {
        return OpenSwathOSWWriter::OSWValue(static_cast<double>(value));
      }

      return makeOSWValueFromString(value.toString());
    }

    OpenSwathOSWWriter::OSWValue oswValue(const double value)
    {
      return OpenSwathOSWWriter::OSWValue(value);
    }

    OpenSwathOSWWriter::OSWValue oswValue(const UInt64 value)
    {
      return OpenSwathOSWWriter::OSWValue(value);
    }

    class ScoreValueList
    {
    public:
      ScoreValueList() = default;

      ScoreValueList(const Feature& feature, const std::string& score_name)
      {
        const DataValue& dv = feature.getMetaValue(score_name);
        if (dv.isEmpty())
        {
          return;
        }

        type_ = dv.valueType();
        if (type_ == DataValue::STRING_LIST)
        {
          string_values_ = dv.toStringList();
        }
        else if (type_ == DataValue::INT_LIST)
        {
          int_values_ = dv.toIntList();
        }
        else if (type_ == DataValue::DOUBLE_LIST)
        {
          double_values_ = dv.toDoubleList();
        }
        else
        {
          single_value_ = oswValue(dv);
        }
      }

      Size size() const
      {
        if (type_ == DataValue::STRING_LIST)
        {
          return string_values_.size();
        }
        if (type_ == DataValue::INT_LIST)
        {
          return int_values_.size();
        }
        if (type_ == DataValue::DOUBLE_LIST)
        {
          return double_values_.size();
        }
        return type_ == DataValue::EMPTY_VALUE ? 0 : 1;
      }

      bool empty() const
      {
        return size() == 0;
      }

      OpenSwathOSWWriter::OSWValue at(Size index) const
      {
        if (index >= size())
        {
          return OpenSwathOSWWriter::OSWValue::null();
        }

        if (type_ == DataValue::STRING_LIST)
        {
          return oswValue(string_values_[index]);
        }
        if (type_ == DataValue::INT_LIST)
        {
          return OpenSwathOSWWriter::OSWValue(static_cast<Int64>(int_values_[index]));
        }
        if (type_ == DataValue::DOUBLE_LIST)
        {
          return oswValue(double_values_[index]);
        }
        return single_value_;
      }

    private:
      DataValue::DataType type_ = DataValue::EMPTY_VALUE;
      StringList string_values_;
      IntList int_values_;
      DoubleList double_values_;
      OpenSwathOSWWriter::OSWValue single_value_;
    };

    OpenSwathOSWWriter::OSWValue oswValueAt(const ScoreValueList& values, Size index)
    {
      return values.at(index);
    }

    OpenSwathOSWWriter::FeatureTransitionRow makeTransitionRow(const OpenSwathOSWWriter::OSWValue& feature_id)
    {
      OpenSwathOSWWriter::FeatureTransitionRow row{};
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

    template <typename RowT>
    void bindRow(sqlite3* db, sqlite3_stmt* stmt, const RowT& row, Size expected_columns)
    {
      if (row.size() != expected_columns)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "OSW row has " + String(row.size()) + " values, expected " + String(expected_columns));
      }

      for (Size i = 0; i < row.size(); ++i)
      {
        int rc = SQLITE_OK;
        if (row[i].isNull())
        {
          rc = sqlite3_bind_null(stmt, static_cast<int>(i + 1));
        }
        else if (row[i].type == OpenSwathOSWWriter::OSWValue::Type::Int64)
        {
          rc = sqlite3_bind_int64(stmt, static_cast<int>(i + 1), static_cast<sqlite3_int64>(row[i].asInt()));
        }
        else if (row[i].type == OpenSwathOSWWriter::OSWValue::Type::Double)
        {
          rc = sqlite3_bind_double(stmt, static_cast<int>(i + 1), row[i].asDouble());
        }
        else if (row[i].type == OpenSwathOSWWriter::OSWValue::Type::Text)
        {
          rc = sqlite3_bind_text(stmt, static_cast<int>(i + 1), row[i].asText().c_str(),
                                 static_cast<int>(row[i].asText().size()), SQLITE_TRANSIENT);
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

    template <typename RowT>
    void writeTableRows(sqlite3* db, const char* table, const std::vector<const char*>& columns,
                        const std::vector<RowT>& rows)
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

    template <typename RowT>
    void appendRows(std::vector<RowT>& target, std::vector<RowT>& source)
    {
      if (target.empty())
      {
        target = std::move(source);
      }
      else
      {
        target.insert(target.end(), std::make_move_iterator(source.begin()), std::make_move_iterator(source.end()));
      }
      source.clear();
    }

    template <typename RowT>
    UInt64 estimateRowVectorMemory(const std::vector<RowT>& rows)
    {
      constexpr UInt64 string_headroom_per_value = 16;
      return static_cast<UInt64>(rows.capacity()) *
             (static_cast<UInt64>(sizeof(RowT)) +
              static_cast<UInt64>(std::tuple_size<RowT>::value) * string_headroom_per_value);
    }

    template <typename RowT>
    void appendInsertSQL(std::stringstream& sql, const char* table, const std::vector<const char*>& columns,
                         const std::vector<RowT>& rows)
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
          sql << oswValueToSQLLiteral(row[i]);
        }
        sql << "); ";
      }
    }
  }

  OpenSwathOSWWriter::OSWValue::Storage::Storage() :
    int_value(0)
  {
  }

  OpenSwathOSWWriter::OSWValue::Storage::~Storage()
  {
  }

  OpenSwathOSWWriter::OSWValue::OSWValue() :
    type(Type::Null),
    storage()
  {
  }

  OpenSwathOSWWriter::OSWValue::~OSWValue()
  {
    destroyText_();
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(std::nullptr_t) :
    OSWValue()
  {
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(const char* value) :
    OSWValue(value == nullptr ? String() : String(value))
  {
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(const String& value) :
    OSWValue(makeOSWValueFromString(value))
  {
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(const std::string& value) :
    OSWValue(makeOSWValueFromString(String(value)))
  {
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(Int64 value) :
    type(Type::Int64),
    storage()
  {
    storage.int_value = value;
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(UInt64 value) :
    type(Type::Int64),
    storage()
  {
    storage.int_value = static_cast<Int64>(value);
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(int value) :
    type(Type::Int64),
    storage()
  {
    storage.int_value = static_cast<Int64>(value);
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(double value) :
    type(std::isnan(value) ? Type::Null : Type::Double),
    storage()
  {
    if (type == Type::Double)
    {
      storage.double_value = value;
    }
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(const OSWValue& rhs) :
    OSWValue()
  {
    copyFrom_(rhs);
  }

  OpenSwathOSWWriter::OSWValue::OSWValue(OSWValue&& rhs) noexcept :
    OSWValue()
  {
    moveFrom_(std::move(rhs));
  }

  OpenSwathOSWWriter::OSWValue& OpenSwathOSWWriter::OSWValue::operator=(const OSWValue& rhs)
  {
    if (this != &rhs)
    {
      destroyText_();
      copyFrom_(rhs);
    }
    return *this;
  }

  OpenSwathOSWWriter::OSWValue& OpenSwathOSWWriter::OSWValue::operator=(OSWValue&& rhs) noexcept
  {
    if (this != &rhs)
    {
      destroyText_();
      moveFrom_(std::move(rhs));
    }
    return *this;
  }

  OpenSwathOSWWriter::OSWValue OpenSwathOSWWriter::OSWValue::null()
  {
    return OSWValue();
  }

  OpenSwathOSWWriter::OSWValue OpenSwathOSWWriter::OSWValue::text(const String& value)
  {
    OSWValue result;
    new (&result.storage.text_value) String(value);
    result.type = Type::Text;
    return result;
  }

  bool OpenSwathOSWWriter::OSWValue::isNull() const
  {
    return type == Type::Null || (type == Type::Double && std::isnan(storage.double_value));
  }

  Int64 OpenSwathOSWWriter::OSWValue::asInt() const
  {
    return storage.int_value;
  }

  double OpenSwathOSWWriter::OSWValue::asDouble() const
  {
    return storage.double_value;
  }

  const String& OpenSwathOSWWriter::OSWValue::asText() const
  {
    return storage.text_value;
  }

  void OpenSwathOSWWriter::OSWValue::destroyText_()
  {
    if (type == Type::Text)
    {
      storage.text_value.~String();
    }
    type = Type::Null;
    storage.int_value = 0;
  }

  void OpenSwathOSWWriter::OSWValue::copyFrom_(const OSWValue& rhs)
  {
    if (rhs.type == Type::Int64)
    {
      storage.int_value = rhs.storage.int_value;
      type = Type::Int64;
    }
    else if (rhs.type == Type::Double)
    {
      storage.double_value = rhs.storage.double_value;
      type = Type::Double;
    }
    else if (rhs.type == Type::Text)
    {
      new (&storage.text_value) String(rhs.storage.text_value);
      type = Type::Text;
    }
    else
    {
      type = Type::Null;
      storage.int_value = 0;
    }
  }

  void OpenSwathOSWWriter::OSWValue::moveFrom_(OSWValue&& rhs)
  {
    if (rhs.type == Type::Int64)
    {
      storage.int_value = rhs.storage.int_value;
      type = Type::Int64;
      rhs.type = Type::Null;
      rhs.storage.int_value = 0;
    }
    else if (rhs.type == Type::Double)
    {
      storage.double_value = rhs.storage.double_value;
      type = Type::Double;
      rhs.type = Type::Null;
      rhs.storage.int_value = 0;
    }
    else if (rhs.type == Type::Text)
    {
      new (&storage.text_value) String(std::move(rhs.storage.text_value));
      type = Type::Text;
      rhs.storage.text_value.~String();
      rhs.type = Type::Null;
      rhs.storage.int_value = 0;
    }
    else
    {
      type = Type::Null;
      storage.int_value = 0;
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

  Size OpenSwathOSWWriter::OSWData::rowCount() const
  {
    return feature_rows.size() + feature_ms1_rows.size() + feature_ms2_rows.size() +
           feature_precursor_rows.size() + feature_transition_rows.size();
  }

  UInt64 OpenSwathOSWWriter::OSWData::estimateMemoryUsage() const
  {
    return estimateRowVectorMemory(feature_rows) +
           estimateRowVectorMemory(feature_ms1_rows) +
           estimateRowVectorMemory(feature_ms2_rows) +
           estimateRowVectorMemory(feature_precursor_rows) +
           estimateRowVectorMemory(feature_transition_rows);
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
    return oswValueToString(oswValue(feature.getMetaValue(score_name)));
  }

  std::vector<String> OpenSwathOSWWriter::getSeparateScore(const Feature& feature, const std::string& score_name) const
  {
    std::vector<String> separated_scores;
    const ScoreValueList score_values(feature, score_name);
    separated_scores.reserve(score_values.size());
    for (Size i = 0; i < score_values.size(); ++i)
    {
      separated_scores.push_back(oswValueToString(score_values.at(i)));
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
    std::vector<FeatureTransitionRow> ms2_transition_rows;
    std::vector<FeatureTransitionRow> uis_transition_rows;
    ms2_transition_rows.reserve(output.size());
    uis_transition_rows.reserve(output.size());

    for (const auto& feature_it : output)
    {
      const OSWValue feature_id = oswValue(Internal::SqliteHelper::clearSignBit(feature_it.getUniqueId()));

      const ScoreValueList masserror_ppm(feature_it, "masserror_ppm");

      const auto& subordinates = feature_it.getSubordinates();
      for (Size i=0; i < subordinates.size(); i++)
      {
        const auto& sub_it = subordinates[i];
        if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2")
        {
          OSWValue total_mi = oswValue(sub_it.getMetaValue("total_mi")); // total_mi is not guaranteed to be set
          OSWValue masserror_ppm_query = oswValueAt(masserror_ppm, i); // masserror_ppm is not guaranteed to be set

          FeatureTransitionRow transition_row = makeTransitionRow(feature_id);
          transition_row[1] = oswValue(sub_it.getMetaValue("native_id"));
          transition_row[2] = oswValue(sub_it.getIntensity());
          transition_row[3] = oswValue(sub_it.getMetaValue("total_xic"));
          transition_row[4] = oswValue(sub_it.getMetaValue("peak_apex_position"));
          transition_row[5] = oswValue(sub_it.getMetaValue("peak_apex_int"));
          transition_row[6] = oswValue(sub_it.getMetaValue("width_at_50"));
          transition_row[7] = masserror_ppm_query;
          transition_row[8] = total_mi;

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
          oswValueToString(oswValue(sub_it.getMetaValue("native_id"))).split(OpenMS::String("Precursor_i"), precursor_id);
          rows.feature_precursor_rows.push_back({
            feature_id,
            precursor_id.size() > 1 ? oswValue(precursor_id[1]) : OSWValue::null(),
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
        oswValue(feature_it.getMetaValue("im_drift")),
        oswValue(norm_rt),
        oswValue(delta_rt),
        oswValue(feature_it.getMetaValue("leftWidth")),
        oswValue(feature_it.getMetaValue("rightWidth")),
        oswValue(feature_it.getMetaValue("im_drift_left")),
        oswValue(feature_it.getMetaValue("im_drift_right"))
      });

      rows.feature_ms2_rows.push_back({
        feature_id,
        oswValue(feature_it.getIntensity()),
        oswValue(feature_it.getMetaValue("total_xic")),
        oswValue(feature_it.getMetaValue("peak_apices_sum")),
        oswValue(feature_it.getMetaValue("im_drift")),
        oswValue(feature_it.getMetaValue("im_drift_left")),
        oswValue(feature_it.getMetaValue("im_drift_right")),
        oswValue(feature_it.getMetaValue("im_delta")),
        oswValue(feature_it.getMetaValue("total_mi")),
        oswValue(feature_it.getMetaValue("var_bseries_score")),
        oswValue(feature_it.getMetaValue("var_dotprod_score")),
        oswValue(feature_it.getMetaValue("var_intensity_score")),
        oswValue(feature_it.getMetaValue("var_isotope_correlation_score")),
        oswValue(feature_it.getMetaValue("var_isotope_overlap_score")),
        oswValue(feature_it.getMetaValue("var_library_corr")),
        oswValue(feature_it.getMetaValue("var_library_dotprod")),
        oswValue(feature_it.getMetaValue("var_library_manhattan")),
        oswValue(feature_it.getMetaValue("var_library_rmsd")),
        oswValue(feature_it.getMetaValue("var_library_rootmeansquare")),
        oswValue(feature_it.getMetaValue("var_library_sangle")),
        oswValue(feature_it.getMetaValue("var_log_sn_score")),
        oswValue(feature_it.getMetaValue("var_manhatt_score")),
        oswValue(feature_it.getMetaValue("var_massdev_score")),
        oswValue(feature_it.getMetaValue("var_massdev_score_weighted")),
        oswValue(feature_it.getMetaValue("var_mi_score")),
        oswValue(feature_it.getMetaValue("var_mi_weighted_score")),
        oswValue(feature_it.getMetaValue("var_mi_ratio_score")),
        oswValue(feature_it.getMetaValue("var_norm_rt_score")),
        oswValue(feature_it.getMetaValue("var_xcorr_coelution")),
        oswValue(feature_it.getMetaValue("var_xcorr_coelution_weighted")),
        oswValue(feature_it.getMetaValue("var_xcorr_shape")),
        oswValue(feature_it.getMetaValue("var_xcorr_shape_weighted")),
        oswValue(feature_it.getMetaValue("var_yseries_score")),
        oswValue(feature_it.getMetaValue("var_elution_model_fit_score")),
        oswValue(feature_it.getMetaValue("var_im_xcorr_shape")),
        oswValue(feature_it.getMetaValue("var_im_xcorr_coelution")),
        oswValue(feature_it.getMetaValue("var_im_delta_score")),
        oswValue(feature_it.getMetaValue("im_log_intensity"))
      });

      bool enable_ms1 = feature_it.metaValueExists("var_ms1_ppm_diff");
      if (enable_ms1) // only write MS1 scores if they are present
      {
        rows.feature_ms1_rows.push_back({
          feature_id,
          oswValue(feature_it.getMetaValue("ms1_area_intensity")),
          oswValue(feature_it.getMetaValue("ms1_apex_intensity")),
          oswValue(feature_it.getMetaValue("im_ms1_drift")),
          oswValue(feature_it.getMetaValue("im_ms1_delta")),
          oswValue(feature_it.getMetaValue("var_ms1_ppm_diff")),
          oswValue(feature_it.getMetaValue("var_im_ms1_delta_score")),
          oswValue(feature_it.getMetaValue("var_ms1_mi_score")),
          oswValue(feature_it.getMetaValue("var_ms1_mi_contrast_score")),
          oswValue(feature_it.getMetaValue("var_ms1_mi_combined_score")),
          oswValue(feature_it.getMetaValue("var_ms1_isotope_correlation")),
          oswValue(feature_it.getMetaValue("var_ms1_isotope_overlap")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_coelution")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_coelution_contrast")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_coelution_combined")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_shape")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_shape_contrast")),
          oswValue(feature_it.getMetaValue("var_ms1_xcorr_shape_combined"))
        });
      }

      if (enable_uis_scoring_)
      {
        ScoreValueList id_target_transition_names(feature_it, "id_target_transition_names");
        ScoreValueList id_target_area_intensity(feature_it, "id_target_area_intensity");
        ScoreValueList id_target_total_area_intensity(feature_it, "id_target_total_area_intensity");
        ScoreValueList id_target_apex_intensity(feature_it, "id_target_apex_intensity");
        ScoreValueList id_target_peak_apex_position(feature_it, "id_target_peak_apex_position");
        ScoreValueList id_target_peak_fwhm(feature_it, "id_target_width_at_50");
        ScoreValueList id_target_total_mi(feature_it, "id_target_total_mi");
        ScoreValueList id_target_intensity_score(feature_it, "id_target_intensity_score");
        ScoreValueList id_target_intensity_ratio_score(feature_it, "id_target_intensity_ratio_score");
        ScoreValueList id_target_log_intensity(feature_it, "id_target_ind_log_intensity");
        ScoreValueList id_target_ind_xcorr_coelution(feature_it, "id_target_ind_xcorr_coelution");
        ScoreValueList id_target_ind_xcorr_shape(feature_it, "id_target_ind_xcorr_shape");
        ScoreValueList id_target_ind_log_sn_score(feature_it, "id_target_ind_log_sn_score");
        ScoreValueList id_target_ind_massdev_score(feature_it, "id_target_ind_massdev_score");
        ScoreValueList id_target_ind_mi_score(feature_it, "id_target_ind_mi_score");
        ScoreValueList id_target_ind_mi_ratio_score(feature_it, "id_target_ind_mi_ratio_score");
        ScoreValueList id_target_ind_isotope_correlation(feature_it, "id_target_ind_isotope_correlation");
        ScoreValueList id_target_ind_isotope_overlap(feature_it, "id_target_ind_isotope_overlap");
        // Ion Mobility scores
        ScoreValueList id_target_ind_im_drift(feature_it, "id_target_ind_im_drift");
        ScoreValueList id_target_ind_im_drift_left(feature_it, "id_target_ind_im_drift_left");
        ScoreValueList id_target_ind_im_drift_right(feature_it, "id_target_ind_im_drift_right");
        ScoreValueList id_target_ind_im_delta(feature_it, "id_target_ind_im_delta");
        ScoreValueList id_target_ind_im_delta_score(feature_it, "id_target_ind_im_delta_score");
        ScoreValueList id_target_ind_im_log_intensity(feature_it, "id_target_ind_im_log_intensity");
        ScoreValueList id_target_ind_im_contrast_coelution(feature_it, "id_target_ind_im_contrast_coelution");
        ScoreValueList id_target_ind_im_contrast_shape(feature_it, "id_target_ind_im_contrast_shape");
        ScoreValueList id_target_ind_im_sum_contrast_coelution(feature_it, "id_target_ind_im_sum_contrast_coelution");
        ScoreValueList id_target_ind_im_sum_contrast_shape(feature_it, "id_target_ind_im_sum_contrast_shape");


        // check if there are compute_peak_shape_metrics scores
        ScoreValueList id_target_ind_start_position_at_5(feature_it, "id_target_ind_start_position_at_5");
        bool enable_compute_peak_shape_metrics = !id_target_ind_start_position_at_5.empty() &&
                                                 oswValueToString(id_target_ind_start_position_at_5.at(0)) != "0";

        // get scores for peak shape metrics will just be empty vector if not present
        const ScoreValueList& start_position_at_5 = id_target_ind_start_position_at_5;
        ScoreValueList end_position_at_5(feature_it, "id_target_ind_end_position_at_5");
        ScoreValueList start_position_at_10(feature_it, "id_target_ind_start_position_at_10");
        ScoreValueList end_position_at_10(feature_it, "id_target_ind_end_position_at_10");
        ScoreValueList start_position_at_50(feature_it, "id_target_ind_start_position_at_50");
        ScoreValueList end_position_at_50(feature_it, "id_target_ind_end_position_at_50");
        ScoreValueList total_width(feature_it, "id_target_ind_total_width");
        ScoreValueList tailing_factor(feature_it, "id_target_ind_tailing_factor");
        ScoreValueList asymmetry_factor(feature_it, "id_target_ind_asymmetry_factor");
        ScoreValueList slope_of_baseline(feature_it, "id_target_ind_slope_of_baseline");
        ScoreValueList baseline_delta_2_height(feature_it, "id_target_ind_baseline_delta_2_height");
        ScoreValueList points_across_baseline(feature_it, "id_target_ind_points_across_baseline");
        ScoreValueList points_across_half_height(feature_it, "id_target_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_target_num_transitions"))
        {
          int id_target_num_transitions = feature_it.getMetaValue("id_target_num_transitions");

          for (int i = 0; i < id_target_num_transitions; ++i)
          {
            FeatureTransitionRow transition_row = makeTransitionRow(feature_id);
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

        ScoreValueList id_decoy_transition_names(feature_it, "id_decoy_transition_names");
        ScoreValueList id_decoy_area_intensity(feature_it, "id_decoy_area_intensity");
        ScoreValueList id_decoy_total_area_intensity(feature_it, "id_decoy_total_area_intensity");
        ScoreValueList id_decoy_apex_intensity(feature_it, "id_decoy_apex_intensity");
        ScoreValueList id_decoy_peak_apex_position(feature_it, "id_decoy_peak_apex_position");
        ScoreValueList id_decoy_peak_fwhm(feature_it, "id_decoy_width_at_50");
        ScoreValueList id_decoy_total_mi(feature_it, "id_decoy_total_mi");
        ScoreValueList id_decoy_intensity_score(feature_it, "id_decoy_intensity_score");
        ScoreValueList id_decoy_intensity_ratio_score(feature_it, "id_decoy_intensity_ratio_score");
        ScoreValueList id_decoy_log_intensity(feature_it, "id_decoy_ind_log_intensity");
        ScoreValueList id_decoy_ind_xcorr_coelution(feature_it, "id_decoy_ind_xcorr_coelution");
        ScoreValueList id_decoy_ind_xcorr_shape(feature_it, "id_decoy_ind_xcorr_shape");
        ScoreValueList id_decoy_ind_log_sn_score(feature_it, "id_decoy_ind_log_sn_score");
        ScoreValueList id_decoy_ind_massdev_score(feature_it, "id_decoy_ind_massdev_score");
        ScoreValueList id_decoy_ind_mi_score(feature_it, "id_decoy_ind_mi_score");
        ScoreValueList id_decoy_ind_mi_ratio_score(feature_it, "id_decoy_ind_mi_ratio_score");
        ScoreValueList id_decoy_ind_isotope_correlation(feature_it, "id_decoy_ind_isotope_correlation");
        ScoreValueList id_decoy_ind_isotope_overlap(feature_it, "id_decoy_ind_isotope_overlap");

        // Ion Mobility scores
        ScoreValueList id_decoy_ind_im_drift(feature_it, "id_decoy_ind_im_drift");
        ScoreValueList id_decoy_ind_im_drift_left(feature_it, "id_decoy_ind_im_drift_left");
        ScoreValueList id_decoy_ind_im_drift_right(feature_it, "id_decoy_ind_im_drift_right");
        ScoreValueList id_decoy_ind_im_delta(feature_it, "id_decoy_ind_im_delta");
        ScoreValueList id_decoy_ind_ind_im_delta_score(feature_it, "id_decoy_ind_im_delta_score");
        ScoreValueList id_decoy_ind_log_intensity(feature_it, "id_decoy_ind_im_log_intensity");
        ScoreValueList id_decoy_ind_im_contrast_coelution(feature_it, "id_decoy_ind_im_contrast_coelution");
        ScoreValueList id_decoy_ind_im_contrast_shape(feature_it, "id_decoy_ind_im_contrast_shape");
        ScoreValueList id_decoy_ind_im_sum_contrast_coelution(feature_it, "id_decoy_ind_im_sum_contrast_coelution");
        ScoreValueList id_decoy_ind_im_sum_contrast_shape(feature_it, "id_decoy_ind_im_sum_contrast_shape");

        // get scores for peak shape metrics will just be empty vector if not present
        ScoreValueList decoy_start_position_at_5(feature_it, "id_decoy_ind_start_position_at_5");
        ScoreValueList decoy_end_position_at_5(feature_it, "id_decoy_ind_end_position_at_5");
        ScoreValueList decoy_start_position_at_10(feature_it, "id_decoy_ind_start_position_at_10");
        ScoreValueList decoy_end_position_at_10(feature_it, "id_decoy_ind_end_position_at_10");
        ScoreValueList decoy_start_position_at_50(feature_it, "id_decoy_ind_start_position_at_50");
        ScoreValueList decoy_end_position_at_50(feature_it, "id_decoy_ind_end_position_at_50");
        ScoreValueList decoy_total_width(feature_it, "id_decoy_ind_total_width");
        ScoreValueList decoy_tailing_factor(feature_it, "id_decoy_ind_tailing_factor");
        ScoreValueList decoy_asymmetry_factor(feature_it, "id_decoy_ind_asymmetry_factor");
        ScoreValueList decoy_slope_of_baseline(feature_it, "id_decoy_ind_slope_of_baseline");
        ScoreValueList decoy_baseline_delta_2_height(feature_it, "id_decoy_ind_baseline_delta_2_height");
        ScoreValueList decoy_points_across_baseline(feature_it, "id_decoy_ind_points_across_baseline");
        ScoreValueList decoy_points_across_half_height(feature_it, "id_decoy_ind_points_across_half_height");

        if (feature_it.metaValueExists("id_decoy_num_transitions"))
        {
          int id_decoy_num_transitions = feature_it.getMetaValue("id_decoy_num_transitions");

          for (int i = 0; i < id_decoy_num_transitions; ++i)
          {
            FeatureTransitionRow transition_row = makeTransitionRow(feature_id);
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
