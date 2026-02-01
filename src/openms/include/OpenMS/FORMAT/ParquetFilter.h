// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Column type for typed parquet filters.
  */
  enum class ColumnType
  {
    INT,
    STRING
  };

  /**
    @brief Single filter condition (column, operator, values).
  */
  struct Condition
  {
    String column;
    String op;
    std::vector<String> values;
    ColumnType type{ColumnType::INT};
  };

  /**
    @brief Simple conjunction/disjunction of conditions.
  */
  struct FilterExpression
  {
    std::vector<Condition> conditions;
    std::vector<String> connectors; // "AND" / "OR"

    bool empty() const
    {
      return conditions.empty();
    }
  };

  /**
    @brief Typed filter builder for parquet-backed datasets.

    This builder provides a structured alternative to string-encoded filters.
    It is intended for use in Parquet-backed readers (e.g., XICParquetFile) to
    construct filter expressions that can be pushed down to Arrow Dataset
    scanners or evaluated in-memory if pushdown is unavailable.

    Conditions are combined left-to-right. Use andNext()/orNext() to control
    how consecutive conditions are connected (default is AND).

    @section ParquetFilter_Example Example
    @code
    using namespace OpenMS;

    ParquetFilter f;
    f.eq("PRECURSOR_ID", 1318)
     .andNext()
     .eq("ANNOTATION", "y3^1");

    // Use f.expression() to pass into a parquet reader, or pass the builder
    // directly when supported.
    @endcode
  */
  class ParquetFilter
  {
  public:
    ParquetFilter& andNext()
    {
      next_connector_ = "AND";
      return *this;
    }

    ParquetFilter& orNext()
    {
      next_connector_ = "OR";
      return *this;
    }

    ParquetFilter& eq(const String& column, Int64 value)
    {
      return addCondition_(column, "=", ColumnType::INT, {String(value)});
    }

    ParquetFilter& ne(const String& column, Int64 value)
    {
      return addCondition_(column, "!=", ColumnType::INT, {String(value)});
    }

    ParquetFilter& lt(const String& column, Int64 value)
    {
      return addCondition_(column, "<", ColumnType::INT, {String(value)});
    }

    ParquetFilter& le(const String& column, Int64 value)
    {
      return addCondition_(column, "<=", ColumnType::INT, {String(value)});
    }

    ParquetFilter& gt(const String& column, Int64 value)
    {
      return addCondition_(column, ">", ColumnType::INT, {String(value)});
    }

    ParquetFilter& ge(const String& column, Int64 value)
    {
      return addCondition_(column, ">=", ColumnType::INT, {String(value)});
    }

    ParquetFilter& in(const String& column, const std::vector<Int64>& values)
    {
      std::vector<String> str_values;
      str_values.reserve(values.size());
      for (const auto& v : values)
      {
        str_values.emplace_back(String(v));
      }
      return addCondition_(column, "IN", ColumnType::INT, str_values);
    }

    ParquetFilter& eq(const String& column, const String& value)
    {
      return addCondition_(column, "=", ColumnType::STRING, {value});
    }

    ParquetFilter& ne(const String& column, const String& value)
    {
      return addCondition_(column, "!=", ColumnType::STRING, {value});
    }

    ParquetFilter& lt(const String& column, const String& value)
    {
      return addCondition_(column, "<", ColumnType::STRING, {value});
    }

    ParquetFilter& le(const String& column, const String& value)
    {
      return addCondition_(column, "<=", ColumnType::STRING, {value});
    }

    ParquetFilter& gt(const String& column, const String& value)
    {
      return addCondition_(column, ">", ColumnType::STRING, {value});
    }

    ParquetFilter& ge(const String& column, const String& value)
    {
      return addCondition_(column, ">=", ColumnType::STRING, {value});
    }

    ParquetFilter& in(const String& column, const std::vector<String>& values)
    {
      return addCondition_(column, "IN", ColumnType::STRING, values);
    }

    const FilterExpression& expression() const
    {
      return expr_;
    }

    bool empty() const
    {
      return expr_.empty();
    }

    static FilterExpression merge(const FilterExpression& lhs,
                                  const FilterExpression& rhs,
                                  const String& connector = "AND")
    {
      if (lhs.conditions.empty())
      {
        return rhs;
      }
      if (rhs.conditions.empty())
      {
        return lhs;
      }

      FilterExpression combined;
      combined.conditions.reserve(lhs.conditions.size() + rhs.conditions.size());
      combined.connectors.reserve(lhs.connectors.size() + rhs.connectors.size() + 1);

      combined.conditions.insert(combined.conditions.end(), lhs.conditions.begin(), lhs.conditions.end());
      combined.connectors.insert(combined.connectors.end(), lhs.connectors.begin(), lhs.connectors.end());
      combined.connectors.push_back(connector);
      combined.conditions.insert(combined.conditions.end(), rhs.conditions.begin(), rhs.conditions.end());
      combined.connectors.insert(combined.connectors.end(), rhs.connectors.begin(), rhs.connectors.end());

      return combined;
    }

  private:
    ParquetFilter& addCondition_(const String& column,
                                 const String& op,
                                 ColumnType type,
                                 const std::vector<String>& values)
    {
      Condition cond;
      cond.column = column;
      cond.op = op;
      cond.values = values;
      cond.type = type;

      if (!expr_.conditions.empty())
      {
        expr_.connectors.push_back(next_connector_);
      }
      expr_.conditions.push_back(std::move(cond));
      next_connector_ = "AND";
      return *this;
    }

    FilterExpression expr_;
    String next_connector_{"AND"};
  };
} // namespace OpenMS
