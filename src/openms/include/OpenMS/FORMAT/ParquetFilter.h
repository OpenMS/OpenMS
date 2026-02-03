// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
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
  struct OPENMS_DLLAPI Condition
  {
    String column;
    String op;
    std::vector<String> values;
    ColumnType type{ColumnType::INT};
  };

  /**
    @brief Simple conjunction/disjunction of conditions.
  */
  struct OPENMS_DLLAPI FilterExpression
  {
    std::vector<Condition> conditions;
    std::vector<String> connectors; // "AND" / "OR"

    bool empty() const;
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
  class OPENMS_DLLAPI ParquetFilter
  {
  public:
    ParquetFilter() = default;
    ParquetFilter(const ParquetFilter&) = default;

    /**
      @brief Combine the next condition with logical AND.
    */
    ParquetFilter& andNext();

    /**
      @brief Combine the next condition with logical OR.
    */
    ParquetFilter& orNext();

    /**
      @brief Add an equality condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& eq(const String& column, Int64 value);

    /**
      @brief Add an inequality condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& ne(const String& column, Int64 value);

    /**
      @brief Add a less-than condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& lt(const String& column, Int64 value);

    /**
      @brief Add a less-than-or-equal condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& le(const String& column, Int64 value);

    /**
      @brief Add a greater-than condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& gt(const String& column, Int64 value);

    /**
      @brief Add a greater-than-or-equal condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilter& ge(const String& column, Int64 value);

    /**
      @brief Add an IN condition for an integer column.
      @param column Column name.
      @param values Integer values to compare.
    */
    ParquetFilter& in(const String& column, const std::vector<Int64>& values);

    /**
      @brief Add an equality condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& eq(const String& column, const String& value);

    /**
      @brief Add an inequality condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& ne(const String& column, const String& value);

    /**
      @brief Add a less-than condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& lt(const String& column, const String& value);

    /**
      @brief Add a less-than-or-equal condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& le(const String& column, const String& value);

    /**
      @brief Add a greater-than condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& gt(const String& column, const String& value);

    /**
      @brief Add a greater-than-or-equal condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilter& ge(const String& column, const String& value);

    /**
      @brief Add an IN condition for a string column.
      @param column Column name.
      @param values String values to compare.
    */
    ParquetFilter& in(const String& column, const std::vector<String>& values);

    /**
      @brief Return the filter expression.
      @return The constructed filter expression.
    */
    const FilterExpression& expression() const;

    /**
      @brief Check if the filter is empty (has no conditions).
      @return True if the filter has no conditions.
    */
    bool empty() const;

    /**
      @brief Merge two filter expressions with a connector.
      @param[in] lhs Left-hand side filter expression.
      @param[in] rhs Right-hand side filter expression.
      @param[in] connector Logical connector ("AND" or "OR").
      @return Combined filter expression.
    */
    static FilterExpression merge(const FilterExpression& lhs,
                                  const FilterExpression& rhs,
                                  const String& connector = "AND");

  private:
    ParquetFilter& addCondition_(const String& column,
                                 const String& op,
                                 ColumnType type,
                                 const std::vector<String>& values);

    FilterExpression expr_;
    String next_connector_{"AND"};
  };

  /**
    @brief Fluent builder for ParquetFilter objects.

    This builder offers a stable, chainable interface for constructing typed
    filters that can be passed to Parquet-backed readers. It wraps a
    ParquetFilter instance and forwards all condition-building operations.

    @section ParquetFilterBuilder_Example Example
    @code
    using namespace OpenMS;

    ParquetFilterBuilder f;
    f.eq("PRECURSOR_ID", 1318)
     .andNext()
     .eq("ANNOTATION", "y3^1");

    // Pass f.filter() into Parquet-backed readers.
    @endcode
  */
  class OPENMS_DLLAPI ParquetFilterBuilder
  {
  public:
    ParquetFilterBuilder() = default;
    ParquetFilterBuilder(const ParquetFilterBuilder&) = default;

    /**
      @brief Combine the next condition with logical AND.
    */
    ParquetFilterBuilder& andNext();

    /**
      @brief Combine the next condition with logical OR.
    */
    ParquetFilterBuilder& orNext();

    /**
      @brief Add an equality condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& eq(const String& column, Int64 value);

    /**
      @brief Add an inequality condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& ne(const String& column, Int64 value);

    /**
      @brief Add a less-than condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& lt(const String& column, Int64 value);

    /**
      @brief Add a less-than-or-equal condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& le(const String& column, Int64 value);

    /**
      @brief Add a greater-than condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& gt(const String& column, Int64 value);

    /**
      @brief Add a greater-than-or-equal condition for an integer column.
      @param column Column name.
      @param value Integer value to compare.
    */
    ParquetFilterBuilder& ge(const String& column, Int64 value);

    /**
      @brief Add an IN condition for an integer column.
      @param column Column name.
      @param values Integer values to compare.
    */
    ParquetFilterBuilder& in(const String& column, const std::vector<Int64>& values);

    /**
      @brief Add an equality condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& eq(const String& column, const String& value);

    /**
      @brief Add an inequality condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& ne(const String& column, const String& value);

    /**
      @brief Add a less-than condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& lt(const String& column, const String& value);

    /**
      @brief Add a less-than-or-equal condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& le(const String& column, const String& value);

    /**
      @brief Add a greater-than condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& gt(const String& column, const String& value);

    /**
      @brief Add a greater-than-or-equal condition for a string column.
      @param column Column name.
      @param value String value to compare.
    */
    ParquetFilterBuilder& ge(const String& column, const String& value);

    /**
      @brief Add an IN condition for a string column.
      @param column Column name.
      @param values String values to compare.
    */
    ParquetFilterBuilder& in(const String& column, const std::vector<String>& values);

    /**
      @brief Return the built filter.
      @return ParquetFilter instance.
    */
    const ParquetFilter& filter() const;

    /**
      @brief Return whether the filter is empty.
    */
    bool empty() const;

  private:
    ParquetFilter filter_;
  };
} // namespace OpenMS
