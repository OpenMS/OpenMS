// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParquetFilter.h>

namespace OpenMS
{
  // FilterExpression
  bool FilterExpression::empty() const
  {
    return conditions.empty();
  }

  // ParquetFilter
  ParquetFilter& ParquetFilter::andNext()
  {
    next_connector_ = "AND";
    return *this;
  }

  ParquetFilter& ParquetFilter::orNext()
  {
    next_connector_ = "OR";
    return *this;
  }

  ParquetFilter& ParquetFilter::eq(const String& column, Int64 value)
  {
    return addCondition_(column, "=", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::ne(const String& column, Int64 value)
  {
    return addCondition_(column, "!=", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::lt(const String& column, Int64 value)
  {
    return addCondition_(column, "<", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::le(const String& column, Int64 value)
  {
    return addCondition_(column, "<=", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::gt(const String& column, Int64 value)
  {
    return addCondition_(column, ">", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::ge(const String& column, Int64 value)
  {
    return addCondition_(column, ">=", ColumnType::INT, {String(value)});
  }

  ParquetFilter& ParquetFilter::in(const String& column, const std::vector<Int64>& values)
  {
    std::vector<String> str_values;
    str_values.reserve(values.size());
    for (const auto& v : values)
    {
      str_values.emplace_back(v);
    }
    return addCondition_(column, "IN", ColumnType::INT, str_values);
  }

  ParquetFilter& ParquetFilter::eq(const String& column, const String& value)
  {
    return addCondition_(column, "=", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::ne(const String& column, const String& value)
  {
    return addCondition_(column, "!=", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::lt(const String& column, const String& value)
  {
    return addCondition_(column, "<", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::le(const String& column, const String& value)
  {
    return addCondition_(column, "<=", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::gt(const String& column, const String& value)
  {
    return addCondition_(column, ">", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::ge(const String& column, const String& value)
  {
    return addCondition_(column, ">=", ColumnType::STRING, {value});
  }

  ParquetFilter& ParquetFilter::in(const String& column, const std::vector<String>& values)
  {
    return addCondition_(column, "IN", ColumnType::STRING, values);
  }

  const FilterExpression& ParquetFilter::expression() const
  {
    return expr_;
  }

  bool ParquetFilter::empty() const
  {
    return expr_.empty();
  }

  FilterExpression ParquetFilter::merge(const FilterExpression& lhs,
                                        const FilterExpression& rhs,
                                        const String& connector)
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

  ParquetFilter& ParquetFilter::addCondition_(const String& column,
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

  // ParquetFilterBuilder
  ParquetFilterBuilder& ParquetFilterBuilder::andNext()
  {
    filter_.andNext();
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::orNext()
  {
    filter_.orNext();
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::eq(const String& column, Int64 value)
  {
    filter_.eq(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::ne(const String& column, Int64 value)
  {
    filter_.ne(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::lt(const String& column, Int64 value)
  {
    filter_.lt(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::le(const String& column, Int64 value)
  {
    filter_.le(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::gt(const String& column, Int64 value)
  {
    filter_.gt(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::ge(const String& column, Int64 value)
  {
    filter_.ge(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::in(const String& column, const std::vector<Int64>& values)
  {
    filter_.in(column, values);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::eq(const String& column, const String& value)
  {
    filter_.eq(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::ne(const String& column, const String& value)
  {
    filter_.ne(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::lt(const String& column, const String& value)
  {
    filter_.lt(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::le(const String& column, const String& value)
  {
    filter_.le(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::gt(const String& column, const String& value)
  {
    filter_.gt(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::ge(const String& column, const String& value)
  {
    filter_.ge(column, value);
    return *this;
  }

  ParquetFilterBuilder& ParquetFilterBuilder::in(const String& column, const std::vector<String>& values)
  {
    filter_.in(column, values);
    return *this;
  }

  const ParquetFilter& ParquetFilterBuilder::filter() const
  {
    return filter_;
  }

  bool ParquetFilterBuilder::empty() const
  {
    return filter_.empty();
  }
} // namespace OpenMS
