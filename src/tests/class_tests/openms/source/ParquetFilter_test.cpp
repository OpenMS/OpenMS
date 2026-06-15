// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ParquetFilter.h>
///////////////////////////

#include <OpenMS/CONCEPT/Types.h>

using namespace OpenMS;
using namespace std;

START_TEST(ParquetFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ParquetFilter* ptr = nullptr;
ParquetFilter* null_ptr = nullptr;

START_SECTION(ParquetFilter())
{
  ptr = new ParquetFilter();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // a fresh filter has no conditions
  TEST_EQUAL(ptr->empty(), true)
  TEST_EQUAL(ptr->expression().empty(), true)
  TEST_EQUAL(ptr->expression().conditions.size(), 0)
  TEST_EQUAL(ptr->expression().connectors.size(), 0)
  delete ptr;
}
END_SECTION

START_SECTION((ParquetFilter& eq(const std::string& column, Int64 value)))
{
  ParquetFilter f;
  f.eq("PRECURSOR_ID", Int64(1318));
  TEST_EQUAL(f.empty(), false)
  const Condition& c = f.expression().conditions[0];
  TEST_STRING_EQUAL(c.column, "PRECURSOR_ID")
  TEST_STRING_EQUAL(c.op, "=")
  TEST_EQUAL(c.type == ColumnType::INT, true)
  TEST_EQUAL(c.values.size(), 1)
  TEST_STRING_EQUAL(c.values[0], "1318")
}
END_SECTION

START_SECTION((ParquetFilter& ne(const std::string& column, Int64 value)))
{
  ParquetFilter f; f.ne("ID", Int64(-5));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "!=")
  TEST_STRING_EQUAL(f.expression().conditions[0].values[0], "-5")
  TEST_EQUAL(f.expression().conditions[0].type == ColumnType::INT, true)
}
END_SECTION

START_SECTION((ParquetFilter& lt(const std::string& column, Int64 value)))
{
  ParquetFilter f; f.lt("ID", Int64(10));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "<")
}
END_SECTION

START_SECTION((ParquetFilter& le(const std::string& column, Int64 value)))
{
  ParquetFilter f; f.le("ID", Int64(10));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "<=")
}
END_SECTION

START_SECTION((ParquetFilter& gt(const std::string& column, Int64 value)))
{
  ParquetFilter f; f.gt("ID", Int64(10));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, ">")
}
END_SECTION

START_SECTION((ParquetFilter& ge(const std::string& column, Int64 value)))
{
  ParquetFilter f; f.ge("ID", Int64(10));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, ">=")
}
END_SECTION

START_SECTION((ParquetFilter& in(const std::string& column, const std::vector<Int64>& values)))
{
  ParquetFilter f; f.in("ID", std::vector<Int64>{1, 2, 3});
  const Condition& c = f.expression().conditions[0];
  TEST_STRING_EQUAL(c.op, "IN")
  TEST_EQUAL(c.type == ColumnType::INT, true)
  TEST_EQUAL(c.values.size(), 3)
  TEST_STRING_EQUAL(c.values[0], "1")
  TEST_STRING_EQUAL(c.values[2], "3")
}
END_SECTION

START_SECTION((ParquetFilter& eq(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.eq("ANNOTATION", std::string("y3^1"));
  const Condition& c = f.expression().conditions[0];
  TEST_STRING_EQUAL(c.op, "=")
  TEST_EQUAL(c.type == ColumnType::STRING, true)
  TEST_STRING_EQUAL(c.values[0], "y3^1")
}
END_SECTION

START_SECTION((ParquetFilter& ne(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.ne("ANNOTATION", std::string("b2"));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "!=")
  TEST_EQUAL(f.expression().conditions[0].type == ColumnType::STRING, true)
}
END_SECTION

START_SECTION((ParquetFilter& lt(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.lt("NAME", std::string("m"));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "<")
  TEST_EQUAL(f.expression().conditions[0].type == ColumnType::STRING, true)
}
END_SECTION

START_SECTION((ParquetFilter& le(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.le("NAME", std::string("m"));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, "<=")
}
END_SECTION

START_SECTION((ParquetFilter& gt(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.gt("NAME", std::string("m"));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, ">")
}
END_SECTION

START_SECTION((ParquetFilter& ge(const std::string& column, const std::string& value)))
{
  ParquetFilter f; f.ge("NAME", std::string("m"));
  TEST_STRING_EQUAL(f.expression().conditions[0].op, ">=")
}
END_SECTION

START_SECTION((ParquetFilter& in(const std::string& column, const std::vector<std::string>& values)))
{
  ParquetFilter f; f.in("NAME", std::vector<std::string>{"a", "b"});
  const Condition& c = f.expression().conditions[0];
  TEST_STRING_EQUAL(c.op, "IN")
  TEST_EQUAL(c.type == ColumnType::STRING, true)
  TEST_EQUAL(c.values.size(), 2)
  TEST_STRING_EQUAL(c.values[1], "b")
}
END_SECTION

START_SECTION((ParquetFilter& andNext()))
{
  // default connector between consecutive conditions is AND
  ParquetFilter f;
  f.eq("A", Int64(1)).eq("B", Int64(2));
  TEST_EQUAL(f.expression().conditions.size(), 2)
  TEST_EQUAL(f.expression().connectors.size(), 1)
  TEST_STRING_EQUAL(f.expression().connectors[0], "AND")

  // an explicit andNext() is equivalent to the default
  ParquetFilter g;
  g.eq("A", Int64(1)).andNext().eq("B", Int64(2));
  TEST_STRING_EQUAL(g.expression().connectors[0], "AND")
}
END_SECTION

START_SECTION((ParquetFilter& orNext()))
{
  ParquetFilter f;
  // orNext() only applies to the immediately-following condition; the
  // connector then resets to AND for the next one.
  f.eq("A", Int64(1)).orNext().eq("B", Int64(2)).eq("C", Int64(3));
  TEST_EQUAL(f.expression().conditions.size(), 3)
  TEST_EQUAL(f.expression().connectors.size(), 2)
  TEST_STRING_EQUAL(f.expression().connectors[0], "OR")
  TEST_STRING_EQUAL(f.expression().connectors[1], "AND")
}
END_SECTION

START_SECTION((const FilterExpression& expression() const))
{
  ParquetFilter f;
  f.eq("A", Int64(1)).orNext().eq("B", std::string("x"));
  const FilterExpression& e = f.expression();
  TEST_EQUAL(e.conditions.size(), 2)
  TEST_STRING_EQUAL(e.conditions[0].column, "A")
  TEST_STRING_EQUAL(e.conditions[1].column, "B")
  TEST_STRING_EQUAL(e.connectors[0], "OR")
}
END_SECTION

START_SECTION((bool empty() const))
{
  ParquetFilter f;
  TEST_EQUAL(f.empty(), true)
  f.eq("A", Int64(1));
  TEST_EQUAL(f.empty(), false)
}
END_SECTION

START_SECTION((static FilterExpression merge(const FilterExpression& lhs, const FilterExpression& rhs, const std::string& connector = "AND")))
{
  ParquetFilter a; a.eq("A", Int64(1)).eq("B", Int64(2)); // 2 conditions, 1 AND connector
  ParquetFilter b; b.eq("C", Int64(3));                   // 1 condition, 0 connectors

  FilterExpression m = ParquetFilter::merge(a.expression(), b.expression(), "OR");
  TEST_EQUAL(m.conditions.size(), 3)
  TEST_EQUAL(m.connectors.size(), 2)
  TEST_STRING_EQUAL(m.connectors[0], "AND") // preserved from lhs
  TEST_STRING_EQUAL(m.connectors[1], "OR")  // the merge connector
  TEST_STRING_EQUAL(m.conditions[2].column, "C")

  // merging with an empty side returns the other side unchanged
  FilterExpression empty_expr;
  FilterExpression m2 = ParquetFilter::merge(empty_expr, b.expression(), "AND");
  TEST_EQUAL(m2.conditions.size(), 1)
  TEST_EQUAL(m2.connectors.size(), 0)
  FilterExpression m3 = ParquetFilter::merge(a.expression(), empty_expr, "AND");
  TEST_EQUAL(m3.conditions.size(), 2)
  TEST_EQUAL(m3.connectors.size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ParquetFilterBuilder (fluent wrapper forwarding to ParquetFilter)
/////////////////////////////////////////////////////////////

START_SECTION((ParquetFilterBuilder& eq(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb;
  fb.eq("PRECURSOR_ID", Int64(1318)).andNext().eq("ANNOTATION", std::string("y3^1"));
  TEST_EQUAL(fb.empty(), false)
  const ParquetFilter& pf = fb.filter();
  TEST_EQUAL(pf.expression().conditions.size(), 2)
  TEST_STRING_EQUAL(pf.expression().conditions[0].column, "PRECURSOR_ID")
  TEST_EQUAL(pf.expression().conditions[0].type == ColumnType::INT, true)
  TEST_STRING_EQUAL(pf.expression().conditions[1].column, "ANNOTATION")
  TEST_EQUAL(pf.expression().conditions[1].type == ColumnType::STRING, true)
  TEST_STRING_EQUAL(pf.expression().connectors[0], "AND")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& ne(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb; fb.ne("ID", Int64(5));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "!=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& lt(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb; fb.lt("ID", Int64(5));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "<")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& le(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb; fb.le("ID", Int64(5));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "<=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& gt(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb; fb.gt("ID", Int64(5));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, ">")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& ge(const std::string& column, Int64 value)))
{
  ParquetFilterBuilder fb; fb.ge("ID", Int64(5));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, ">=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& in(const std::string& column, const std::vector<Int64>& values)))
{
  ParquetFilterBuilder fb; fb.in("ID", std::vector<Int64>{7, 8});
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "IN")
  TEST_EQUAL(fb.filter().expression().conditions[0].values.size(), 2)
}
END_SECTION

START_SECTION((ParquetFilterBuilder& eq(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.eq("NAME", std::string("x"));
  TEST_EQUAL(fb.filter().expression().conditions[0].type == ColumnType::STRING, true)
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& ne(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.ne("NAME", std::string("x"));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "!=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& lt(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.lt("NAME", std::string("x"));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "<")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& le(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.le("NAME", std::string("x"));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "<=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& gt(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.gt("NAME", std::string("x"));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, ">")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& ge(const std::string& column, const std::string& value)))
{
  ParquetFilterBuilder fb; fb.ge("NAME", std::string("x"));
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, ">=")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& in(const std::string& column, const std::vector<std::string>& values)))
{
  ParquetFilterBuilder fb; fb.in("NAME", std::vector<std::string>{"a", "b", "c"});
  TEST_STRING_EQUAL(fb.filter().expression().conditions[0].op, "IN")
  TEST_EQUAL(fb.filter().expression().conditions[0].type == ColumnType::STRING, true)
  TEST_EQUAL(fb.filter().expression().conditions[0].values.size(), 3)
}
END_SECTION

START_SECTION((ParquetFilterBuilder& andNext()))
{
  ParquetFilterBuilder fb;
  fb.eq("A", Int64(1)).andNext().eq("B", Int64(2));
  TEST_STRING_EQUAL(fb.filter().expression().connectors[0], "AND")
}
END_SECTION

START_SECTION((ParquetFilterBuilder& orNext()))
{
  ParquetFilterBuilder fb;
  fb.eq("A", Int64(1)).orNext().eq("B", Int64(2));
  TEST_STRING_EQUAL(fb.filter().expression().connectors[0], "OR")
}
END_SECTION

START_SECTION((const ParquetFilter& filter() const))
{
  ParquetFilterBuilder fb; fb.eq("A", Int64(1));
  TEST_EQUAL(fb.filter().expression().conditions.size(), 1)
}
END_SECTION

START_SECTION((bool empty() const))
{
  ParquetFilterBuilder fb;
  TEST_EQUAL(fb.empty(), true)
  fb.eq("A", Int64(1));
  TEST_EQUAL(fb.empty(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
