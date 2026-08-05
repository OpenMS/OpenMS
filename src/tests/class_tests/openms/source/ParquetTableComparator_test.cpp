// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/ParquetTableComparator.h>

#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>

#include <limits>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// A minimal three-column table: key (string), value (float32), tags (list<string>).
  std::shared_ptr<arrow::Table> makeTable(const std::vector<std::string>& keys,
                                          const std::vector<double>& values,
                                          const std::vector<std::vector<std::string>>& tags)
  {
    arrow::StringBuilder key_b;
    for (const auto& k : keys) { (void)key_b.Append(k); }
    std::shared_ptr<arrow::Array> key_a;
    (void)key_b.Finish(&key_a);

    arrow::FloatBuilder value_b;
    for (double v : values) { (void)value_b.Append(static_cast<float>(v)); }
    std::shared_ptr<arrow::Array> value_a;
    (void)value_b.Finish(&value_a);

    auto tag_values = std::make_shared<arrow::StringBuilder>();
    arrow::ListBuilder tag_b(arrow::default_memory_pool(), tag_values);
    for (const auto& row : tags)
    {
      (void)tag_b.Append();
      for (const auto& t : row) { (void)tag_values->Append(t); }
    }
    std::shared_ptr<arrow::Array> tag_a;
    (void)tag_b.Finish(&tag_a);

    auto schema = arrow::schema({arrow::field("key", arrow::utf8(), false),
                                 arrow::field("value", arrow::float32(), true),
                                 arrow::field("tags", arrow::list(arrow::utf8()), true)});
    return arrow::Table::Make(schema, {key_a, value_a, tag_a});
  }

  /// key (string) + value (float64), for NaN / infinity cases.
  std::shared_ptr<arrow::Table> makeDoubleTable(const std::vector<std::string>& keys,
                                                const std::vector<double>& values)
  {
    arrow::StringBuilder key_b;
    for (const auto& k : keys) { (void)key_b.Append(k); }
    std::shared_ptr<arrow::Array> key_a;
    (void)key_b.Finish(&key_a);

    arrow::DoubleBuilder value_b;
    for (double v : values) { (void)value_b.Append(v); }
    std::shared_ptr<arrow::Array> value_a;
    (void)value_b.Finish(&value_a);

    return arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false),
                                             arrow::field("value", arrow::float64(), true)}),
                              {key_a, value_a});
  }

  /// key (string) + value (int64), for the wide-integer case.
  std::shared_ptr<arrow::Table> makeIntTable(const std::vector<std::string>& keys,
                                             const std::vector<int64_t>& values)
  {
    arrow::StringBuilder key_b;
    for (const auto& k : keys) { (void)key_b.Append(k); }
    std::shared_ptr<arrow::Array> key_a;
    (void)key_b.Finish(&key_a);

    arrow::Int64Builder value_b;
    for (int64_t v : values) { (void)value_b.Append(v); }
    std::shared_ptr<arrow::Array> value_a;
    (void)value_b.Finish(&value_a);

    return arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false),
                                             arrow::field("value", arrow::int64(), true)}),
                              {key_a, value_a});
  }

  /// key (string) + nums (list<double>), for multiset pairing of numeric lists.
  std::shared_ptr<arrow::Table> makeNumListTable(const std::vector<std::string>& keys,
                                                 const std::vector<std::vector<double>>& nums)
  {
    arrow::StringBuilder key_b;
    for (const auto& k : keys) { (void)key_b.Append(k); }
    std::shared_ptr<arrow::Array> key_a;
    (void)key_b.Finish(&key_a);

    auto num_values = std::make_shared<arrow::DoubleBuilder>();
    arrow::ListBuilder num_b(arrow::default_memory_pool(), num_values);
    for (const auto& row : nums)
    {
      (void)num_b.Append();
      for (double v : row) { (void)num_values->Append(v); }
    }
    std::shared_ptr<arrow::Array> num_a;
    (void)num_b.Finish(&num_a);

    return arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false),
                                             arrow::field("nums", arrow::list(arrow::float64()), true)}),
                              {key_a, num_a});
  }

  /// Write @p table to @p path. The caller supplies the path via NEW_TMP_FILE_EXT, which is
  /// line-based and so must be expanded at the call site rather than inside this helper.
  bool writeTo(const std::shared_ptr<arrow::Table>& table, const std::string& path)
  {
    return ArrowIOHelpers::writeTableToParquet(table, path);
  }
}

START_TEST(ParquetTableComparator, "$Id$")

/////////////////////////////////////////////////////////////

ParquetDiffSettings pk_settings;
pk_settings.primary_key = {"key"};

START_SECTION((static ParquetDiffResult compare(const std::string& file_1, const std::string& file_2, const ParquetDiffSettings& settings)))
{
  // identical content -> equal
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p2"}, {10.0, 20.0}, {{"x"}, {"y"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p2"}, {10.0, 20.0}, {{"x"}, {"y"}}), b))
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, true)
  TEST_EQUAL(r.rows_compared, 2)
  TEST_EQUAL(r.schema_errors.empty(), true)

  // row order must not matter: the whole point of keying on the primary key
  std::string c;
  NEW_TMP_FILE_EXT(c, "parquet")
  TEST_TRUE(writeTo(makeTable({"p2", "p1"}, {20.0, 10.0}, {{"y"}, {"x"}}), c))
  r = ParquetTableComparator::compare(a, c, pk_settings);
  TEST_EQUAL(r.equal, true)
  TEST_EQUAL(r.rows_compared, 2)

}
END_SECTION

START_SECTION([EXTRA] numeric tolerance is satisfied by either ratio or absdiff)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {100.0}, {{"x"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {101.0}, {{"x"}}), b))

  // exact comparison fails
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.value_errors.size(), 1)

  // ratio 101/100 = 1.01 accepted
  ParquetDiffSettings ratio_settings = pk_settings;
  ratio_settings.acceptable_ratio = 1.02;
  r = ParquetTableComparator::compare(a, b, ratio_settings);
  TEST_EQUAL(r.equal, true)

  // absdiff alone is also sufficient
  ParquetDiffSettings abs_settings = pk_settings;
  abs_settings.acceptable_absdiff = 2.0;
  r = ParquetTableComparator::compare(a, b, abs_settings);
  TEST_EQUAL(r.equal, true)

}
END_SECTION

START_SECTION([EXTRA] zero vs. epsilon is only reachable through absdiff)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {0.0}, {{"x"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1e-6}, {{"x"}}), b))

  // no ratio can bridge zero to non-zero
  ParquetDiffSettings huge_ratio = pk_settings;
  huge_ratio.acceptable_ratio = 1e9;
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, huge_ratio);
  TEST_EQUAL(r.equal, false)

  ParquetDiffSettings abs_settings = pk_settings;
  abs_settings.acceptable_absdiff = 1e-5;
  r = ParquetTableComparator::compare(a, b, abs_settings);
  TEST_EQUAL(r.equal, true)

}
END_SECTION

START_SECTION([EXTRA] rows present in only one file are reported as key differences)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p2"}, {1.0, 2.0}, {{"x"}, {"y"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x"}}), b))

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.key_errors.size(), 1)
  TEST_EQUAL(r.rows_compared, 1)
  TEST_EQUAL(r.rows_1, 2)
  TEST_EQUAL(r.rows_2, 1)

}
END_SECTION

START_SECTION([EXTRA] a duplicate primary key is itself an error)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p1"}, {1.0, 2.0}, {{"x"}, {"y"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x"}}), b))

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_NOT_EQUAL(r.key_errors.size(), 0)

}
END_SECTION

START_SECTION([EXTRA] a list-valued cell compares element-wise or as a multiset on request)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x", "y"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"y", "x"}}), b))

  // sequence semantics by default: the reordering is a difference
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)

  ParquetDiffSettings unordered = pk_settings;
  unordered.unordered_lists = true;
  r = ParquetTableComparator::compare(a, b, unordered);
  TEST_EQUAL(r.equal, true)

}
END_SECTION

START_SECTION([EXTRA] a column present in only one file is a schema difference)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x"}}), a))

  arrow::StringBuilder key_b;
  (void)key_b.Append("p1");
  std::shared_ptr<arrow::Array> key_a;
  (void)key_b.Finish(&key_a);
  auto narrow = arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false)}),
                                   {key_a});
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(narrow, b))

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.schema_errors.size(), 2) // 'value' and 'tags' only in file 1

}
END_SECTION

START_SECTION([EXTRA] ignored columns are excluded from value comparison)
{
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {99.0}, {{"x"}}), b))

  ParquetDiffSettings ignore_value = pk_settings;
  ignore_value.ignore_columns = {"value"};
  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, ignore_value);
  TEST_EQUAL(r.equal, true)

}
END_SECTION

START_SECTION((static std::vector<std::string> qpxPrimaryKey(const std::string& view)))
{
  const auto pg = ParquetTableComparator::qpxPrimaryKey("pg");
  TEST_EQUAL(pg.size(), 3)
  TEST_STRING_EQUAL(pg[0], "anchor_protein")
  TEST_STRING_EQUAL(pg[1], "grouped_runs")
  TEST_STRING_EQUAL(pg[2], "label")

  const auto feature = ParquetTableComparator::qpxPrimaryKey("feature");
  TEST_EQUAL(feature.size(), 4)

  TEST_EXCEPTION(Exception::IllegalArgument, ParquetTableComparator::qpxPrimaryKey("nosuchview"))
}
END_SECTION

START_SECTION((static std::string viewFromFileType(const std::string& file_type)))
{
  TEST_STRING_EQUAL(ParquetTableComparator::viewFromFileType("pg_file"), "pg")
  TEST_STRING_EQUAL(ParquetTableComparator::viewFromFileType("feature_file"), "feature")
  TEST_STRING_EQUAL(ParquetTableComparator::viewFromFileType("psm_file"), "psm")
  TEST_STRING_EQUAL(ParquetTableComparator::viewFromFileType("something_else"), "")
}
END_SECTION

START_SECTION((static ParquetDiffResult validate(const std::string& file, const std::string& view, const ParquetDiffSettings& settings)))
{
  // a table that is plainly not a QPX pg view: every required column is missing
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1"}, {1.0}, {{"x"}}), a))

  ParquetDiffSettings settings;
  settings.primary_key = {"key"}; // the file has no QPX key; check schema reporting only
  settings.max_reported = 0;      // unlimited
  const ParquetDiffResult r = ParquetTableComparator::validate(a, "pg", settings);
  TEST_EQUAL(r.equal, false)
  TEST_NOT_EQUAL(r.schema_errors.size(), 0)

  TEST_EXCEPTION(Exception::IllegalArgument,
                 ParquetTableComparator::validate(a, "nosuchview", settings))

}
END_SECTION

START_SECTION([EXTRA] schema_only reports schema drift without building keys or comparing values)
{
  // Differs in BOTH schema and values, and carries a duplicate primary key: under schema_only
  // only the schema difference may be reported.
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p1"}, {1.0, 2.0}, {{"x"}, {"y"}}), a))

  arrow::StringBuilder key_b;
  (void)key_b.Append("p1");
  (void)key_b.Append("p1");
  std::shared_ptr<arrow::Array> key_a;
  (void)key_b.Finish(&key_a);
  auto narrow = arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false)}),
                                   {key_a});
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(narrow, b))

  ParquetDiffSettings schema_only = pk_settings;
  schema_only.schema_only = true;
  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, schema_only);
  TEST_EQUAL(r.equal, false)
  TEST_NOT_EQUAL(r.schema_errors.size(), 0)
  TEST_EQUAL(r.key_errors.size(), 0)     // no key building, so no duplicate-key report
  TEST_EQUAL(r.value_errors.size(), 0)
  TEST_EQUAL(r.rows_compared, 0)
}
END_SECTION

START_SECTION([EXTRA] max_reported caps each diagnostic and the remainder is counted as suppressed)
{
  // Three rows that all differ in 'value'.
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p2", "p3"}, {1.0, 2.0, 3.0}, {{"x"}, {"y"}, {"z"}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeTable({"p1", "p2", "p3"}, {9.0, 8.0, 7.0}, {{"x"}, {"y"}, {"z"}}), b))

  ParquetDiffSettings uncapped = pk_settings;
  uncapped.max_reported = 0;   // unlimited
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, uncapped);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.value_errors.size(), 3)
  TEST_EQUAL(r.suppressed, 0)

  ParquetDiffSettings capped = pk_settings;
  capped.max_reported = 1;
  r = ParquetTableComparator::compare(a, b, capped);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.value_errors.size(), 1)
  TEST_EQUAL(r.suppressed, 2)  // the two it did not list
}
END_SECTION

START_SECTION([EXTRA] two NaNs compare equal and a NaN never matches a number)
{
  const double nan_v = std::numeric_limits<double>::quiet_NaN();
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeDoubleTable({"p1"}, {nan_v}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeDoubleTable({"p1"}, {nan_v}), b))

  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, true)

  // ... but NaN against a real number is a difference no tolerance can accept
  std::string c;
  NEW_TMP_FILE_EXT(c, "parquet")
  TEST_TRUE(writeTo(makeDoubleTable({"p1"}, {1.0}), c))
  ParquetDiffSettings loose = pk_settings;
  loose.acceptable_absdiff = 1e9;
  r = ParquetTableComparator::compare(a, c, loose);
  TEST_EQUAL(r.equal, false)
}
END_SECTION

START_SECTION([EXTRA] a 64-bit integer difference survives the widening to double)
{
  // 2^53 and 2^53+1 are the same double; the difference must still be seen exactly.
  const int64_t big = 9007199254740992LL;   // 2^53
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeIntTable({"p1"}, {big}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeIntTable({"p1"}, {big + 1}), b))

  // exact comparison must report the difference rather than lose it to rounding
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)

  // and a tolerance of 1 must accept it
  ParquetDiffSettings tol = pk_settings;
  tol.acceptable_absdiff = 1.0;
  r = ParquetTableComparator::compare(a, b, tol);
  TEST_EQUAL(r.equal, true)
}
END_SECTION

START_SECTION([EXTRA] unordered numeric lists are paired by value not by canonical text)
{
  // Lexically "10" < "20" but "20.1" < "9.99"; pairing by text would compare 10 against 20.1.
  std::string a;
  NEW_TMP_FILE_EXT(a, "parquet")
  TEST_TRUE(writeTo(makeNumListTable({"p1"}, {{10.0, 20.0}}), a))
  std::string b;
  NEW_TMP_FILE_EXT(b, "parquet")
  TEST_TRUE(writeTo(makeNumListTable({"p1"}, {{20.1, 9.99}}), b))

  ParquetDiffSettings unordered = pk_settings;
  unordered.unordered_lists = true;
  unordered.acceptable_absdiff = 0.11;
  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, unordered);
  TEST_EQUAL(r.equal, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
