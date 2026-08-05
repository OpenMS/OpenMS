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

  /// Write @p table to a fresh temporary file and return its path.
  std::string writeTemp(const std::shared_ptr<arrow::Table>& table, const std::string& stem)
  {
    const std::string path = File::getTempDirectory() + "/" + stem + "_"
                             + File::getUniqueName() + ".parquet";
    ArrowIOHelpers::writeTableToParquet(table, path);
    return path;
  }
}

START_TEST(ParquetTableComparator, "$Id$")

/////////////////////////////////////////////////////////////

ParquetDiffSettings pk_settings;
pk_settings.primary_key = {"key"};

START_SECTION((static ParquetDiffResult compare(const std::string& file_1, const std::string& file_2, const ParquetDiffSettings& settings)))
{
  // identical content -> equal
  const std::string a = writeTemp(makeTable({"p1", "p2"}, {10.0, 20.0}, {{"x"}, {"y"}}), "pd_a");
  const std::string b = writeTemp(makeTable({"p1", "p2"}, {10.0, 20.0}, {{"x"}, {"y"}}), "pd_b");
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, true)
  TEST_EQUAL(r.rows_compared, 2)
  TEST_EQUAL(r.schema_errors.empty(), true)

  // row order must not matter: the whole point of keying on the primary key
  const std::string c = writeTemp(makeTable({"p2", "p1"}, {20.0, 10.0}, {{"y"}, {"x"}}), "pd_c");
  r = ParquetTableComparator::compare(a, c, pk_settings);
  TEST_EQUAL(r.equal, true)
  TEST_EQUAL(r.rows_compared, 2)

  File::remove(a);
  File::remove(b);
  File::remove(c);
}
END_SECTION

START_SECTION([EXTRA] numeric tolerance is satisfied by either ratio or absdiff)
{
  const std::string a = writeTemp(makeTable({"p1"}, {100.0}, {{"x"}}), "pd_tol_a");
  const std::string b = writeTemp(makeTable({"p1"}, {101.0}, {{"x"}}), "pd_tol_b");

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

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] zero vs. epsilon is only reachable through absdiff)
{
  const std::string a = writeTemp(makeTable({"p1"}, {0.0}, {{"x"}}), "pd_zero_a");
  const std::string b = writeTemp(makeTable({"p1"}, {1e-6}, {{"x"}}), "pd_zero_b");

  // no ratio can bridge zero to non-zero
  ParquetDiffSettings huge_ratio = pk_settings;
  huge_ratio.acceptable_ratio = 1e9;
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, huge_ratio);
  TEST_EQUAL(r.equal, false)

  ParquetDiffSettings abs_settings = pk_settings;
  abs_settings.acceptable_absdiff = 1e-5;
  r = ParquetTableComparator::compare(a, b, abs_settings);
  TEST_EQUAL(r.equal, true)

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] rows present in only one file are reported as key differences)
{
  const std::string a = writeTemp(makeTable({"p1", "p2"}, {1.0, 2.0}, {{"x"}, {"y"}}), "pd_k_a");
  const std::string b = writeTemp(makeTable({"p1"}, {1.0}, {{"x"}}), "pd_k_b");

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.key_errors.size(), 1)
  TEST_EQUAL(r.rows_compared, 1)
  TEST_EQUAL(r.rows_1, 2)
  TEST_EQUAL(r.rows_2, 1)

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] a duplicate primary key is itself an error)
{
  const std::string a = writeTemp(makeTable({"p1", "p1"}, {1.0, 2.0}, {{"x"}, {"y"}}), "pd_dup_a");
  const std::string b = writeTemp(makeTable({"p1"}, {1.0}, {{"x"}}), "pd_dup_b");

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_NOT_EQUAL(r.key_errors.size(), 0)

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] a list-valued cell compares element-wise or as a multiset on request)
{
  const std::string a = writeTemp(makeTable({"p1"}, {1.0}, {{"x", "y"}}), "pd_l_a");
  const std::string b = writeTemp(makeTable({"p1"}, {1.0}, {{"y", "x"}}), "pd_l_b");

  // sequence semantics by default: the reordering is a difference
  ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)

  ParquetDiffSettings unordered = pk_settings;
  unordered.unordered_lists = true;
  r = ParquetTableComparator::compare(a, b, unordered);
  TEST_EQUAL(r.equal, true)

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] a column present in only one file is a schema difference)
{
  const std::string a = writeTemp(makeTable({"p1"}, {1.0}, {{"x"}}), "pd_s_a");

  arrow::StringBuilder key_b;
  (void)key_b.Append("p1");
  std::shared_ptr<arrow::Array> key_a;
  (void)key_b.Finish(&key_a);
  auto narrow = arrow::Table::Make(arrow::schema({arrow::field("key", arrow::utf8(), false)}),
                                   {key_a});
  const std::string b = writeTemp(narrow, "pd_s_b");

  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, pk_settings);
  TEST_EQUAL(r.equal, false)
  TEST_EQUAL(r.schema_errors.size(), 2) // 'value' and 'tags' only in file 1

  File::remove(a);
  File::remove(b);
}
END_SECTION

START_SECTION([EXTRA] ignored columns are excluded from value comparison)
{
  const std::string a = writeTemp(makeTable({"p1"}, {1.0}, {{"x"}}), "pd_i_a");
  const std::string b = writeTemp(makeTable({"p1"}, {99.0}, {{"x"}}), "pd_i_b");

  ParquetDiffSettings ignore_value = pk_settings;
  ignore_value.ignore_columns = {"value"};
  const ParquetDiffResult r = ParquetTableComparator::compare(a, b, ignore_value);
  TEST_EQUAL(r.equal, true)

  File::remove(a);
  File::remove(b);
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
  const std::string a = writeTemp(makeTable({"p1"}, {1.0}, {{"x"}}), "pd_v_a");

  ParquetDiffSettings settings;
  settings.primary_key = {"key"}; // the file has no QPX key; check schema reporting only
  settings.max_reported = 0;      // unlimited
  const ParquetDiffResult r = ParquetTableComparator::validate(a, "pg", settings);
  TEST_EQUAL(r.equal, false)
  TEST_NOT_EQUAL(r.schema_errors.size(), 0)

  TEST_EXCEPTION(Exception::IllegalArgument,
                 ParquetTableComparator::validate(a, "nosuchview", settings))

  File::remove(a);
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
