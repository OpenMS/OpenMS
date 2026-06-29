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
#include <OpenMS/FORMAT/ParquetFile.h>
///////////////////////////

#include <arrow/api.h>
#include <arrow/io/file.h>

using namespace OpenMS;
using namespace std;

namespace
{
  // small single-value Array builders (build cannot fail for valid appends)
  std::shared_ptr<arrow::Array> mkI32(int32_t v)  { arrow::Int32Builder b;   (void)b.Append(v); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }
  std::shared_ptr<arrow::Array> mkI64(int64_t v)  { arrow::Int64Builder b;   (void)b.Append(v); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }
  std::shared_ptr<arrow::Array> mkI64Null()       { arrow::Int64Builder b;   (void)b.AppendNull(); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }
  std::shared_ptr<arrow::Array> mkDbl(double v)   { arrow::DoubleBuilder b;  (void)b.Append(v); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }
  std::shared_ptr<arrow::Array> mkBool(bool v)    { arrow::BooleanBuilder b; (void)b.Append(v); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }
  std::shared_ptr<arrow::Array> mkStr(const std::string& v) { arrow::StringBuilder b; (void)b.Append(v); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a); return a; }

  // a 2-row table: name(utf8), id(int64), val(float64)
  std::shared_ptr<arrow::Table> mkTable()
  {
    arrow::StringBuilder nb; (void)nb.Append("alpha"); (void)nb.Append("beta");
    arrow::Int64Builder  ib; (void)ib.Append(static_cast<int64_t>(10)); (void)ib.Append(static_cast<int64_t>(20));
    arrow::DoubleBuilder db; (void)db.Append(1.5); (void)db.Append(2.5);
    std::shared_ptr<arrow::Array> na, ia, da;
    (void)nb.Finish(&na); (void)ib.Finish(&ia); (void)db.Finish(&da);
    auto schema = arrow::schema({arrow::field("name", arrow::utf8()),
                                 arrow::field("id", arrow::int64()),
                                 arrow::field("val", arrow::float64())});
    return arrow::Table::Make(schema, {na, ia, da});
  }
}

START_TEST(ParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void writeTable(const std::shared_ptr<arrow::Table>& table, const std::string& filename, int64_t row_group_size = 262144)))
{
  auto table = mkTable();
  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  ParquetFile::writeTable(table, fn);
  TEST_EQUAL(File::exists(fn), true)
  // the written file reports the right number of rows via the metadata reader
  TEST_EQUAL(ParquetFile::rowCount(fn), 2)
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> readTable(const std::string& filename)))
{
  auto table = mkTable();
  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  ParquetFile::writeTable(table, fn);

  auto rt = ParquetFile::readTable(fn);
  TEST_NOT_EQUAL(rt, nullptr)
  TEST_EQUAL(rt->num_rows(), 2)
  TEST_EQUAL(rt->num_columns(), 3)

  auto cname = ParquetFile::getColumn(rt, "name");
  TEST_EQUAL(ParquetFile::getString(cname, 0), "alpha")
  TEST_EQUAL(ParquetFile::getString(cname, 1), "beta")
  auto cid = ParquetFile::getColumn(rt, "id");
  TEST_EQUAL(ParquetFile::getInt64(cid, 0, -1, false), 10)
  TEST_EQUAL(ParquetFile::getInt64(cid, 1, -1, false), 20)
  auto cval = ParquetFile::getColumn(rt, "val");
  TEST_REAL_SIMILAR(ParquetFile::getDouble(cval, 0, -1.0, false), 1.5)
  TEST_REAL_SIMILAR(ParquetFile::getDouble(cval, 1, -1.0, false), 2.5)
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Table> readTable(const std::shared_ptr<arrow::io::RandomAccessFile>& infile)))
{
  auto table = mkTable();
  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  ParquetFile::writeTable(table, fn);

  auto maybe_file = arrow::io::ReadableFile::Open(fn);
  TEST_EQUAL(maybe_file.ok(), true)
  std::shared_ptr<arrow::io::RandomAccessFile> raf = *maybe_file;
  auto rt = ParquetFile::readTable(raf);
  TEST_NOT_EQUAL(rt, nullptr)
  TEST_EQUAL(rt->num_rows(), 2)
  TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(rt, "name"), 0), "alpha")
}
END_SECTION

START_SECTION((static int64_t rowCount(const std::string& filename)))
{
  auto table = mkTable();
  std::string fn;
  NEW_TMP_FILE(fn)
  fn += ".parquet";
  ParquetFile::writeTable(table, fn);
  TEST_EQUAL(ParquetFile::rowCount(fn), 2)
  // a non-existent file reports 0 rows (no throw)
  TEST_EQUAL(ParquetFile::rowCount("/no/such/file_xyz.parquet"), 0)
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Array> getColumn(const std::shared_ptr<arrow::Table>& table, const std::string& name)))
{
  auto table = mkTable();
  auto col = ParquetFile::getColumn(table, "id");
  TEST_NOT_EQUAL(col, nullptr)
  TEST_EQUAL(ParquetFile::getInt64(col, 0, -1, false), 10)
  // a required column that is absent throws
  TEST_EXCEPTION(Exception::MissingInformation, ParquetFile::getColumn(table, "does_not_exist"))
}
END_SECTION

START_SECTION((static std::shared_ptr<arrow::Array> getOptionalColumn(const std::shared_ptr<arrow::Table>& table, const std::string& name)))
{
  auto table = mkTable();
  auto col = ParquetFile::getOptionalColumn(table, "name");
  TEST_NOT_EQUAL(col, nullptr)
  TEST_EQUAL(ParquetFile::getString(col, 0), "alpha")
  // an absent optional column returns nullptr (no throw)
  auto missing = ParquetFile::getOptionalColumn(table, "does_not_exist");
  TEST_EQUAL(missing == nullptr, true)
}
END_SECTION

START_SECTION((static int64_t getInt64(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_value, bool allow_null)))
{
  // plain int64
  TEST_EQUAL(ParquetFile::getInt64(mkI64(123), 0, -1, false), 123)
  // type coercion from a narrower integer column (Int32 -> Int64)
  TEST_EQUAL(ParquetFile::getInt64(mkI32(7), 0, -1, false), 7)
  // null + allow_null -> default; null + !allow_null -> throw
  TEST_EQUAL(ParquetFile::getInt64(mkI64Null(), 0, 42, true), 42)
  TEST_EXCEPTION(Exception::MissingInformation, ParquetFile::getInt64(mkI64Null(), 0, 0, false))
  // a null array pointer with allow_null returns the default
  std::shared_ptr<arrow::Array> none;
  TEST_EQUAL(ParquetFile::getInt64(none, 0, 99, true), 99)
}
END_SECTION

START_SECTION((static double getDouble(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_value, bool allow_null)))
{
  TEST_REAL_SIMILAR(ParquetFile::getDouble(mkDbl(3.25), 0, -1.0, false), 3.25)
  // integer columns are coerced to double
  TEST_REAL_SIMILAR(ParquetFile::getDouble(mkI64(5), 0, -1.0, false), 5.0)
  TEST_REAL_SIMILAR(ParquetFile::getDouble(mkI64Null(), 0, -2.5, true), -2.5)
  TEST_EXCEPTION(Exception::MissingInformation, ParquetFile::getDouble(mkI64Null(), 0, 0.0, false))
}
END_SECTION

START_SECTION((static bool getBool(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_value, bool allow_null)))
{
  TEST_EQUAL(ParquetFile::getBool(mkBool(true), 0, false, false), true)
  TEST_EQUAL(ParquetFile::getBool(mkBool(false), 0, true, false), false)
  // integer columns coerce: non-zero -> true, zero -> false
  TEST_EQUAL(ParquetFile::getBool(mkI64(5), 0, false, false), true)
  TEST_EQUAL(ParquetFile::getBool(mkI64(0), 0, true, false), false)
  TEST_EQUAL(ParquetFile::getBool(mkI64Null(), 0, true, true), true)
}
END_SECTION

START_SECTION((static std::string getString(const std::shared_ptr<arrow::Array>& array, int64_t row)))
{
  TEST_EQUAL(ParquetFile::getString(mkStr("hello"), 0), "hello")
  // null value / null array -> empty string
  arrow::StringBuilder b; (void)b.AppendNull(); std::shared_ptr<arrow::Array> a; (void)b.Finish(&a);
  TEST_EQUAL(ParquetFile::getString(a, 0), "")
  std::shared_ptr<arrow::Array> none;
  TEST_EQUAL(ParquetFile::getString(none, 0), "")
}
END_SECTION

START_SECTION((static std::vector<std::string> getStringList(const std::shared_ptr<arrow::Array>& array, int64_t row)))
{
  // a String column is split on ';'
  auto v = ParquetFile::getStringList(mkStr("a;b;c"), 0);
  TEST_EQUAL(v.size(), 3)
  TEST_EQUAL(v[0], "a")
  TEST_EQUAL(v[1], "b")
  TEST_EQUAL(v[2], "c")
  // null array -> empty list
  std::shared_ptr<arrow::Array> none;
  TEST_EQUAL(ParquetFile::getStringList(none, 0).empty(), true)
}
END_SECTION

START_SECTION((static std::string jsonEscape(const std::string& input)))
{
  TEST_EQUAL(ParquetFile::jsonEscape("plain"), "plain")
  TEST_EQUAL(ParquetFile::jsonEscape("a\"b"), "a\\\"b")   // double quote escaped
  TEST_EQUAL(ParquetFile::jsonEscape("a\\b"), "a\\\\b")   // backslash escaped
}
END_SECTION

START_SECTION((static void appendOrThrow(const arrow::Status& status, const std::string& column)))
{
  // an OK status passes through silently
  ParquetFile::appendOrThrow(arrow::Status::OK(), "col");
  TEST_EQUAL(true, true) // reached -> no throw
  // a non-OK status raises Exception::InvalidValue
  TEST_EXCEPTION(Exception::InvalidValue,
                 ParquetFile::appendOrThrow(arrow::Status::Invalid("boom"), "col"))
}
END_SECTION

START_SECTION((template <typename BuilderT> static std::shared_ptr<arrow::Array> finishArray(BuilderT& builder, const std::string& name)))
{
  arrow::Int64Builder b;
  TEST_EQUAL(b.Append(static_cast<int64_t>(77)).ok(), true)
  auto arr = ParquetFile::finishArray(b, "id");
  TEST_NOT_EQUAL(arr, nullptr)
  TEST_EQUAL(arr->length(), 1)
  TEST_EQUAL(ParquetFile::getInt64(arr, 0, -1, false), 77)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
