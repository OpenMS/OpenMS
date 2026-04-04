#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <memory>
#include <arrow/status.h>

using namespace OpenMS;

// Tests Arrow/Parquet I/O — the format OpenMS actually uses.

static ::arrow::Status RunParquetRoundtrip(const std::string& in_path,
                                           const std::string& out_path)
{
  // Build a small table
  arrow::Int8Builder int8builder;

  int8_t days_raw[5] = {1, 12, 17, 23, 28};
  ARROW_RETURN_NOT_OK(int8builder.AppendValues(days_raw, 5));
  std::shared_ptr<arrow::Array> days;
  ARROW_ASSIGN_OR_RAISE(days, int8builder.Finish());

  int8_t months_raw[5] = {1, 3, 5, 7, 1};
  ARROW_RETURN_NOT_OK(int8builder.AppendValues(months_raw, 5));
  std::shared_ptr<arrow::Array> months;
  ARROW_ASSIGN_OR_RAISE(months, int8builder.Finish());

  arrow::Int16Builder int16builder;
  int16_t years_raw[5] = {1990, 2000, 1995, 2000, 1995};
  ARROW_RETURN_NOT_OK(int16builder.AppendValues(years_raw, 5));
  std::shared_ptr<arrow::Array> years;
  ARROW_ASSIGN_OR_RAISE(years, int16builder.Finish());

  std::vector<std::shared_ptr<arrow::Array>> columns = {days, months, years};

  auto field_day = arrow::field("Day", arrow::int8());
  auto field_month = arrow::field("Month", arrow::int8());
  auto field_year = arrow::field("Year", arrow::int16());
  auto schema = arrow::schema({field_day, field_month, field_year});

  auto table = arrow::Table::Make(schema, columns);

  // Write Parquet
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open(in_path));
  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 5));

  // Read Parquet
  std::shared_ptr<arrow::io::ReadableFile> infile;
  ARROW_ASSIGN_OR_RAISE(infile, arrow::io::ReadableFile::Open(in_path));
  std::unique_ptr<parquet::arrow::FileReader> reader;
  PARQUET_ASSIGN_OR_THROW(reader,
                          parquet::arrow::OpenFile(infile, arrow::default_memory_pool()));

  std::shared_ptr<arrow::Table> read_table;
  PARQUET_THROW_NOT_OK(reader->ReadTable(&read_table));

  // Write back out to verify roundtrip
  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open(out_path));
  PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
      *read_table, arrow::default_memory_pool(), outfile, 5));

  return ::arrow::Status::OK();
}

START_TEST(Arrow, "$Id$")

START_SECTION(RunParquetRoundtrip())
{
  String in_parquet, out_parquet;
  NEW_TMP_FILE(in_parquet);
  NEW_TMP_FILE(out_parquet);
  TEST_EQUAL(RunParquetRoundtrip(in_parquet, out_parquet).ok(), true)
}
END_SECTION

END_TEST
