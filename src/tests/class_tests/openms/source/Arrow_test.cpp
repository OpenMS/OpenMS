#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#include <iostream>
#include <vector>
#include <memory>
#include <arrow/status.h>

using namespace OpenMS;

// Adapted from https://arrow.apache.org/docs/cpp/tutorials/io_tutorial.html
// Tests IPC and Parquet I/O — the formats OpenMS actually uses.

static ::arrow::Status GenInitialFile()
{
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

  // Write IPC file
  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_in.arrow"));
  ARROW_ASSIGN_OR_RAISE(auto ipc_writer, arrow::ipc::MakeFileWriter(outfile, schema));
  ARROW_RETURN_NOT_OK(ipc_writer->WriteTable(*table));
  ARROW_RETURN_NOT_OK(ipc_writer->Close());

  // Write Parquet file
  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_in.parquet"));
  PARQUET_THROW_NOT_OK(
      parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 5));

  return ::arrow::Status::OK();
}

static ::arrow::Status RunMain()
{
  ARROW_RETURN_NOT_OK(GenInitialFile());

  // --- IPC read/write ---
  std::shared_ptr<arrow::io::ReadableFile> infile;
  ARROW_ASSIGN_OR_RAISE(infile, arrow::io::ReadableFile::Open(
                                    "test_in.arrow", arrow::default_memory_pool()));
  ARROW_ASSIGN_OR_RAISE(auto ipc_reader, arrow::ipc::RecordBatchFileReader::Open(infile));

  std::shared_ptr<arrow::RecordBatch> rbatch;
  ARROW_ASSIGN_OR_RAISE(rbatch, ipc_reader->ReadRecordBatch(0));

  std::shared_ptr<arrow::io::FileOutputStream> outfile;
  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_out.arrow"));
  ARROW_ASSIGN_OR_RAISE(auto ipc_writer,
                        arrow::ipc::MakeFileWriter(outfile, rbatch->schema()));
  ARROW_RETURN_NOT_OK(ipc_writer->WriteRecordBatch(*rbatch));
  ARROW_RETURN_NOT_OK(ipc_writer->Close());

  // --- Parquet read/write ---
  ARROW_ASSIGN_OR_RAISE(infile, arrow::io::ReadableFile::Open("test_in.parquet"));
  std::unique_ptr<parquet::arrow::FileReader> reader;
  PARQUET_ASSIGN_OR_THROW(reader,
                          parquet::arrow::OpenFile(infile, arrow::default_memory_pool()));

  std::shared_ptr<arrow::Table> parquet_table;
  PARQUET_THROW_NOT_OK(reader->ReadTable(&parquet_table));

  ARROW_ASSIGN_OR_RAISE(outfile, arrow::io::FileOutputStream::Open("test_out.parquet"));
  PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(
      *parquet_table, arrow::default_memory_pool(), outfile, 5));

  return ::arrow::Status::OK();
}

START_TEST(Arrow, "$Id$")

START_SECTION(GenInitialFile())
{
  TEST_EQUAL(GenInitialFile().ok(), true)
}
END_SECTION

START_SECTION(RunMain())
{
  TEST_EQUAL(RunMain().ok(), true)
}
END_SECTION

END_TEST
