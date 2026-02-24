// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// minizip-ng headers
#include <mz.h>
#include <mz_strm.h>
#include <mz_strm_os.h>
#include <mz_zip.h>
#include <mz_zip_rw.h>

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/io/interfaces.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <arrow/buffer.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>
#endif

#include <string>
#include <vector>

#ifdef WITH_PARQUET
/// An arrow::io::RandomAccessFile backed by a region of an uncompressed ZIP entry.
/// For STORE entries, byte N of the logical file = byte (data_offset_ + N) in the
/// physical archive file.  This lets Arrow/Parquet read directly from the archive
/// without extracting.
class ZipEntryRandomAccessFile : public arrow::io::RandomAccessFile
{
public:
  ZipEntryRandomAccessFile(std::shared_ptr<arrow::io::RandomAccessFile> archive_file,
                           int64_t data_offset, int64_t entry_size)
    : archive_file_(std::move(archive_file)),
      data_offset_(data_offset),
      entry_size_(entry_size),
      position_(0),
      closed_(false) {}

  arrow::Status Close() override { closed_ = true; return arrow::Status::OK(); }
  bool closed() const override { return closed_; }
  arrow::Result<int64_t> Tell() const override { return position_; }
  arrow::Result<int64_t> GetSize() override { return entry_size_; }

  arrow::Status Seek(int64_t position) override
  {
    if (position < 0 || position > entry_size_)
      return arrow::Status::IOError("Seek out of bounds");
    position_ = position;
    return arrow::Status::OK();
  }

  arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
  {
    int64_t bytes_available = entry_size_ - position_;
    int64_t to_read = std::min(nbytes, bytes_available);
    if (to_read <= 0) return 0;
    ARROW_RETURN_NOT_OK(archive_file_->Seek(data_offset_ + position_));
    ARROW_ASSIGN_OR_RAISE(int64_t bytes_read, archive_file_->Read(to_read, out));
    position_ += bytes_read;
    return bytes_read;
  }

  arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
  {
    int64_t bytes_available = entry_size_ - position_;
    int64_t to_read = std::min(nbytes, bytes_available);
    if (to_read <= 0) return std::make_shared<arrow::Buffer>("");
    ARROW_ASSIGN_OR_RAISE(auto buf, arrow::AllocateResizableBuffer(to_read));
    ARROW_ASSIGN_OR_RAISE(int64_t bytes_read, Read(to_read, buf->mutable_data()));
    ARROW_RETURN_NOT_OK(buf->Resize(bytes_read));
    return buf;
  }

  arrow::Result<int64_t> ReadAt(int64_t position, int64_t nbytes, void* out) override
  {
    ARROW_RETURN_NOT_OK(Seek(position));
    return Read(nbytes, out);
  }

  arrow::Result<std::shared_ptr<arrow::Buffer>> ReadAt(int64_t position, int64_t nbytes) override
  {
    ARROW_RETURN_NOT_OK(Seek(position));
    return Read(nbytes);
  }

private:
  std::shared_ptr<arrow::io::RandomAccessFile> archive_file_;
  int64_t data_offset_;
  int64_t entry_size_;
  int64_t position_;
  bool closed_;
};

/// Helper: get data offset and size for a named STORE entry in a ZIP archive.
/// Uses the low-level minizip-ng API to parse the local header and determine
/// where the raw data begins in the file.
static std::pair<int64_t, int64_t> get_zip_entry_offset_and_size(
  const std::string& zip_path, const std::string& entry_name)
{
  void* file_stream = mz_stream_os_create();
  mz_stream_open(file_stream, zip_path.c_str(), MZ_OPEN_MODE_READ);
  void* zip = mz_zip_create();
  mz_zip_open(zip, file_stream, MZ_OPEN_MODE_READ);

  mz_zip_locate_entry(zip, entry_name.c_str(), 0);
  mz_zip_entry_read_open(zip, 0, nullptr);

  // Get archive stream and record data start position
  void* archive_stream = nullptr;
  mz_zip_get_stream(zip, &archive_stream);
  int64_t data_offset = mz_stream_tell(archive_stream);

  // Get entry size
  mz_zip_file* info = nullptr;
  mz_zip_entry_get_info(zip, &info);
  int64_t entry_size = info->uncompressed_size;

  mz_zip_entry_close(zip);
  mz_zip_close(zip);
  mz_zip_delete(&zip);
  mz_stream_close(file_stream);
  mz_stream_delete(&file_stream);

  return {data_offset, entry_size};
}
#endif // WITH_PARQUET

using namespace OpenMS;

START_TEST(MinizipNG, "$Id$")

/////////////////////////////////////////////////////////////
// Test: Create an uncompressed ZIP64 archive with multiple
// files, then random-access each file by name and verify content.
/////////////////////////////////////////////////////////////

// Test data: three "files" with known content
const std::string file1_name = "spectrum_001.mzML";
const std::string file1_data = "<?xml version=\"1.0\"?><mzML>spectrum data 1 with some padding to make it non-trivial</mzML>";

const std::string file2_name = "spectrum_002.mzML";
const std::string file2_data = "<?xml version=\"1.0\"?><mzML>spectrum data 2 - different content here for testing</mzML>";

const std::string file3_name = "subdir/metadata.json";
const std::string file3_data = R"({"experiment":"test","num_spectra":2,"version":"1.0"})";

std::string zip_path;

START_SECTION(create_uncompressed_zip64_archive)
{
  NEW_TMP_FILE(zip_path)
  zip_path += ".zip";

  // Create writer
  void* writer = mz_zip_writer_create();
  TEST_NOT_EQUAL(writer, nullptr)

  // Configure: uncompressed (STORE) + force ZIP64
  mz_zip_writer_set_compress_method(writer, MZ_COMPRESS_METHOD_STORE);
  mz_zip_writer_set_compress_level(writer, 0);

  // Open archive file
  int32_t err = mz_zip_writer_open_file(writer, zip_path.c_str(), 0 /* disk_size */, 0 /* append */);
  TEST_EQUAL(err, MZ_OK)

  // Add file 1
  {
    mz_zip_file file_info = {};
    file_info.filename = file1_name.c_str();
    file_info.uncompressed_size = static_cast<int64_t>(file1_data.size());
    file_info.compression_method = MZ_COMPRESS_METHOD_STORE;
    file_info.zip64 = MZ_ZIP64_FORCE;

    err = mz_zip_writer_entry_open(writer, &file_info);
    TEST_EQUAL(err, MZ_OK)

    int32_t written = mz_zip_writer_entry_write(writer, file1_data.data(), static_cast<int32_t>(file1_data.size()));
    TEST_EQUAL(written, static_cast<int32_t>(file1_data.size()))

    err = mz_zip_writer_entry_close(writer);
    TEST_EQUAL(err, MZ_OK)
  }

  // Add file 2
  {
    mz_zip_file file_info = {};
    file_info.filename = file2_name.c_str();
    file_info.uncompressed_size = static_cast<int64_t>(file2_data.size());
    file_info.compression_method = MZ_COMPRESS_METHOD_STORE;
    file_info.zip64 = MZ_ZIP64_FORCE;

    err = mz_zip_writer_entry_open(writer, &file_info);
    TEST_EQUAL(err, MZ_OK)

    int32_t written = mz_zip_writer_entry_write(writer, file2_data.data(), static_cast<int32_t>(file2_data.size()));
    TEST_EQUAL(written, static_cast<int32_t>(file2_data.size()))

    err = mz_zip_writer_entry_close(writer);
    TEST_EQUAL(err, MZ_OK)
  }

  // Add file 3 (in a subdirectory)
  {
    mz_zip_file file_info = {};
    file_info.filename = file3_name.c_str();
    file_info.uncompressed_size = static_cast<int64_t>(file3_data.size());
    file_info.compression_method = MZ_COMPRESS_METHOD_STORE;
    file_info.zip64 = MZ_ZIP64_FORCE;

    err = mz_zip_writer_entry_open(writer, &file_info);
    TEST_EQUAL(err, MZ_OK)

    int32_t written = mz_zip_writer_entry_write(writer, file3_data.data(), static_cast<int32_t>(file3_data.size()));
    TEST_EQUAL(written, static_cast<int32_t>(file3_data.size()))

    err = mz_zip_writer_entry_close(writer);
    TEST_EQUAL(err, MZ_OK)
  }

  // Close and clean up writer
  err = mz_zip_writer_close(writer);
  TEST_EQUAL(err, MZ_OK)
  mz_zip_writer_delete(&writer);
}
END_SECTION

START_SECTION(random_access_by_filename)
{
  // Open the archive for reading
  void* reader = mz_zip_reader_create();
  TEST_NOT_EQUAL(reader, nullptr)

  int32_t err = mz_zip_reader_open_file(reader, zip_path.c_str());
  TEST_EQUAL(err, MZ_OK)

  // Random access: read file 3 first (out of order)
  {
    err = mz_zip_reader_locate_entry(reader, file3_name.c_str(), 0 /* case sensitive */);
    TEST_EQUAL(err, MZ_OK)

    err = mz_zip_reader_entry_open(reader);
    TEST_EQUAL(err, MZ_OK)

    // Check entry info
    mz_zip_file* file_info = nullptr;
    err = mz_zip_reader_entry_get_info(reader, &file_info);
    TEST_EQUAL(err, MZ_OK)
    TEST_EQUAL(std::string(file_info->filename), file3_name)
    TEST_EQUAL(file_info->compression_method, static_cast<uint16_t>(MZ_COMPRESS_METHOD_STORE))

    // Read content
    std::vector<char> buf(file3_data.size());
    int32_t bytes_read = mz_zip_reader_entry_read(reader, buf.data(), static_cast<int32_t>(buf.size()));
    TEST_EQUAL(bytes_read, static_cast<int32_t>(file3_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file3_data)

    err = mz_zip_reader_entry_close(reader);
    TEST_EQUAL(err, MZ_OK)
  }

  // Random access: read file 1
  {
    err = mz_zip_reader_locate_entry(reader, file1_name.c_str(), 0);
    TEST_EQUAL(err, MZ_OK)

    err = mz_zip_reader_entry_open(reader);
    TEST_EQUAL(err, MZ_OK)

    std::vector<char> buf(file1_data.size());
    int32_t bytes_read = mz_zip_reader_entry_read(reader, buf.data(), static_cast<int32_t>(buf.size()));
    TEST_EQUAL(bytes_read, static_cast<int32_t>(file1_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file1_data)

    err = mz_zip_reader_entry_close(reader);
    TEST_EQUAL(err, MZ_OK)
  }

  // Random access: read file 2
  {
    err = mz_zip_reader_locate_entry(reader, file2_name.c_str(), 0);
    TEST_EQUAL(err, MZ_OK)

    err = mz_zip_reader_entry_open(reader);
    TEST_EQUAL(err, MZ_OK)

    std::vector<char> buf(file2_data.size());
    int32_t bytes_read = mz_zip_reader_entry_read(reader, buf.data(), static_cast<int32_t>(buf.size()));
    TEST_EQUAL(bytes_read, static_cast<int32_t>(file2_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file2_data)

    err = mz_zip_reader_entry_close(reader);
    TEST_EQUAL(err, MZ_OK)
  }

  // Verify: locating a non-existent entry fails
  {
    err = mz_zip_reader_locate_entry(reader, "nonexistent.txt", 0);
    TEST_NOT_EQUAL(err, MZ_OK)
  }

  err = mz_zip_reader_close(reader);
  TEST_EQUAL(err, MZ_OK)
  mz_zip_reader_delete(&reader);
}
END_SECTION

START_SECTION(seek_within_uncompressed_entry)
{
  // For uncompressed (STORE) entries, we can seek within the entry data
  // by working with the raw file stream. After opening an entry and reading
  // the local header, the stream position is at the data start. We record
  // this offset and then seek directly on the file stream.

  // Open archive using low-level API
  void* file_stream = mz_stream_os_create();
  TEST_NOT_EQUAL(file_stream, nullptr)

  int32_t err = mz_stream_open(file_stream, zip_path.c_str(), MZ_OPEN_MODE_READ);
  TEST_EQUAL(err, MZ_OK)

  void* zip = mz_zip_create();
  TEST_NOT_EQUAL(zip, nullptr)

  err = mz_zip_open(zip, file_stream, MZ_OPEN_MODE_READ);
  TEST_EQUAL(err, MZ_OK)

  // Locate file1
  err = mz_zip_locate_entry(zip, file1_name.c_str(), 0);
  TEST_EQUAL(err, MZ_OK)

  // Open entry for reading - this parses the local header and positions
  // the stream at the start of the entry data
  err = mz_zip_entry_read_open(zip, 0 /* raw */, nullptr /* password */);
  TEST_EQUAL(err, MZ_OK)

  // Get the underlying archive stream (positioned at data start)
  void* archive_stream = nullptr;
  err = mz_zip_get_stream(zip, &archive_stream);
  TEST_EQUAL(err, MZ_OK)

  // Record the data start position by reading current stream position
  int64_t data_start = mz_stream_tell(archive_stream);
  TEST_NOT_EQUAL(data_start, -1)

  // Close the entry (we'll read directly from the stream)
  mz_zip_entry_close(zip);

  // Now seek and read directly from the file stream.
  // For STORE entries, the data is stored verbatim in the archive,
  // so we can read it directly at (data_start + offset).

  // Read 10 bytes from offset 0
  std::vector<char> buf(10);
  err = mz_stream_seek(archive_stream, data_start + 0, MZ_SEEK_SET);
  TEST_EQUAL(err, MZ_OK)
  int32_t bytes_read = mz_stream_read(archive_stream, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(0, 10))

  // Seek to byte offset 20 from data start
  err = mz_stream_seek(archive_stream, data_start + 20, MZ_SEEK_SET);
  TEST_EQUAL(err, MZ_OK)
  bytes_read = mz_stream_read(archive_stream, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(20, 10))

  // Seek backwards to offset 5
  err = mz_stream_seek(archive_stream, data_start + 5, MZ_SEEK_SET);
  TEST_EQUAL(err, MZ_OK)
  bytes_read = mz_stream_read(archive_stream, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(5, 10))

  // Seek to near the end and read last 5 bytes
  int64_t near_end = static_cast<int64_t>(file1_data.size()) - 5;
  err = mz_stream_seek(archive_stream, data_start + near_end, MZ_SEEK_SET);
  TEST_EQUAL(err, MZ_OK)
  std::vector<char> tail_buf(5);
  bytes_read = mz_stream_read(archive_stream, tail_buf.data(), 5);
  TEST_EQUAL(bytes_read, 5)
  TEST_EQUAL(std::string(tail_buf.data(), 5), file1_data.substr(file1_data.size() - 5, 5))

  // Also test seeking into file2 (different entry) directly
  // First locate file2 and get its data offset
  err = mz_zip_locate_entry(zip, file2_name.c_str(), 0);
  TEST_EQUAL(err, MZ_OK)
  err = mz_zip_entry_read_open(zip, 0, nullptr);
  TEST_EQUAL(err, MZ_OK)
  int64_t data2_start = mz_stream_tell(archive_stream);
  TEST_NOT_EQUAL(data2_start, -1)
  mz_zip_entry_close(zip);

  // Read 10 bytes from the middle of file2
  int64_t mid2 = static_cast<int64_t>(file2_data.size()) / 2;
  err = mz_stream_seek(archive_stream, data2_start + mid2, MZ_SEEK_SET);
  TEST_EQUAL(err, MZ_OK)
  bytes_read = mz_stream_read(archive_stream, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file2_data.substr(mid2, 10))

  // Clean up
  mz_zip_close(zip);
  mz_zip_delete(&zip);
  mz_stream_close(file_stream);
  mz_stream_delete(&file_stream);
}
END_SECTION

START_SECTION(verify_entries_are_uncompressed_zip64)
{
  // Re-open and iterate all entries to verify STORE method and ZIP64 flag
  void* reader = mz_zip_reader_create();
  int32_t err = mz_zip_reader_open_file(reader, zip_path.c_str());
  TEST_EQUAL(err, MZ_OK)

  int entry_count = 0;
  err = mz_zip_reader_goto_first_entry(reader);
  while (err == MZ_OK)
  {
    mz_zip_file* file_info = nullptr;
    err = mz_zip_reader_entry_get_info(reader, &file_info);
    TEST_EQUAL(err, MZ_OK)

    // Verify uncompressed
    TEST_EQUAL(file_info->compression_method, static_cast<uint16_t>(MZ_COMPRESS_METHOD_STORE))

    // Verify ZIP64 (zip64 field should be non-zero when ZIP64 extensions are used)
    TEST_NOT_EQUAL(file_info->zip64, static_cast<uint16_t>(MZ_ZIP64_DISABLE))

    entry_count++;
    err = mz_zip_reader_goto_next_entry(reader);
  }

  TEST_EQUAL(entry_count, 3)

  mz_zip_reader_close(reader);
  mz_zip_reader_delete(&reader);
}
END_SECTION

#ifdef WITH_PARQUET
START_SECTION(parquet_roundtrip_through_zip_archive)
{
  // End-to-end test: write two Parquet tables to a ZIP archive (uncompressed),
  // then read each one back directly from the archive using Arrow's Parquet
  // reader with a ZipEntryRandomAccessFile adapter — no extraction needed.

  std::string parquet_zip_path;
  NEW_TMP_FILE(parquet_zip_path)
  parquet_zip_path += ".zip";

  // --- Build two small Arrow tables ---
  auto schema1 = arrow::schema({
    arrow::field("mz", arrow::float64()),
    arrow::field("intensity", arrow::float32())
  });

  arrow::DoubleBuilder mz_builder;
  arrow::FloatBuilder int_builder;
  TEST_TRUE(mz_builder.AppendValues({100.5, 200.3, 300.1, 400.7, 500.9}).ok());
  TEST_TRUE(int_builder.AppendValues({1000.0f, 2000.0f, 3000.0f, 4000.0f, 5000.0f}).ok());
  std::shared_ptr<arrow::Array> mz_arr, int_arr;
  TEST_TRUE(mz_builder.Finish(&mz_arr).ok());
  TEST_TRUE(int_builder.Finish(&int_arr).ok());
  auto table1 = arrow::Table::Make(schema1, {mz_arr, int_arr});

  auto schema2 = arrow::schema({
    arrow::field("protein", arrow::utf8()),
    arrow::field("score", arrow::float64())
  });

  arrow::StringBuilder str_builder;
  arrow::DoubleBuilder score_builder;
  TEST_TRUE(str_builder.AppendValues({"BSA", "Lysozyme", "Myoglobin"}).ok());
  TEST_TRUE(score_builder.AppendValues({99.5, 87.2, 76.1}).ok());
  std::shared_ptr<arrow::Array> str_arr, score_arr;
  TEST_TRUE(str_builder.Finish(&str_arr).ok());
  TEST_TRUE(score_builder.Finish(&score_arr).ok());
  auto table2 = arrow::Table::Make(schema2, {str_arr, score_arr});

  // --- Write each table to a temporary Parquet file, then add to ZIP ---
  std::string pq1_path, pq2_path;
  NEW_TMP_FILE(pq1_path) pq1_path += ".parquet";
  NEW_TMP_FILE(pq2_path) pq2_path += ".parquet";

  {
    auto outfile = arrow::io::FileOutputStream::Open(pq1_path).ValueOrDie();
    PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(*table1, arrow::default_memory_pool(), outfile, 1024));
  }
  {
    auto outfile = arrow::io::FileOutputStream::Open(pq2_path).ValueOrDie();
    PARQUET_THROW_NOT_OK(parquet::arrow::WriteTable(*table2, arrow::default_memory_pool(), outfile, 1024));
  }

  // --- Create ZIP archive with both Parquet files ---
  void* writer = mz_zip_writer_create();
  mz_zip_writer_set_compress_method(writer, MZ_COMPRESS_METHOD_STORE);
  mz_zip_writer_set_compress_level(writer, 0);

  int32_t err = mz_zip_writer_open_file(writer, parquet_zip_path.c_str(), 0, 0);
  TEST_EQUAL(err, MZ_OK)

  // Add both parquet files to archive
  err = mz_zip_writer_add_file(writer, pq1_path.c_str(), "spectra.parquet");
  TEST_EQUAL(err, MZ_OK)
  err = mz_zip_writer_add_file(writer, pq2_path.c_str(), "proteins.parquet");
  TEST_EQUAL(err, MZ_OK)

  err = mz_zip_writer_close(writer);
  TEST_EQUAL(err, MZ_OK)
  mz_zip_writer_delete(&writer);

  // --- Read spectra.parquet directly from the ZIP archive ---
  {
    auto [data_offset, entry_size] = get_zip_entry_offset_and_size(parquet_zip_path, "spectra.parquet");
    TEST_NOT_EQUAL(entry_size, 0)

    auto archive_file = arrow::io::ReadableFile::Open(parquet_zip_path).ValueOrDie();
    auto zip_file = std::make_shared<ZipEntryRandomAccessFile>(archive_file, data_offset, entry_size);

    std::unique_ptr<parquet::arrow::FileReader> pq_reader;
    PARQUET_ASSIGN_OR_THROW(pq_reader,
      parquet::arrow::OpenFile(zip_file, arrow::default_memory_pool()));

    std::shared_ptr<arrow::Table> read_table;
    PARQUET_THROW_NOT_OK(pq_reader->ReadTable(&read_table));

    TEST_EQUAL(read_table->num_rows(), 5)
    TEST_EQUAL(read_table->num_columns(), 2)
    TEST_EQUAL(read_table->schema()->field(0)->name(), "mz")
    TEST_EQUAL(read_table->schema()->field(1)->name(), "intensity")

    // Verify data values
    auto mz_col = std::static_pointer_cast<arrow::DoubleArray>(read_table->column(0)->chunk(0));
    TEST_REAL_SIMILAR(mz_col->Value(0), 100.5)
    TEST_REAL_SIMILAR(mz_col->Value(4), 500.9)
  }

  // --- Read proteins.parquet directly from the same ZIP archive ---
  {
    auto [data_offset, entry_size] = get_zip_entry_offset_and_size(parquet_zip_path, "proteins.parquet");
    TEST_NOT_EQUAL(entry_size, 0)

    auto archive_file = arrow::io::ReadableFile::Open(parquet_zip_path).ValueOrDie();
    auto zip_file = std::make_shared<ZipEntryRandomAccessFile>(archive_file, data_offset, entry_size);

    std::unique_ptr<parquet::arrow::FileReader> pq_reader;
    PARQUET_ASSIGN_OR_THROW(pq_reader,
      parquet::arrow::OpenFile(zip_file, arrow::default_memory_pool()));

    std::shared_ptr<arrow::Table> read_table;
    PARQUET_THROW_NOT_OK(pq_reader->ReadTable(&read_table));

    TEST_EQUAL(read_table->num_rows(), 3)
    TEST_EQUAL(read_table->num_columns(), 2)
    TEST_EQUAL(read_table->schema()->field(0)->name(), "protein")
    TEST_EQUAL(read_table->schema()->field(1)->name(), "score")

    // Verify data values
    auto protein_col = std::static_pointer_cast<arrow::StringArray>(read_table->column(0)->chunk(0));
    TEST_EQUAL(protein_col->GetString(0), "BSA")
    TEST_EQUAL(protein_col->GetString(2), "Myoglobin")

    auto score_col = std::static_pointer_cast<arrow::DoubleArray>(read_table->column(1)->chunk(0));
    TEST_REAL_SIMILAR(score_col->Value(1), 87.2)
  }
}
END_SECTION
#endif // WITH_PARQUET

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
