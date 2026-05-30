// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// libzip header
#include <zip.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/io/interfaces.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <arrow/buffer.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

/// An arrow::io::RandomAccessFile backed by a STORE entry inside a ZIP archive.
/// Uses libzip's zip_fseek / zip_fread / zip_ftell for true random access
/// on uncompressed entries (internally mapped to fseek on the physical file).
///
/// Owns both the zip_t archive and the zip_file_t entry handle.
/// Not thread-safe — each instance should be used by a single thread.
/// To read multiple entries from the same archive, open a separate
/// ZipEntryRandomAccessFile per entry (each gets its own zip_t handle
/// since libzip doesn't support concurrent entry reads from a single
/// archive handle).
///
/// libzip handles ZIP64 extensions transparently — no explicit flags needed.
class ZipEntryRandomAccessFile : public arrow::io::RandomAccessFile
{
public:
  ZipEntryRandomAccessFile(zip_t* archive, zip_file_t* entry, int64_t entry_size)
    : archive_(archive),
      entry_(entry),
      entry_size_(entry_size),
      closed_(false) {}

  ~ZipEntryRandomAccessFile() override
  {
    if (!closed_)
    {
      if (entry_) zip_fclose(entry_);
      if (archive_) zip_close(archive_);
    }
  }

  arrow::Status Close() override
  {
    if (!closed_)
    {
      if (entry_) { zip_fclose(entry_); entry_ = nullptr; }
      if (archive_) { zip_close(archive_); archive_ = nullptr; }
      closed_ = true;
    }
    return arrow::Status::OK();
  }

  bool closed() const override { return closed_; }

  arrow::Result<int64_t> Tell() const override
  {
    if (closed_) return arrow::Status::Invalid("File is closed");
    zip_int64_t pos = zip_ftell(entry_);
    if (pos < 0) return arrow::Status::IOError("zip_ftell failed");
    return static_cast<int64_t>(pos);
  }

  arrow::Result<int64_t> GetSize() override { return entry_size_; }

  arrow::Status Seek(int64_t position) override
  {
    if (closed_) return arrow::Status::Invalid("File is closed");
    if (position < 0 || position > entry_size_)
      return arrow::Status::IOError("Seek out of bounds");
    if (zip_fseek(entry_, position, SEEK_SET) < 0)
      return arrow::Status::IOError("zip_fseek failed");
    return arrow::Status::OK();
  }

  arrow::Result<int64_t> Read(int64_t nbytes, void* out) override
  {
    if (closed_) return arrow::Status::Invalid("File is closed");
    zip_int64_t bytes_read = zip_fread(entry_, out, static_cast<zip_uint64_t>(nbytes));
    if (bytes_read < 0) return arrow::Status::IOError("zip_fread failed");
    return static_cast<int64_t>(bytes_read);
  }

  arrow::Result<std::shared_ptr<arrow::Buffer>> Read(int64_t nbytes) override
  {
    ARROW_ASSIGN_OR_RAISE(auto buf, arrow::AllocateResizableBuffer(nbytes));
    ARROW_ASSIGN_OR_RAISE(int64_t bytes_read, Read(nbytes, buf->mutable_data()));
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
  zip_t* archive_;
  zip_file_t* entry_;
  int64_t entry_size_;
  bool closed_;
};

/// Helper: open a ZIP archive and a named STORE entry, returning a
/// ZipEntryRandomAccessFile that Arrow/Parquet can read from directly.
static std::shared_ptr<ZipEntryRandomAccessFile> open_zip_entry_for_arrow(
  const std::string& zip_path, const std::string& entry_name)
{
  int err = 0;
  zip_t* archive = zip_open(zip_path.c_str(), ZIP_RDONLY, &err);
  if (!archive) throw std::runtime_error("Failed to open ZIP archive: " + zip_path);

  zip_stat_t st;
  zip_stat_init(&st);
  if (zip_stat(archive, entry_name.c_str(), 0, &st) != 0)
  {
    zip_close(archive);
    throw std::runtime_error("Entry not found in ZIP: " + entry_name);
  }

  zip_file_t* entry = zip_fopen(archive, entry_name.c_str(), 0);
  if (!entry)
  {
    zip_close(archive);
    throw std::runtime_error("Failed to open ZIP entry: " + entry_name);
  }

  return std::make_shared<ZipEntryRandomAccessFile>(archive, entry, static_cast<int64_t>(st.size));
}

using namespace OpenMS;

START_TEST(Libzip, "$Id$")

/////////////////////////////////////////////////////////////
// Test: Create an uncompressed ZIP archive with multiple
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

START_SECTION(create_uncompressed_zip_archive)
{
  NEW_TMP_FILE(zip_path)
  zip_path += ".zip";

  int err = 0;
  zip_t* archive = zip_open(zip_path.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &err);
  TEST_NOT_EQUAL(archive, nullptr)

  // Add file 1 (from buffer, uncompressed)
  {
    zip_source_t* src = zip_source_buffer(archive, file1_data.data(), file1_data.size(), 0);
    TEST_NOT_EQUAL(src, nullptr)
    zip_int64_t idx = zip_file_add(archive, file1_name.c_str(), src, ZIP_FL_OVERWRITE);
    TEST_NOT_EQUAL(idx, -1)
    TEST_EQUAL(zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE, 0), 0)
  }

  // Add file 2
  {
    zip_source_t* src = zip_source_buffer(archive, file2_data.data(), file2_data.size(), 0);
    TEST_NOT_EQUAL(src, nullptr)
    zip_int64_t idx = zip_file_add(archive, file2_name.c_str(), src, ZIP_FL_OVERWRITE);
    TEST_NOT_EQUAL(idx, -1)
    TEST_EQUAL(zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE, 0), 0)
  }

  // Add file 3 (in a subdirectory)
  {
    zip_source_t* src = zip_source_buffer(archive, file3_data.data(), file3_data.size(), 0);
    TEST_NOT_EQUAL(src, nullptr)
    zip_int64_t idx = zip_file_add(archive, file3_name.c_str(), src, ZIP_FL_OVERWRITE);
    TEST_NOT_EQUAL(idx, -1)
    TEST_EQUAL(zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE, 0), 0)
  }

  // Close writes the archive to disk
  TEST_EQUAL(zip_close(archive), 0)
}
END_SECTION

START_SECTION(random_access_by_filename)
{
  int err = 0;
  zip_t* archive = zip_open(zip_path.c_str(), ZIP_RDONLY, &err);
  TEST_NOT_EQUAL(archive, nullptr)

  // Random access: read file 3 first (out of order)
  {
    zip_stat_t st;
    zip_stat_init(&st);
    err = zip_stat(archive, file3_name.c_str(), 0, &st);
    TEST_EQUAL(err, 0)
    TEST_EQUAL(std::string(st.name), file3_name)
    TEST_EQUAL(st.comp_method, ZIP_CM_STORE)

    zip_file_t* f = zip_fopen(archive, file3_name.c_str(), 0);
    TEST_NOT_EQUAL(f, nullptr)

    std::vector<char> buf(file3_data.size());
    zip_int64_t bytes_read = zip_fread(f, buf.data(), buf.size());
    TEST_EQUAL(bytes_read, static_cast<zip_int64_t>(file3_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file3_data)

    zip_fclose(f);
  }

  // Random access: read file 1
  {
    zip_file_t* f = zip_fopen(archive, file1_name.c_str(), 0);
    TEST_NOT_EQUAL(f, nullptr)

    std::vector<char> buf(file1_data.size());
    zip_int64_t bytes_read = zip_fread(f, buf.data(), buf.size());
    TEST_EQUAL(bytes_read, static_cast<zip_int64_t>(file1_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file1_data)

    zip_fclose(f);
  }

  // Random access: read file 2
  {
    zip_file_t* f = zip_fopen(archive, file2_name.c_str(), 0);
    TEST_NOT_EQUAL(f, nullptr)

    std::vector<char> buf(file2_data.size());
    zip_int64_t bytes_read = zip_fread(f, buf.data(), buf.size());
    TEST_EQUAL(bytes_read, static_cast<zip_int64_t>(file2_data.size()))
    TEST_EQUAL(std::string(buf.data(), bytes_read), file2_data)

    zip_fclose(f);
  }

  // Verify: locating a non-existent entry fails
  {
    zip_stat_t st;
    zip_stat_init(&st);
    err = zip_stat(archive, "nonexistent.txt", 0, &st);
    TEST_NOT_EQUAL(err, 0)
  }

  zip_close(archive);
}
END_SECTION

START_SECTION(seek_within_uncompressed_entry)
{
  // For uncompressed (STORE) entries, zip_fseek performs a real fseek
  // on the underlying file, giving true random access.

  int err = 0;
  zip_t* archive = zip_open(zip_path.c_str(), ZIP_RDONLY, &err);
  TEST_NOT_EQUAL(archive, nullptr)

  // Open file1 entry
  zip_file_t* f = zip_fopen(archive, file1_name.c_str(), 0);
  TEST_NOT_EQUAL(f, nullptr)

  // Seekability is implicitly verified by the zip_fseek/zip_ftell assertions below.
  // (zip_file_is_seekable() would be a direct check but requires libzip >= 1.9.0,
  // while Ubuntu 24.04 ships 1.7.3. zip_fseek/zip_ftell are available since 1.2.0.)

  // Read 10 bytes from offset 0
  std::vector<char> buf(10);
  zip_int64_t bytes_read = zip_fread(f, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(0, 10))

  // Seek to byte offset 20
  TEST_EQUAL(zip_fseek(f, 20, SEEK_SET), 0)
  TEST_EQUAL(zip_ftell(f), 20)
  bytes_read = zip_fread(f, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(20, 10))

  // Seek backwards to offset 5
  TEST_EQUAL(zip_fseek(f, 5, SEEK_SET), 0)
  TEST_EQUAL(zip_ftell(f), 5)
  bytes_read = zip_fread(f, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file1_data.substr(5, 10))

  // Seek to near the end and read last 5 bytes
  int64_t near_end = static_cast<int64_t>(file1_data.size()) - 5;
  TEST_EQUAL(zip_fseek(f, near_end, SEEK_SET), 0)
  std::vector<char> tail_buf(5);
  bytes_read = zip_fread(f, tail_buf.data(), 5);
  TEST_EQUAL(bytes_read, 5)
  TEST_EQUAL(std::string(tail_buf.data(), 5), file1_data.substr(file1_data.size() - 5, 5))

  // Seek from end
  TEST_EQUAL(zip_fseek(f, -5, SEEK_END), 0)
  bytes_read = zip_fread(f, tail_buf.data(), 5);
  TEST_EQUAL(bytes_read, 5)
  TEST_EQUAL(std::string(tail_buf.data(), 5), file1_data.substr(file1_data.size() - 5, 5))

  zip_fclose(f);

  // Also test seeking within file2 (different entry)
  f = zip_fopen(archive, file2_name.c_str(), 0);
  TEST_NOT_EQUAL(f, nullptr)

  int64_t mid2 = static_cast<int64_t>(file2_data.size()) / 2;
  TEST_EQUAL(zip_fseek(f, mid2, SEEK_SET), 0)
  bytes_read = zip_fread(f, buf.data(), 10);
  TEST_EQUAL(bytes_read, 10)
  TEST_EQUAL(std::string(buf.data(), 10), file2_data.substr(mid2, 10))

  zip_fclose(f);
  zip_close(archive);
}
END_SECTION

START_SECTION(verify_entries_are_uncompressed)
{
  // Note: libzip handles ZIP64 extensions transparently for large archives
  // (>4 GiB or >65535 entries). Unlike minizip-ng, no explicit ZIP64 flags
  // are needed, so we only verify the STORE compression method here.
  int err = 0;
  zip_t* archive = zip_open(zip_path.c_str(), ZIP_RDONLY, &err);
  TEST_NOT_EQUAL(archive, nullptr)

  zip_int64_t num_entries = zip_get_num_entries(archive, 0);
  TEST_EQUAL(num_entries, 3)

  for (zip_int64_t i = 0; i < num_entries; ++i)
  {
    zip_stat_t st;
    zip_stat_init(&st);
    err = zip_stat_index(archive, static_cast<zip_uint64_t>(i), 0, &st);
    TEST_EQUAL(err, 0)

    // Verify uncompressed (STORE)
    TEST_EQUAL(st.comp_method, ZIP_CM_STORE)
  }

  zip_close(archive);
}
END_SECTION

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

  // --- Create ZIP archive with both Parquet files (STORE, no compression) ---
  {
    int err = 0;
    zip_t* archive = zip_open(parquet_zip_path.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &err);
    TEST_NOT_EQUAL(archive, nullptr)

    // Add spectra.parquet from file on disk
    zip_source_t* src1 = zip_source_file(archive, pq1_path.c_str(), 0, -1);
    TEST_NOT_EQUAL(src1, nullptr)
    zip_int64_t idx1 = zip_file_add(archive, "spectra.parquet", src1, ZIP_FL_OVERWRITE);
    TEST_NOT_EQUAL(idx1, -1)
    TEST_EQUAL(zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx1), ZIP_CM_STORE, 0), 0)

    // Add proteins.parquet from file on disk
    zip_source_t* src2 = zip_source_file(archive, pq2_path.c_str(), 0, -1);
    TEST_NOT_EQUAL(src2, nullptr)
    zip_int64_t idx2 = zip_file_add(archive, "proteins.parquet", src2, ZIP_FL_OVERWRITE);
    TEST_NOT_EQUAL(idx2, -1)
    TEST_EQUAL(zip_set_file_compression(archive, static_cast<zip_uint64_t>(idx2), ZIP_CM_STORE, 0), 0)

    TEST_EQUAL(zip_close(archive), 0)
  }

  // --- Read spectra.parquet directly from the ZIP archive ---
  {
    auto zip_file = open_zip_entry_for_arrow(parquet_zip_path, "spectra.parquet");

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
    auto zip_file = open_zip_entry_for_arrow(parquet_zip_path, "proteins.parquet");

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
