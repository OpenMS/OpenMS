// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ParquetFile.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#ifdef WITH_PARQUET
#include <parquet/file_reader.h>
#include <filesystem>
#include <fstream>
#include <vector>
#include <limits>
#if __has_include(<zip.h>)
#include <zip.h>
#define OPENMS_HAVE_LIBZIP 1
#endif
#endif

namespace OpenMS
{

#ifdef WITH_PARQUET

  // ---- Arrow builder helpers ------------------------------------------------

  void ParquetFile::appendOrThrow(const arrow::Status& status, const std::string& column)
  {
    if (!status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to append value for " + column, status.ToString());
    }
  }

  void ParquetFile::throw_finish_error_(const std::string& name, const std::string& error)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Failed to finish array for " + name, error);
  }

  // ---- Parquet file I/O -----------------------------------------------------

  void ParquetFile::writeTable(const std::shared_ptr<arrow::Table>& table, const String& filename)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
    if (!outfile_result.ok())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    auto outfile = outfile_result.ValueOrDie();
    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024);
    if (!status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to write parquet table", status.ToString());
    }
  }

  std::shared_ptr<arrow::Table> ParquetFile::readTable(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    if (!infile_result.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to open parquet file", filename);
    }
    std::shared_ptr<arrow::io::ReadableFile> infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to create parquet reader", filename);
    }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to read parquet table", filename);
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to combine parquet chunks", filename);
    }

    return *combined;
  }

  // ---- Column accessors -----------------------------------------------------

  std::shared_ptr<arrow::Array> ParquetFile::getColumn(const std::shared_ptr<arrow::Table>& table,
                                                       const std::string& name)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing required column '" + name + "'");
    }
    if (column->num_chunks() == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Column has no chunks", name);
    }
    return column->chunk(0);
  }

  std::shared_ptr<arrow::Array> ParquetFile::getOptionalColumn(const std::shared_ptr<arrow::Table>& table,
                                                               const std::string& name)
  {
    auto column = table->GetColumnByName(name);
    if (!column)
    {
      return nullptr;
    }
    if (column->num_chunks() == 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Optional column has no chunks", name);
    }
    return column->chunk(0);
  }

  // ---- Type-safe value accessors --------------------------------------------

  int64_t ParquetFile::getInt64(const std::shared_ptr<arrow::Array>& array, int64_t row,
                                int64_t default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Required integer value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::INT64:
        return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
      case arrow::Type::INT32:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int32Array>(array)->Value(row));
      case arrow::Type::INT16:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int16Array>(array)->Value(row));
      case arrow::Type::INT8:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::Int8Array>(array)->Value(row));
      case arrow::Type::UINT64:
      {
        // Guard narrowing from uint64_t to int64_t to avoid silent wrap-around
        const uint64_t v = std::static_pointer_cast<arrow::UInt64Array>(array)->Value(row);
        if (v > static_cast<uint64_t>(std::numeric_limits<int64_t>::max()))
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Unsigned integer value too large to fit into int64_t", String(std::to_string(v)));
        }
        return static_cast<int64_t>(v);
      }
      case arrow::Type::UINT32:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt32Array>(array)->Value(row));
      case arrow::Type::UINT16:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt16Array>(array)->Value(row));
      case arrow::Type::UINT8:
        return static_cast<int64_t>(std::static_pointer_cast<arrow::UInt8Array>(array)->Value(row));
      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported integer column type", "");
    }
  }

  double ParquetFile::getDouble(const std::shared_ptr<arrow::Array>& array, int64_t row,
                                double default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Required numeric value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::DOUBLE:
        return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
      case arrow::Type::FLOAT:
        return static_cast<double>(std::static_pointer_cast<arrow::FloatArray>(array)->Value(row));
      case arrow::Type::INT64:
      case arrow::Type::INT32:
      case arrow::Type::INT16:
      case arrow::Type::INT8:
      case arrow::Type::UINT64:
      case arrow::Type::UINT32:
      case arrow::Type::UINT16:
      case arrow::Type::UINT8:
        return static_cast<double>(getInt64(array, row, 0, false));
      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported numeric column type", "");
    }
  }

  bool ParquetFile::getBool(const std::shared_ptr<arrow::Array>& array, int64_t row,
                            bool default_value, bool allow_null)
  {
    if (!array)
    {
      return default_value;
    }
    if (array->IsNull(row))
    {
      if (allow_null) return default_value;
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Required boolean value is null");
    }
    switch (array->type_id())
    {
      case arrow::Type::BOOL:
        return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
      case arrow::Type::INT8:
      case arrow::Type::INT16:
      case arrow::Type::INT32:
      case arrow::Type::INT64:
      case arrow::Type::UINT8:
      case arrow::Type::UINT16:
      case arrow::Type::UINT32:
      case arrow::Type::UINT64:
        return getInt64(array, row, 0, false) != 0;
      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported boolean column type", "");
    }
  }

  std::string ParquetFile::getString(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    if (!array)
    {
      return "";
    }
    if (array->IsNull(row)) return "";
    switch (array->type_id())
    {
      case arrow::Type::STRING:
        return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
      case arrow::Type::LARGE_STRING:
        return std::static_pointer_cast<arrow::LargeStringArray>(array)->GetString(row);
      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Unsupported string column type", "");
    }
  }

  std::vector<std::string> ParquetFile::getStringList(const std::shared_ptr<arrow::Array>& array, int64_t row)
  {
    std::vector<std::string> values;
    if (!array) return values;
    if (array->IsNull(row)) return values;

    // semicolon-delimited string
    if (array->type_id() == arrow::Type::STRING || array->type_id() == arrow::Type::LARGE_STRING)
    {
      String raw = getString(array, row);
      if (!raw.empty())
      {
        std::vector<String> parts;
        raw.split(';', parts);
        values.reserve(parts.size());
        for (const auto& part : parts)
        {
          values.push_back(part);
        }
      }
      return values;
    }

    // list of strings
    if (array->type_id() == arrow::Type::LIST)
    {
      auto list_array = std::static_pointer_cast<arrow::ListArray>(array);
      auto values_array = list_array->values();
      auto start = list_array->value_offset(row);
      auto end = list_array->value_offset(row + 1);
      for (int64_t i = start; i < end; ++i)
      {
        values.push_back(getString(values_array, i));
      }
      return values;
    }

    if (array->type_id() == arrow::Type::LARGE_LIST)
    {
      auto list_array = std::static_pointer_cast<arrow::LargeListArray>(array);
      auto values_array = list_array->values();
      auto start = list_array->value_offset(row);
      auto end = list_array->value_offset(row + 1);
      for (int64_t i = start; i < end; ++i)
      {
        values.push_back(getString(values_array, i));
      }
      return values;
    }

    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Unsupported list column type", "");
  }

#endif // WITH_PARQUET

  // ---- Parquet archive utilities --------------------------------------------

  void ParquetFile::zipDirectory(const String& directory_path, const String& output_zip)
  {
#ifdef WITH_PARQUET
  const std::filesystem::path dirpath = std::filesystem::u8path(std::string(directory_path));
  const std::filesystem::path outpath = std::filesystem::u8path(std::string(output_zip));
    const String output_zip_abs = File::absolutePath(output_zip);
    if (File::exists(output_zip_abs))
    {
      File::remove(output_zip_abs);
    }

#if defined(OPENMS_HAVE_LIBZIP)
    int err = 0;
    zip_t* za = zip_open(outpath.string().c_str(), ZIP_CREATE | ZIP_TRUNCATE, &err);
    if (!za)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to create zip archive", output_zip_abs);
    }

    for (auto it = std::filesystem::recursive_directory_iterator(dirpath); it != std::filesystem::recursive_directory_iterator(); ++it)
    {
      if (it->is_directory()) continue;
      const auto full = it->path();
      // compute relative path inside the zip
      std::string rel = std::filesystem::relative(full, dirpath).generic_string();
      zip_source_t* zs = zip_source_file(za, full.string().c_str(), 0, 0);
      if (!zs)
      {
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to add file to zip", full.string());
      }
      zip_int64_t idx = zip_file_add(za, rel.c_str(), zs, ZIP_FL_ENC_UTF_8);
      if (idx < 0)
      {
        zip_source_free(zs);
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to add file to zip", full.string());
      }
      // Create an uncompressed archive
      if (zip_set_file_compression(za, static_cast<zip_uint64_t>(idx), ZIP_CM_STORE, 0) < 0)
      {
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to set stored (no compression) for file in zip", full.string());
      }
    }

    if (zip_close(za) < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to finalize zip archive", output_zip_abs);
    }
#else
    // libzip not available -> zipDirectory not implemented
    (void)directory_path;
    (void)output_zip;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
#else
    (void)directory_path;
    (void)output_zip;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }

  String ParquetFile::unzipDirectory(const String& input_path,
                                     std::unique_ptr<File::TempDir>& temp_dir)
  {
#ifdef WITH_PARQUET
    if (File::isDirectory(input_path))
    {
      return input_path;
    }

    if (!File::readable(input_path))
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, input_path);
    }

    temp_dir = std::make_unique<File::TempDir>();
    const String unpack_dir = temp_dir->getPath() + "/parquet_unpacked";
    File::makeDir(unpack_dir);

#if defined(OPENMS_HAVE_LIBZIP)
    int err = 0;
    zip_t* za = zip_open(input_path.c_str(), ZIP_RDONLY, &err);
    if (!za)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to open zip archive", input_path);
    }

    zip_int64_t num = zip_get_num_entries(za, 0);
    const std::filesystem::path base_path = std::filesystem::u8path(std::string(unpack_dir)).lexically_normal();
    for (zip_uint64_t i = 0; i < static_cast<zip_uint64_t>(num); ++i)
    {
      const char* name = zip_get_name(za, i, 0);
      if (!name) continue;
      std::string entry_name(name);

      // Protect against path traversal: do not allow absolute paths or '..' to escape the base unpack dir.
      std::filesystem::path entry_path = std::filesystem::u8path(entry_name);
      if (entry_path.is_absolute())
      {
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Zip entry has absolute path", entry_name);
      }

      std::filesystem::path outpath = (base_path / entry_path).lexically_normal();
      const std::string base_str = base_path.string();
      const std::string out_str = outpath.string();
      if (out_str.size() < base_str.size() || out_str.compare(0, base_str.size(), base_str) != 0)
      {
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Zip entry would extract outside target directory", entry_name);
      }

      // create parent dirs or directory entries
      if (!entry_name.empty() && entry_name.back() == '/')
      {
        std::filesystem::create_directories(outpath);
        continue;
      }
      if (outpath.has_parent_path()) std::filesystem::create_directories(outpath.parent_path());

      zip_file_t* zf = zip_fopen_index(za, i, 0);
      if (!zf)
      {
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read zip entry", entry_name);
      }

      std::ofstream ofs(outpath.string(), std::ios::binary);
      if (!ofs.is_open())
      {
        zip_fclose(zf);
        zip_close(za);
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
      }

      constexpr size_t BUF_SIZE = 1 << 16;
      std::vector<char> buffer(BUF_SIZE);
      zip_int64_t nread = 0;
      while ((nread = zip_fread(zf, buffer.data(), buffer.size())) > 0)
      {
        ofs.write(buffer.data(), static_cast<std::streamsize>(nread));
        if (ofs.fail())
        {
          zip_fclose(zf);
          zip_close(za);
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, outpath.string());
        }
      }
      if (nread < 0)
      {
        zip_fclose(zf);
        zip_close(za);
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to read zip entry", entry_name);
      }

      ofs.close();
      zip_fclose(zf);
    }

    zip_close(za);

    return unpack_dir;
#else
    // libzip not available -> unzipDirectory not implemented
    (void)input_path;
    (void)temp_dir;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
#else
    (void)input_path;
    (void)temp_dir;
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
  }

#ifdef WITH_PARQUET
  std::string ParquetFile::jsonEscape(const String& input)
  {
    std::string out;
    out.reserve(input.size());
    for (Size i = 0; i < input.size(); ++i)
    {
      const char c = input[i];
      switch (c)
      {
        case '\\': out += "\\\\"; break;
        case '"': out += "\\\""; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
          if (static_cast<unsigned char>(c) < 0x20)
          {
            const char hex[] = "0123456789abcdef";
            out += "\\u00";
            out += hex[(c >> 4) & 0x0f];
            out += hex[c & 0x0f];
          }
          else
          {
            out += c;
          }
          break;
      }
    }
    return out;
  }

  int64_t ParquetFile::rowCount(const String& filename)
  {
    if (!File::exists(filename)) return 0;
    std::unique_ptr<parquet::ParquetFileReader> reader = parquet::ParquetFileReader::OpenFile(std::string(filename), false);
    return reader->metadata()->num_rows();
  }
#endif

} // namespace OpenMS
