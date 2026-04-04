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

#include <parquet/file_reader.h>
#include <filesystem>
#include <fstream>
#include <vector>
#include <limits>
#if __has_include(<zip.h>)
#include <zip.h>
#define OPENMS_HAVE_LIBZIP 1
#endif

namespace OpenMS
{

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

  void ParquetFile::writeTable(const std::shared_ptr<arrow::Table>& table, const String& filename,
                               int64_t row_group_size)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
    if (!outfile_result.ok())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    const auto& outfile = outfile_result.ValueOrDie();
    // Use a larger default row_group_size than 1024 to improve compression and reduce metadata overhead.
    // Default is configurable by callers via the row_group_size parameter.
    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, static_cast<int>(row_group_size));
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
    const std::shared_ptr<arrow::io::ReadableFile>& infile = *infile_result;

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

    // Only combine chunks when necessary: if every column already has a single chunk,
    // avoid CombineChunks() which copies data and doubles memory usage for large tables.
    bool need_combine = false;
    for (int i = 0; i < static_cast<int>(table->num_columns()); ++i)
    {
      const auto& col = table->column(i);
      if (col->num_chunks() > 1)
      {
        need_combine = true;
        break;
      }
    }
    if (!need_combine)
    {
      return table;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to combine parquet chunks", filename);
    }

    return *combined;
  }

  std::shared_ptr<arrow::Table> ParquetFile::readTable(const std::shared_ptr<arrow::io::RandomAccessFile>& infile)
  {
    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    if (!reader_result.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to create parquet reader from RandomAccessFile", "");
    }
    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());

    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    if (!read_status.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to read parquet table from RandomAccessFile", "");
    }

    bool need_combine = false;
    for (int i = 0; i < static_cast<int>(table->num_columns()); ++i)
    {
      const auto& col = table->column(i);
      if (col->num_chunks() > 1)
      {
        need_combine = true;
        break;
      }
    }
    if (!need_combine)
    {
      return table;
    }

    auto combined = table->CombineChunks(arrow::default_memory_pool());
    if (!combined.ok())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Failed to combine parquet chunks from RandomAccessFile", "");
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

  // ---- Parquet archive utilities --------------------------------------------


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

} // namespace OpenMS
