// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/table.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>

namespace OpenMS
{
namespace ArrowIOHelpers
{

String generateUuidV4()
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint32_t> dist;
  uint8_t bytes[16];
  for (int i = 0; i < 4; ++i)
  {
    uint32_t r = dist(gen);
    std::memcpy(bytes + i * 4, &r, 4);
  }
  bytes[6] = (bytes[6] & 0x0F) | 0x40; // version 4
  bytes[8] = (bytes[8] & 0x3F) | 0x80; // variant 1
  char buf[37];
  std::snprintf(buf, sizeof(buf),
    "%02x%02x%02x%02x-%02x%02x-%02x%02x-%02x%02x-%02x%02x%02x%02x%02x%02x",
    bytes[0], bytes[1], bytes[2], bytes[3],
    bytes[4], bytes[5], bytes[6], bytes[7],
    bytes[8], bytes[9], bytes[10], bytes[11],
    bytes[12], bytes[13], bytes[14], bytes[15]);
  return String(buf);
}

namespace
{
  arrow::Compression::type toArrowCompression(ParquetWriteConfig::Compression c)
  {
    switch (c)
    {
      case ParquetWriteConfig::Compression::NONE:   return arrow::Compression::UNCOMPRESSED;
      case ParquetWriteConfig::Compression::SNAPPY: return arrow::Compression::SNAPPY;
      case ParquetWriteConfig::Compression::GZIP:   return arrow::Compression::GZIP;
      case ParquetWriteConfig::Compression::LZ4:    return arrow::Compression::LZ4;
      case ParquetWriteConfig::Compression::ZSTD:   return arrow::Compression::ZSTD;
    }
    return arrow::Compression::ZSTD;
  }
}

bool writeTableToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const String& filename,
  const ParquetWriteConfig& config)
{
  if (!table)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: null table passed to writeTableToParquet (" << filename << ")" << std::endl;
    return false;
  }

  auto file_result = arrow::io::FileOutputStream::Open(filename);
  if (!file_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to open " << filename
                     << " for writing: " << file_result.status().ToString() << std::endl;
    return false;
  }
  const auto& outfile = file_result.ValueOrDie();

  auto builder = parquet::WriterProperties::Builder();
  builder.compression(toArrowCompression(config.compression));
  if (config.compression == ParquetWriteConfig::Compression::ZSTD ||
      config.compression == ParquetWriteConfig::Compression::GZIP)
  {
    builder.compression_level(config.compression_level);
  }
  builder.data_pagesize(config.data_page_size);
  if (config.write_statistics) { builder.enable_statistics(); }
  else                         { builder.disable_statistics(); }

  auto writer_properties = builder.build();
  auto arrow_properties = parquet::ArrowWriterProperties::Builder().store_schema()->build();

  auto status = parquet::arrow::WriteTable(
    *table,
    arrow::default_memory_pool(),
    outfile,
    config.row_group_size,
    writer_properties,
    arrow_properties);

  if (!status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to write " << filename
                     << ": " << status.ToString() << std::endl;
    return false;
  }

  auto close_status = outfile->Close();
  if (!close_status.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to close " << filename
                     << ": " << close_status.ToString() << std::endl;
    return false;
  }

  return true;
}

bool concatenateAndWriteToParquet(
  const std::vector<std::shared_ptr<arrow::Table>>& tables,
  const String& filename,
  const ParquetWriteConfig& config)
{
  if (tables.empty()) { return true; }

  auto concat_result = arrow::ConcatenateTables(tables);
  if (!concat_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to concatenate tables for " << filename
                     << ": " << concat_result.status().ToString() << std::endl;
    return false;
  }

  return writeTableToParquet(concat_result.ValueOrDie(), filename, config);
}

// ---------------------------------------------------------------------------
// Read helpers
// ---------------------------------------------------------------------------

std::shared_ptr<arrow::Array> getColumn(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& name,
  bool required)
{
  auto column = table->GetColumnByName(name);
  if (!column)
  {
    if (required)
    {
      OPENMS_LOG_ERROR << "ArrowIOHelpers: Missing required column '" << name << "'" << std::endl;
    }
    return nullptr;
  }
  if (column->num_chunks() == 0)
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Column '" << name << "' has no chunks" << std::endl;
    return nullptr;
  }
  if (column->num_chunks() == 1)
  {
    return column->chunk(0);
  }
  auto combined = arrow::Concatenate(column->chunks(), arrow::default_memory_pool());
  if (!combined.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to combine chunks for column '" << name << "'" << std::endl;
    return nullptr;
  }
  return *combined;
}

String getStringValue(const std::shared_ptr<arrow::Array>& array, int64_t row)
{
  if (!array || array->IsNull(row)) return "";
  return std::static_pointer_cast<arrow::StringArray>(array)->GetString(row);
}

double getDoubleValue(const std::shared_ptr<arrow::Array>& array, int64_t row, double default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::DoubleArray>(array)->Value(row);
}

float getFloatValue(const std::shared_ptr<arrow::Array>& array, int64_t row, float default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::FloatArray>(array)->Value(row);
}

int32_t getInt32Value(const std::shared_ptr<arrow::Array>& array, int64_t row, int32_t default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::Int32Array>(array)->Value(row);
}

int64_t getInt64Value(const std::shared_ptr<arrow::Array>& array, int64_t row, int64_t default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::Int64Array>(array)->Value(row);
}

bool getBoolValue(const std::shared_ptr<arrow::Array>& array, int64_t row, bool default_val)
{
  if (!array || array->IsNull(row)) return default_val;
  return std::static_pointer_cast<arrow::BooleanArray>(array)->Value(row);
}

bool isNull(const std::shared_ptr<arrow::Array>& array, int64_t row)
{
  return !array || array->IsNull(row);
}

void readMetaValues(
  const std::shared_ptr<arrow::Array>& array,
  int64_t row,
  MetaInfoInterface& target,
  const std::unordered_set<std::string>& excluded_keys)
{
  if (!array || array->IsNull(row)) return;
  auto list_arr = std::static_pointer_cast<arrow::ListArray>(array);
  auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->value_slice(row));
  if (!struct_arr || struct_arr->length() == 0) return;

  auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
  auto value_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));
  auto type_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(2));

  for (int64_t i = 0; i < struct_arr->length(); ++i)
  {
    std::string name = name_arr->GetString(i);
    if (excluded_keys.count(name)) continue;

    std::string value_str = value_arr->GetString(i);
    std::string type_str = type_arr->GetString(i);

    if (type_str == "int")
    {
      try { target.setMetaValue(name, static_cast<int>(std::stol(value_str))); }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "double" || type_str == "float")
    {
      try { target.setMetaValue(name, std::stod(value_str)); }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "int_list")
    {
      try
      {
        String s(value_str);
        if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<Int>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "double_list")
    {
      try
      {
        String s(value_str);
        if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<double>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "string_list")
    {
      try
      {
        String s(value_str);
        if (s.hasPrefix("[") && s.hasSuffix("]")) { s = s.substr(1, s.size() - 2); }
        auto sl = ListUtils::create<String>(s);
        for (auto& e : sl) { e = e.trim(); }
        target.setMetaValue(name, DataValue(sl));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else
    {
      target.setMetaValue(name, value_str);
    }
  }
}

} // namespace ArrowIOHelpers
} // namespace OpenMS
