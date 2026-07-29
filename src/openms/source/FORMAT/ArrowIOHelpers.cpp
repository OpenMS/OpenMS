// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/SpectrumNativeIDParser.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/table.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <cstdint>

namespace OpenMS
{
namespace ArrowIOHelpers
{

std::string generateUuidV4()
{
  return UniqueIdGenerator::getUUID();
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

  /// QPX `compression_format` token for @p c, or "" if QPX does not define one.
  /// The spec's vocabulary is zstd|snappy|gzip|lzo|none — LZ4 has no token.
  std::string qpxCompressionName(ParquetWriteConfig::Compression c)
  {
    switch (c)
    {
      case ParquetWriteConfig::Compression::NONE:   return "none";
      case ParquetWriteConfig::Compression::SNAPPY: return "snappy";
      case ParquetWriteConfig::Compression::GZIP:   return "gzip";
      case ParquetWriteConfig::Compression::ZSTD:   return "zstd";
      case ParquetWriteConfig::Compression::LZ4:    return "";
    }
    return "";
  }
}

std::shared_ptr<const arrow::KeyValueMetadata> qpxFileMetadata(
  const std::string& file_type,
  const ParquetWriteConfig& config,
  const std::map<std::string, std::string>& extra)
{
  const std::string compression = qpxCompressionName(config.compression);
  if (compression.empty())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: QPX does not define a compression_format token for LZ4; "
                        "use ZSTD (default), SNAPPY, GZIP or NONE for QPX output." << std::endl;
    return nullptr;
  }

  std::vector<std::string> keys{
    "qpx_version", "file_type", "creator", "software_provider", "creation_date",
    "compression_format", "uuid"};
  std::vector<std::string> values{
    "1.0",
    file_type,
    "OpenMS",
    "OpenMS " + VersionInfo::getVersion(),
    DateTime::nowUTC().toString("yyyy-MM-ddThh:mm:ssZ"),
    compression,
    generateUuidV4()};

  for (const auto& [k, v] : extra)
  {
    keys.push_back(k);
    values.push_back(v);
  }

  return std::make_shared<arrow::KeyValueMetadata>(std::move(keys), std::move(values));
}

std::string qpxScanFormat(const std::string& native_id)
{
  if (StringUtils::hasPrefix(native_id, "index=")) { return "index"; }
  if (SpectrumNativeIDParser::isNativeID(native_id)) { return "scan"; }
  return "";
}

std::string qpxRunFileName(const std::string& ms_run_path)
{
  if (ms_run_path.empty()) { return ""; }
  return File::stemName(ms_run_path);
}

std::string qpxScanFormat(const std::vector<std::string>& native_ids)
{
  std::string resolved;
  for (const auto& id : native_ids)
  {
    const std::string fmt = qpxScanFormat(id);
    if (fmt.empty()) { continue; } // unrecognized IDs carry no evidence either way
    if (resolved.empty()) { resolved = fmt; }
    else if (resolved != fmt)
    {
      // Two conventions in one export: neither token would describe the whole file, so
      // omit scan_format rather than assert a wrong one.
      OPENMS_LOG_WARN << "ArrowIOHelpers: spectrum native IDs mix '" << resolved << "' and '"
                      << fmt << "' conventions; omitting QPX scan_format." << std::endl;
      return "";
    }
  }
  return resolved;
}

bool writeTableToParquet(
  const std::shared_ptr<arrow::Table>& table,
  const std::string& filename,
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
  const std::string& filename,
  const ParquetWriteConfig& config,
  const std::shared_ptr<const arrow::KeyValueMetadata>& metadata)
{
  if (tables.empty()) { return true; }

  auto concat_result = arrow::ConcatenateTables(tables);
  if (!concat_result.ok())
  {
    OPENMS_LOG_ERROR << "ArrowIOHelpers: Failed to concatenate tables for " << filename
                     << ": " << concat_result.status().ToString() << std::endl;
    return false;
  }

  auto table = concat_result.ValueOrDie();
  // The concatenated table inherits the first input's schema metadata, which carries that
  // input's uuid/creation_date. Stamp one identity for the merged file instead.
  if (metadata) { table = table->ReplaceSchemaMetadata(metadata); }
  return writeTableToParquet(table, filename, config);
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

std::string getStringValue(const std::shared_ptr<arrow::Array>& array, int64_t row)
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
  const int64_t off = list_arr->value_offset(row);
  const int64_t len = list_arr->value_length(row);
  if (len == 0) return;

  // Index the shared child struct via value_offset/length instead of value_slice(row),
  // which allocates a sliced array on every call (once per metavalue list per row).
  auto struct_arr = std::static_pointer_cast<arrow::StructArray>(list_arr->values());
  auto name_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(0));
  auto value_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(1));
  auto type_arr = std::static_pointer_cast<arrow::StringArray>(struct_arr->field(2));

  for (int64_t i = off; i < off + len; ++i)
  {
    std::string name = name_arr->GetString(i);
    if (excluded_keys.contains(name)) continue;

    std::string value_str = value_arr->GetString(i);
    std::string type_str = type_arr->GetString(i);

    if (type_str == "int")
    {
      try { target.setMetaValue(name, static_cast<SignedSize>(std::stoll(value_str))); }
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
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<Int>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "double_list")
    {
      try
      {
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        target.setMetaValue(name, DataValue(ListUtils::create<double>(s)));
      }
      catch (...) { target.setMetaValue(name, value_str); }
    }
    else if (type_str == "string_list")
    {
      try
      {
        std::string s(value_str);
        if (StringUtils::hasPrefix(s, "[") && StringUtils::hasSuffix(s, "]")) { s = s.substr(1, s.size() - 2); }
        auto sl = ListUtils::create<std::string>(s);
        for (auto& e : sl) { StringUtils::trim(e); }
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
