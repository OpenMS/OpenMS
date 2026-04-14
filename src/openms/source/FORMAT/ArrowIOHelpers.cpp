// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <OpenMS/CONCEPT/LogStream.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <arrow/table.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

namespace OpenMS
{
namespace ArrowIOHelpers
{

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

} // namespace ArrowIOHelpers
} // namespace OpenMS
