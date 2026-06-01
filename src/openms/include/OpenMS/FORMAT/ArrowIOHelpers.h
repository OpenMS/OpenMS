// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>  // for ParquetWriteConfig

#include <cstdint>
#include <memory>
#include <unordered_set>
#include <vector>

// Forward declarations
namespace arrow
{
  class Array;
  class Table;
}

namespace OpenMS
{
class MetaInfoInterface;

/**
  @brief Public helpers for writing and concatenating Arrow tables to Parquet files

  TOPP tools link against libOpenMS (which exports these helpers) but not directly
  against Arrow/Parquet. These wrappers keep all Arrow/Parquet API calls inside
  libOpenMS so downstream binaries don't need to import Arrow symbols.

  @ingroup FileIO
*/
namespace ArrowIOHelpers
{
  /**
    @brief Generate a lowercase hyphenated RFC 4122 version-4 UUID string

    Used by QPX Parquet exporters when attaching file metadata.

    @return UUID string, e.g. "550e8400-e29b-41d4-a716-446655440000"
  */
  OPENMS_DLLAPI String generateUuidV4();

  /**
    @brief Write an Arrow table to a Parquet file

    @param[in] table The Arrow table to write (must not be null)
    @param[in] filename Output file path
    @param[in] config Parquet writer configuration (compression, row group size, ...)
    @return true on success, false on error (errors are logged)
  */
  OPENMS_DLLAPI bool writeTableToParquet(
    const std::shared_ptr<arrow::Table>& table,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Concatenate a vector of Arrow tables and write the result to a Parquet file

    All tables must share the same schema. An empty input vector is a no-op
    (returns true without writing).

    @param[in] tables Vector of Arrow tables to concatenate (must share schema)
    @param[in] filename Output file path
    @param[in] config Parquet writer configuration
    @return true on success (or if @p tables is empty), false on error
  */
  OPENMS_DLLAPI bool concatenateAndWriteToParquet(
    const std::vector<std::shared_ptr<arrow::Table>>& tables,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ---------------------------------------------------------------------------
  // Read helpers
  // ---------------------------------------------------------------------------

  /**
    @brief Fetch a named column from a table, combining chunks if needed

    Returns nullptr if the column is missing or contains no chunks. When
    @p required is true, missing columns are logged as errors.
  */
  OPENMS_DLLAPI std::shared_ptr<arrow::Array> getColumn(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    bool required = true);

  /// Read a string at @p row, or "" if null/out-of-bounds
  OPENMS_DLLAPI String getStringValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row);

  /// Read a double at @p row, or @p default_val if null
  OPENMS_DLLAPI double getDoubleValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    double default_val = 0.0);

  /// Read a float at @p row, or @p default_val if null
  OPENMS_DLLAPI float getFloatValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    float default_val = 0.0f);

  /// Read an int32 at @p row, or @p default_val if null
  OPENMS_DLLAPI int32_t getInt32Value(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    int32_t default_val = 0);

  /// Read an int64 at @p row, or @p default_val if null
  OPENMS_DLLAPI int64_t getInt64Value(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    int64_t default_val = 0);

  /// Read a bool at @p row, or @p default_val if null
  OPENMS_DLLAPI bool getBoolValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    bool default_val = false);

  /// Whether @p array is null at @p row (or unset)
  OPENMS_DLLAPI bool isNull(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row);

  /**
    @brief Read metavalues from a list<struct{name,value,value_type}> column

    Decodes typed entries (int, double/float, *_list, string) and assigns
    them to @p target. Keys in @p excluded_keys are skipped.
  */
  OPENMS_DLLAPI void readMetaValues(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    MetaInfoInterface& target,
    const std::unordered_set<std::string>& excluded_keys = {});
}

} // namespace OpenMS
