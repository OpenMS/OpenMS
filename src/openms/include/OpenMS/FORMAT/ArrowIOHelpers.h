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

#include <memory>
#include <vector>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{

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
}

} // namespace OpenMS
