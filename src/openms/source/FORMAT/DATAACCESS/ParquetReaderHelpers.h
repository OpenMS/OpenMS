// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/ParquetFilter.h>

#include <arrow/api.h>
#include <arrow/compute/api.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Shared internal Arrow/Parquet helpers for OpenSWATH readers.

      These helpers centralize schema access, value extraction, and filter
      execution for the XIC/XIM parquet readers while keeping Arrow types out
      of the public OpenMS headers.
    */
    namespace ParquetReaderHelpers
    {
      /**
        @brief Controls how filter expressions are pruned against a schema.
      */
      struct FilterPruningOptions
      {
        /// Columns that are known to exist but cannot be evaluated by the helper.
        std::unordered_set<std::string> unsupported_columns;
        /// Throw if unsupported columns remain after schema-aware pruning.
        bool reject_unsupported_after_pruning{false};
        /// Throw if pruning removes all remaining filter conditions.
        bool reject_empty_after_pruning{false};
      };

      /** @name Table and Schema Access
      */
      //@{
      /**
        @brief Read a parquet table from a single file.

        @param[in] filename Path to the parquet file.

        @return The chunk-combined table.
      */
      std::shared_ptr<arrow::Table> readParquetTable(const std::string& filename);

      /**
        @brief Read and concatenate parquet tables from multiple files.

        @param[in] filenames Paths to parquet files with compatible schemas.

        @return The concatenated, chunk-combined table.
      */
      std::shared_ptr<arrow::Table> readParquetTable(const std::vector<std::string>& filenames);

      /**
        @brief Read the schema from a single parquet file.

        @param[in] filename Path to the parquet file.

        @return The parquet schema.
      */
      std::shared_ptr<arrow::Schema> readParquetSchema(const std::string& filename);

      /**
        @brief Read and validate schemas across multiple parquet files.

        @param[in] filenames Paths to parquet files with expected identical schemas.

        @return The shared schema.
      */
      std::shared_ptr<arrow::Schema> readParquetSchemaAllFiles(const std::vector<std::string>& filenames);

      /**
        @brief Read a subset of columns from a single parquet file.

        @param[in] filename Path to the parquet file.
        @param[in] columns Requested column names.

        @return The chunk-combined table containing the requested columns.
      */
      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::string& filename,
                                                            const std::vector<std::string>& columns);

      /**
        @brief Read and concatenate a subset of columns from multiple parquet files.

        @param[in] filenames Paths to parquet files with compatible schemas.
        @param[in] columns Requested column names.

        @return The concatenated, chunk-combined table containing the requested columns.
      */
      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::vector<std::string>& filenames,
                                                            const std::vector<std::string>& columns);
      //@}

      /** @name Column Access and Value Extraction
      */
      //@{
      /**
        @brief Return a single logical Arrow array for a named table column.

        If the table column is chunked, all chunks are concatenated first.

        @param[in] table Source table.
        @param[in] name Column name.
        @param[in] required Throw if the column is missing.

        @return The resolved Arrow array or `nullptr` for missing optional columns.
      */
      std::shared_ptr<arrow::Array> getColumn(const std::shared_ptr<arrow::Table>& table,
                                              const std::string& name,
                                              bool required = true);

      /**
        @brief Read an optional integer-like value from an Arrow array.

        @param[in] array Source Arrow array.
        @param[in] row Zero-based row index.
        @param[out] value Extracted integer value.

        @return `true` if the value exists, `false` for null cells.
      */
      bool getOptionalInt(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value);

      /**
        @brief Read an optional floating-point value from an Arrow array.

        @param[in] array Source Arrow array.
        @param[in] row Zero-based row index.
        @param[out] value Extracted floating-point value.

        @return `true` if the value exists, `false` for null cells.
      */
      bool getOptionalDouble(const std::shared_ptr<arrow::Array>& array, int64_t row, double& value);

      /**
        @brief Read an optional string value from an Arrow array.

        @param[in] array Source Arrow array.
        @param[in] row Zero-based row index.
        @param[out] value Extracted string value.

        @return `true` if the value exists, `false` for null cells.
      */
      bool getOptionalString(const std::shared_ptr<arrow::Array>& array, int64_t row, std::string& value);

      /**
        @brief Return a binary cell as a `std::string` view copy.

        @param[in] array Source Arrow array.
        @param[in] row Zero-based row index.

        @return Binary payload for the requested cell, or an empty string for null cells.
      */
      std::string getBinaryView(const std::shared_ptr<arrow::Array>& array, int64_t row);
      //@}

      /**
        @brief Parse a textual filter expression into a typed helper expression.

        @param[in] filter User-provided filter expression.
        @param[in] column_types Mapping from uppercase column names to expected scalar types.

        @return Parsed filter expression.
      */
      FilterExpression parseFilter(const std::string& filter,
                                   const std::unordered_map<std::string, ColumnType>& column_types);

#ifdef WITH_ARROW_DATASET
      /**
        @brief Read parquet files through Arrow Dataset with filter pushdown.

        @param[in] filenames Paths to parquet files.
        @param[in] parsed Parsed filter expression.
        @param[in] filter_context Original filter string for diagnostics.
        @param[in] pruning_options Filter pruning policy.

        @return The filtered, chunk-combined table.
      */
      std::shared_ptr<arrow::Table> readParquetTableDataset(const std::vector<std::string>& filenames,
                                                            const FilterExpression& parsed,
                                                            const std::string& filter_context,
                                                            const FilterPruningOptions& pruning_options);
#endif

      /**
        @brief Apply an Arrow compute filter to an in-memory table.

        @param[in] table Source table.
        @param[in] parsed Parsed filter expression.
        @param[in] filter_context Original filter string for diagnostics.
        @param[in] pruning_options Filter pruning policy.

        @return The filtered, chunk-combined table.
      */
      std::shared_ptr<arrow::Table> applyArrowFilter(const std::shared_ptr<arrow::Table>& table,
                                                     const FilterExpression& parsed,
                                                     const std::string& filter_context,
                                                     const FilterPruningOptions& pruning_options);
    } // namespace ParquetReaderHelpers
  } // namespace Internal
} // namespace OpenMS
