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
    namespace ParquetReaderHelpers
    {
      struct FilterPruningOptions
      {
        std::unordered_set<std::string> unsupported_columns;
        bool reject_unsupported_after_pruning{false};
        bool reject_empty_after_pruning{false};
      };

      std::shared_ptr<arrow::Table> readParquetTable(const std::string& filename);
      std::shared_ptr<arrow::Table> readParquetTable(const std::vector<std::string>& filenames);

      std::shared_ptr<arrow::Schema> readParquetSchema(const std::string& filename);
      std::shared_ptr<arrow::Schema> readParquetSchemaAllFiles(const std::vector<std::string>& filenames);

      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::string& filename,
                                                            const std::vector<std::string>& columns);
      std::shared_ptr<arrow::Table> readParquetTableColumns(const std::vector<std::string>& filenames,
                                                            const std::vector<std::string>& columns);

      std::shared_ptr<arrow::Array> getColumn(const std::shared_ptr<arrow::Table>& table,
                                              const std::string& name,
                                              bool required = true);

      bool getOptionalInt(const std::shared_ptr<arrow::Array>& array, int64_t row, Int64& value);
      bool getOptionalDouble(const std::shared_ptr<arrow::Array>& array, int64_t row, double& value);
      bool getOptionalString(const std::shared_ptr<arrow::Array>& array, int64_t row, std::string& value);
      std::string getBinaryView(const std::shared_ptr<arrow::Array>& array, int64_t row);

      FilterExpression parseFilter(const std::string& filter,
                                   const std::unordered_map<std::string, ColumnType>& column_types);

#ifdef WITH_ARROW_DATASET
      std::shared_ptr<arrow::Table> readParquetTableDataset(const std::vector<std::string>& filenames,
                                                            const FilterExpression& parsed,
                                                            const std::string& filter_context,
                                                            const FilterPruningOptions& pruning_options);
#endif

      std::shared_ptr<arrow::Table> applyArrowFilter(const std::shared_ptr<arrow::Table>& table,
                                                     const FilterExpression& parsed,
                                                     const std::string& filter_context,
                                                     const FilterPruningOptions& pruning_options);
    } // namespace ParquetReaderHelpers
  } // namespace Internal
} // namespace OpenMS
