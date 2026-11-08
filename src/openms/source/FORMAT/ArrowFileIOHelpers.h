// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace arrow
{
  class StringBuilder;
  class StructBuilder;
  class Table;
}

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers shared by the Arrow/Parquet file classes (FeatureMapArrowIO,
      ConsensusMapArrowIO, ProteinIdentificationArrowIO, ConsensusMapArrowExport, QPXFile).

      Writing and reading the Parquet files of these classes (with their file-level metadata),
      the list<struct{name,value,type}> metavalue columns and the JSON encoding of
      DataProcessing and metavalues in the file metadata. Not installed; Arrow types stay out of
      the public headers.

      @p context is the name of the calling class; it prefixes the log messages.
    */
    namespace ArrowFileIO
    {
      /// Append all metavalues of @p mii (except @p excluded_keys) to the builders of a
      /// list<struct{name,value,type}> column.
      void appendMetaValues(const MetaInfoInterface& mii,
                            std::shared_ptr<arrow::StringBuilder>& name_b,
                            std::shared_ptr<arrow::StringBuilder>& value_b,
                            std::shared_ptr<arrow::StringBuilder>& type_b,
                            std::shared_ptr<arrow::StructBuilder>& struct_b,
                            const std::unordered_set<std::string>& excluded_keys);

      /// Escape a string for JSON (quotes, backslashes, control characters).
      std::string escapeJsonString(const std::string& s);

      /// Undo escapeJsonString().
      std::string unescapeJsonString(const std::string& s);

      /// Skip JSON whitespace at @p pos.
      void skipWhitespace(const std::string& s, size_t& pos);

      /// Parse the JSON string starting at @p pos (which must point at the opening quote) and
      /// advance @p pos past it.
      std::string parseJsonString(const std::string& s, size_t& pos);

      /// Serialize DataProcessing entries to a JSON array.
      std::string serializeDataProcessing(const std::vector<DataProcessing>& dps);

      /// Parse the JSON array written by serializeDataProcessing(). Unknown processing actions
      /// are logged (prefixed with @p context) and skipped.
      std::vector<DataProcessing> deserializeDataProcessing(const std::string& json, const char* context);

      /// Serialize all metavalues of @p mii to a JSON array of {name, value, type} objects.
      std::string serializeMetaValues(const MetaInfoInterface& mii);

      /// Parse the JSON written by serializeMetaValues() into @p target.
      void deserializeMetaValues(const std::string& json, MetaInfoInterface& target);

      /**
        @brief Write @p table to @p filename with the file-level metadata of these formats
        (creator, @p file_type, creation date, a fresh UUID, plus @p extra_metadata).

        A partially written file is removed. Errors are logged; returns false on failure.
      */
      bool writeParquetTable(std::shared_ptr<arrow::Table> table,
                             const std::string& filename,
                             const std::string& file_type,
                             const ParquetWriteConfig& config,
                             const char* context,
                             const std::unordered_map<std::string, std::string>& extra_metadata = {});

      /// Read a whole Parquet file into one table with combined chunks. Errors are logged;
      /// returns nullptr on failure.
      std::shared_ptr<arrow::Table> readParquetTable(const std::string& filename, const char* context);

      /**
        @brief The name of metavalue @p index, cached in @p cache (indexed by metavalue index).

        Saves the lock of MetaInfoRegistry::getName() per value in export loops. @p cache must not
        be shared between threads.
      */
      inline const std::string& cachedMetaName(UInt index, std::vector<std::string>& cache)
      {
        if (index >= cache.size()) { cache.resize(static_cast<size_t>(index) + 1); }
        std::string& name = cache[index];
        if (name.empty()) { name = MetaInfoInterface::metaRegistry().getName(index); }
        return name;
      }
    } // namespace ArrowFileIO
  } // namespace Internal
} // namespace OpenMS
