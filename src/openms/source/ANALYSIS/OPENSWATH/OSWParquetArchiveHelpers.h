// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <arrow/api.h>

#include <memory>
#include <string>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Private helpers for reading the Parquet files of an .oswpq archive (or directory). Not installed.
    */
    namespace OSWParquetArchiveHelpers
    {
      /**
        @brief Reads the Parquet table stored as @p entry of @p oswpq_dir

        Tries to open the entry directly from the archive using a RandomAccessFile. If that fails
        (e.g., compressed entry or libzip not available), falls back to extracting the entry to a
        temporary file and reading from disk.

        @param[in] oswpq_dir Path of the .oswpq archive or directory
        @param[in] entry Name of the entry, e.g. "library/precursors.parquet"
        @param[in,out] temp_dir Temporary directory used for extracted entries (created on demand)
        @param[in] missing_message Start of the exception message if the entry is missing, e.g. "Missing required parquet file"
        @param[in] function Name of the calling function, reported by the exception

        @exception Exception::MissingInformation if the entry does not exist
      */
      inline std::shared_ptr<arrow::Table> openTableFromEntry(const std::string& oswpq_dir, const std::string& entry, std::unique_ptr<TempDir>& temp_dir,
                                                              const char* missing_message, const char* function)
      {
        auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
        if (ra_res.ok())
        {
          const auto& raf = ra_res.ValueOrDie();
          return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
        }
        // Fallback: extract to temp file and read
        const std::string path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
        if (!File::exists(path))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, function,
                                              std::string(missing_message) + " '" + entry + "' in '" + oswpq_dir + "'");
        }
        return ParquetFile::readTable(path);
      }
    } // namespace OSWParquetArchiveHelpers
  } // namespace Internal
} // namespace OpenMS
