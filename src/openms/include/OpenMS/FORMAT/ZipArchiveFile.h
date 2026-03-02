// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{

  /**
    @brief Small libzip-based helpers for working with ZIP archives.

    These helpers were extracted from ParquetFile to provide a reusable
    place for ZIP-related utilities. The implementation uses libzip when
    available; otherwise methods will throw Exception::NotImplemented.

    The API is intentionally small: zip a directory, unzip an archive to a
    temp directory, and add/replace an entry from a filesystem path.
  */
  class OPENMS_DLLAPI ZipArchiveFile
  {
  public:
    /// Create a store-only zip archive from a directory (no additional compression)
    static void zipDirectory(const String& directory_path, const String& output_zip);

    /// Unpack a zip archive into a temporary directory and return the usable path
    static String unzipDirectory(const String& input_path, std::unique_ptr<File::TempDir>& temp_dir);

    /// Add or replace an entry inside an existing zip archive from a file on disk
    /// If the archive does not exist it will be created.
    static void addOrReplaceFromFile(const String& archive_path, const String& entry_name, const String& source_file_path);

    /// List entries in a zip archive (returns empty list if not available)
    static std::vector<String> listEntries(const String& archive_path);
  };

} // namespace OpenMS
