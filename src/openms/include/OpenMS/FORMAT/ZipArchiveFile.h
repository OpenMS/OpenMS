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

    The implementation uses libzip when available; otherwise methods will throw Exception::NotImplemented.

    Zip a directory, unzip an archive to a
    temp directory, and add/replace an entry from a filesystem path.
  */
  class OPENMS_DLLAPI ZipArchiveFile
  {
  public:
    /**
      @brief Create a store-only zip archive from a directory (no additional compression)

      The created archive will store files without additional compression (ZIP_CM_STORE).
      If libzip is not available this method will throw Exception::NotImplemented.

      @param[in] directory_path Path to the directory whose contents will be archived.
      @param[in] output_zip Path to the output zip file to create (existing file will be overwritten).
      @throws Exception::FileNotWritable if the archive cannot be written.
      @throws Exception::NotImplemented if libzip support is unavailable.
    */
    static void zipDirectory(const String& directory_path, const String& output_zip);

    /**
      @brief Unpack a zip archive into a temporary directory and return the usable path

      Extracts the archive into a newly created File::TempDir and returns the path to
      the unpacked directory. This function protects against path traversal and
      absolute paths in archive entries. If the provided path already points to a
      directory it is returned unchanged and no TempDir is created.

      @param[in] input_path Path to the zip archive (or a directory).
      @param[out] temp_dir A unique_ptr which will be set to the owned TempDir when an archive is extracted. If input_path is a directory this will remain unchanged.
      @return The path where the archive was unpacked (or input_path if already a directory).
      @throws Exception::FileNotFound if the input archive is not readable.
      @throws Exception::InvalidValue on archive extraction errors.
      @throws Exception::NotImplemented if libzip support is unavailable.
    */
    static String unzipDirectory(const String& input_path, std::unique_ptr<File::TempDir>& temp_dir);

    /**
      @brief Add or replace an entry inside an existing zip archive from a file on disk

      Streams the contents of @p source_file_path into @p archive_path as the entry
      @p entry_name using a file-backed zip source. The archive is created if it
      does not exist. This avoids loading the entire file into memory and is suitable
      for large parquet blobs.

      @param[in] archive_path Path to the zip archive to modify (created if missing).
      @param[in] entry_name Relative entry name inside the archive (use '/' separators).
      @param[in] source_file_path Path to the source file on disk to add or replace.
      @throws Exception::FileNotFound if the source file is not readable.
      @throws Exception::InvalidValue on libzip errors.
      @throws Exception::NotImplemented if libzip support is unavailable.
    */
    static void addOrReplaceFromFile(const String& archive_path, const String& entry_name, const String& source_file_path);

    /**
      @brief List entries in a zip archive (returns empty list if not available)

      @param[in] archive_path Path to the zip archive to inspect.
      @return A vector of entry names contained in the archive. Returns an empty
              vector if the archive cannot be read or libzip is unavailable.
    */
    static std::vector<String> listEntries(const String& archive_path);

    /**
      @brief Write a small JSON sidecar index for the archive listing entries and sizes.

        For zip archives this helper will write the sidecar JSON as an entry inside
        the archive (named 'basename.idx.json', e.g. 'results.oswpq.idx.json'
        for an archive named 'results.oswpq'). This keeps the archive self-contained
        and portable. When the provided path is a directory (an unpacked layout)
        the helper will still write an external 'archive_path.idx.json' file.

        The JSON object maps entry name -> size in bytes and can be used by
        readers to implement fast index-based access.

        @param[in] archive_path Path to the zip archive (or directory) to index.

        If libzip is unavailable this method will throw Exception::NotImplemented.
    */
    static void writeSidecarIndex(const String& archive_path);

    /**
    @brief Extract a single entry from a zip archive into a temporary file and return its path

      The function will create the provided @p temp_dir (if null) and extract the
      requested entry into that directory, preserving any subdirectory structure
      present in the entry name. This is a lightweight alternative to extracting
      the whole archive when only a few files are needed.

      @param[in] archive_path Path to the zip archive
      @param[in] entry_name Relative name of the entry to extract
      @param[out] temp_dir Unique pointer to a TempDir which will be created if null and will own the extracted file(s)
      @return The absolute path to the extracted file on disk
      @throws Exception::FileNotFound if the archive or entry is not found
      @throws Exception::InvalidValue on extraction errors
    */
    static String extractEntryToTempFile(const String& archive_path, const String& entry_name, std::unique_ptr<File::TempDir>& temp_dir);

  };

} // namespace OpenMS
