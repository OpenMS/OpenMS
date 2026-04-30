// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <set>

namespace OpenMS
{
  /**
    @brief Handles temporary files, ensuring cleanup on destruction.

    This class provides a simple interface to register temporary files that should be removed when the manager is destroyed. 
    It also maintains a persistent registry of temporary files to enable cleanup of files from prior runs that for example crashed before cleanup.
    By doing this, we can assure that temporary files do not accumulate on the endusers drive over time, even in the case of unexpected program termination.

    // an funktion, wo die exceptions passieren
    @throw Exception::InvalidParameter If an empty file path is provided to managed-file APIs.
    @throw Exception::FileNotReadable If an existing registry file cannot be opened for reading.
    @throw Exception::FileNotWritable If the registry file cannot be written.
    @throw Exception::ParseError If the registry file content cannot be read/parsing fails.

    @ingroup System
  */
  class OPENMS_DLLAPI TempFileManager
  {
    public:
      /**
        @brief Constructor with explicit registry id

        Each usage is expected to provide its own id to isolate temporary
        file registries between different tools/components.

        @param[in] registry_id Logical id used to derive registry file name
      */
      explicit TempFileManager(const String& registry_id);

      /**
        @brief Move constructor

        @param[in,out] temp_file_manager Source manager
      */
      TempFileManager(TempFileManager&& temp_file_manager) noexcept;

      // Deleted copy constructor to enforce unique ownership semantics
      TempFileManager(const TempFileManager& temp_file_manager) = delete;

      /**
        @brief Destructor

        Removes all stored temporary files known to this manager and the persisted registry
      */
      virtual ~TempFileManager() noexcept;

      /**
        @brief Move assignment operator

        @param[in,out] temp_file_manager Source manager
        @return Reference to @p *this
      */
      TempFileManager& operator=(TempFileManager&& temp_file_manager) noexcept;
      /// Disable the copy constructor
      TempFileManager& operator=(const TempFileManager& temp_file_manager) = delete;

      /**
        @brief Adds a temporary file to be managed
        @param[in] file_path The path to the temporary file
      */
      void addFile(const String& file_path);

      /**
        @brief Releases a file from management without removing it from disk
        @param[in] file_path The path to the temporary file
      */
      void releaseFile(const String& file_path);

      /**
        @brief Removes a managed file immediately and unregisters it.

        @param[in] file_path The path to the temporary file
        @return True if the file was removed successfully or did not exist, false otherwise
      */
      bool removeFileNow(const String& file_path) noexcept;

      /**
        @brief Removes all currently managed files and clears the persisted registry for this scope
      */
      void cleanupNow() noexcept;

    private:
      /**
        @brief Sanitizes a registry id so it can be used as filename component.

        @param[in] registry_id Input registry id string
        @return Sanitized registry id
      */
      static String sanitizeRegistryId_(const String& registry_id);

      /**
        @brief Builds the registry file path from registry id and scope options.

        @param[in] registry_id Logical id used for file naming
        @return Absolute or relative registry file path
      */
      static String getRegistryFilePath_(const String& registry_id);

      /**
        @brief Loads the registry file and populates the set of temporary files.
        @param[in,out] files The set to populate with file paths from the registry
      */
      void loadRegistryFile_(std::set<String>& files) const;

      /**
        @brief Writes the complete registry set to disk.

        @param[in] files The file set to persist
      */
      void writeRegistryFile_(const std::set<String>& files) const;

      /**
        @brief Adds or removes one file entry from the persisted registry.

        @param[in] file_path The file path to modify
        @param[in] action Action to apply to the registry entry
      */
      enum class RegistryAction
      {
        ADD,
        REMOVE
      };

      void updateRegistryEntry_(const String& file_path, RegistryAction action);

      const String registry_file_path_;
      std::set<String> files_;
  };
} // namespace OpenMS
