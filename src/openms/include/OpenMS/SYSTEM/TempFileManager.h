// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/InterProcessFileLock.h>

#include <memory>
#include <mutex>
#include <utility>
#include <set>

namespace OpenMS
{
  /**
    @brief Manages temporary files via a lock-protected on-disk registry.

    Designed for typical single-user desktop/tool usage (e.g. TOPPView):
    one process owns one registry slot at a time, with occasional crash recovery.

    The implementation provides practical safeguards already used by OpenMS:
    inter-process lock files, per-instance state locking, and persistent registry
    cleanup across runs. It is not a hard security boundary; if used in new or
    adversarial environments, behavior should be validated with scenario-specific tests.

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
        @param[in] lock_timeout_ms Timeout in milliseconds for inter-process lock acquisition
      */
      explicit TempFileManager(const String& registry_id, UInt lock_timeout_ms = 3000);
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

        @throw Exception::InvalidParameter If @p file_path is empty.
        @throw Exception::FileNotReadable If an existing registry file cannot be opened for reading.
        @throw Exception::FileNotWritable If the registry file cannot be written.
        @throw Exception::ParseError If reading/parsing the registry file fails.
      */
      void addFile(const String& file_path);

      /**
        @brief Releases a file from management without removing it from disk

        @param[in] file_path The path to the temporary file

        @throw Exception::InvalidParameter If @p file_path is empty.
        @throw Exception::FileNotReadable If an existing registry file cannot be opened for reading.
        @throw Exception::FileNotWritable If the registry file cannot be written.
        @throw Exception::ParseError If reading/parsing the registry file fails.
      */
      void releaseFile(const String& file_path);

      /**
        @brief Removes a managed file immediately and unregisters it.

        @param[in] file_path The path to the temporary file

        @return True if the file was removed successfully or did not exist, false otherwise
      */
      bool removeFileNow(const String& file_path) noexcept;

    private:
      TempFileManager(std::pair<String, std::unique_ptr<InterProcessFileLock>>&& selected_slot, UInt lock_timeout_ms);

      /**
        @brief Removes all currently managed files and clears the persisted registry for this scope
      */
      void cleanup_() noexcept;

      /**
        @brief Sanitizes a registry id so it can be used as filename component.

        @param[in] registry_id Input registry id string
        @return Sanitized registry id
      */
      static String sanitizeRegistryId_(const String& registry_id);

      /**
        @brief Builds the registry file path for a specific instance slot.

        @param[in] registry_id Logical id used for file naming
        @param[in] instance_index 1-based instance slot index

        @return Absolute or relative registry file path
      */
      static String getRegistryFilePathForInstance_(const String& registry_id, Size instance_index);

      /**
        @brief Builds the lock file path for a specific instance slot.

        @param[in] registry_id Logical id used for file naming
        @param[in] instance_index 1-based instance slot index

        @return Absolute or relative lock file path
      */
      static String getLockFilePathForInstance_(const String& registry_id, Size instance_index);

      /**
        @brief Selects a free instance registry path and cleans up stale crashed registries.

        @param[in] registry_id Logical id used for file naming
        @param[in] lock_timeout_ms Timeout in milliseconds for lock acquisition
        
        @return Registry file path for this manager instance
      */
      static String selectRegistryFilePath_(const String& registry_id, UInt lock_timeout_ms, std::unique_ptr<InterProcessFileLock>& registry_lock);

      // Internal helper to atomically select a slot path and acquire its persistent lock
      static std::pair<String, std::unique_ptr<InterProcessFileLock>> selectRegistrySlot_(const String& registry_id, UInt lock_timeout_ms);

      /**
        @brief Cleanup helper for a specific registry path.

        @param[in] registry_file_path Path to a registry file
      */
      static void cleanupRegistryFile_(const String& registry_file_path) noexcept;

      /**
        @brief Loads the registry file and populates the set of temporary files.

        @param[in,out] files The set to populate with file paths from the registry

        @throw Exception::FileNotReadable If an existing registry file cannot be opened for reading.
        @throw Exception::ParseError If reading/parsing the registry file fails.
      */
      void loadRegistryFile_(std::set<String>& files) const;

      /**
        @brief Writes the complete registry set to disk.

        @param[in] files The file set to persist

        @throw Exception::FileNotWritable If the registry file cannot be written.
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
      const String lock_file_path_;
      UInt lock_timeout_ms_;
      std::unique_ptr<InterProcessFileLock> registry_lock_;
      mutable std::mutex state_mutex_;
      std::set<String> files_;
  };
} // namespace OpenMS
