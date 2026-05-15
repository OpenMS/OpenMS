// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <memory>

namespace boost::interprocess
{
  class file_lock;
}

namespace OpenMS
{
  /**
    @brief Lightweight RAII wrapper for OS-backed file locks.

    This class wraps `boost::interprocess::file_lock` for simple lock-file based
    coordination. It provides blocking, non-blocking, and timed acquisition
    methods and releases the lock in the destructor.

    It is intended as a practical coordination helper (e.g. for temporary-file
    registries), not as a standalone security mechanism. Behavior depends on
    OS/file-system locking semantics and should be validated in each concrete use.

    @ingroup System
  */
  class OPENMS_DLLAPI InterProcessFileLock
  {
    public:
      /**
        @brief Default constructor.

        Creates an unlocked instance without acquiring any inter-process lock.
        Call lock(), tryLock(), or timedLock() to acquire a lock for a target file.

        @note This constructor does not lock any file.
        @see lock(), tryLock(), timedLock()
      */
      InterProcessFileLock() noexcept;

      /**
        @brief Acquire an inter-process lock for the given file.

        @param[in] target_file_path Path to the lock file
        @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts. Default is 100 ms.

        This constructor does not throw; failures are logged and the instance remains unlocked.
      */
      explicit InterProcessFileLock(const String& target_file_path, UInt timeout_ms = 100) noexcept;

      /**
        @brief Acquire (or re-acquire) an inter-process lock for a file.

        @param[in] target_file_path Path to the lock file

        Blocks until the lock becomes available.

        @return True if lock acquisition succeeded, false otherwise
      */
      bool lock(const String& target_file_path) noexcept;

      /**
        @brief Try to acquire an inter-process lock for a file without waiting.

        @param[in] target_file_path Path to the lock file

        @return True if lock acquisition succeeded, false otherwise
      */
      bool tryLock(const String& target_file_path) noexcept;

      /**
        @brief Try to acquire an inter-process lock for a file within a timeout.

        @param[in] target_file_path Path to the lock file
        @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts.
          A value of 0 is normalized to 1 ms.

        @return True if lock acquisition succeeded, false otherwise
      */
      bool timedLock(const String& target_file_path, UInt timeout_ms = 100) noexcept;

      /**
        @brief Returns whether this instance currently holds a lock.
      */
      bool isLocked() const noexcept { return locked_; }

      /**
        @brief Release the held inter-process file lock.
        Unlocks the file in a non-throwing cleanup path.
      */
      ~InterProcessFileLock() noexcept;

      // Non-copyable to prevent multiple instances managing the same lock
      InterProcessFileLock(const InterProcessFileLock&) = delete;
      InterProcessFileLock& operator=(const InterProcessFileLock&) = delete;

    private:
      /**
       @brief Internal implementation of lock acquisition.
       This method is called by the constructor and public lock methods.

       @param[in] target_file_path Path to the lock file
       @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts.
       @param[in] wait_forever If true, acquire lock by blocking until available.
       @param[in] non_blocking If true, attempt exactly once and return immediately.
        
       @return True if lock acquisition succeeded, false otherwise
       */
      bool lockImpl_(const String& target_file_path, UInt timeout_ms, bool wait_forever, bool non_blocking);

      String lock_file_path_;
      std::unique_ptr<boost::interprocess::file_lock> lock_;
      bool locked_ {false};
  };
} // namespace OpenMS
