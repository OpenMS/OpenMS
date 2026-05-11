// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>

#include <memory>

namespace boost::interprocess
{
  class file_lock;
}

namespace OpenMS
{
  /**
    @brief RAII wrapper for inter-process file locking

    The lock is acquired in the constructor and released in the destructor.
    Acquisition waits until the lock becomes available.
  */
  class OPENMS_DLLAPI InterProcessFileLock
  {
    public:
    
      // Default constructor; call lock() to acquire a lock later
      InterProcessFileLock() noexcept;

      /**
        @brief Acquire an inter-process lock for the given file.

        @param[in] target_file_path Path to the lock file
        @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts. Default is 3000 ms.

        This constructor does not throw; failures are logged and the instance remains unlocked.
      */
      explicit InterProcessFileLock(const String& target_file_path, uint timeout_ms = 3000) noexcept;

      /**
        @brief Acquire (or re-acquire) an inter-process lock for a file.

        @param[in] target_file_path Path to the lock file
        @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts. Default is 3000 ms.
        @return True if lock acquisition succeeded, false otherwise
      */
      bool lock(const String& target_file_path, uint timeout_ms = 3000) noexcept;

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
       @brief Internal implementation of lock acquisition with timeout.
       This method is called by the constructor and the public lock() method.

       @param[in] target_file_path Path to the lock file
       @param[in] timeout_ms Timeout in milliseconds for lock acquisition attempts. Default is 3000 ms.
       
       @return True if lock acquisition succeeded, false otherwise
       */
      bool lockImpl_(const String& target_file_path, uint timeout_ms);

      String lock_file_path_;
      std::unique_ptr<boost::interprocess::file_lock> lock_;
      bool locked_ {false};
  };
} // namespace OpenMS
