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
    @brief RAII wrapper for inter-process file locking.

    The lock is acquired in the constructor and released in the destructor.
    Acquisition waits until the lock becomes available.
  */
  class OPENMS_DLLAPI InterProcessFileLock
  {
    public:
      /**
        @brief Acquire an inter-process lock for the given file.

        Creates the lock file if required and acquires the OS-level file lock.

        @param[in] target_file_path Path to the lock file

        @throws Exception::FileNotWritable if the lock file cannot be created/opened
        @throws Exception::FailedAPICall if lock acquisition fails or times out
      */
      explicit InterProcessFileLock(const String& target_file_path);

      /**
        @brief Release the held inter-process file lock.

        Unlocks the file in a non-throwing cleanup path.
      */
      ~InterProcessFileLock() noexcept;

      // Non-copyable to prevent multiple instances managing the same lock
      InterProcessFileLock(const InterProcessFileLock&) = delete;
      InterProcessFileLock& operator=(const InterProcessFileLock&) = delete;

    private:
      String lock_file_path_;
      std::unique_ptr<boost::interprocess::file_lock> lock_;
      bool locked_ {false};
  };
} // namespace OpenMS
