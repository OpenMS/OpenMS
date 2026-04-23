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
      explicit InterProcessFileLock(const String& target_file_path);
      ~InterProcessFileLock() noexcept;

      InterProcessFileLock(const InterProcessFileLock&) = delete;
      InterProcessFileLock& operator=(const InterProcessFileLock&) = delete;

    private:
      String lock_file_path_;
      std::unique_ptr<boost::interprocess::file_lock> lock_;
      bool locked_ {false};
  };
} // namespace OpenMS
