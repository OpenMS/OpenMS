// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>

#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/sync/file_lock.hpp>

#include <chrono>
#include <fstream>
#include <thread>

#include <OpenMS/SYSTEM/InterProcessFileLock.h>

namespace OpenMS
{
  InterProcessFileLock::InterProcessFileLock(const String& target_file_path) :
    lock_file_path_(target_file_path)
  {
    // Ensure the lock file exists before trying to lock it. This also creates the file if it doesn't exist, which is necessary for locking
    std::ofstream lock_file_stream(lock_file_path_.c_str(), std::ios::out | std::ios::app);
    if (!lock_file_stream.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, lock_file_path_);
    }
    lock_file_stream.close();

    try
    {
      // Create the file lock object
      lock_ = std::make_unique<boost::interprocess::file_lock>(lock_file_path_.c_str());

      // Polling acquisition with timeout prevents indefinite waits
      // If a process crashes, OS-level file locks are released automatically
      const auto timeout = std::chrono::seconds(3);
      const auto retry_interval = std::chrono::milliseconds(50);
      const auto start = std::chrono::steady_clock::now();

      // Try to acquire the lock, with timeout
      while (!lock_->try_lock())
      {
        if (std::chrono::steady_clock::now() - start >= timeout)
        {
          throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Timeout while waiting to acquire lock file '" + lock_file_path_ + "'.");
        }
        std::this_thread::sleep_for(retry_interval);
      }
      locked_ = true;
    }
    catch (const Exception::FailedAPICall& e)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Failed to acquire lock file '" + lock_file_path_ + "': " + String(e.what()));
    }
    catch (const Exception::BaseException& e)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Failed to acquire lock file '" + lock_file_path_ + "': " +
                                       String(e.getName()) + " (" + String(e.what()) + ")");
    }
    catch (const boost::interprocess::interprocess_exception& e)
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Failed to acquire lock file '" + lock_file_path_ + "': " + e.what());
    }
  }

  InterProcessFileLock::~InterProcessFileLock() noexcept
  {
    // Ensure we only try to unlock if we actually hold the lock
    if (!locked_ || lock_ == nullptr)
    {
      return;
    }

    // unlock
    try
    {
      lock_->unlock();
      locked_ = false;
    }
    catch (const boost::interprocess::interprocess_exception&)
    {
    }
    catch (const Exception::BaseException&)
    {
    }
    catch (const std::exception&)
    {
    }
  }
} // namespace OpenMS
