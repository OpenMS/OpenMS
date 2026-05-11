// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/sync/file_lock.hpp>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <thread>

#include <OpenMS/SYSTEM/InterProcessFileLock.h>

namespace OpenMS
{
  InterProcessFileLock::InterProcessFileLock() noexcept = default;

  InterProcessFileLock::InterProcessFileLock(const String& target_file_path, UInt timeout_ms) noexcept
  {
    lock(target_file_path, timeout_ms);
  }

  bool InterProcessFileLock::lock(const String& target_file_path, UInt timeout_ms) noexcept
  {
    try
    {
      return lockImpl_(target_file_path, timeout_ms);
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock acquisition failed: "
                      << e.getName() << ": " << e.what() << "\n";
    }
    catch (const boost::interprocess::interprocess_exception& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock acquisition failed: "
                      << e.what() << "\n";
    }
    catch (const std::exception& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock acquisition failed: "
                      << e.what() << "\n";
    }
    catch (...)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock acquisition failed: unknown exception\n";
    }
    return false;
  }

  bool InterProcessFileLock::lockImpl_(const String& target_file_path, UInt timeout_ms)
  {
    // If we already hold a lock, release it before acquiring a new one
    if (locked_ && lock_ != nullptr)
    {
      try
      {
        lock_->unlock();
      }
      catch (const boost::interprocess::interprocess_exception& e)
      {
        OPENMS_LOG_WARN << "InterProcessFileLock failed to unlock current lock file '"
                        << lock_file_path_ << "' before relocking: " << e.what() << "\n";
        return false;
      }
      catch (const std::exception& e)
      {
        OPENMS_LOG_WARN << "InterProcessFileLock failed to unlock current lock file '"
                        << lock_file_path_ << "' before relocking: " << e.what() << "\n";
        return false;
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "InterProcessFileLock failed to unlock current lock file '"
                        << lock_file_path_ << "' before relocking: unknown exception\n";
        return false;
      }
      locked_ = false;
      lock_.reset();
    }

    lock_file_path_ = target_file_path;

    // Ensure the lock file exists before trying to lock it
    // boost::interprocess::file_lock requires an existing file and does not create one implicitly
    std::ofstream lock_file_stream(lock_file_path_.c_str(), std::ios::out | std::ios::app);
    if (!lock_file_stream.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, lock_file_path_);
    }
    lock_file_stream.close();

    // Create the file lock object
    lock_ = std::make_unique<boost::interprocess::file_lock>(lock_file_path_.c_str());

    // Polling acquisition with timeout prevents indefinite waits
    // If a process crashes, OS-level file locks are released automatically
    const UInt effective_timeout_ms = std::max(timeout_ms, static_cast<UInt>(1));
    const auto timeout = std::chrono::milliseconds(effective_timeout_ms);
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
    return true;
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
    catch (const boost::interprocess::interprocess_exception& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock destructor failed to unlock '"
                      << lock_file_path_ << "': " << e.what() << "\n";
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock destructor failed to unlock '"
                      << lock_file_path_ << "': " << e.getName() << ": " << e.what() << "\n";
    }
    catch (const std::exception& e)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock destructor failed to unlock '"
                      << lock_file_path_ << "': " << e.what() << "\n";
    }
    catch (...)
    {
      OPENMS_LOG_WARN << "InterProcessFileLock destructor failed to unlock '"
                      << lock_file_path_ << "': unknown exception\n";
    }
  }
} // namespace OpenMS
