// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <sys/stat.h>
#include <ctime>
#include <cerrno>
#include <cstring>
#include <vector>

using namespace std;

namespace OpenMS
{
  FileWatcher::FileWatcher() :
    watched_files_(),
    pending_timers_(),
    next_timer_id_(0),
    delay_in_seconds_(1.0),
    callback_(),
    monitor_thread_(),
    timer_thread_(),
    stop_monitoring_(false),
    mutex_()
  {
    // Start monitoring thread
    monitor_thread_ = std::make_unique<std::thread>(&FileWatcher::monitorFiles_, this);

    // Start timer processing thread
    timer_thread_ = std::make_unique<std::thread>(&FileWatcher::processTimers_, this);
  }

  FileWatcher::~FileWatcher()
  {
    // Signal threads to stop
    stop_monitoring_ = true;

    // Wait for threads to finish
    if (monitor_thread_ && monitor_thread_->joinable())
    {
      monitor_thread_->join();
    }
    if (timer_thread_ && timer_thread_->joinable())
    {
      timer_thread_->join();
    }
  }

  void FileWatcher::setDelayInSeconds(double delay)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    delay_in_seconds_ = delay;
  }

  bool FileWatcher::addFile(const String& path)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    // Get initial modification time
    struct stat file_stat;
    if (stat(path.c_str(), &file_stat) == 0)
    {
      watched_files_[path.c_str()] = file_stat.st_mtime;
      return true;
    }
    else
    {
      // stat() failed - capture errno immediately before it can be overwritten
      int err = errno;
      OPENMS_LOG_ERROR << "FileWatcher: Failed to add file '" << path
                       << "': " << std::strerror(err)
                       << " (errno=" << err << ")" << std::endl;
      return false;
    }
  }

  void FileWatcher::removeFile(const String& path)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    watched_files_.erase(path.c_str());
  }

  void FileWatcher::setCallback(FileChangeCallback callback)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    callback_ = callback;
  }

  void FileWatcher::monitorFiles_()
  {
    while (!stop_monitoring_)
    {
      // Copy watched files map to reduce lock contention
      // (stat() calls can be slow on network filesystems)
      std::map<std::string, std::time_t> files_snapshot;
      int delay_ms;

      {
        std::lock_guard<std::mutex> lock(mutex_);
        files_snapshot = watched_files_;
        delay_ms = static_cast<int>(delay_in_seconds_ * 1000);
      }

      // Check each watched file for changes (without holding lock)
      for (auto& entry : files_snapshot)
      {
        const std::string& file_path = entry.first;
        std::time_t snapshot_mtime = entry.second;

        struct stat file_stat;
        if (stat(file_path.c_str(), &file_stat) == 0)
        {
          // Check if file was modified
          if (file_stat.st_mtime != snapshot_mtime)
          {
            // File changed - update under lock
            std::lock_guard<std::mutex> lock(mutex_);

            // Check if file is still being watched (might have been removed)
            auto it = watched_files_.find(file_path);
            if (it == watched_files_.end())
            {
              continue; // File was removed while we were checking
            }

            // Update last modification time
            it->second = file_stat.st_mtime;

            // Look for existing timer for this file
            int existing_timer_id = -1;
            for (auto& timer_entry : pending_timers_)
            {
              if (timer_entry.second.file_path == file_path && timer_entry.second.active)
              {
                existing_timer_id = timer_entry.first;
                break;
              }
            }

            // Create new timer or reset existing one
            if (existing_timer_id == -1)
            {
              // Create new timer
              PendingTimer timer;
              timer.file_path = file_path;
              timer.trigger_time = std::chrono::steady_clock::now() +
                                  std::chrono::milliseconds(delay_ms);
              timer.active = true;
              pending_timers_[next_timer_id_++] = timer;
            }
            else
            {
              // Reset existing timer
              pending_timers_[existing_timer_id].trigger_time =
                  std::chrono::steady_clock::now() +
                  std::chrono::milliseconds(delay_ms);
            }
          }
        }
      }

      // Sleep for a short interval (100ms) to avoid excessive CPU usage
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  }

  void FileWatcher::processTimers_()
  {
    while (!stop_monitoring_)
    {
      // Collect expired timers to invoke outside the lock
      std::vector<std::string> expired_files;
      FileChangeCallback callback_copy;

      {
        std::lock_guard<std::mutex> lock(mutex_);

        auto now = std::chrono::steady_clock::now();

        // Copy callback under lock
        callback_copy = callback_;

        // Check for expired timers and collect them
        for (auto it = pending_timers_.begin(); it != pending_timers_.end(); )
        {
          if (it->second.active && now >= it->second.trigger_time)
          {
            // Timer expired - collect file path
            expired_files.push_back(it->second.file_path);

            // Remove timer
            it = pending_timers_.erase(it);
          }
          else
          {
            ++it;
          }
        }
      }
      // mutex_ is now released

      // Invoke callbacks outside the lock to avoid deadlock
      if (callback_copy)
      {
        for (const auto& file_path : expired_files)
        {
          callback_copy(file_path);
        }
      }

      // Sleep for a short interval (50ms) to check timers frequently
      std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
  }

  //OPENMS_DLLAPI FileWatcher myFileWatcher_instance; // required, such that the moc file get generated during building OpenMS.dll, not later during OpenMS_GUI.dll as DLL flags are wrong then

} // namespace OpenMS
