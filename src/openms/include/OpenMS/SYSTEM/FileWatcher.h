// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#pragma once

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

//STL
#include <map>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <memory>

namespace OpenMS
{
  class String;

  /**
      @brief Watcher that monitors file changes.

      This class monitors files for changes with a delayed callback mechanism.
      The delay prevents multiple callbacks for large files that are written in chunks.

      @ingroup System
  */
  class OPENMS_DLLAPI FileWatcher
  {
public:
    /// Callback type for file change notifications
    using FileChangeCallback = std::function<void(const String&)>;

    /// Constructor
    FileWatcher();

    /// Destructor
    ~FileWatcher();

    /// Sets the delay in seconds (default: 1s)
    void setDelayInSeconds(double delay);

    /// Adds a file to the watcher
    /// @return true if the file was successfully added, false otherwise
    bool addFile(const String& path);

    /// Removes a file from the watcher
    void removeFile(const String& path);

    /// Sets the callback function to be called when a file changes
    void setCallback(FileChangeCallback callback);

protected:
    /// Structure to track pending file change notifications
    struct PendingTimer
    {
      std::string file_path;
      std::chrono::steady_clock::time_point trigger_time;
      bool active;
    };

    /// Monitors file changes (platform-specific implementation)
    void monitorFiles_();

    /// Processes pending timers
    void processTimers_();

    /// Map of watched files to their last modification time
    std::map<std::string, std::time_t> watched_files_;

    /// Pending timers for delayed notifications
    std::map<int, PendingTimer> pending_timers_;
    int next_timer_id_;

    /// Delay (seconds)
    double delay_in_seconds_;

    /// Callback for file changes
    FileChangeCallback callback_;

    /// Monitoring thread
    std::unique_ptr<std::thread> monitor_thread_;

    /// Timer processing thread
    std::unique_ptr<std::thread> timer_thread_;

    /// Flag to stop threads
    std::atomic<bool> stop_monitoring_;

    /// Mutex for thread safety
    std::mutex mutex_;
  };

  // OPENMS_DLLAPI extern FileWatcher myFileWatcher_instance;
}

