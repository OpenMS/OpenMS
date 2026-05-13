// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/InterProcessFileLock.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/TempFileManager.h>

#include <cassert>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <mutex>
#include <string>
#include <utility>

namespace OpenMS
{
  namespace
  {
    // Process-local slot ownership to disambiguate multiple managers in one process
    std::set<String> active_registry_paths;
    std::mutex active_registry_paths_mutex;
  }

  TempFileManager::TempFileManager(std::pair<String, std::unique_ptr<InterProcessFileLock>>&& selected_slot, UInt lock_timeout_ms) :
    registry_file_path_(std::move(selected_slot.first)),
    lock_file_path_(registry_file_path_ + ".lock"),
    lock_timeout_ms_(lock_timeout_ms),
    registry_lock_(std::move(selected_slot.second))
  {
    loadRegistryFile_(files_);
  }

  TempFileManager::TempFileManager(const String& registry_id, UInt lock_timeout_ms) :
    TempFileManager(selectRegistrySlot_(registry_id, lock_timeout_ms), lock_timeout_ms)
  {
  }

  std::pair<String, std::unique_ptr<InterProcessFileLock>> TempFileManager::selectRegistrySlot_(const String& registry_id, UInt lock_timeout_ms)
  {
    std::unique_ptr<InterProcessFileLock> selected_lock;
    String selected_path = selectRegistryFilePath_(registry_id, lock_timeout_ms, selected_lock);
    if (selected_lock == nullptr || !selected_lock->isLocked())
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Could not acquire lock for temp file registry lockfile for registry id '" + registry_id + "'.");
    }
    return {selected_path, std::move(selected_lock)};
  }

  TempFileManager::TempFileManager(TempFileManager&& temp_file_manager) noexcept :
    registry_file_path_(temp_file_manager.registry_file_path_),
    lock_file_path_(temp_file_manager.lock_file_path_),
    lock_timeout_ms_(temp_file_manager.lock_timeout_ms_),
    registry_lock_(std::move(temp_file_manager.registry_lock_)),
    files_(std::move(temp_file_manager.files_))
  {
  }

  TempFileManager::~TempFileManager() noexcept
  {
    cleanup_();
  }

  TempFileManager& TempFileManager::operator=(TempFileManager&& temp_file_manager) noexcept
  {
    if (this == &temp_file_manager)
    {
      return *this;
    }

    // Remove currently owned resources before taking over new state
    cleanup_();

    const_cast<String&>(registry_file_path_) = std::move(temp_file_manager.registry_file_path_);
    const_cast<String&>(lock_file_path_) = std::move(temp_file_manager.lock_file_path_);
    lock_timeout_ms_ = temp_file_manager.lock_timeout_ms_;
    registry_lock_ = std::move(temp_file_manager.registry_lock_);
    files_ = std::move(temp_file_manager.files_);
    return *this;
  }

  void TempFileManager::addFile(const String& file_path)
  {
    std::lock_guard<std::mutex> state_lock(state_mutex_);

    // Validate input early to avoid storing unusable registry entries
    if (file_path.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "TempFileManager::addFile() received an empty file path.");
    }

    // Track file in memory and only persist newly inserted entries
    const auto insert_result = files_.insert(file_path);
    if (insert_result.second)
    {
      updateRegistryEntry_(file_path, RegistryAction::ADD);
    }
  }

  void TempFileManager::releaseFile(const String& file_path)
  {
    std::lock_guard<std::mutex> state_lock(state_mutex_);

    // Release management without deleting the file
    if (file_path.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "TempFileManager::releaseFile() received an empty file path.");
    }

    files_.erase(file_path);
    updateRegistryEntry_(file_path, RegistryAction::REMOVE);
  }

  bool TempFileManager::removeFileNow(const String& file_path) noexcept
  {
    std::lock_guard<std::mutex> state_lock(state_mutex_);

    try
    {
      if (File::exists(file_path) && !File::remove(file_path))
      {
        return false;
      }
      files_.erase(file_path);
      updateRegistryEntry_(file_path, RegistryAction::REMOVE);
      return true;
    }
    catch (...)
    {
      return false;
    }
  }

  void TempFileManager::cleanup_() noexcept
  {
    std::lock_guard<std::mutex> state_lock(state_mutex_);

    if (registry_file_path_.empty())
    {
      files_.clear();
      return;
    }

    try
    {
      cleanupRegistryFile_(registry_file_path_);
      files_.clear();
      registry_lock_.reset();
      File::remove(lock_file_path_);
      {
        std::lock_guard<std::mutex> lock(active_registry_paths_mutex);
        active_registry_paths.erase(registry_file_path_);
      }
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "TempFileManager cleanup encountered errors: "
                      << e.getName() << ": " << e.what() << "\n";
    }
    catch (const std::exception& e)
    {
      OPENMS_LOG_WARN << "TempFileManager cleanup encountered errors: "
                      << e.what() << "\n";
    }
  }

  String TempFileManager::sanitizeRegistryId_(const String& registry_id)
  {
    String sanitized = registry_id;

    // Sanitize the registry id by replacing invalid filename characters with underscores
    for (char& c : sanitized)
    {
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.'))
      {
        c = '_';
      }
    }

    // Fallback
    if (sanitized.empty())
    {
      return "openms_temp_file_registry";
    }

    return sanitized;
  }

  String TempFileManager::getRegistryFilePathForInstance_(const String& registry_id, Size instance_index)
  {
    String temp_dir = File::getTempDirectory();
    while (temp_dir.hasSuffix("/") || temp_dir.hasSuffix("\\"))
    {
      temp_dir = temp_dir.prefix(temp_dir.size() - 1);
    }

    String basename = sanitizeRegistryId_(registry_id);
    const String registry_filename = basename + "_" + String(instance_index) + ".list";

    if (temp_dir.empty())
    {
      return registry_filename;
    }

    return temp_dir + "/" + registry_filename;
  }

  String TempFileManager::getLockFilePathForInstance_(const String& registry_id, Size instance_index)
  {
    return getRegistryFilePathForInstance_(registry_id, instance_index) + ".lock";
  }

  void TempFileManager::cleanupRegistryFile_(const String& registry_file_path) noexcept
  {
    if (registry_file_path.empty())
    {
      return;
    }

    std::set<String> files_to_remove;
    std::ifstream registry_stream(registry_file_path.c_str());
    if (registry_stream.is_open())
    {
      std::string line;
      while (TextFile::getLine(registry_stream, line))
      {
        if (!line.empty())
        {
          files_to_remove.insert(line);
        }
      }
    }

    std::set<String> remaining_files = files_to_remove;
    for (const String& file : files_to_remove)
    {
      if (file.empty() || file == registry_file_path)
      {
        remaining_files.erase(file);
        continue;
      }
      bool removed = false;
      if (!File::exists(file))
      {
        removed = true;
      }
      else
      {
        removed = File::remove(file);
      }

      if (removed)
      {
        remaining_files.erase(file);
      }
    }

    if (remaining_files.empty())
    {
      File::remove(registry_file_path);
      return;
    }

    // Persist survivors conservatively so crashed cleanup can resume next run
    std::ofstream registry_out(registry_file_path.c_str(), std::ios::out | std::ios::trunc);
    if (!registry_out.is_open())
    {
      return;
    }
    for (const String& file : remaining_files)
    {
      registry_out << file << '\n';
    }
  }

  String TempFileManager::selectRegistryFilePath_(const String& registry_id, UInt lock_timeout_ms, std::unique_ptr<InterProcessFileLock>& registry_lock)
  {
    const Size max_instances = 15;

    // pass 1: cleanup stale registries (existing lock + acquirable lock)
    for (Size i = 1; i <= max_instances; ++i)
    {
      const String candidate_path = getRegistryFilePathForInstance_(registry_id, i);
      const String candidate_lock_path = getLockFilePathForInstance_(registry_id, i);
      // skip paths already owned by other TempFileManager instances in this process
      {
        std::lock_guard<std::mutex> lock(active_registry_paths_mutex);
        if (active_registry_paths.find(candidate_path) != active_registry_paths.end())
        {
          continue;
        }
      }
      const bool lock_exists = File::exists(candidate_lock_path);
      const bool registry_exists = File::exists(candidate_path);

      // If neither exists, there is nothing stale to clean up.
      if (!lock_exists && !registry_exists)
      {
        continue;
      }

      {
        InterProcessFileLock probe_lock(candidate_lock_path, lock_timeout_ms);
        if (!probe_lock.isLocked())
        {
          continue;
        }

        if (registry_exists)
        {
          cleanupRegistryFile_(candidate_path);
        }
      }

      if (!File::remove(candidate_lock_path))
      {
        OPENMS_LOG_WARN << "TempFileManager stale lock file cleanup failed for '"
                        << candidate_lock_path << "'.\n";
      }
    }

    // pass 2: pick lowest free slot by acquiring the slot lock first
    for (Size i = 1; i <= max_instances; ++i)
    {
      const String candidate_path = getRegistryFilePathForInstance_(registry_id, i);
      const String candidate_lock_path = getLockFilePathForInstance_(registry_id, i);

      {
        std::lock_guard<std::mutex> lock(active_registry_paths_mutex);
        if (active_registry_paths.find(candidate_path) != active_registry_paths.end())
        {
          continue;
        }
      }

      // Attempt to acquire the persistent inter-process lock for this slot
      std::unique_ptr<InterProcessFileLock> candidate_registry_lock =
        std::make_unique<InterProcessFileLock>(candidate_lock_path, lock_timeout_ms);
      if (!candidate_registry_lock->isLocked())
      {
        continue;
      }

      // We own the slot lock for this scope; now reserve registry path in-process
      {
        std::lock_guard<std::mutex> lock(active_registry_paths_mutex);
        if (active_registry_paths.find(candidate_path) != active_registry_paths.end())
        {
          continue;
        }

        // Reserve the slot immediately so another instance in the same process...
        // ...cannot pick the same index before constructor finishes
        std::ofstream registry_stream(candidate_path.c_str(), std::ios::out | std::ios::app);
        if (!registry_stream.is_open())
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, candidate_path);
        }
        registry_stream.close();
        active_registry_paths.insert(candidate_path);
      }

      registry_lock = std::move(candidate_registry_lock);
      return candidate_path;
    }

    throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   "Could not find a free temp registry slot in range 1..15 for registry id '" + registry_id + "'.");
  }

  void TempFileManager::loadRegistryFile_(std::set<String>& files) const
  {
    // Merge entries from the persistent registry into the provided set
    if (registry_file_path_.empty())
    {
      return;
    }

    std::ifstream registry_stream(registry_file_path_.c_str());
    if (!registry_stream.is_open())
    {
      // Missing registry is valid on first run; existing-but-unreadable is an error
      if (File::exists(registry_file_path_))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
      }
      return;
    }
    std::string line;
    while (TextFile::getLine(registry_stream, line))
    {
      // Ignore empty lines and add valid paths to the set
      if (!line.empty())
      {
        files.insert(line);
      }
    }

    // Distinguish parse/read errors from normal EOF
    if (registry_stream.bad() && !registry_stream.eof())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_,
                                  "Error while reading temporary file registry.");
    }
  }

  void TempFileManager::writeRegistryFile_(const std::set<String>& files) const
  {
    assert(!registry_file_path_.empty());

    // Write atomically via temporary file + rename to avoid truncation-loss on crashes
    const String tmp_registry_path = registry_file_path_ + ".tmp." + File::getUniqueName(false);
    std::ofstream registry_stream(tmp_registry_path.c_str(), std::ios::out | std::ios::trunc);
    if (!registry_stream.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tmp_registry_path);
    }

    for (const String& file : files)
    {
      registry_stream << file << '\n';
    }

    registry_stream.close();
    if (!File::rename(tmp_registry_path, registry_file_path_, true, false))
    {
      File::remove(tmp_registry_path);
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
    }
  }

  void TempFileManager::updateRegistryEntry_(const String& file_path, RegistryAction action)
  {
    if (registry_lock_ == nullptr || !registry_lock_->isLocked())
    {
      throw Exception::FailedAPICall(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "TempFileManager registry lock is not held for '" + registry_file_path_ + "'.");
    }

    if (action == RegistryAction::ADD)
    {
      std::ofstream registry_stream(registry_file_path_.c_str(), std::ios::out | std::ios::app);
      if (!registry_stream.is_open())
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
      }
      registry_stream << file_path << '\n';
      return;
    }

    std::set<String> persisted_files;
    loadRegistryFile_(persisted_files);

    if (action == RegistryAction::REMOVE)
    {
      persisted_files.erase(file_path);
      writeRegistryFile_(persisted_files);
    }
  }

} // namespace OpenMS
