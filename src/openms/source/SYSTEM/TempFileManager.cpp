// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Tilman Aurich $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/TempFileManager.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/InterProcessFileLock.h>

#include <cctype>
#include <exception>
#include <fstream>
#include <string>
#include <utility>

namespace OpenMS
{

  TempFileManager::TempFileManager(const String& registry_id) :
    // Use a per-usage registry path derived from registry id
    registry_file_path_(getRegistryFilePath_(registry_id))
  {
    // Restore persisted entries to support crash-recovery cleanup
    InterProcessFileLock lock(registry_file_path_);
    loadRegistryFile_(files_);
  }

  TempFileManager::TempFileManager(TempFileManager&& temp_file_manager) noexcept :
    registry_file_path_(std::move(temp_file_manager.registry_file_path_)),
    files_(std::move(temp_file_manager.files_))
  {
    // Leave moved-from object in inert state
    temp_file_manager.registry_file_path_.clear();
    temp_file_manager.files_.clear();
  }

  TempFileManager::~TempFileManager() noexcept
  {
    try
    {
      // Destructors must not throw
      cleanupNow();
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "TempFileManager cleanup failed during destruction: "
                      << e.getName() << ": " << e.what() << "\n";
    }
    // Handle any other exceptions to avoid termination, log them as well
    catch (const std::exception& e)
    {
      OPENMS_LOG_WARN << "TempFileManager cleanup failed during destruction: "
                      << e.what() << "\n";
    }
  }

  TempFileManager& TempFileManager::operator=(TempFileManager&& temp_file_manager) noexcept
  {
    if (this == &temp_file_manager)
    {
      return *this;
    }

    // Remove currently owned resources before taking over new state
    cleanupNow();

    registry_file_path_ = std::move(temp_file_manager.registry_file_path_);
    files_ = std::move(temp_file_manager.files_);
    temp_file_manager.registry_file_path_.clear();
    temp_file_manager.files_.clear();
    return *this;
  }

  void TempFileManager::addFile(const String& file_path)
  {
    // Validate input early to avoid storing unusable registry entries
    if (file_path.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "TempFileManager::addFile() received an empty file path.");
    }

    // Track file in memory ...
    files_.insert(file_path);
    // ... and persist it for crash-recovery cleanup
    updateRegistryEntry_(file_path, true);
  }

  void TempFileManager::releaseFile(const String& file_path)
  {
    // Release management without deleting the file
    if (file_path.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "TempFileManager::releaseFile() received an empty file path.");
    }

    files_.erase(file_path);
    updateRegistryEntry_(file_path, false);
  }

  bool TempFileManager::removeFileNow(const String& file_path) noexcept
  {
    // Indicate success via return value
    if (file_path.empty())
    {
      return false;
    }

    bool success = true;
    try
    {
      if (File::exists(file_path) && !File::remove(file_path))
      {
        success = false;
      }
    }
    catch (const Exception::BaseException&)
    {
      success = false;
    }
    catch (const std::exception&)
    {
      success = false;
    }

    try
    {
      files_.erase(file_path);
      updateRegistryEntry_(file_path, false);
    }
    catch (const Exception::BaseException&)
    {
      success = false;
    }
    catch (const std::exception&)
    {
      success = false;
    }

    return success;
  }

  void TempFileManager::cleanupNow() noexcept
  {
    try
    {
      // Avoid unnecessary work if there is no registry to clean up
      if (registry_file_path_.empty())
      {
        files_.clear();
        return;
      }

      // lock the registry to prevent concurrent modifications while we clean up
      InterProcessFileLock lock(registry_file_path_);

      // Merge in-memory and persisted entries
      std::set<String> files_to_remove = files_;
      loadRegistryFile_(files_to_remove);

      for (const String& file : files_to_remove)
      {
        if (file.empty() || file == registry_file_path_)
        {
          continue;
        }

        File::remove(file);
      }

      files_.clear();
      writeRegistryFile_(std::set<String>());
      File::remove(registry_file_path_);
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
    // If the registry id is empty, return a default name to avoid creating files with empty names
    if (registry_id.empty())
    {
      return "openms_temp_file_registry";
    }

    // Sanitize the registry id by replacing invalid filename characters with underscores
    String sanitized;
    sanitized.reserve(registry_id.size());
    for (char c : registry_id)
    {
      if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')
      {
        sanitized.push_back(c);
      }
      else
      {
        sanitized.push_back('_');
      }
    }

    if (sanitized.empty())
    {
      return "openms_temp_file_registry";
    }

    return sanitized;
  }

  String TempFileManager::getRegistryFilePath_(const String& registry_id)
  {
    // Resolve the OS temp directory and normalize trailing separators
    String temp_dir = File::getTempDirectory();
    while (temp_dir.hasSuffix("/") || temp_dir.hasSuffix("\\"))
    {
      temp_dir = temp_dir.prefix(temp_dir.size() - 1);
    }

    // Build file name from caller-provided registry id
    String basename = sanitizeRegistryId_(registry_id);
    const String registry_filename = basename + ".list";

    if (temp_dir.empty())
    {
      return registry_filename;
    }

    return temp_dir + "/" + registry_filename;
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
      // Missing registry is valid on first run; existing-but-unreadable is an error.
      if (File::exists(registry_file_path_))
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
      }
      return;
    }

    String line;
    while (std::getline(registry_stream, line))
    {
      // Normalize CRLF line endings when files are edited or shared across platforms
      if (!line.empty() && line.back() == '\r')
      {
        line.pop_back();
      }
      // Ignore empty lines and add valid paths to the set
      if (!line.empty())
      {
        files.insert(line);
      }
    }

    // Distinguish parse/read errors from normal EOF
    if (registry_stream.bad())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_,
                                  "Error while reading temporary file registry.");
    }
  }

  void TempFileManager::writeRegistryFile_(const std::set<String>& files) const
  {
    if (registry_file_path_.empty())
    {
      return;
    }

    // Write the current set of files to the registry, overwriting any existing content
    std::ofstream registry_stream(registry_file_path_.c_str(), std::ios::out | std::ios::trunc);
    if (!registry_stream.is_open())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
    }

    for (const String& file : files)
    {
      registry_stream << file << '\n';
    }

    if (!registry_stream.good())
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, registry_file_path_);
    }
  }

  void TempFileManager::updateRegistryEntry_(const String& file_path, bool add_file)
  {
    InterProcessFileLock lock(registry_file_path_);

    std::set<String> persisted_files;
    loadRegistryFile_(persisted_files);

    if (add_file)
    {
      persisted_files.insert(file_path);
    }
    else
    {
      persisted_files.erase(file_path);
    }

    writeRegistryFile_(persisted_files);
  }

} // namespace OpenMS
