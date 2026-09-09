// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/TempFiles.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SystemSettings.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h>
#include <rpc.h>
#pragma comment(lib, "Rpcrt4.lib")
#else
#include <unistd.h> // for mkdtemp on macOS
#endif

namespace OpenMS
{
namespace
{
  std::string createUniqueDir_(const std::string& prefix)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    // Ensure the parent directory chain exists first (old fs::create_directories
    // behavior created the whole chain; mkdtemp/CreateDirectoryW only create the leaf).
    std::filesystem::path parent_path = to_path(prefix).parent_path();
    if (!parent_path.empty())
    {
      std::error_code ec;
      std::filesystem::create_directories(parent_path, ec);
      if (ec)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prefix, ec.message());
      }
    }

    for (int attempt = 0; attempt < 100; ++attempt)
    {
      UUID uuid;
      RPC_STATUS create_status = UuidCreate(&uuid);
      if (create_status != RPC_S_OK && create_status != RPC_S_UUID_LOCAL_ONLY)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prefix, "UuidCreate failed with status " + std::to_string(create_status));
      }

      RPC_WSTR wstr = nullptr;
      RPC_STATUS str_status = UuidToStringW(&uuid, &wstr);
      if (str_status != RPC_S_OK || wstr == nullptr)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prefix, "UuidToStringW failed with status " + std::to_string(str_status));
      }
      std::wstring wsuffix(reinterpret_cast<wchar_t*>(wstr));
      RpcStringFreeW(&wstr);
      std::string suffix;
      suffix.reserve(wsuffix.size());
      for (wchar_t wc : wsuffix)
      {
        suffix.push_back(static_cast<char>(wc));
      }

      std::string candidate = prefix + "_" + suffix;
      std::filesystem::path wcandidate = to_path(candidate);

      if (CreateDirectoryW(wcandidate.native().c_str(), NULL))
      {
        return candidate + "/";
      }
      DWORD err = GetLastError();
      if (err != ERROR_ALREADY_EXISTS)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, candidate, "GetLastError() = " + std::to_string(err));
      }
      // else: collision, loop and try again with a new UUID
    }
    throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prefix, "exceeded 100 attempts");
#else
    // Ensure the parent directory chain exists first (same reasoning as above).
    std::filesystem::path parent_path = to_path(prefix).parent_path();
    if (!parent_path.empty())
    {
      std::error_code ec;
      std::filesystem::create_directories(parent_path, ec);
      if (ec)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prefix, ec.message());
      }
    }

    std::string tmpl = prefix + "_XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (::mkdtemp(buf.data()) == nullptr)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, tmpl, std::strerror(errno));
    }
    return std::string(buf.data()) + "/";
#endif
  }
}

  TempDir::TempDir(bool keep_dir)
    : keep_dir_(keep_dir)
  {
    std::string prefix = SystemSettings::getTempDirectory() + "/" +File::getUniqueName();
    temp_dir_ = createUniqueDir_(prefix);
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'\n";
  };

  TempDir::TempDir(const std::string& base_dir, bool keep_dir)
    : keep_dir_(keep_dir)
  {
    std::string prefix = base_dir;
    if (!prefix.empty() && !StringUtils::hasSuffix(prefix,"/"))
    {
      prefix += "/";
    }
    prefix += "OpenMSTempDir_" + File::getUniqueName();
    temp_dir_ = createUniqueDir_(prefix);
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'\n";
  };

  TempDir::~TempDir()
  {
    if (keep_dir_)
    {
      OPENMS_LOG_DEBUG << "Keeping temporary files in directory '" << temp_dir_ << '\n';
      return;
    }

    File::removeDirRecursively(temp_dir_);
  };

  const std::string& TempDir::getPath() const
  {
    return temp_dir_;
  }

  std::string TempFiles::getTemporaryFile(const std::string& alternative_file)
  {
    // take no action
    if (!alternative_file.empty())
    {
      return alternative_file;
    }
    // create temporary (and schedule for deletion)
    return registry_.newFile();
  }


  TempFiles::Registry_::Registry_()
    : filenames_()
  {
  }

  std::string TempFiles::Registry_::newFile()
  {
    std::string s = SystemSettings::getTempDirectory(); StringUtils::ensureLastChar(s, '/'); s += File::getUniqueName();
    std::lock_guard<std::mutex> _(mtx_);
    filenames_.push_back(s);
    // do NOT return filenames_.back() by ref, since another thread might resize the vector and invalidate the reference!
    return s; // uses RVO, so its efficient
  }

  TempFiles::Registry_::~Registry_()
  {
    std::lock_guard<std::mutex> _(mtx_);
    for (const std::string& filename : filenames_)
    {
      if (File::exists(filename) && !File::remove(filename))
      {
        std::cerr << "Warning: unable to remove temporary file '" << filename << "'" << std::endl;
      }
    }
  }

  TempFiles::Registry_ TempFiles::registry_;

  // --------------------------------------------------------------------------
  // Deprecated File forwarder (kept for one release; see TempFiles)
  // --------------------------------------------------------------------------
  std::string File::getTemporaryFile(const std::string& alternative_file)
  {
    return TempFiles::getTemporaryFile(alternative_file);
  }

} // namespace OpenMS
