// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/VersionInfo.h>

#include <filesystem>
#include <iostream>
#include <string>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <link.h>
#endif

// Consumer program of an OpenMS installation whose prefix was renamed after
// `cmake --install` (relocated/CMakeLists.txt). Deliberately free of the OpenMS
// shared data: the installation this is built against leaves out the share
// component, so nothing here may reach for <prefix>/share/OpenMS.

using namespace OpenMS;

namespace
{
  /**
    @brief Path of the OpenMS core library this process has loaded, empty when the loader reports none.

    Asked of the loader rather than of the build system, because which library is
    resolved at run time is the question here: the RPATH of this executable is only
    one of the places the loader looks, and a same-version OpenMS elsewhere on the
    search path would otherwise satisfy this program and hide a broken relocation.
    That is not hypothetical -- the install RPATH of a macOS build names the build
    tree of this very build (CMAKE_INSTALL_RPATH in the top-level CMakeLists.txt),
    which holds an OpenMS of exactly this version.

    Identified by module name rather than by the address of an exported function:
    for a function defined in a shared library, taking its address in an executable
    yields the address of the PLT stub inside the executable, so dladdr() would
    report the executable instead of the library.
  */
  std::string loadedOpenMSLibrary();

#if !defined(_WIN32)
  /**
    @brief Whether @p name is a path to the OpenMS core library.

    Matched on everything before the first '.' of the file name, so that the
    versioned forms ("libOpenMS.so.3.6") count and the libraries of the other
    layers ("libOpenMS_CLI.so") do not.
  */
  bool isOpenMSCoreLibrary(const char* name)
  {
    if (name == nullptr)
    {
      return false;
    }
    const std::string file = std::filesystem::path(name).filename().string();
    return file.substr(0, file.find('.')) == "libOpenMS";
  }
#endif

#if defined(_WIN32)
  std::string loadedOpenMSLibrary()
  {
    // "OpenMS" for an MSVC build, "libOpenMS" for the MinGW naming convention
    HMODULE module = GetModuleHandleA("OpenMS");
    if (module == nullptr)
    {
      module = GetModuleHandleA("libOpenMS");
    }
    if (module == nullptr)
    {
      return {};
    }
    char path[MAX_PATH] = {};
    const DWORD length = GetModuleFileNameA(module, path, MAX_PATH);
    if (length == 0 || length >= MAX_PATH)
    {
      return {};
    }
    return std::string(path, length);
  }
#elif defined(__APPLE__)
  std::string loadedOpenMSLibrary()
  {
    for (uint32_t i = 0; i < _dyld_image_count(); ++i)
    {
      const char* name = _dyld_get_image_name(i);
      if (isOpenMSCoreLibrary(name))
      {
        return name;
      }
    }
    return {};
  }
#else
  /// Callback of the dl_iterate_phdr() walk below; stores the first match in @p data.
  int collectOpenMSLibrary(struct dl_phdr_info* info, size_t /* size */, void* data)
  {
    if (!isOpenMSCoreLibrary(info->dlpi_name))
    {
      return 0;
    }
    *static_cast<std::string*>(data) = info->dlpi_name;
    return 1; // anything but 0 ends the walk
  }

  std::string loadedOpenMSLibrary()
  {
    std::string path;
    dl_iterate_phdr(collectOpenMSLibrary, &path);
    return path;
  }
#endif

  /// Directory of @p path, resolved as far as it exists; empty when it cannot be resolved.
  std::filesystem::path resolvedDirectory(const std::filesystem::path& path)
  {
    std::error_code ec;
    const std::filesystem::path resolved = std::filesystem::weakly_canonical(path, ec);
    return ec ? std::filesystem::path() : resolved;
  }
}

/// Checks that this program runs on the OpenMS installation the package describes.
int main()
{
  // Compiled against the headers of the moved prefix and linked against the
  // libOpenMS installed there.
  const std::string version = VersionInfo::getVersion();

  // The library that was loaded has to be the one the package describes, not
  // another OpenMS that happened to be earlier on the loader's search path.
  // Both halves of that: it comes out of the directory the package reports ...
  const std::string module = loadedOpenMSLibrary();
  if (module.empty())
  {
    std::cerr << "cannot determine which OpenMS library this process loaded\n";
    return 1;
  }
  const std::filesystem::path module_dir = resolvedDirectory(std::filesystem::path(module).parent_path());
  const std::filesystem::path expected_dir = resolvedDirectory(OPENMS_EXPECTED_LIB_DIR);
  if (module_dir.empty() || expected_dir.empty())
  {
    std::cerr << "cannot resolve the directory of " << module << " or of " << OPENMS_EXPECTED_LIB_DIR << "\n";
    return 1;
  }
  if (module_dir != expected_dir)
  {
    std::cerr << "loaded the OpenMS library " << module << ", but the package reports its libraries in "
              << OPENMS_EXPECTED_LIB_DIR << "\n";
    return 1;
  }
  std::cout << "Loaded OpenMS " << version << " from " << module << std::endl;

  // ... and it is the version the package reports. Compared through the parsed
  // struct rather than as a string: a pre-release build appends "-pre-<identifier>"
  // to the library's version, while the package reports the plain major.minor.patch
  // in OPENMS_EXPECTED_VERSION.
  const VersionInfo::VersionDetails expected = VersionInfo::VersionDetails::create(OPENMS_EXPECTED_VERSION);
  if (expected == VersionInfo::VersionDetails::EMPTY)
  {
    std::cerr << "cannot parse the version the package reported: " << OPENMS_EXPECTED_VERSION << "\n";
    return 1;
  }
  const VersionInfo::VersionDetails loaded = VersionInfo::getVersionStruct();
  if (loaded.version_major != expected.version_major
      || loaded.version_minor != expected.version_minor
      || loaded.version_patch != expected.version_patch)
  {
    std::cerr << "the package reports OpenMS " << OPENMS_EXPECTED_VERSION
              << ", but the library loaded from it is " << version << "\n";
    return 1;
  }

  std::cout << "All good and well!" << std::endl;
  return 0;
}
