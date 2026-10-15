// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/VersionInfo.h>

#include <filesystem>
#include <initializer_list>
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
    @brief Path of the loaded module whose file name is @p file_name, empty when the loader reports none.

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
    report the executable instead of the library. The names are not guessed here
    either -- relocated/CMakeLists.txt passes what $<TARGET_FILE_NAME:...> of the
    imported targets says, so a SOVERSION, an OUTPUT_NAME or a debug postfix on the
    OpenMS libraries changes what is looked for along with what is built.
  */
  std::string loadedModulePath(const std::string& file_name);

#if defined(_WIN32)
  std::string loadedModulePath(const std::string& file_name)
  {
    const HMODULE module = GetModuleHandleA(file_name.c_str());
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
  std::string loadedModulePath(const std::string& file_name)
  {
    for (uint32_t i = 0; i < _dyld_image_count(); ++i)
    {
      const char* name = _dyld_get_image_name(i);
      if (name != nullptr && std::filesystem::path(name).filename().string() == file_name)
      {
        return name;
      }
    }
    return {};
  }
#else
  /// What collectModule() is looking for, and where it puts the answer.
  struct ModuleQuery
  {
    const std::string* file_name = nullptr;
    std::string path;
  };

  /// Callback of the dl_iterate_phdr() walk below; stores the first match in @p data.
  int collectModule(struct dl_phdr_info* info, size_t /* size */, void* data)
  {
    ModuleQuery& query = *static_cast<ModuleQuery*>(data);
    if (info->dlpi_name == nullptr
        || std::filesystem::path(info->dlpi_name).filename().string() != *query.file_name)
    {
      return 0;
    }
    query.path = info->dlpi_name;
    return 1; // anything but 0 ends the walk
  }

  std::string loadedModulePath(const std::string& file_name)
  {
    ModuleQuery query;
    query.file_name = &file_name;
    dl_iterate_phdr(collectModule, &query);
    return query.path;
  }
#endif

  /// @p path resolved as far as it exists; empty when it cannot be resolved.
  std::filesystem::path resolvedPath(const std::filesystem::path& path)
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
  // libraries installed there.
  const std::string version = VersionInfo::getVersion();

  const std::filesystem::path expected_dir = resolvedPath(OPENMS_EXPECTED_LIB_DIR);
  if (expected_dir.empty())
  {
    std::cerr << "cannot resolve the library directory the package reports: " << OPENMS_EXPECTED_LIB_DIR << "\n";
    return 1;
  }

  // The libraries that were loaded have to be the ones the package describes, not
  // another OpenMS that happened to be earlier on the loader's search path. Both
  // libraries of the core layer, because the executable links OpenMS::OpenMS and
  // gets OpenMS::OpenSwathAlgo with it: each is resolved by the loader in its own
  // right, so an RPATH that sends one of them back to the build tree is exactly
  // the regression this test exists to catch.
  for (const std::string& module_name : {std::string(OPENMS_CORE_MODULE_NAME),
                                         std::string(OPENMS_OPENSWATHALGO_MODULE_NAME)})
  {
    const std::string module = loadedModulePath(module_name);
    if (module.empty())
    {
      std::cerr << "the loader does not report a module named " << module_name << " in this process\n";
      return 1;
    }
    const std::filesystem::path module_dir = resolvedPath(std::filesystem::path(module).parent_path());
    if (module_dir.empty())
    {
      std::cerr << "cannot resolve the directory of " << module << "\n";
      return 1;
    }
    if (module_dir != expected_dir)
    {
      std::cerr << "loaded " << module << ", but the package reports its libraries in "
                << OPENMS_EXPECTED_LIB_DIR << "\n";
      return 1;
    }
    std::cout << "Loaded " << module << std::endl;
  }

  // And the library is the version the package reports. Compared through the parsed
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

  std::cout << "OpenMS " << version << ": all good and well!" << std::endl;
  return 0;
}
