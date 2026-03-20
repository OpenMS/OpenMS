// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/build_config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstddef>
#include <string>
#include <fstream>
#include <cstring>

#ifdef _WIN32
  #include <windows.h>
  #include <winternl.h>
#elif defined(__APPLE__)
  #include <sys/types.h>
  #include <sys/sysctl.h>
#endif

#ifdef _OPENMP
  #include "omp.h"
#endif

namespace OpenMS
{
  namespace Internal
  {

    enum class OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX, SIZE_OF_OPENMS_OS};
    inline const std::string OpenMS_OSNames[] = {"unknown", "MacOS", "Windows", "Linux"};
    enum class OpenMS_Architecture {ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT, SIZE_OF_OPENMS_ARCHITECTURE};
    inline const std::string OpenMS_ArchNames[] = {"unknown", "32 bit", "64 bit"};

    /// @brief Helper: get OS version string using platform APIs
    inline std::string getOSVersionString_()
    {
#if defined(_WIN32)
      // Use RtlGetVersion to get the real version (not compatibility-shimmed)
      typedef NTSTATUS (WINAPI* RtlGetVersionPtr)(PRTL_OSVERSIONINFOW);
      HMODULE hMod = GetModuleHandleW(L"ntdll.dll");
      if (hMod)
      {
        auto fxPtr = reinterpret_cast<RtlGetVersionPtr>(GetProcAddress(hMod, "RtlGetVersion"));
        if (fxPtr)
        {
          RTL_OSVERSIONINFOW rovi{};
          rovi.dwOSVersionInfoSize = sizeof(rovi);
          if (fxPtr(&rovi) == 0)
          {
            return std::to_string(rovi.dwMajorVersion) + "." + std::to_string(rovi.dwMinorVersion)
                   + "." + std::to_string(rovi.dwBuildNumber);
          }
        }
      }
      return "unknown";
#elif defined(__APPLE__)
      char version[64] = {};
      size_t len = sizeof(version);
      if (sysctlbyname("kern.osproductversion", version, &len, nullptr, 0) == 0)
      {
        return std::string(version);
      }
      return "unknown";
#elif defined(__unix__)
      // Read VERSION_ID from /etc/os-release
      std::ifstream ifs("/etc/os-release");
      std::string line;
      while (std::getline(ifs, line))
      {
        if (line.compare(0, 11, "VERSION_ID=") == 0)
        {
          std::string val = line.substr(11);
          // Strip surrounding quotes if present
          if (val.size() >= 2 && val.front() == '"' && val.back() == '"')
          {
            val = val.substr(1, val.size() - 2);
          }
          return val;
        }
      }
      return "unknown";
#else
      return "unknown";
#endif
    }

    class OPENMS_DLLAPI OpenMSOSInfo
    {
      OpenMS_OS os_;
      String os_version_;
      OpenMS_Architecture arch_;

    public:
      OpenMSOSInfo() :
          os_(OpenMS_OS::OS_UNKNOWN),
          os_version_("unknown"),
          arch_(OpenMS_Architecture::ARCH_UNKNOWN)
      {}

      /// @brief Get the current operating system (Windows, MacOS, Linux)
      String getOSAsString() const
      {
        return OpenMS_OSNames[static_cast<size_t>(os_)];
      }

      /// @brief Get the current architecture (32-bit or 64-bit)
      String getArchAsString() const
      {
        return OpenMS_ArchNames[static_cast<size_t>(arch_)];
      }

      /// @brief Get the OS version (e.g. 10.15 for macOS or 10 for Windows)
      String getOSVersionAsString() const
      {
        return os_version_;
      }

      /// @brief Get Architecture of this binary (simply by looking at size of a pointer, i.e. size_t).
      static String getBinaryArchitecture()
      {
        size_t bytes = sizeof(size_t);
        switch (bytes)
        {
          case 4:
            return OpenMS_ArchNames[static_cast<size_t>(OpenMS_Architecture::ARCH_32BIT)];
          case 8:
            return OpenMS_ArchNames[static_cast<size_t>(OpenMS_Architecture::ARCH_64BIT)];
          default:
            return OpenMS_ArchNames[static_cast<size_t>(OpenMS_Architecture::ARCH_UNKNOWN)];
        }
      }

      /// @brief Obtain a list of SIMD extensions which are currently in use (i.e. used by the compiler during optimization, as well as for SIMDe code within OpenMS)
      static String getActiveSIMDExtensions();

      /// @brief Constructs and returns an OpenMSOSInfo object
      static OpenMSOSInfo getOSInfo()
      {
        OpenMSOSInfo info;
        #if defined(WIN32)  // Windows
        info.os_ = OpenMS_OS::OS_WINDOWS;
        #elif (defined(__MACH__) && defined(__APPLE__)) // MacOS
        info.os_ = OpenMS_OS::OS_MACOS;
        #elif (defined(__unix__)) //Linux/FreeBSD TODO make a difference?
        info.os_ = OpenMS_OS::OS_LINUX;
        #endif // else stays unknown

        // Get OS version using platform APIs
        info.os_version_ = getOSVersionString_();

        // identify architecture by pointer size
        if (sizeof(void*) == 4)
        {
          info.arch_ = OpenMS_Architecture::ARCH_32BIT;
        }
        else
        {
          info.arch_ = OpenMS_Architecture::ARCH_64BIT;
        }

        return info;
      }
    };

    /// @brief Struct with some static methods to get informations on the build configuration
    struct OpenMSBuildInfo

    {
    public:

      /// @brief Checks if OpenMP was enabled during build, based on the _OPENMP macro
      static bool isOpenMPEnabled()
      {
        #ifdef _OPENMP
        return true;
        #else
        return false;
        #endif
      }

      /// @brief Get the build type used during building the OpenMS library
      static String getBuildType()
      {
        return OPENMS_BUILD_TYPE;
      }

      /// @brief Get the maximum number of threads that OpenMP will use (including hyperthreads)
      /// Note: This could also be limited by the OMP_NUM_THREADS environment variable
      /// Returns 1 if OpenMP was disabled.
      static Size getOpenMPMaxNumThreads()
      {
        #ifdef _OPENMP
        return omp_get_max_threads();
        #else
        return 1;
        #endif
      }
      /// @brief Set the number of threads that OpenMP will use (including hyperthreads)
      /// Note: Can be initialized by the OMP_NUM_THREADS environment variable. This function can overwrite this at runtime.
      static void setOpenMPNumThreads(Int num_threads)
      {
        #ifdef _OPENMP
        omp_set_num_threads(num_threads);
        #endif
        (void)num_threads; // avoid 'unreferenced formal parameter' C4100 on Windows
      }
    };

  } // NS Internal
} // NS OpenMS
