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

#ifdef _OPENMP
  #include "omp.h"
#endif

// Platform-specific includes for OS version
#ifdef _WIN32
  #include <windows.h>
#elif defined(__APPLE__)
  #include <sys/sysctl.h>
#else
  #include <sys/utsname.h>
#endif

namespace OpenMS
{
  namespace Internal
  {

    enum OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX};
    std::string OpenMS_OSNames[] = {"unknown", "MacOS", "Windows", "Linux"};
    enum OpenMS_Architecture {ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT};
    std::string OpenMS_ArchNames[] = {"unknown", "32 bit", "64 bit"};

    class OPENMS_DLLAPI OpenMSOSInfo
    {
      OpenMS_OS os_;
      String os_version_;
      OpenMS_Architecture arch_;

    public:
      OpenMSOSInfo() :
          os_(OS_UNKNOWN),
          os_version_("unknown"),
          arch_(ARCH_UNKNOWN)
      {}

      /// @brief Get the current operating system (Windows, MacOS, Linux)
      String getOSAsString() const
      {
        return OpenMS_OSNames[os_];
      }

      /// @brief Get the current architecture (32-bit or 64-bit)
      String getArchAsString() const
      {
        return OpenMS_ArchNames[arch_];
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
            return OpenMS_ArchNames[ARCH_32BIT];
          case 8:
            return OpenMS_ArchNames[ARCH_64BIT];
          default:
            return OpenMS_ArchNames[ARCH_UNKNOWN];
        }
      }

      /// @brief Obtain a list of SIMD extensions which are currently in use (i.e. used by the compiler during optimization, as well as for SIMDe code within OpenMS)
      static String getActiveSIMDExtensions();

      /// @brief Constructs and returns an OpenMSOSInfo object
      static OpenMSOSInfo getOSInfo()
      {
        OpenMSOSInfo info;
        #if defined(WIN32)  // Windows
        info.os_ = OS_WINDOWS;
        #elif (defined(__MACH__) && defined(__APPLE__)) // MacOS
        info.os_ = OS_MACOS;
        #elif (defined(__unix__)) //Linux/FreeBSD TODO make a difference?
        info.os_ = OS_LINUX;
        #endif // else stays unknown

        // Get OS version using platform-specific APIs
        #ifdef _WIN32
        // Windows version - use RtlGetVersion for accurate Windows 10+ reporting
        typedef LONG (WINAPI* RtlGetVersionPtr)(PRTL_OSVERSIONINFOW);
        HMODULE hNtdll = GetModuleHandleA("ntdll.dll");
        RtlGetVersionPtr RtlGetVersion = nullptr;
        bool version_obtained = false;

        if (hNtdll)
        {
          RtlGetVersion = (RtlGetVersionPtr)GetProcAddress(hNtdll, "RtlGetVersion");
          if (RtlGetVersion)
          {
            RTL_OSVERSIONINFOW osvi;
            ZeroMemory(&osvi, sizeof(RTL_OSVERSIONINFOW));
            osvi.dwOSVersionInfoSize = sizeof(RTL_OSVERSIONINFOW);
            if (RtlGetVersion(&osvi) == 0) // STATUS_SUCCESS
            {
              info.os_version_ = std::to_string(osvi.dwMajorVersion) + "." + std::to_string(osvi.dwMinorVersion);
              version_obtained = true;
            }
          }
        }

        // Fallback to GetVersionEx if RtlGetVersion unavailable
        if (!version_obtained)
        {
          OSVERSIONINFOEX osvi;
          ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));
          osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
          #pragma warning(push)
          #pragma warning(disable: 4996)
          if (GetVersionEx((OSVERSIONINFO*)&osvi))
          {
            info.os_version_ = std::to_string(osvi.dwMajorVersion) + "." + std::to_string(osvi.dwMinorVersion);
          }
          else
          {
            info.os_version_ = "unknown";
          }
          #pragma warning(pop)
        }
        #elif defined(__APPLE__)
        // macOS version
        char version[256];
        size_t size = sizeof(version);
        if (sysctlbyname("kern.osproductversion", version, &size, NULL, 0) == 0)
        {
          info.os_version_ = version;
        }
        else
        {
          info.os_version_ = "unknown";
        }
        #else
        // Linux version
        struct utsname uts;
        if (uname(&uts) == 0)
        {
          info.os_version_ = uts.release;
        }
        else
        {
          info.os_version_ = "unknown";
        }
        #endif

        // Identify architecture using sizeof(size_t)
        size_t word_size = sizeof(size_t);
        if (word_size == 4)
        {
          info.arch_ = ARCH_32BIT;
        }
        else if (word_size == 8)
        {
          info.arch_ = ARCH_64BIT;
        }
        else
        {
          info.arch_ = ARCH_UNKNOWN;
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
