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

    /// Operating systems @c OpenMSOSInfo can report. @c SIZE_OF_OPENMS_OS is a sentinel for array sizing.
    enum class OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX, SIZE_OF_OPENMS_OS};
    /// String labels paired with @ref OpenMS_OS; index by @c static_cast<size_t>(enum_value).
    inline const std::string OpenMS_OSNames[] = {"unknown", "MacOS", "Windows", "Linux"};
    /// Pointer-width architectures @c OpenMSOSInfo can report. @c SIZE_OF_OPENMS_ARCHITECTURE is a sentinel for array sizing.
    enum class OpenMS_Architecture {ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT, SIZE_OF_OPENMS_ARCHITECTURE};
    /// String labels paired with @ref OpenMS_Architecture; index by @c static_cast<size_t>(enum_value).
    inline const std::string OpenMS_ArchNames[] = {"unknown", "32 bit", "64 bit"};

    /**
      @brief Return the running OS version using a platform-specific probe.

      Windows uses @c RtlGetVersion from @c ntdll.dll to bypass the compatibility shim
      and return the real version as @c "MAJOR.MINOR.BUILD".
      macOS uses @c sysctlbyname("kern.osproductversion", ...).
      Linux (and other @c __unix__) parses @c /etc/os-release and returns the
      @c VERSION_ID value with surrounding quotes stripped.

      Returns @c "unknown" if the platform query fails or the platform is unsupported.
    */
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

    /**
      @brief Snapshot of the runtime operating system: name, version, and pointer-width architecture.

      Populated via the static factory @ref getOSInfo, which probes the running platform
      with preprocessor defines (@c WIN32, @c __MACH__ @c && @c __APPLE__, @c __unix__).
      All accessors return string labels suitable for logging or display; the underlying
      enum values are kept private and only exposed via the @c *AsString accessors.

      Distinct from the build-time view of the binary: @ref getBinaryArchitecture and
      @ref getActiveSIMDExtensions return what the compiler emitted, not what the
      running kernel reports.

      @ingroup System
    */
    class OPENMS_DLLAPI OpenMSOSInfo
    {
      OpenMS_OS os_;
      String os_version_;
      OpenMS_Architecture arch_;

    public:
      /// Construct an empty instance with all fields set to "unknown". Use @ref getOSInfo to probe the runtime platform.
      OpenMSOSInfo() :
          os_(OpenMS_OS::OS_UNKNOWN),
          os_version_("unknown"),
          arch_(OpenMS_Architecture::ARCH_UNKNOWN)
      {}

      /// Return the running operating system name (@c "Windows", @c "MacOS", @c "Linux", or @c "unknown").
      String getOSAsString() const
      {
        return OpenMS_OSNames[static_cast<size_t>(os_)];
      }

      /// Return the running architecture (@c "32 bit", @c "64 bit", or @c "unknown") as detected by @ref getOSInfo from the pointer width.
      String getArchAsString() const
      {
        return OpenMS_ArchNames[static_cast<size_t>(arch_)];
      }

      /// Return the OS version captured by @ref getOSInfo (e.g. @c "10.15" on macOS or @c "10.0.19045" on Windows); @c "unknown" if the platform query failed.
      String getOSVersionAsString() const
      {
        return os_version_;
      }

      /**
        @brief Return the pointer width of this binary, derived from @c sizeof(size_t).

        @c "32 bit" if @c sizeof(size_t) is 4, @c "64 bit" if 8, otherwise @c "unknown".
        This reflects how the binary was compiled, independent of the running OS.
      */
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

      /**
        @brief Return a comma-separated list of SIMD instruction sets compiled into this binary.

        Reflects the SIMDe preprocessor probes that were defined at build time: any of
        @c neon, @c SSE, @c SSE2, @c SSE3, @c SSE4.1, @c SSE4.2, @c AVX, @c AVX2, @c FMA.
        Returns an empty string if none of those is set.
      */
      static String getActiveSIMDExtensions();

      /**
        @brief Probe the running platform and return a populated @ref OpenMSOSInfo.

        Sets @c os_ from compile-time platform macros (@c WIN32 / @c __MACH__ @c && @c __APPLE__ /
        @c __unix__), @c os_version_ via @ref getOSVersionString_, and @c arch_ from
        @c sizeof(void*) (@c 4 → 32-bit, otherwise 64-bit). Fields stay at their
        @c OS_UNKNOWN / @c ARCH_UNKNOWN defaults if the running platform isn't recognised.
      */
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

    /**
      @brief Build-configuration accessors: build type, OpenMP availability, and OpenMP thread count.

      Static utility methods only; the struct carries no state. All values reflect how the
      OpenMS library was compiled, not the running platform — see @ref OpenMSOSInfo for the
      runtime view.

      @ingroup System
    */
    struct OpenMSBuildInfo

    {
    public:

      /// Return @c true when the library was compiled with OpenMP support (@c _OPENMP defined), @c false otherwise.
      static bool isOpenMPEnabled()
      {
        #ifdef _OPENMP
        return true;
        #else
        return false;
        #endif
      }

      /// Return the CMake build type baked into the library (e.g. @c "Release", @c "Debug"), captured from the @c OPENMS_BUILD_TYPE configure-time macro.
      static String getBuildType()
      {
        return OPENMS_BUILD_TYPE;
      }

      /**
        @brief Return the maximum number of OpenMP threads (including hyperthreads).

        Wraps @c omp_get_max_threads(). Returns @c 1 if the library was compiled without
        OpenMP. The value can be capped by the @c OMP_NUM_THREADS environment variable
        and/or by @ref setOpenMPNumThreads.
      */
      static Size getOpenMPMaxNumThreads()
      {
        #ifdef _OPENMP
        return omp_get_max_threads();
        #else
        return 1;
        #endif
      }
      /**
        @brief Override the OpenMP thread count for subsequent parallel regions.

        Wraps @c omp_set_num_threads. No-op when the library was compiled without OpenMP.
        Takes precedence over the initial @c OMP_NUM_THREADS environment variable.
      */
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
