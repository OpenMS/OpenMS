// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/build_config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <QSysInfo>
#include <QString>
#ifdef _OPENMP
  #include "omp.h"
#endif

namespace OpenMS
{
  namespace Internal
  {

    enum OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX};
    std::string OpenMS_OSNames[] = {"unknown", "MacOS", "Windows", "Linux"};
    enum OpenMS_Architecture {ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT};
    std::string OpenMS_ArchNames[] = {"unknown", "32 bit", "64 bit"};

    class OpenMSOSInfo
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

        // returns something meaningful for basically all important platforms
        info.os_version_ = QSysInfo::productVersion();

        // identify architecture
        if (QSysInfo::WordSize == 32)
        {
          info.arch_ = ARCH_32BIT;
        }
        else
        {
          info.arch_ = ARCH_64BIT;
        }

        return info;
      }
    };


    /// @brief Struct with some static methods to get infos on the build configuration
    struct OpenMSBuildInfo
    {

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
    };

  } // NS Internal
} // NS OpenMS