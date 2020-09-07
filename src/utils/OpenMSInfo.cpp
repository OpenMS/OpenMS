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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#include <QSysInfo>
#include <QDir>

#if defined(Q_WS_MAC)
  #include <CoreServices/CoreServices.h>
#endif

#include <iostream>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    enum OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX};
    std::string OpenMS_OSNames[] = {"unknown", "MacOS", "Windows", "Linux"};
    enum OpenMS_Architecture {ARCH_UNKNOWN, ARCH_32BIT, ARCH_64BIT};
    std::string OpenMS_ArchNames[] = {"unknown", "32 bit", "64 bit"};

#if WIN32
    OpenMS_Architecture getArchOnWin();
    String getWinOSVersion();
#endif

    struct OpenMSOSInfo
    {
      OpenMSOSInfo() :
        os(OS_UNKNOWN),
        os_version("unknown"),
        arch(ARCH_UNKNOWN)
      {}

      OpenMS_OS os;
      String os_version;
      OpenMS_Architecture arch;

      String getOSAsString() const
      {
        return OpenMS_OSNames[os];
      }

      String getArchAsString() const
      {
        return OpenMS_ArchNames[arch];
      }

    };

    static OpenMSOSInfo getOSInfo()
    {
      OpenMSOSInfo info;
#if defined(WIN32)  // Windows
      info.os = OS_WINDOWS;
      info.arch = getArchOnWin();
      info.os_version = getWinOSVersion();
#elif (defined(__MACH__) && defined(__APPLE__)) // MacOS
      info.os = OS_MACOS;

// check if we can use QSysInfo
#if defined(Q_WS_MAC)

      // query gestalt for detailed osx version information
      // NOTE: Gestalt will be deprecated at some point in the future where we
      //       have to look for a better solution
      SInt32 majorVersion,minorVersion,bugFixVersion;
      Gestalt(gestaltSystemVersionMajor, &majorVersion);
      Gestalt(gestaltSystemVersionMinor, &minorVersion);
      Gestalt(gestaltSystemVersionBugFix, &bugFixVersion);
      info.os_version = String(majorVersion) + "." + String(minorVersion) + String(".") + String(bugFixVersion);

      // identify architecture
      if (QSysInfo::WordSize == 32)
      {
        info.arch = ARCH_32BIT;
      }
      else
      {
        info.arch = ARCH_64BIT;
      }
#endif

#else //Linux
      info.os = OS_LINUX;
      //TODO
#endif

      return info;
    }

//********************
//  Windows specific API calls
//********************
#ifdef WIN32
#include <windows.h>
#include <stdio.h>

    typedef BOOL (WINAPI * LPFN_ISWOW64PROCESS)(HANDLE, PBOOL);

    LPFN_ISWOW64PROCESS fnIsWow64Process;

    OpenMS_Architecture getArchOnWin()
    {
#ifdef OPENMS_64BIT_ARCHITECTURE
      return ARCH_64BIT;

#else
      BOOL bIsWow64 = FALSE;

      //IsWow64Process is not available on all supported versions of Windows.
      //Use GetModuleHandle to get a handle to the DLL that contains the function
      //and GetProcAddress to get a pointer to the function if available.

      fnIsWow64Process = (LPFN_ISWOW64PROCESS) GetProcAddress(
        GetModuleHandle(TEXT("kernel32")), "IsWow64Process");

      if (NULL != fnIsWow64Process)
      {
        if (!fnIsWow64Process(GetCurrentProcess(), &bIsWow64))
        {
          return ARCH_UNKNOWN;
        }
      }
      if (bIsWow64)
      {
        return ARCH_64BIT;
      }
      else
      {
        return ARCH_32BIT;
      }
#endif
    }

    String getWinOSVersion()
    {
      OSVERSIONINFO osvi;
      ZeroMemory(&osvi, sizeof(OSVERSIONINFO));
      osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
      GetVersionEx(&osvi);
      return String(osvi.dwMajorVersion) + "." + String(osvi.dwMinorVersion);
    }

#endif // WIN32 API functions

    /*
      @brief Get Architecture of this binary (simply by looking at size of a pointer, i.e. size_t).
    */
    OpenMS_Architecture getBinaryArchitecture()
    {
      size_t bytes = sizeof(size_t);
      switch (bytes)
      {
        case 4:
          return ARCH_32BIT;
        case 8:
          return ARCH_64BIT;
        default:
          return ARCH_UNKNOWN;
      }
    }

  } // NS Internal
} // NS OpenMS

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int main(int /*argc*/, const char ** /*argv*/)
{
  cout << "OpenMS Version:" << "\n";
  cout << "==================" << "\n";
  cout << "Version      : " << VersionInfo::getVersion() << "\n";
  cout << "Build time   : " << VersionInfo::getTime() << "\n";
  cout << "Git sha1     : " << VersionInfo::getRevision() << "\n";
  cout << "Git branch   : " << VersionInfo::getBranch() << "\n";
  cout << "\n";
  cout << "Installation information:" << "\n";
  cout << "==================" << "\n";
  cout << "Data path    : " << File::getOpenMSDataPath() << "\n";
  cout << "Temp path    : " << File::getTempDirectory() << "\n";
  cout << "Userdata path: " << File::getUserDirectory() << "\n";

  cout << "\n";
  cout << "Build information:" << "\n";
  cout << "==================" << "\n";
  cout << "Source path  : " << OPENMS_SOURCE_PATH << "\n";
  cout << "Binary path  : " << OPENMS_BINARY_PATH << "\n";
  cout << "Binary arch  : " << Internal::OpenMS_ArchNames[Internal::getBinaryArchitecture()] << "\n";
  cout << "\n";

  OpenMS::Internal::OpenMSOSInfo info = OpenMS::Internal::getOSInfo();
  // experimental: OS information
  cout << "OS Information:" << "\n";
  cout << "==================" << "\n";
  cout << "Name: " << info.getOSAsString() << "\n";
  cout << "Version: " << info.os_version << "\n";
  cout << "Architecture: " << info.getArchAsString() << "\n";
  cout << "\n";


  return 0;
}

/// @endcond
