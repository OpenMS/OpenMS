// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QSysInfo>
#include <QDir>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
  namespace Internal
  {

    enum OpenMS_OS {OS_UNKNOWN, OS_MACOS, OS_WINDOWS, OS_LINUX};
    std::string OpenMS_OSNames[] = {"unkown","MacOS","Windows","Linux"};
    enum OpenMS_Architecture {ARCH_UNKOWN, ARCH_32BIT, ARCH_64BIT};
    std::string OpenMS_ArchNames[] = {"unkown","32bit","64bit"};

#if WIN32
    OpenMS_Architecture getArchOnWin();
    String getWinOSVersion();
#endif

    struct OpenMSOSInfo
    {
      OpenMSOSInfo()
        : os(OS_UNKNOWN),
          os_version("unkown"),
          arch(ARCH_UNKOWN)
      {}

      OpenMS_OS os;
      String os_version;
      OpenMS_Architecture arch;

      String getOSAsString()
      {
        return OpenMS_OSNames[os];
      }
      String getArchAsString()
      {
        return OpenMS_ArchNames[arch];
      }
    };

    static OpenMSOSInfo getOSInfo()
    {
      OpenMSOSInfo info;
#if defined (WIN32) // Windows
      info.os = OS_WINDOWS;
      info.arch = getArchOnWin();
      info.os_version = getWinOSVersion();
#elif (defined(__MACH__) && defined(__APPLE__)) // MacOS
      info.os = OS_MACOS;

// check if we can use QSysInfo 
#if defined(Q_WS_MAC)
			// identify 
			switch (QSysInfo::MacintoshVersion)
			{
				case QSysInfo::MV_10_2: 
				{
					info.os_version = "10.2";
					break;
				}
				case QSysInfo::MV_10_3: 
				{
					info.os_version = "10.3";
					break;
				}
				case QSysInfo::MV_10_4: 
				{
					info.os_version = "10.4";
					break;
				}
				case QSysInfo::MV_10_5:
				{
					info.os_version = "10.5";
					break;
				}
				case QSysInfo::MV_10_6:
				{
					info.os_version = "10.6";
					break;
				}
				default:
				{
					info.os_version = "unsupported";
				}				
			}

			// identify architecture
			if(QSysInfo::WordSize == 32) 
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

    typedef BOOL (WINAPI *LPFN_ISWOW64PROCESS) (HANDLE, PBOOL);

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
          GetModuleHandle(TEXT("kernel32")),"IsWow64Process");

      if(NULL != fnIsWow64Process)
      {
          if (!fnIsWow64Process(GetCurrentProcess(),&bIsWow64))
          {
             return ARCH_UNKOWN;
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

  } // NS Internal
} // NS OpenMS


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int main( int /*argc*/, const char** /*argv*/ )
{
	cout << "OpenMS Version:" << endl;
	cout << "==================" << endl;
	cout << "Version      : " << VersionInfo::getVersion() << endl;
	cout << "Build time   : " << VersionInfo::getTime() << endl;
	cout << "SVN revision : " << VersionInfo::getRevision() << endl;
	cout << endl;
	cout << "Installation information:" << endl;
	cout << "==================" << endl;
	cout << "Data path    : " << OPENMS_DATA_PATH << endl;
  cout << "Temp path    : " << String(QDir::tempPath ()) << endl;
  cout << "Userdata path: " << String(QDir::homePath ()) << endl;

	cout << endl;
	cout << "Build information:" << endl;
	cout << "==================" << endl;
	cout << "Source path  : " << OPENMS_SOURCE_PATH << endl;
	cout << "Binary path  : " << OPENMS_BINARY_PATH << endl;
	cout << endl;

  OpenMS::Internal::OpenMSOSInfo info = OpenMS::Internal::getOSInfo();
  // experimental: OS information
  cout << "OS Information:" << endl;
	cout << "==================" << endl;
  cout << "Name: " << info.getOSAsString() << endl;
  cout << "Version: " << info.os_version << endl;
  cout << "Architecture: " << info.getArchAsString() << endl;
	cout << endl;


	return 0;
}

/// @endcond



