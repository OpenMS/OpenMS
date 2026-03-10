// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/SYSTEM/BuildInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#ifdef __linux__
#include <sys/utsname.h>
#endif

#include <simde/simde-arch.h>

using namespace std;

namespace OpenMS
{

  String Internal::OpenMSOSInfo::getActiveSIMDExtensions()
  {
    StringList ret;
#ifdef SIMDE_ARCH_ARM_NEON
    ret.push_back("neon");
#endif
#ifdef SIMDE_ARCH_X86_SSE
    ret.push_back("SSE");
#endif
#ifdef SIMDE_ARCH_X86_SSE2
    ret.push_back("SSE2");
#endif
#ifdef SIMDE_ARCH_X86_SSE3
    ret.push_back("SSE3");
#endif
#ifdef SIMDE_ARCH_X86_SSE4_1
    ret.push_back("SSE4.1");
#endif
#ifdef SIMDE_ARCH_X86_SSE4_2
    ret.push_back("SSE4.2");
#endif
#ifdef SIMDE_ARCH_X86_AVX
    ret.push_back("AVX");
#endif
#ifdef SIMDE_ARCH_X86_AVX2
    ret.push_back("AVX2");
#endif
#ifdef SIMDE_ARCH_X86_FMA
    ret.push_back("FMA");
#endif

    return ListUtils::concatenate(ret, ", ");
  }

  String Internal::OpenMSOSInfo::getOSVersion()
  {
#ifdef OPENMS_WINDOWSPLATFORM
    // Windows version detection
    OSVERSIONINFO osvi;
    ZeroMemory(&osvi, sizeof(OSVERSIONINFO));
    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
    
    if (GetVersionEx(&osvi))
    {
      return String(std::to_string(osvi.dwMajorVersion) + "." + std::to_string(osvi.dwMinorVersion));
    }
    return "unknown";
#elif defined(__APPLE__)
    // macOS version detection
    size_t size = 0;
    if (sysctlbyname("kern.osproductversion", nullptr, &size, nullptr, 0) == 0)
    {
      std::string version(size, '\0');
      if (sysctlbyname("kern.osproductversion", &version[0], &size, nullptr, 0) == 0)
      {
        return String(version.c_str());
      }
    }
    return "unknown";
#elif defined(__linux__)
    // Linux version detection
    struct utsname uts;
    if (uname(&uts) == 0)
    {
      return String(uts.release);
    }
    return "unknown";
#else
    return "unknown";
#endif
  }

} // namespace OpenMS