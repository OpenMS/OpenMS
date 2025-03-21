// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/SysInfo.h>

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef OPENMS_WINDOWSPLATFORM
  #include "windows.h"
  #include "psapi.h"
#elif __APPLE__
  #include <mach/mach.h>
  #include <mach/mach_init.h>
#else
  #define OMS_USELINUXMEMORYPLATFORM
  #include <cstdio>
  #include <unistd.h>
  #ifdef OPENMS_HAS_SYS_RESOURCE_H
    #include <sys/resource.h> // for rusage
  #endif
#endif

namespace OpenMS
{

  std::string bytesToHumanReadable(UInt64 bytes)
  {
    std::array units {"byte", "KiB", "MiB", "GiB", "TiB", "PiB"};

    const int divisor = 1024;
    double db = double(bytes);
    for (const auto u : units)
    {
      if (db < divisor)
      {
        std::stringstream ss;
        ss << std::setprecision(4) /* 4 digits overall, i.e. 1000 or 1.456 */ << db << ' ' << u;
        return ss.str();
      }
      db /= divisor;
    }
    // wow ... you made it here...
    return std::string("Congrats. That's a lot of bytes: ") + std::to_string(bytes);
  }

#ifdef OMS_USELINUXMEMORYPLATFORM
  // see http://stackoverflow.com/questions/1558402/memory-usage-of-current-process-in-c
  typedef struct {
      long size,resident,share,text,lib,data,dt;
  } statm_t;

  bool read_off_memory_status_linux(statm_t& result)
  {
    const char* statm_path = "/proc/self/statm";

    FILE *f = fopen(statm_path,"r");
    if (!f)
    {
      return false;
    }

    // get 'data size (heap + stack)'  (residence size (vmRSS) is usually too
    // small and not changing, total memory (vmSize) is changing but usually
    // too large)

    // From the proc(5) man-page:
    //
    //    /proc/[pid]/statm
    //           Provides information about memory usage, measured in pages.
    //           The columns are:
    //
    //               size       total program size
    //                          (same as VmSize in /proc/[pid]/status)
    //               resident   resident set size
    //                          (same as VmRSS in /proc/[pid]/status)
    //               share      shared pages (from shared mappings)
    //               text       text (code)
    //               lib        library (unused in Linux 2.6)
    //               data       data + stack
    //               dt         dirty pages (unused in Linux 2.6)


    if (7 != fscanf(f,"%ld %ld %ld %ld %ld %ld %ld",
              &result.size,&result.resident,&result.share,&result.text,&result.lib,&result.data,&result.dt))
    {
      fclose(f);
      return false;
    }
    fclose(f);
    return true;
  }
#endif

  bool SysInfo::getProcessMemoryConsumption(size_t& mem_virtual)
  {
    mem_virtual = 0;
#ifdef OPENMS_WINDOWSPLATFORM
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc)))
    {
      return false;
    }
    mem_virtual = pmc.WorkingSetSize / 1024; // byte to KB
#elif __APPLE__
    struct task_basic_info_64 t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_64_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
                                  TASK_BASIC_INFO_64, (task_info_t)&t_info,
                                  &t_info_count))
    {
      return false;
    }
    mem_virtual = t_info.resident_size / 1024; // byte to KB
#else // Linux
    statm_t mem;
    if (!read_off_memory_status_linux(mem))
    {
      return false;
    }
    mem_virtual = (size_t)mem.resident * (size_t)sysconf(_SC_PAGESIZE) / 1024; // byte to KB
#endif
    return true;
  }

  bool SysInfo::getProcessPeakMemoryConsumption(size_t& mem_virtual)
  {
    mem_virtual = 0;
#ifdef OPENMS_WINDOWSPLATFORM
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc)))
    {
      return false;
    }
    mem_virtual = pmc.PeakWorkingSetSize / 1024; // byte to KB
    return true;
#elif __APPLE__
    rusage ru;
    if (getrusage(0, &ru) == 0) // success;
    {
      mem_virtual = ru.ru_maxrss / 1024; // reported in bytes (whereas Linux is KB!). Convert to KB
      return true;
    }
    return false;
#else // Linux
    #ifdef OPENMS_HAS_SYS_RESOURCE_H
    rusage ru;
    if (getrusage(0, &ru) == 0) // success;
    {
      mem_virtual = ru.ru_maxrss; // in KB already
      return true;
    }
    #endif
    return false;
#endif
  }

  SysInfo::MemUsage::MemUsage()
    : mem_before(0), mem_before_peak(0), mem_after(0), mem_after_peak(0)
  {
    before();
  }

  void SysInfo::MemUsage::reset()
  {
    mem_before = mem_before_peak = mem_after = mem_after_peak = 0;
  }

  void SysInfo::MemUsage::before()
  {
    SysInfo::getProcessMemoryConsumption(mem_before);
    SysInfo::getProcessPeakMemoryConsumption(mem_before_peak);
  }

  void SysInfo::MemUsage::after()
  {
    SysInfo::getProcessMemoryConsumption(mem_after);
    SysInfo::getProcessPeakMemoryConsumption(mem_after_peak);
  }

  String SysInfo::MemUsage::delta(const String& event)
  {
    if (mem_after == 0)
    {
      after(); // collect data if missing; do not test using mem_after_peak, since it might be unsupported on the platform
    }
    String s = String("Memory usage (") + event + "): ";
    s += diff_str_(mem_before, mem_after) + " (working set delta)";
    if (mem_after_peak > 0)
    { // only if supported
      s+= ", " + diff_str_(mem_before_peak, mem_after_peak) + " (peak working set delta)";
    }
    return s;
  }

  String SysInfo::MemUsage::usage()
  {
    if (mem_after == 0)
    {
      after(); // collect data if missing; do not test using mem_after_peak, since it might be unsupported on the platform
    }
    String s("Memory usage: ");
    s += diff_str_(0, mem_after) + " (working set)";
    if (mem_after_peak > 0)
    { // only if supported
      s += ", " + diff_str_(0, mem_after_peak) + " (peak working set)";
    }
    return s;
  }

  String SysInfo::MemUsage::diff_str_(size_t mem_before, size_t mem_after)
  {
    String s;
    if (mem_after < mem_before)
    {
      s += String("-");
    }
    s = String(std::abs(((long long)mem_after - (long long)mem_before) / 1024)) + " MB";
    return s;
  }

} // namespace OpenMS
