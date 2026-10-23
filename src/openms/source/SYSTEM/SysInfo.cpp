// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/SysInfo.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <thread>
#include <vector>

#ifdef OPENMS_WINDOWSPLATFORM
  #include "windows.h"
  #include "psapi.h"
  #include <process.h>  // for _getpid()
#elif __APPLE__
  #include <mach/mach.h>
  #include <mach/mach_host.h>
  #include <mach/mach_init.h>
  #include <sys/sysctl.h> // sysctlbyname()
  #include <unistd.h>   // for getpid()
#else
  #define OMS_USELINUXMEMORYPLATFORM
  #include <cstdio>
  #include <unistd.h>
  #ifdef OPENMS_HAS_SYS_RESOURCE_H
    #include <sys/resource.h> // for rusage
  #endif
  #ifdef __linux__
    #include <sched.h> // sched_getaffinity()
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

  bool read_available_memory_linux(size_t& mem_available)
  {
    mem_available = 0;

    std::ifstream meminfo("/proc/meminfo");
    std::string line;
    while (std::getline(meminfo, line))
    {
      std::istringstream iss(line);
      std::string key;
      UInt64 value = 0;
      if (!(iss >> key >> value))
      {
        continue;
      }

      if (key == "MemAvailable:")
      {
        if (value > static_cast<UInt64>(std::numeric_limits<size_t>::max()))
        {
          mem_available = std::numeric_limits<size_t>::max();
        }
        else
        {
          mem_available = static_cast<size_t>(value);
        }
        return true;
      }
    }

    const long available_pages = sysconf(_SC_AVPHYS_PAGES);
    const long page_size = sysconf(_SC_PAGESIZE);
    if (available_pages <= 0 || page_size <= 0)
    {
      return false;
    }

    const UInt64 available_kb = static_cast<UInt64>(available_pages) *
                                static_cast<UInt64>(page_size) / 1024ull;
    if (available_kb > static_cast<UInt64>(std::numeric_limits<size_t>::max()))
    {
      mem_available = std::numeric_limits<size_t>::max();
    }
    else
    {
      mem_available = static_cast<size_t>(available_kb);
    }
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

  bool SysInfo::getFreeSystemMemory(size_t& mem_available)
  {
    mem_available = 0;
#ifdef OPENMS_WINDOWSPLATFORM
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (!GlobalMemoryStatusEx(&status))
    {
      return false;
    }
    mem_available = static_cast<size_t>(status.ullAvailPhys / 1024ull);
#elif __APPLE__
    vm_size_t page_size = 0;
    if (KERN_SUCCESS != host_page_size(mach_host_self(), &page_size))
    {
      return false;
    }

    vm_statistics64_data_t vm_stat;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    if (KERN_SUCCESS != host_statistics64(mach_host_self(),
                                          HOST_VM_INFO64,
                                          reinterpret_cast<host_info64_t>(&vm_stat),
                                          &count))
    {
      return false;
    }

    const UInt64 available_pages = static_cast<UInt64>(vm_stat.free_count) +
                                   static_cast<UInt64>(vm_stat.inactive_count);
    mem_available = static_cast<size_t>(available_pages * static_cast<UInt64>(page_size) / 1024ull);
#else // Linux
    if (!read_available_memory_linux(mem_available))
    {
      return false;
    }
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

  Int64 SysInfo::getProcessId()
  {
#ifdef OPENMS_WINDOWSPLATFORM
    return static_cast<Int64>(_getpid());
#else
    return static_cast<Int64>(getpid());
#endif
  }

  UInt SysInfo::getAvailablePhysicalCores([[maybe_unused]] const std::vector<int>& cpus)
  {
    const UInt hardware_threads = std::max(1u, std::thread::hardware_concurrency());
#ifdef OPENMS_WINDOWSPLATFORM
    // processor cores (of the processor group of this thread) with a logical processor in the affinity mask
    DWORD_PTR process_mask = 0, system_mask = 0;
    DWORD length = 0;
    GetLogicalProcessorInformation(nullptr, &length);
    std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> info(length / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
    if (!GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask) || info.empty()
        || !GetLogicalProcessorInformation(info.data(), &length))
    {
      return hardware_threads;
    }
    UInt cores = 0;
    for (const auto& i : info)
    {
      if (i.Relationship == RelationProcessorCore && (i.ProcessorMask & process_mask) != 0) { ++cores; }
    }
    return cores > 0 ? cores : hardware_threads;
#elif __APPLE__
    int cores = 0;
    size_t size = sizeof(cores);
    return sysctlbyname("hw.physicalcpu", &cores, &size, nullptr, 0) == 0 && cores > 0 ? UInt(cores) : hardware_threads;
#elif defined(__linux__)
    // @p cpus, or the CPUs of the affinity mask (which can be larger than CPU_SETSIZE)
    std::vector<int> counted = cpus;
    for (int n = CPU_SETSIZE; counted.empty() && n <= (1 << 20); n *= 2)
    {
      const std::unique_ptr<cpu_set_t, void (*)(cpu_set_t*)> set(CPU_ALLOC(n), [](cpu_set_t* s) { CPU_FREE(s); });
      if (set == nullptr) { break; }
      const size_t set_size = CPU_ALLOC_SIZE(n);
      if (sched_getaffinity(0, set_size, set.get()) == 0)
      {
        for (int cpu = 0; cpu < n; ++cpu)
        {
          if (CPU_ISSET_S(cpu, set_size, set.get())) { counted.push_back(cpu); }
        }
        break;
      }
      if (errno != EINVAL) { break; } // EINVAL: the mask is larger
    }
    // one per core: the hardware threads of a core share their list of siblings
    std::set<std::string> cores;
    for (const int cpu : counted)
    {
      std::ifstream topology("/sys/devices/system/cpu/cpu" + std::to_string(cpu) + "/topology/thread_siblings_list");
      std::string siblings;
      cores.insert(std::getline(topology, siblings) && !siblings.empty() ? siblings : "cpu" + std::to_string(cpu));
    }
    UInt available = cores.empty() ? hardware_threads : UInt(cores.size());
    // cgroup v2: the CPU quota (cpu.max: "<quota> <period>" or "max <period>") of the cgroup of the process and of
    // the cgroups above it, as far as a cgroup2 mount shows them (/proc/self/mountinfo: its root within the hierarchy
    // and its mount point, with octal escapes of blanks and backslashes)
    std::string path; // of the cgroup within the hierarchy
    std::ifstream cgroups("/proc/self/cgroup");
    for (std::string line; path.empty() && std::getline(cgroups, line);)
    {
      if (line.rfind("0::", 0) == 0) { path = line.substr(3); }
    }
    const auto unescape = [](const std::string& field)
    {
      std::string text;
      for (Size i = 0; i < field.size(); ++i)
      {
        const auto octal = [&field](Size at) { return at < field.size() && field[at] >= '0' && field[at] <= '7'; };
        if (field[i] == '\\' && octal(i + 1) && octal(i + 2) && octal(i + 3))
        {
          text += static_cast<char>((field[i + 1] - '0') * 64 + (field[i + 2] - '0') * 8 + (field[i + 3] - '0'));
          i += 3;
        }
        else { text += field[i]; }
      }
      return text;
    };
    std::ifstream mounts("/proc/self/mountinfo");
    for (std::string line; !path.empty() && std::getline(mounts, line);)
    {
      const size_t separator = line.find(" - ");
      if (separator == std::string::npos || line.compare(separator + 3, 8, "cgroup2 ") != 0) { continue; }
      std::istringstream fields(line.substr(0, separator));
      std::string id, parent, device, root, point;
      if (!(fields >> id >> parent >> device >> root >> point)) { continue; }
      root = unescape(root);
      point = unescape(point);
      if (root == "/") { root.clear(); }
      if (path.compare(0, root.size(), root) != 0 || (path.size() > root.size() && path[root.size()] != '/')) { continue; }
      for (std::string below = path.substr(root.size());; below.erase(below.rfind('/')))
      {
        std::ifstream max(point + below + "/cpu.max");
        std::string quota;
        long long period = 0;
        if (max >> quota >> period && quota != "max" && period > 0)
        {
          const long long q = std::atoll(quota.c_str());
          if (q > 0) { available = std::min<UInt>(available, UInt(std::max(1ll, (q + period - 1) / period))); }
        }
        if (below.find('/') == std::string::npos) { break; } // the root of the mount
      }
    }
    return available;
#else
    return hardware_threads;
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

  std::string SysInfo::MemUsage::delta(const std::string& event)
  {
    if (mem_after == 0)
    {
      after(); // collect data if missing; do not test using mem_after_peak, since it might be unsupported on the platform
    }
    std::string s =std::string("Memory usage (") + event + "): ";
    s += diff_str_(mem_before, mem_after) + " (working set delta)";
    if (mem_after_peak > 0)
    { // only if supported
      s+= ", " + diff_str_(mem_before_peak, mem_after_peak) + " (peak working set delta)";
    }
    return s;
  }

  std::string SysInfo::MemUsage::usage()
  {
    if (mem_after == 0)
    {
      after(); // collect data if missing; do not test using mem_after_peak, since it might be unsupported on the platform
    }
    std::string s("Memory usage: ");
    s += diff_str_(0, mem_after) + " (working set)";
    if (mem_after_peak > 0)
    { // only if supported
      s += ", " + diff_str_(0, mem_after_peak) + " (peak working set)";
    }
    return s;
  }

  std::string SysInfo::MemUsage::diff_str_(size_t mem_before, size_t mem_after)
  {
    std::string s;
    if (mem_after < mem_before)
    {
      s +=std::string("-");
    }
    s =StringUtils::toStr(std::abs(((long long)mem_after - (long long)mem_before) / 1024)) + " MB";
    return s;
  }

} // namespace OpenMS
