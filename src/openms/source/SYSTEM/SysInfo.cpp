// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/SysInfo.h>

#ifdef OPENMS_WINDOWSPLATFORM
#include "windows.h"
#include "psapi.h"
#elif __APPLE__
#include <mach/mach.h>
#include <mach/mach_init.h>
#else
#include <cstdio>
#include <unistd.h>
#include <stdlib.h>
#define OMS_USELINUXMEMORYPLATFORM
#endif

namespace OpenMS
{

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
    mem_virtual = pmc.PrivateUsage / 1024; // byte to KB
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

} // namespace OpenMS
