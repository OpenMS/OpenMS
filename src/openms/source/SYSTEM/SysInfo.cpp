// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#elif _APPLE__
#include <mach/mach.h>
#else
#include <cstdio>
#include <unistd.h>
#endif

namespace OpenMS
{
	bool SysInfo::getProcessMemoryConsumption(size_t& mem_virtual) 
	{
		mem_virtual = 0;
#ifdef OPENMS_WINDOWSPLATFORM
		PROCESS_MEMORY_COUNTERS_EX pmc;
		if (!GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc)))
		{
			return false;
		}
		mem_virtual = pmc.PrivateUsage / 1024;  // byte to KB
#elif _APPLE__
		struct task_basic_info t_info;
		mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

		if (KERN_SUCCESS != task_info(mach_task_self(),
			TASK_BASIC_INFO, (task_info_t)&t_info,
			&t_info_count))
		{
			return false;
		}
		mem_virtual = t_info.resident_size / 1024;  // byte to KB
		//mem_resident = t_info.virtual_size / 1024;  // byte to KB
#else // Linux
		
        long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
	{
		return false;
	}
	char buf[1024];
	fread(buf, 1, 1024, fp);
	//printf("%s", buf);
        // get 'data size (heap + stack)'  (residence size (vmRSS) is usually too small and not changing, total memory (vmSize) is changing but usually too large)
	if ( sscanf( buf, "%*s%*s%*s%*s%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
                return false;

	}
	fclose( fp );
	mem_virtual = (size_t)rss * (size_t)sysconf( _SC_PAGESIZE) / 1024;

#endif
		return true;
	}

} // namespace OpenMS
