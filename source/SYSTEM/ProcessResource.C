// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/SYSTEM/ProcessResource.h>

#ifdef OPENMS_WINDOWSPLATFORM
  #define _WIN32_WINNT 0x0500
  #include <Windows.h>
#else
  #include <sys/resource.h>
	#include <iostream>
#endif


namespace OpenMS
{
	
 	void ProcessResource::LimitCPUTime(const Int& seconds)
	{
		#ifdef OPENMS_WINDOWSPLATFORM
			// Create a named job object
			HANDLE hjob = CreateJobObject(NULL, "self");
			
			JOBOBJECT_EXTENDED_LIMIT_INFORMATION jobli = { 0 }; 
			
			// in 100 of nanoseconds
			
			LONGLONG ss = seconds;
			jobli.BasicLimitInformation.PerJobUserTimeLimit.QuadPart = ss * 10000000;

			
			jobli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_JOB_TIME;
			SetInformationJobObject(hjob, JobObjectExtendedLimitInformation, &jobli, sizeof(jobli));	
			
			// Put our own process in the job 
			AssignProcessToJobObject(hjob, GetCurrentProcess());
			
			// Closing the job does not kill our process or the job
			CloseHandle(hjob);			
		#else // Linux
			struct rlimit r_limit;	
			r_limit.rlim_cur = seconds;
			//r_limit.rlim_max = 15;
			if (setrlimit(RLIMIT_CPU, &r_limit) != 0)
		  {
		    std::cerr << "Error in ProcessResource::LimitCPUTime(): could not set limit!" << std::endl;
		  };
		#endif
	}
		
	
}
