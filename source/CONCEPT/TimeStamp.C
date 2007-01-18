// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/TimeStamp.h>

#ifdef OPENMS_HAS_SYS_TIME_H
#include <sys/time.h>
#endif 
#ifdef OPENMS_HAS_TIME_H
#include <time.h>
#endif 

#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
#include <windows.h>
#include <sys/timeb.h>
#endif

using namespace std;

namespace OpenMS 
{

	PreciseTime::PreciseTime()
		throw()
		: secs_(0),
			usecs_(0)
	{
		#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
			LARGE_INTEGER t;
			QueryPerformanceFrequency(&t);
			ticks_ = (long)t.QuadPart;
		#endif
	}

	PreciseTime::PreciseTime(const PreciseTime& time)
		throw()
		:	secs_(time.secs_),
			usecs_(time.usecs_)
	{
	}

#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
	long PreciseTime::ticks_;
#endif

	TimeStamp::TimeStamp()
		throw()
		:	time_()
	{
	}


	PreciseTime PreciseTime::now() 
		throw()
	{
		#ifdef OPENMS_COMPILER_MSVC
			#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
				LARGE_INTEGER tvl;
				QueryPerformanceCounter(&tvl);
				long sec = tvl.QuadPart / ticks_;
				long usec = (tvl.QuadPart - sec * ticks_) * 1000000 / ticks_;
				return PreciseTime(sec, usec);
			#else
				struct _timeb tv;
				_ftime(&tv);
				return PreciseTime(tv.time, tv.millitm * 1000);
			#endif
		#else
			// get the current time via the system call
			// gettimeofday()
			struct timeval tv;
			gettimeofday(&tv, 0);
			return PreciseTime(tv.tv_sec, tv.tv_usec);
		#endif
	}

	const PreciseTime PreciseTime::ZERO;

  ostream& operator << (ostream& os, const PreciseTime& time)
		throw()
	{
		int usecs = time.getMicroSeconds();
		time_t secs = (time_t)time.getSeconds();
		static char buf1[128];
		static char buf2[128];
		strftime(buf1, 127, "%Y%m%d%H%M%S", localtime(&secs));
		sprintf(buf2, "%.06d", usecs);

		return os << buf1 << buf2;
	}

  ostream& operator << (ostream& os, const TimeStamp& stamp)
		throw()
	{
		return os << stamp.getTime();
	}
}
