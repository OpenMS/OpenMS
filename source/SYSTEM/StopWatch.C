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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/StopWatch.h>

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h>
#endif
#ifdef OPENMS_HAS_TIME_H
#include <ctime>
#endif
#ifdef OPENMS_HAS_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef OPENMS_HAS_SYS_TIMES_H
#include <sys/times.h>
#endif
#ifdef OPENMS_HAS_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef OPENMS_WINDOWSPLATFORM
#include <windows.h>
#include <sys/timeb.h>
#endif


namespace OpenMS 
{

	PointerSizeInt StopWatch::cpu_speed_ = 0L;

	#ifdef OPENMS_WINDOWSPLATFORM
		PointerSizeInt StopWatch::clock_speed_ = 0L;
	#endif

	StopWatch::StopWatch()
		:	is_running_(false),
			last_secs_(0),
			last_usecs_(0),
			last_user_time_(0),
			last_system_time_(0),
			current_secs_(0),
			current_usecs_(0),
			current_user_time_(0),
			current_system_time_(0)
	{
		#ifdef OPENMS_HAS_SYSCONF
		if (cpu_speed_ == 0L)
		{
			cpu_speed_ = sysconf(_SC_CLK_TCK);
		}
		#endif
	
		#ifdef OPENMS_WINDOWSPLATFORM
			if (cpu_speed_ == 0L)
			{
				LARGE_INTEGER ticks;
				if (QueryPerformanceFrequency(&ticks))
				{
					cpu_speed_ = (PointerSizeInt) ticks.QuadPart;
				}
				else 
				{
					cpu_speed_ = 1L;
				}
				clock_speed_ = CLOCKS_PER_SEC;
			}
		#endif
	}

	StopWatch::StopWatch(const StopWatch& stop_watch)
		:	is_running_(stop_watch.is_running_),
			last_secs_(stop_watch.last_secs_),
			last_usecs_(stop_watch.last_usecs_),
			last_user_time_(stop_watch.last_user_time_),
			last_system_time_(stop_watch.last_system_time_),
			current_secs_(stop_watch.current_secs_),
			current_usecs_(stop_watch.current_usecs_),
			current_user_time_(stop_watch.current_user_time_),
			current_system_time_(stop_watch.current_system_time_)
	{
	}

	StopWatch::~StopWatch()
	{
	}

	void StopWatch::clear()
	{
		is_running_ = false;
		current_secs_ = 0L;
		current_usecs_ = 0L;
		current_user_time_ = 0L;
		current_system_time_ = (clock_t)0;
	}

	bool StopWatch::start()
	{
		if (is_running_ == true)
		{ 
			/* tried to start a running stop_watch */
			return false;
		}

		#ifdef OPENMS_WINDOWSPLATFORM
			LARGE_INTEGER tms;
			FILETIME kt,ut,ct,et;
			
			QueryPerformanceCounter(&tms);
			HANDLE my_id = GetCurrentProcess();
			GetProcessTimes(my_id, &ct, &et, &kt, &ut);
			ULARGE_INTEGER kernel_time; 
			kernel_time.HighPart = kt.dwHighDateTime;
			kernel_time.LowPart = kt.dwLowDateTime;
			ULARGE_INTEGER user_time; 
			user_time.HighPart = ut.dwHighDateTime;
			user_time.LowPart = ut.dwLowDateTime;

			last_secs_  = tms.QuadPart / cpu_speed_;			
			last_usecs_ = (PointerSizeInt)((DoubleReal)(tms.QuadPart - (last_secs_*cpu_speed_)) / (DoubleReal)(cpu_speed_) * 1000000.0);

			last_user_time_ = (clock_t) (user_time.QuadPart / 10);
			last_system_time_ = (clock_t) (kernel_time.QuadPart / 10);

		#else

			struct tms tms_buffer; 
			struct timeval timeval_buffer;	
			struct timezone timezone_buffer; 

			gettimeofday(&timeval_buffer, &timezone_buffer);
			times(&tms_buffer);
		
			last_secs_ = timeval_buffer.tv_sec;
			last_usecs_ = timeval_buffer.tv_usec;
			last_user_time_ = tms_buffer.tms_utime;
			last_system_time_ = tms_buffer.tms_stime;
		#endif

		is_running_ = true;

		return true;
	}

	bool StopWatch::stop()
	{
		if (is_running_ == false)
		{ /* tried to stop a stopped stop_watch */
			return false;
		}
		#ifdef OPENMS_WINDOWSPLATFORM
			LARGE_INTEGER tms;
	
			QueryPerformanceCounter(&tms);
			FILETIME kt,ut,ct,et;
			
			HANDLE my_id=GetCurrentProcess();
			GetProcessTimes(my_id,&ct,&et,&kt, &ut);
			
			ULARGE_INTEGER kernel_time; 
			kernel_time.HighPart = kt.dwHighDateTime;
			kernel_time.LowPart = kt.dwLowDateTime;
			ULARGE_INTEGER user_time; 
			user_time.HighPart = ut.dwHighDateTime;
			user_time.LowPart = ut.dwLowDateTime;

			PointerSizeInt secs_to_add = tms.QuadPart/cpu_speed_;
			current_secs_ += secs_to_add - last_secs_;
			PointerSizeInt usecs_to_add = (PointerSizeInt)((DoubleReal)(tms.QuadPart - secs_to_add*cpu_speed_) /(DoubleReal)(cpu_speed_) * 1000000.0);
			current_usecs_ += usecs_to_add - last_usecs_;
			
			current_user_time_ += (clock_t) (user_time.QuadPart / 10 - last_user_time_);
			current_system_time_ += (clock_t) (kernel_time.QuadPart / 10 - last_system_time_);
		#else
			struct tms tms_buffer;
			struct timeval timeval_buffer;
			struct timezone timezone_buffer;

			gettimeofday(&timeval_buffer, &timezone_buffer);
			times(&tms_buffer);

			current_secs_ += timeval_buffer.tv_sec - last_secs_;
			current_usecs_ += timeval_buffer.tv_usec - last_usecs_;

			current_user_time_ += tms_buffer.tms_utime - last_user_time_;
			current_system_time_ += tms_buffer.tms_stime - last_system_time_;
		#endif
		
		is_running_ = false;

		return true;
	}

	void StopWatch::reset()
	{
		if (is_running_ == false)
		{
			clear();
		} 
		else 
		{
			stop();
			clear();
			start();
		}
	}


	//getClockTime returns the current amount of real (clock) time			
	//accumulated by this stop_watch.  If the stop_watch is stopped, this is just		
	//the total accumulated time.  If the stop_watch is running, this is the		
	//accumulated time + the time since the stop_watch was last started.				
	DoubleReal StopWatch::getClockTime() const
	{
		PointerSizeInt elapsed_seconds;
		PointerSizeInt elapsed_useconds;

		if (is_running_ == false)
		{ 
			/* stop_watch is currently off, so just return accumulated time */
			elapsed_seconds = current_secs_;
			elapsed_useconds = current_usecs_;
		} 
		else 
		{ 
			/* stop_watch is currently running, so add the elapsed time since */
			/* the stop_watch was last started to the accumulated time        */
			#ifdef OPENMS_WINDOWSPLATFORM
				LARGE_INTEGER tms;
				if (QueryPerformanceCounter(&tms))
				{
					PointerSizeInt secs_to_add = tms.QuadPart / cpu_speed_;
					elapsed_seconds = current_secs_ + secs_to_add - last_secs_;
					PointerSizeInt usecs_to_add = (PointerSizeInt)((DoubleReal)(tms.QuadPart - secs_to_add * cpu_speed_) /(DoubleReal)(cpu_speed_) * 1000000.0);
					elapsed_useconds  = current_usecs_ + usecs_to_add - last_usecs_;
				}
			#else
				struct timeval timeval_buffer;
				struct timezone timezone_buffer;

				gettimeofday(&timeval_buffer, &timezone_buffer);

				elapsed_seconds = current_secs_ + timeval_buffer.tv_sec - last_secs_;
				elapsed_useconds = current_usecs_ + timeval_buffer.tv_usec - last_usecs_;
			#endif
		}

		/* Adjust for the fact that the useconds may be negative. */
		/* If they are, take away 1 second and add 1 million      */
		/* microseconds until they are positive.                  */
		while (elapsed_useconds < 0L)
		{
			elapsed_useconds += 1000000L;
			elapsed_seconds--;
		}

		/* convert into floating point number of seconds */
		return (DoubleReal)((DoubleReal)elapsed_seconds + (DoubleReal)elapsed_useconds / 1000000.0);
	}

	//getUserTime reports the current amount of user cpu time
	//accumulated by this StopWatch.  If the stop_watch is currently off,
	//this is just the accumulated time.  If the StopWatch is running, this
	//is the accumulated time plust the time since the stop_watch was last started.
	DoubleReal StopWatch::getUserTime() const
	{
		DoubleReal temp_value;

		#ifdef OPENMS_WINDOWSPLATFORM
			FILETIME kt,ut,ct,et;
		#else
			struct tms tms_buffer;	
		#endif
		if (is_running_ == false)
		{ 
			/* stop_watch is off, just return accumulated time */
			temp_value = (DoubleReal)current_user_time_;
		}	
		else 
		{
			/* stop_watch is on, add current running time to accumulated time */
			#ifdef OPENMS_WINDOWSPLATFORM
				HANDLE my_id=GetCurrentProcess();
				GetProcessTimes(my_id,&ct,&et,&kt,&ut);
				
				ULARGE_INTEGER kernel_time; 
				kernel_time.HighPart = kt.dwHighDateTime;
				kernel_time.LowPart = kt.dwLowDateTime;
				ULARGE_INTEGER user_time; 
				user_time.HighPart = ut.dwHighDateTime;
				user_time.LowPart = ut.dwLowDateTime;
				
				temp_value = (DoubleReal)(current_user_time_ + user_time.QuadPart / 10.0 - last_user_time_);
			#else
				times(&tms_buffer);
				temp_value = (DoubleReal)(current_user_time_ + tms_buffer.tms_utime - last_user_time_);
			#endif
		}

		#ifdef OPENMS_WINDOWSPLATFORM
			return (DoubleReal)(temp_value / 1000000.0);
		#else		
			/* convert from clock ticks to seconds using the */
			/* cpu-speed value obtained in the constructor   */
			return (DoubleReal)(temp_value / (DoubleReal)cpu_speed_);
		#endif	
	}

	// system_time reports the current amount of system cpu time 
	// accumulated by this StopWatch.  If the stop_watch is currently off,
	// this is just the accumulated time.  If the StopWatch is running, this
	// is the accumulated time plus  the time since the stop_watch was last started
	DoubleReal StopWatch::getSystemTime() const
	{
		DoubleReal temp_value = 0.0;

		#ifdef OPENMS_WINDOWSPLATFORM
			//struct tms tms_buffer;
			FILETIME kt,ut,ct,et;
		#endif												
		
		if (is_running_ == false)
		{ 
			/* stop_watch is off, just return accumulated time */
			temp_value = (DoubleReal)current_system_time_;
		} 
		else 
		{ 
			/* stop_watch is on, return accumulated plus current */
			#ifdef OPENMS_WINDOWSPLATFORM
				//times(&tms_buffer);
				HANDLE my_id=GetCurrentProcess();
				GetProcessTimes(my_id,&ct,&et,&kt,&ut);
			
				ULARGE_INTEGER kernel_time; 
				kernel_time.HighPart = kt.dwHighDateTime;
				kernel_time.LowPart = kt.dwLowDateTime;
				ULARGE_INTEGER user_time; 
				user_time.HighPart = ut.dwHighDateTime;
				user_time.LowPart = ut.dwLowDateTime;
				temp_value = (DoubleReal)((DoubleReal)current_system_time_ + kernel_time.QuadPart / 10.0 - (DoubleReal)last_system_time_);
			#endif
		}

		/* convert from clock ticks to seconds using the */
		/* cpu-speed value obtained by the constructor   */
		#ifndef OPENMS_WINDOWSPLATFORM
			return (DoubleReal)(temp_value / 1000000.0);
		#else 
			return 0.0;
		#endif
	}

	StopWatch& StopWatch::operator = (const StopWatch& stop_watch)
	{
		if (this == &stop_watch)
		{
			return *this;
		}

		is_running_ = stop_watch.is_running_;
		last_secs_ = stop_watch.last_secs_;	
		last_usecs_ = stop_watch.last_usecs_;		
		last_user_time_ = stop_watch.last_user_time_;   
		last_system_time_ = stop_watch.last_system_time_; 
		current_secs_ = stop_watch.current_secs_;		
		current_usecs_ = stop_watch.current_usecs_;		
		current_user_time_ = stop_watch.current_user_time_;		
		current_system_time_ = stop_watch.current_system_time_;

		return *this;
	}

	bool StopWatch::operator == (const StopWatch& stop_watch) const
	{
		return (last_secs_ == stop_watch.last_secs_
									&& last_usecs_ == stop_watch.last_usecs_
									&& last_user_time_ == stop_watch.last_user_time_
									&& last_system_time_ == stop_watch.last_system_time_
									&& current_secs_ == stop_watch.current_secs_
									&& current_usecs_ == stop_watch.current_usecs_
									&& current_user_time_ == stop_watch.current_user_time_
									&& current_system_time_ == stop_watch.current_system_time_);
	}

} // namespace OpenMS
