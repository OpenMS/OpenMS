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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_SYSTEM_STOPWATCH_H
#define OPENMS_SYSTEM_STOPWATCH_H

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>

#ifdef OPENMS_HAS_SYS_TIME_H
#include <sys/time.h>
#endif

#ifdef OPENMS_HAS_TIME_H
#include <ctime>
#endif

#include <iostream>

namespace OpenMS
{
	/**
		@brief StopWatch Class.
		
		This class is used to determine the current process time.
		
		@ingroup System
	*/
	class OPENMS_DLLAPI StopWatch
	{
		public:
		
		/**	@name	Constructors and Destructors 
		*/
		//@{

		/**	Default constructor.
				Create a new stop watch. The stop watch is stopped.
		*/
		StopWatch();			        

		/**	Copy constructor.
				Create a new stop watch from an existing stop watch.
		*/
		StopWatch(const StopWatch& stop_watch);			        

		/**	Destructor.
				Destructs a stop watch object.
		*/
		virtual ~StopWatch();

		//@}
		/**	Starting, Stopping and Resetting the stop watch 
		*/
		//@{

		/**	Clear and stop the stop watch.
				This sets the stop watch to zero and stops it when running.
				@see	reset
		*/
		void clear();			
	
		/** Start the stop watch.
				The stop watch is started. If the stop watch is already running, <b>false</b>
				is returned.
				@return bool <b>false</b> if the stop watch was already running, <b>true</b> otherwise
		*/
		bool start();			

		/** Stop the stop watch.
				The stop watch is stopped. If the stop watch was not running, <b>false</b>
				is returned.
				@return bool <b>false</b> if the was not running, <b>true</b> otherwise
		*/
		bool stop();

		/** Clear the stop watch without stopping.
				The stop watch is cleared, but not stopped (if running).
				@see clear
		*/
		void reset();

		//@}
			
		/**	@name Readout of the StopWatch 
		*/
		//@{

		/**	Get clock time.
				Return the accumulated clock (real) time in seconds.
		*/
		DoubleReal getClockTime() const;

		/**	Get user time.
				Return the accumulated user time in seconds.
		*/
		DoubleReal getUserTime() const;		

		/**	Get user time.
				Return the accumulated system time in seconds.
		*/
		DoubleReal getSystemTime() const;

		/**	Get CPU time.
				Return the accumulated CPU time in seconds.
				CPU time is the sum of user time and system time.
		*/
		inline DoubleReal getCPUTime() const 
		{ 
			return (getUserTime() + getSystemTime());
		}


		//@}

		/**	@name	Assignment
		*/
		//@{

		/**	Assignment operator.
				Assigns a stop watch from another. The two stop watch will then run 
				synchronously.
				@return StopWatch <tt>*this</tt>
		*/
		StopWatch& operator = (const StopWatch& stop_watch);

		//@}

		/**	@name	Predicates 
		*/
		//@{

		/**	Return true if the stop watch is running.
				@return bool <b>true</b> if the stop watch is running, <b>false</b> otherwise
		*/
		bool isRunning() const 
		{
			return is_running_;
		}

		/**	Equality operator.
				Return <b>true</b> if two stop watchs are equal, i.e. they contain exactly 
				the same time intervals for clock, user and system time and have the
				same running status.
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> on equality, <b>false</b> otherwise
		*/
		bool operator == (const StopWatch& stop_watch) const;

		/**	Inequality operator.
				Return <b>false</b> if two stop watchs differ in any way, i.e. they differ
				in either the clock, user, or system time or have a different 
				running status.
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> on inequality, <b>false</b> otherwise
		*/	
		inline bool operator != (const StopWatch& stop_watch) const 
		{
			return !(*this == stop_watch);
		}


		/**	Lesser than operator.
				Return true, if the stop watch is in all timings lesser than the
				stop watch to be compared with (clock, user and system time).
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> if all times are lesser
		*/
		inline bool operator < (const StopWatch& stop_watch) const 
		{
			return (getCPUTime() < stop_watch.getCPUTime());
		}


		/**	Lesser or equal operator.
				Return true, if the stop watch is in all timings lesser or equal than the
				stop watch to be compared with (clock, user and system time).
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> if all times are lesser or equal
		*/
		inline bool operator <= (const StopWatch& stop_watch) const 
		{
			return !(stop_watch < *this);
		}


		/**	Greater or equal operator.
				Return true, if the stop watch is in all timings greater or equal than the
				stop watch to be compared with (clock, user and system time).
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> if all times are greater or equal
		*/
		inline bool operator >= (const StopWatch& stop_watch) const 
		{
			return !(*this < stop_watch);
		}


		/**	Greater operator.
				Return true, if the stop watch is in all timings greater than the
				stop watch to be compared with (clock, user and system time).
				@param stop_watch the stop watch to compare with
				@return bool <b>true</b> if all times are greater 
		*/
		inline bool operator > (const StopWatch& stop_watch) const 
		{
			return (stop_watch < *this);
		}


		//@}


		private:

		static PointerSizeInt cpu_speed_;

		#ifdef OPENMS_WINDOWSPLATFORM
			static PointerSizeInt clock_speed_;
		#endif

		// state of stop watch, either true(on) or false(off) 
		bool is_running_;

		// clock seconds value when the stop watch was last started 
		PointerSizeInt last_secs_;	

		// clock useconds value when the stop watch was last started 
		PointerSizeInt last_usecs_;		

		// user time when the stop watch was last started 
		clock_t last_user_time_;   

		// system time when the stop watch was last started 
		clock_t last_system_time_; 
		 
		// current accumulated clock seconds 
		PointerSizeInt current_secs_;		

		// current accumulated clock useconds 
		PointerSizeInt current_usecs_;		
		
		// current accumulated user time 
		clock_t current_user_time_;		

		// current accumulated user time 
		clock_t current_system_time_;
	};

}

#endif // OPENMS_SYSTEM_STOPWATCH_H
