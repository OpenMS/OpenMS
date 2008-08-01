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
// $Maintainer: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_TIMESTAMP_H
#define OPENMS_CONCEPT_TIMESTAMP_H

#include <OpenMS/config.h>

#include <iostream>
#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
#include <windows.h>
#include <sys/timeb.h>
#endif

namespace OpenMS 
{

	/**	
		@brief Time class.
		
		Used to store a point of time.
			
		@ingroup Concept
	*/
	class PreciseTime
	{

		public:

		/**	@name Constructors and Destructors.
		*/
		//@{

		/**	Default constructor.
				Initialize with zero.
		*/
		PreciseTime()
			;
			
		/**	Copy constructor
		*/
		PreciseTime(const PreciseTime& time)
			;

		/**	Detailed constructor
		*/
		PreciseTime(long secs, long usecs)
			;
			
		/**	Destructor
		*/
		virtual ~PreciseTime()
			;

		//@}
		/**	@name Constants.
		*/
		//@{

		/**	Zero object.
		*/
		static const PreciseTime ZERO;

		//@}
		/**	@name Assignment
		*/
		//@{

		/** Assignment method
		*/
		void set(long secs, long usecs) 
			;

		/** Assignment method
		*/
		void set(const PreciseTime& time) 
			;

		/**	Assignment operator
		*/
		const PreciseTime& operator = (const PreciseTime& time) 
			;

		/**	Clear method
		*/
		virtual void clear() 
			;

		//@}
		/**	@name Predicates
		*/
		//@{

		/**	Greater than operator.
		*/
		bool operator < (const PreciseTime& time) const 
			;

		/**	Lesser than operator.
		*/
		bool operator > (const PreciseTime& time) const 
			;

		/**	Equality operator.
		*/
		bool operator == (const PreciseTime& time) const 
			;

		//@}
		/**	@name Accessors 
		*/
		//@{

		/**	Return the seconds since Jan. 1, 1970.
		*/
		long getSeconds() const 
			;

		/**	Return the microseconds.
		*/
		long getMicroSeconds() const 
			;

		/**	Return the current time.
				@return PreciseTime the current time in seconds since Jan. 1, 1970
		*/
		static PreciseTime now() 
			;

		//@}

		protected:

		long secs_;
		long usecs_;

		#ifdef OPENMS_HAS_WINDOWS_PERFORMANCE_COUNTER
			static long ticks_;
		#endif
	};

	/**	
		@brief Time stamp class.
		
		This class implements a so-called time stamp. It is used to 
		store modification or creation times of objects.
			
		@ingroup Concept
	*/
	class TimeStamp
	{
		public:

		/**	@name	Constructors and Destructors
		*/
		//@{

		/** Default constructor
		*/
		TimeStamp()
			;

		/** Destructor
		*/
		virtual ~TimeStamp()
			;

		//@}
		/**	@name	Predicates
		*/
		//@{

		/**	Check the time stamp.
		*/
		bool isNewerThan(const PreciseTime& time) const 
			;
		
		/**	Check the time stamp.
		*/
		bool isOlderThan(const PreciseTime& time) const 
			;

		/**	Check the time stamp.
		*/
		bool isNewerThan(const TimeStamp& stamp) const 
			;
		
		/**	Check the time stamp.
		*/
		bool isOlderThan(const TimeStamp& stamp) const 
			;

		/** Equality operator
		*/
		bool operator == (const TimeStamp& stamp) const 
			;

		/**	Lesser than operator.
		*/
		bool operator < (const TimeStamp& stamp) const 
			;

		/**	Greater than operator.
		*/
		bool operator > (const TimeStamp& stamp) const 
			;

		
		//@}
		/**	@name Accessors
		*/	
		//@{

		/**	Update the time stamp.
				Store the value of <tt>time</tt> in the internal time stamp.
				If <tt>time</tt> is 0, use the current time (as given by
				 \link PreciseTime::now PreciseTime::now \endlink ).
				@param time the new time stamp (default =  \link PreciseTime::now PreciseTime::now \endlink )														
		*/
		virtual void stamp(const PreciseTime& time = PreciseTime::ZERO) 
			;

		/**	Return the time of last modification
				@return the time stamp
		*/
		const PreciseTime& getTime() const 
			;

		//@}
		/**	@name	Assignment
		*/
		//@{

		/**	Assignment operator
		*/
		const PreciseTime& operator = (const PreciseTime& time) 
			;

		/**	Clear method
		*/
		virtual void clear() 
			;

		//@}

		protected:

		/**	The time stamp.
		*/
		PreciseTime time_;
	};

	/**	Global stream operators for PreciseTime and TimeStamp
	*/
	//@{

	/**	Print the contents of a PreciseTime object to a stream.
	*/
	std::ostream& operator << (std::ostream& os, const PreciseTime& time)
		;

	/**	Print the contents of a TimeStamp object to a stream.
	*/
	std::ostream& operator << (std::ostream& os, const TimeStamp& stamp)
		;

	//@}

	inline
	PreciseTime::PreciseTime(long secs, long usecs)
		
		:	secs_(secs),
			usecs_(usecs)
	{
	}

	inline
	PreciseTime::~PreciseTime()
		
	{
	}

	inline
	void PreciseTime::set(const PreciseTime& time)
		
	{
		secs_ = time.secs_;
		usecs_ = time.usecs_;
	}

	inline
	void PreciseTime::set(long secs, long usecs)
		
	{
		secs_ = secs;
		usecs_ = usecs;
	}

	inline
	const PreciseTime& PreciseTime::operator = (const PreciseTime& time)
		
	{
		set(time);
		return *this;
	}

	inline
	void PreciseTime::clear() 
		
	{
		secs_ = 0;
		usecs_ = 0;
	}

	inline
	bool PreciseTime::operator < (const PreciseTime& time) const
		
	{
		return ((secs_ < time.secs_) || ((secs_ == time.secs_) && (usecs_ < time.usecs_)));
	}

	inline
	bool PreciseTime::operator > (const PreciseTime& time) const 
		
	{
		return ((secs_ > time.secs_) || ((secs_ == time.secs_) && (usecs_ > time.usecs_)));
	}

	inline
	bool PreciseTime::operator == (const PreciseTime& time) const 
		
	{
		return ((secs_ == time.secs_) && (usecs_ == time.usecs_));
	}

	inline
	long PreciseTime::getSeconds() const 
		
	{
		return secs_;
	}

	inline
	long PreciseTime::getMicroSeconds() const 
		
	{
		return usecs_;
	}

	inline
	bool TimeStamp::isOlderThan(const PreciseTime& time) const 
		
	{
		return (time_ < time);
	}

	inline
	bool TimeStamp::isNewerThan(const PreciseTime& time) const 
		
	{
		return (time_ > time);
	}

	inline
	bool TimeStamp::isOlderThan(const TimeStamp& stamp) const 
		
	{
		return (time_ < stamp.time_);
	}

	inline
	bool TimeStamp::isNewerThan(const TimeStamp& stamp) const 
		
	{
		return (time_ > stamp.time_);
	}

	inline
	void TimeStamp::clear()
		
	{
		time_.clear();
	}

	inline
	TimeStamp::~TimeStamp()
		
	{
		clear();
	}

	inline
	void TimeStamp::stamp(const PreciseTime& time) 
		
	{
		// in the default case, stamp with the current 
		// time
		if (time == PreciseTime::ZERO)
		{
			time_ = PreciseTime::now();
		}
		else 
		{
			time_ = time;
		}
	}

	inline
	const PreciseTime& TimeStamp::getTime() const 
		
	{
		return time_;
	}

	inline
	bool TimeStamp::operator == (const TimeStamp& stamp) const 
		
	{
		return time_ == stamp.time_;
	}

	inline
	bool TimeStamp::operator < (const TimeStamp& stamp) const
		
	{
		return time_ < stamp.time_;
	}

	inline
	bool TimeStamp::operator > (const TimeStamp& stamp) const 
		
	{
		return time_ > stamp.time_;
	}

} // namespace OpenMS

#endif // OPENMS_CONCEPT_TIMESTAMP_H
