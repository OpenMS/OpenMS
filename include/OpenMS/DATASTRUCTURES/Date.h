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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATE_H
#define OPENMS_DATASTRUCTURES_DATE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
	/**	
		@brief Date Class.
	
		This class implements date handling.
		Import and export to/from both string and integers is possible.
		
		@ingroup Datastructures
	*/
	class Date
	{
		public:
		
			/**
				@brief Default constructor
				
				Fills the object with an undefined date: 00/00/0000
			*/
			Date();
			/// Copy constructor
			Date(const Date& date);
			/// Desctructor
			~Date();
			
			/// Assignment operator
			Date& operator= (const Date& source);
			
			/// Equality operator
			bool operator == (const Date& rhs) const;

			/// Equality operator
			bool operator != (const Date& rhs) const;
			
			/**
				@brief sets data from a string
				
				Reads both english, german and iso/ansi date formats: 'mm/dd/yyyy', 'dd.mm.yyyy' or 'yyyy-mm-dd'
			*/
			void set(const String& date) throw (Exception::ParseError);
				
			/**
				@brief sets data from three integers
				
				Give the numbers in the following order: month, day and year.
			*/
			void set(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw (Exception::ParseError);
		
			/// sets to date to today
			void today();

			/// returns the current date and time as a String in the format (YYYY-MM-DD HH:MM:SS)
			static std::string now();
	
			/**
				@brief Fills the string @p date with the iso/ansi date
				
				Uses the iso/ansi date format: 'yyyy-mm-dd'
			*/
			void get(String& date) const;
			
			/**
				@brief Fills the arguments with the date
			 	
			 	Give the numbers in the following order: month, day and year.
			*/
			void get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const;
			
			///Sets the undefined date: 00/00/0000
			void clear();
			
			/// return true if the given year @p year is a leap year
			bool isLeapYear(UnsignedInt year) const;
			
		protected:
			UnsignedInt day_;
			UnsignedInt month_;
			UnsignedInt year_;
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATE_H
