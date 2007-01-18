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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATETIME_H
#define OPENMS_DATASTRUCTURES_DATETIME_H

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
	/**	
		@brief DateTime Class.
	
		This class implements date handling.
		Import and export to/from both string and integers is possible.
		
		@ingroup Datastructures
	*/
	class DateTime : Date
	{
		public:
		
			/**
				@brief Default constructor
				
				Fills the object with an undefined date: 00/00/0000
			*/
			DateTime();
			/// Copy constructor
			DateTime(const DateTime& date);
			/// Desctructor
			~DateTime();
			
			/// Assignment operator
			DateTime& operator= (const DateTime& source);
			
			/// Equality operator
			bool operator == (const DateTime& rhs) const;

			/// Equality operator
			bool operator != (const DateTime& rhs) const;
			
			/**
				@brief sets date from a string
				
				Reads both english, german and iso/ansi date formats: 'mm/dd/yyyy', 'dd.mm.yyyy' or 'yyyy-mm-dd'
			*/
			void setDate(const String& date) throw (Exception::ParseError);
				
			/**
				@brief sets time from a string
				
				Reads time format: 'HH:MM:SS'
			*/
			void setTime(const String& date) throw (Exception::ParseError);
				
			/**
				@brief sets data from three integers
				
				Give the numbers in the following order: month, day and year.
			*/
			void setDate(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw (Exception::ParseError);
		
			/**
				@brief sets time from three integers
				
				Give the numbers in the following order: hour, minute and second.
			*/
			void setTime(UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw (Exception::ParseError);
		
			/**
				@brief sets data from six integers
				
				Give the numbers in the following order: month, day, year, hour, minute, second.
			*/
			void set(UnsignedInt month, UnsignedInt day, UnsignedInt year, UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw (Exception::ParseError);
		
			/**
				@brief Fills the arguments with the date and the time
			 	
			 	Give the numbers in the following order: month, day and year, hour minute, second.
			*/
			void get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year, UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const;

			/**
				@brief Fills the arguments with the date
			 	
			 	Give the numbers in the following order: month, day and year.
			*/
			void getDate(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const;

			/**
				@brief Fills the arguments with the date
			 	
				The format of the string will be yyyy-mm-dd
			*/
			void getDate(String& date) const;
			
			/**
				@brief Fills the arguments with the time
			 	
				The arguments are all UnsignedInts and the order is hour minute second
			*/
			void getTime(UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const;
			
			/**
				@brief Fills the arguments with the time
			 	
				The format of the string will be hh:mm:ss
			*/
			void getTime(String& time) const;
			
			/// sets instance to actual date and time
			void now();
			
			///Sets the undefined date: 00/00/0000 00:00:00
			void clear();
			
			/**
				@brief Fills the argument with the date and time
			 	
			 	The format of the string will be yyyy-mm-dd hh:mm:ss
			*/
			void get(String& date) const;
			
			/**
				@brief Sets date and time
			 	
			 	The format of the string is yyyy-mm-dd hh:mm:ss
			*/
			void set(const String& date) throw (Exception::ParseError);
			
		protected:
			UnsignedInt hour_;
			UnsignedInt minute_; 
			UnsignedInt second_;
	};
	
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATETIME_H
