// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATETIME_H
#define OPENMS_DATASTRUCTURES_DATETIME_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QDateTime>

namespace OpenMS
{
	class String;

	/**	
		@brief DateTime Class.
	
		This class implements date handling.
		Import and export to/from both string and integers is possible.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI DateTime 
		: public QDateTime
	{
		public:
		
			/**
				@brief Default constructor
				
				Fills the object with an undefined date: 00/00/0000
			*/
			DateTime();
			/// Copy constructor
			DateTime(const DateTime& date);
			/// Copy constructor from Qt base class
			DateTime(const QDateTime& date);
			
			/// Assignment operator
			DateTime& operator= (const DateTime& source);
			
			/**
				@brief sets date from a string
				
				Reads both english, german and iso/ansi date formats: 'MM/dd/yyyy', 'dd.MM.yyyy' or 'yyyy-MM-dd'

				@exception Exception::ParseError
			*/
			void setDate(const String& date);
				
			/**
				@brief sets time from a string
				
				Reads time format: 'hh:mm:ss'

				@exception Exception::ParseError
			*/
			void setTime(const String& date);
				
			/**
				@brief sets data from three integers
				
				Give the numbers in the following order: month, day and year.

				@exception Exception::ParseError
			*/
			void setDate(UInt month, UInt day, UInt year);
		
			/**
				@brief sets time from three integers
				
				Give the numbers in the following order: hour, minute and second.

				@exception Exception::ParseError
			*/
			void setTime(UInt hour, UInt minute, UInt second);
		
			/**
				@brief sets data from six integers
				
				Give the numbers in the following order: month, day, year, hour, minute, second.

				@exception Exception::ParseError
			*/
			void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second);
		
			/**
				@brief Fills the arguments with the date and the time
			 	
			 	Give the numbers in the following order: month, day and year, hour minute, second.
			*/
			void get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const;

			/**
				@brief Fills the arguments with the date
			 	
			 	Give the numbers in the following order: month, day and year.
			*/
			void getDate(UInt& month, UInt& day, UInt& year) const;

			/**
				@brief Returns the date as string
			 	
				The format of the string is yyyy-MM-dd
			*/
			String getDate() const;
			
			/**
				@brief Fills the arguments with the time
			 	
				The arguments are all UInts and the order is hour minute second
			*/
			void getTime(UInt& hour, UInt& minute, UInt& second) const;
			
			/**
				@brief Returns the time as string
			 	
				The format of the string is hh:mm:ss
			*/
			String getTime() const;
			
			/// Returns the current date and time
			static DateTime now();
			
			///Sets the undefined date: 00/00/0000 00:00:00
			void clear();
			
			/**
				@brief Returns a string representation of the date and time
			 	
			 	The format of the string will be yyyy-MM-dd hh:mm:ss
			*/
			String get() const;
			
			/**
				@brief Sets date and time
			 	
			 	The following formats are supported:
			 	- MM/dd/yyyy hh:mm:ss
			 	- dd.MM.yyyy hh:mm:ss
			 	- yyyy-MM-dd hh:mm:ss
			 	- yyyy-MM-ddThh:mm:ss (ISO 8601 format)
			 	- yyyy-MM-ddZ (ISO 8601 format)
			 	- yyyy-MM-dd+hh:mm (ISO 8601 format)

				@exception Exception::ParseError
			*/
			void set(const String& date);
						
		protected:
	};
	
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATETIME_H
