// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATE_H
#define OPENMS_DATASTRUCTURES_DATE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QDate>

namespace OpenMS
{
	/**	
		@brief Date Class.
	
		This class implements date handling.
		Import and export to/from both string and integers is possible.
		
		@ingroup Datastructures
	*/
	class OPENMS_DLLAPI Date : public QDate 
	{
		public:
		
			/**
				@brief Default constructor
				
				Fills the object with an undefined date: 00/00/0000
			*/
			Date();
			/// Copy constructor
			Date(const Date& date);
			/// Copy constructor from Qt base class
			Date(const QDate& date);

			/// Assignment operator
			Date& operator= (const Date& source);
			
			/**
				@brief sets data from a string
				
				The following date formats are supoported:
				- mm/dd/yyyy
				- dd.mm.yyyy
				- yyyy-mm-dd
				
				@exception Exception::ParseError is thrown if the date is given in the wrong format
			*/
			void set(const String& date);
				
			/**
				@brief sets data from three integers
				
				@exception Exception::ParseError is thrown if an invalid date is given
			*/
			void set(UInt month, UInt day, UInt year);
		
			/// Returns the current date
			static Date today();

			/**
				@brief Returns a string representation of the date
				
				Uses the iso/ansi date format: 'yyyy-mm-dd'
			*/
			String get() const;
			
			/**
				@brief Fills the arguments with the date
			 	
			 	Give the numbers in the following order: month, day and year.
			*/
			void get(UInt& month, UInt& day, UInt& year) const;
			
			///Sets the undefined date: 00/00/0000
			void clear();
			
		protected:
	};
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATE_H
