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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_STRINGLIST_H
#define OPENMS_DATASTRUCTURES_STRINGLIST_H

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
	/**
		@brief String list
		
		This class is based on std::vector<String> but adds some methods for convenience.
		
		@ingroup Datastructures
	*/
	class StringList
		: public std::vector<String>
	{
		public:

			///@name Constructors and assignment operators
			//@{
			/// Default constructor
			StringList();
			/// Copy constructor
			StringList(const StringList& rhs);
			/// Constructor from vector<String>
			StringList(const std::vector<String>& rhs);
			/// Constructor from vector<string>
			StringList(const std::vector<std::string>& rhs);
			///  Assignment operator
			StringList& operator=(const StringList& rhs);
			///  Assignment operator from vector<String>
			StringList& operator=(const std::vector<String>& rhs);
			///  Assignment operator vector<string>
			StringList& operator=(const std::vector<std::string>& rhs);
			//@}
			
			///Operator for appending entries with less code
			template<typename StringType>
			StringList& operator<<(const StringType& string)
			{
				this->push_back(string);
				return *this;
			}

			/// Returns a list that is created by splitting the given comma-separated string (String are not trimmed!)
			static StringList create(const String& list);
			/// Returns if a string is contains in the list
			bool contains(const String& s) const;
			/// Transforms all contained strings to upper case
			void toUpper();
			/// Transforms all contained strings to lower case
			void toLower();
			
			/// output stream operator
			friend std::ostream& operator<<(std::ostream& os, const StringList& p);
			
	};
	
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_STRINGLIST_H
