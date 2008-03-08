// -*- Mode: C++; tab-width: 2; -*-
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
	class StringList:
		public std::vector<String>
	{
		public:
			//Operator for appending entries with less code
			template<typename StringType>
			StringList& operator<<(const StringType& string)
			{
				this->push_back(string);
				return *this;
			}
			/// Returns a list that contains "yes" and "no"
			static StringList getYesNoList();
	};
	
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_STRINGLIST_H
