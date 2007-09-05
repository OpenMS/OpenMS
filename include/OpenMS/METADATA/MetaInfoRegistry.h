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

#ifndef OPENMS_METADATA_METAINFOREGISTRY_H
#define OPENMS_METADATA_METAINFOREGISTRY_H

#include <map>
#include <string>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

	/**
		@brief Registry which assigns unique integer indices to strings.
		
		When registering a new name an index >= 1024 is assigned.
		Indices from 1 to 1023 are reserved for fast access and will never change:<BR>
		1 - isotopic_range<BR>
		2 - cluster_id<BR>
		3 - label<BR>
		4 - icon<BR>
		5 - color<BR>
		
		@todo Automatically register new names (Marc)
		
		@ingroup Metadata
	*/
	class MetaInfoRegistry
	{
    public:	
			///default constructor
			MetaInfoRegistry();
	
			///copy constructor
			MetaInfoRegistry(const MetaInfoRegistry& rhs);
	
			///destructor
			~MetaInfoRegistry();
	
			///assignment operator
			MetaInfoRegistry& operator = (const MetaInfoRegistry& rhs);
	
			/**
				Registers a string, stores it's description and unit, and returns the corresponding index.
				If the string is already registered, it retuns the index of the string.
			*/
			UInt registerName(const String& name, const String& description, const String& unit="");
	
			///Returns the corresponding integer to a string
			UInt getIndex(const String& name) const throw(Exception::InvalidValue);
	
			///Returns the corresponding name to an index
			String getName(UInt index) const throw(Exception::InvalidValue);
	
	
			/// retuns the description of an index
			String getDescription(UInt index) const throw(Exception::InvalidValue);
			/// retuns the description of a name
			String getDescription(const String& name) const throw(Exception::InvalidValue);
	
			/// retuns the unit of an index
			String getUnit(UInt index) const throw(Exception::InvalidValue);
			/// retuns the unit of a name
			String getUnit(const String& name) const throw(Exception::InvalidValue);

		private:
			/// internal counter, that stores the next index to assign
			UInt next_index_;
			/// map from name to index
			std::map<String,UInt> name_to_index_;
			/// map from index to name
			std::map<UInt,String> index_to_name_;
			/// map from index to description
			std::map<UInt,String> index_to_description_;
			/// map from index to unit
			std::map<UInt,String> index_to_unit_;

	};

} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFOREGISTRY_H
