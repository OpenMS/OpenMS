// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_METAINFOREGISTRY_H
#define OPENMS_METADATA_METAINFOREGISTRY_H

#include <map>
#include <string>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( push )
	#pragma warning( disable : 4251 ) // disable MSVC dll-interface warning
#endif

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
		6 - RT<BR>
		7 - MZ<BR>
		8 - predicted_RT<BR>
		9 - predicted_RT_p_value<BR>
		10 - spectrum_reference<BR>
		11 - ID<BR>
		12 - low_quality<BR>
		13 - charge<BR>

		@ingroup Metadata
	*/
	class OPENMS_DLLAPI MetaInfoRegistry
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
				Registers a string, stores its description and unit, and returns the corresponding index.
				If the string is already registered, it returns the index of the string.

			  @note This method is const, because getIndex(..) const must be able to call this method if the requested
			  string is not registered yet. Therefor all changed fields are declared mutable.
			*/
			UInt registerName(const String& name, const String& description, const String& unit="") const;
	
			/**
			  @brief Sets the description (String), corresponding to an index

			  @exception Exception::InvalidValue is thrown for unregistered indices
			*/
			void setDescription(UInt index, const String& description);
			
			/**
			  @brief Sets the description (String), corresponding to a name

			  @exception Exception::InvalidValue is thrown for unregistered names
			*/
			void setDescription(const String& name, const String& description);
			
			/**
			  @brief Sets the unit (String), corresponding to an index

			  @exception Exception::InvalidValue is thrown for unregistered indices
			*/
			void setUnit(UInt index, const String& unit);
			
			/**
			  @brief Sets the unit (String), corresponding to a name

			  @exception Exception::InvalidValue is thrown for unregistered names
			*/
			void setUnit(const String& name, const String& unit);
			
			/**
			  Returns the corresponding integer to a string. If the string is not registered yet, it 
			  registers the string (with empty description and empty unit) and returns the corresponding index.
			*/
			UInt getIndex(const String& name) const;
	
			/**
			  @brief Returns the corresponding name to an index

			  @exception Exception::InvalidValue is thrown for unregistered indices
			*/
			String getName(UInt index) const;
	
			/**
			  @brief returns the description of an index

			  @exception Exception::InvalidValue is thrown for unregistered indices
			*/
			String getDescription(UInt index) const;
			/**
			  @brief returns the description of a name

			  @exception Exception::InvalidValue is thrown for unregistered names
			*/
			String getDescription(const String& name) const;
	
			/**
			  @brief returns the unit of an index

			  @exception Exception::InvalidValue is thrown for unregistered indices
			*/
			String getUnit(UInt index) const;
			/**
			  @brief returns the unit of a name

			  @exception Exception::InvalidValue is thrown for unregistered names
			*/
			String getUnit(const String& name) const;

		private:
			/// internal counter, that stores the next index to assign
			mutable UInt next_index_;
			/// map from name to index
			mutable std::map<String,UInt> name_to_index_;
			/// map from index to name
			mutable std::map<UInt,String> index_to_name_;
			/// map from index to description
			mutable std::map<UInt,String> index_to_description_;
			/// map from index to unit
			mutable std::map<UInt,String> index_to_unit_;

	};

} // namespace OpenMS

#ifdef OPENMS_COMPILER_MSVC
	#pragma warning( pop ) 
#endif

#endif // OPENMS_METADATA_METAINFOREGISTRY_H
