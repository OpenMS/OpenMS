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
#include <OpenMS/FORMAT/Serialization.h>

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
		
		@ingroup Metadata, Serialization
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
			UnsignedInt registerName(const std::string& name, const std::string& description, const std::string& unit="");
	
			///Returns the corresponding integer to a string
			UnsignedInt getIndex(const std::string& name) const throw(Exception::InvalidValue);
	
			///Returns the corresponding name to an index
			std::string getName(UnsignedInt index) const throw(Exception::InvalidValue);
	
	
			/// retuns the description of an index
			std::string getDescription(UnsignedInt index) const throw(Exception::InvalidValue);
			/// retuns the description of a name
			std::string getDescription(const std::string& name) const throw(Exception::InvalidValue);
	
			/// retuns the unit of an index
			std::string getUnit(UnsignedInt index) const throw(Exception::InvalidValue);
			/// retuns the unit of a name
			std::string getUnit(const std::string& name) const throw(Exception::InvalidValue);

		private:
			/// internal counter, that stores the next index to assign
			UnsignedInt next_index_;
			/// map from name to index
			std::map<std::string,UnsignedInt> name_to_index_;
			/// map from index to name
			std::map<UnsignedInt,std::string> index_to_name_;
			/// map from index to description
			std::map<UnsignedInt,std::string> index_to_description_;
			/// map from index to unit
			std::map<UnsignedInt,std::string> index_to_unit_;

		///@name Serialization
		//@{
	 private:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{ 
			ar & boost::serialization::make_nvp("name_to_index_",name_to_index_);
			ar & boost::serialization::make_nvp("index_to_name_",index_to_name_);
			ar & boost::serialization::make_nvp("index_to_description_",index_to_description_);
			ar & boost::serialization::make_nvp("index_to_unit_",index_to_unit_);
    }
		//@}

		/// Serialization
		friend class boost::serialization::access;

	};

} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFOREGISTRY_H
