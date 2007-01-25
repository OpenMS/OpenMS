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

#ifndef OPENMS_METADATA_METAINFO_H
#define OPENMS_METADATA_METAINFO_H

#include <map>
#include <string>
#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{

	/**
		@brief A Type-Name-Value tuple class.
		
		MetaInfo maps an index ( an integer corresponding to a string ) to DataValue objects.
		The mapping of strings to the index is perfomed by the MetaInfoRegistry,
		which can be accessed by the method registry().
		
		There are two versions of nearly all members. One which operates with a string name and another
		one which operates on an index. The index version is always faster, as it does not need to look
		up the index corresponding to the the string in the MetaInfoRegistry.
		
		If you wish to add one MetaInfo member to a class, consider deriving that class from 
		MetaInfoInterface, instead of simply adding MetaInfo as member. MetaInfoInterface implements
		a full interface to a MetaInfo member.
		
		@ingroup Metadata, Serialization
	*/
	class MetaInfo
	{
    public:
		///constructor
		MetaInfo();

		///copy constructor
		MetaInfo(const MetaInfo& rhs);

		///destructor
		~MetaInfo();

		///assignment operator
		MetaInfo& operator = (const MetaInfo& rhs);

    /// Equality operator
    bool operator== (const MetaInfo& rhs) const;
    /// Equality operator
    bool operator!= (const MetaInfo& rhs) const;

		/// returns the value corresponding to a string
		const DataValue& getValue(const std::string& name) const;
		/// returns the value corresponding to an index
		const DataValue& getValue(UnsignedInt index) const;

		/// returns if this MetaInfo is set
		bool exists(const std::string& name) const;
		/// returns if this MetaInfo is set
		bool exists(UnsignedInt index) const;

		/// sets the value (string) corresponding to a name
		void setValue(const std::string& name, const std::string& value);
		/// sets the value (string) corresponding to an index
		void setValue(UnsignedInt index, const std::string& value);
		/// sets the value (integer) corresponding to a name
		void setValue(const std::string& name, SignedInt value);
		/// sets the value (integer) corresponding to an index
		void setValue(UnsignedInt index, SignedInt value);
		/// sets the value (double) corresponding to a name
		void setValue(const std::string& name, double value);
		/// sets the value (double) corresponding to an index
		void setValue(UnsignedInt index, double value);
		/// sets the DataValue corresponding to a name
		void setValue(const std::string& name, const DataValue& value);
		///  sets the DataValue corresponding to an index
		void setValue(UnsignedInt index, const DataValue& value);
		
		/// Removes the DataValue corresponding to @p name if it exists
		void removeValue(const std::string& name);
		/// Removes the DataValue corresponding to @p index if it exists
		void removeValue(UnsignedInt index);		
		
		/// retuns a reference to the MetaInfoRegistry
		static MetaInfoRegistry& registry();
    
    /// fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<std::string>& keys) const;

		/// fills the given vector with a list of all keys for which a value is set
    void getKeys(std::vector<UnsignedInt>& keys) const;

    /// returns if the MetaInfo is empty
    bool empty() const;
    
    /// removes all meta values
    void clear();
    
		private:
		/// static MetaInfoRegistry
		static MetaInfoRegistry registry_;
		/// the actual mapping of index to the DataValue
		std::map<UnsignedInt,DataValue> index_to_value_;

		///@name Serialization
		//@{
	 private:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{
      ar & boost::serialization::make_nvp("metainfo_registry",registry_);
      ar & boost::serialization::make_nvp("index_to_value",index_to_value_);
		}
		//@}

		/// Serialization
		friend class boost::serialization::access;


	};

} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFO_H
