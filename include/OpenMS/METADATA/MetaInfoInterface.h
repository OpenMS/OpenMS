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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_METAINFOINTERFACE_H
#define OPENMS_METADATA_METAINFOINTERFACE_H

#include <map>
#include <string>

#include <iostream>

#include <OpenMS/METADATA/MetaInfo.h>

namespace OpenMS
{
	class String;
	
	/**
		@brief Interface for classes that can store arbitrary meta information
		(Type-Name-Value tuples).
		
		MetaInfoInterface is a base class for all classes that use one MetaInfo
		object as member.  If you want to add meta information to a class, let it
		publically inherit the MetaInfoInterface.  Meta information is an array of
		Type-Name-Value tupels.
		
		@ingroup Metadata
	*/

	class OPENMS_DLLAPI MetaInfoInterface
	{
    public:

			///constructor
			MetaInfoInterface();
			///copy constructor
			MetaInfoInterface(const MetaInfoInterface& rhs);
			///destructor
			~MetaInfoInterface();
	
			///assignment operator
			MetaInfoInterface& operator = (const MetaInfoInterface& rhs);

      /// Equality operator
      bool operator== (const MetaInfoInterface& rhs) const;
      /// Equality operator
      bool operator!= (const MetaInfoInterface& rhs) const;
	
      /// returns the value corresponding to a string
			const DataValue& getMetaValue(const String& name) const;
			/// returns the value corresponding to an index
			const DataValue& getMetaValue(UInt index) const;
	
			/// returns if this MetaInfo is set
			bool metaValueExists(const String& name) const;
			/// returns if this MetaInfo is set
			bool metaValueExists(UInt index) const;
	
			/// sets the DataValue corresponding to a name
			void setMetaValue(const String& name, const DataValue& value);
			///  sets the DataValue corresponding to an index
			void setMetaValue(UInt index, const DataValue& value);
			
			/// Removes the DataValue corresponding to @p name if it exists
			void removeMetaValue(const String& name);
			/// Removes the DataValue corresponding to @p index if it exists
			void removeMetaValue(UInt index);		
			
			/// retuns a reference to the MetaInfoRegistry
		  static MetaInfoRegistry& metaRegistry();
	
	    /// fills the given vector with a list of all keys for which a value is set
	    void getKeys(std::vector<String>& keys) const;

	    /// fills the given vector with a list of all keys for which a value is set
	    void getKeys(std::vector<UInt>& keys) const;
	    
	    /// returns if the MetaInfo is empty
	    bool isMetaEmpty() const;

	    /// removes all meta values
	    void clearMetaInfo();

		protected:
			/// creates the MetaInfo object if it does not exist
			inline void createIfNotExists_();
			/// pointer to the MetaInfo object. 0 by default
			MetaInfo* meta_;
	};

} // namespace OpenMS

#endif // OPENMS_METADATA_METAINFOINTERFACE_H
