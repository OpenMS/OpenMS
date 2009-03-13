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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IDTAGGER_H
#define OPENMS_METADATA_IDTAGGER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Tags OpenMS file containers with a DocumentID
    
    Intented usage is from within a TOPP tool. An instance of this class
		is present in TOPPBase and can be used by all derived TOPP tools
		to assign a unique ID which is fetched from an ID pool in ./share/OpenMS/IDPool/.
				
		@ingroup Metadata
  */
  class OPENMS_DLLAPI IDTagger
  {
		private:
			/// Constructor (declared away)
      IDTagger();
		
			/// retrieve an ID from the pool
			bool getID_(String& id, Int& free, bool idcount_only) const;

		public:
			/// Constructor
      IDTagger(String toolname);
			/// Copy constructor
      IDTagger(const IDTagger& source);
      /// Destructor
      ~IDTagger();
      
      /// Assignment operator
      IDTagger& operator = (const IDTagger& source);
      
      /// Equality operator
      bool operator == (const IDTagger& source) const;
      /// Equality operator
      bool operator != (const IDTagger& source) const;
			
			/**
				@brief Tags any structure which is derived from DocumentIdentifier with a unique tag
			
				@exception Exception::DepletedIDPool is thrown if the ID pool is empty
			*/
			bool tag(DocumentIdentifier& map) const;
			
			/**
				@brief return the number of available IDs in the pool.

				Returns 0 if pool is depleted or not available and a positive integer otherwise.
			*/
			bool countFreeIDs(Int& free) const;

    protected:
			/// name of the calling TOPP tool
			String toolname_;
  };
 
} // namespace OpenMS

#endif // OPENMS_METADATA_IDTAGGER_H

