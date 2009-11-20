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

#ifndef OPENMS_METADATA_DOCUMENTIDTAGGER_H
#define OPENMS_METADATA_DOCUMENTIDTAGGER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Tags OpenMS file containers with a DocumentID
    
    Intended usage is from within a TOPP tool. An instance of this class
		is present in TOPPBase and can be used by all derived TOPP tools
		to assign a unique ID which is fetched from an ID pool in ./share/OpenMS/IDPool/.
				
		@ingroup Metadata
  */
  class OPENMS_DLLAPI DocumentIDTagger
  {
		private:
			/// Constructor (declared away)
      DocumentIDTagger();
		
			/**
				@brief retrieve an ID from the pool

				Uses boost filelocks to safely retrieve an ID from an ID pool.
				
				@param id Unique identifier returned from ID pool
				@param free Number of available identifiers in ID pool (before this query)
				@param idcount_only Only count available identifiers, do NOT retrieve one 
							 (the id string will nevertheless be filled)

				Return true if all file operations could be executed successfully (this does not imply there was an ID left over - check free>0)
			*/
			bool getID_(String& id, Int& free, bool idcount_only) const;

		public:
			/// Constructor
      DocumentIDTagger(String toolname);
			/// Copy constructor
      DocumentIDTagger(const DocumentIDTagger& source);
      /// Destructor
      ~DocumentIDTagger();
      
      /// Assignment operator
      DocumentIDTagger& operator = (const DocumentIDTagger& source);
      
      /// Equality operator
      bool operator == (const DocumentIDTagger& source) const;
      /// Equality operator
      bool operator != (const DocumentIDTagger& source) const;


			/**
				@brief Return the file used as ID pool

				The default ID pool file is in /share/OpenMS/IDPool/IDPool.txt
				A custom file can be set by setIDPoolFile()
			*/
			String getPoolFile() const;

			/// Set the file used as ID pool
			void setPoolFile(const String& file);

			/**
				@brief Tags any structure which is derived from DocumentIdentifier with a unique tag
			
				Tags any structure which is derived from DocumentIdentifier with a unique tag
				Returns true if ID could be assigned, otherwise an Exception::DepletedIDPool is thrown

				@param map Some class (derived from a DocumentIdentifier class) which needs a unique id
				@exception Exception::DepletedIDPool when no identifier (for whatever reason) could be aquired
			*/
			bool tag(DocumentIdentifier& map) const;
			
			/**
				@brief return the number of available IDs in the pool.

				Retrieve the number of available IDs in the pool.
				Returns true of count was successful, false otherwise (locking error, file creation error ...)
				
				@param free Number of available identifiers. You should worry if it's 0!
			*/
			bool countFreeIDs(Int& free) const;

    protected:
			/// name of the calling TOPP tool
			String toolname_;

			/// location of the ID pool
			String pool_file_;
  };
 
} // namespace OpenMS

#endif // OPENMS_METADATA_DOCUMENTIDTAGGER_H

