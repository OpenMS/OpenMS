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

#ifndef OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
#define OPENMS_FORMAT_CONTROLLEDVOCABULARY_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <map>

namespace OpenMS
{
	/**
		@brief Representation of a controlled vocabulary.
		
		This represenation only contains the information used for parsing and validation.
	*/
	class ControlledVocabulary
	{
		friend std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv);
			
		public:
			/// Representation of a CV term
			struct CVTerm
			{
				String name;									///< Text name
				String id;										///< Identifier
				std::vector<String> parents;	///< List of parent IDs
				bool obsolete; 								///< Flag that indicates of the term is obsolete
				
				///Default constructor
				CVTerm()
					: name(),
						id(),
						parents(),
						obsolete(false)
				{
				}
			};
			
			/// Constructor
			ControlledVocabulary();
			
			///Destructor
			virtual ~ControlledVocabulary();
			
			/// Returns the CV name (set in the load method)
			const String& name() const;
			
			///load the CV from a OBO file
			void loadFromOBO(const String& name, const String& filename) throw (Exception::FileNotFound, Exception::ParseError);
			
			/// Returns true if the term is in the CV. Returns false otherwise.
			bool exists(const String& id) const;
			
			///Returns a term specified by ID
			const CVTerm& getTerm(const String& id) const throw (Exception::InvalidValue);

			///Returns if @p child is a child of @p parent
			bool isChildOf(const String& child, const String& parent) const  throw (Exception::InvalidValue);
			
		protected:
			///Map from ID to CVTerm
			std::map<String, CVTerm> terms_;
			///Name set in the load method
			String name_;			
	};
	
	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const ControlledVocabulary& cv);
	
	
} // namespace OpenMS

#endif // OPENMS_FORMAT_CONTROLLEDVOCABULARY_H
