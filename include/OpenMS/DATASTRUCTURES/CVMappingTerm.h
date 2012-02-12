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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
	/**
		@brief Representation of controlled vocabulary term
		
		This class simply stores CV terms read from e.g. an OBO-file
		
		@ingroup Datastructures
	*/
			///Represenation of a CV term used by CVMappings
			class OPENMS_DLLAPI CVMappingTerm
			{
				public:
				
				/// Defaults constructor
				CVMappingTerm();

				/// Copy constructor
				CVMappingTerm(const CVMappingTerm& rhs);

				/// Destructor
				virtual ~CVMappingTerm();

				/// Assignment operator
				CVMappingTerm& operator = (const CVMappingTerm& rhs);

				/** @name Accessors
				*/
				//@{
				/// sets the accession string of the term
				void setAccession(const String& accession);

				/// returns the accession string of the term
				const String& getAccession() const;

				/// sets whether the term name should be used, instead of the accession
				void setUseTermName(bool use_term_name);

				/// returns whether the term name should be used, instead of the accession
				bool getUseTermName() const;

				/// sets whether the term itself can be used (or only its children)
				void setUseTerm(bool use_term);

				/// returns true if the term can be used, false if only children are allowed
				bool getUseTerm() const;

				/// sets the name of the term
				void setTermName(const String& term_name);

				/// returns the name of the term
				const String& getTermName() const;

				/// sets whether this term can be repeated
				void setIsRepeatable(bool is_repeatable);

				/// returns true if this term can be repeated, false otherwise
				bool getIsRepeatable() const;

				/// sets whether children of this term are allowed
				void setAllowChildren(bool allow_children);

				/// returns true if the children of this term are allowed to be used
				bool getAllowChildren() const;

				/// sets the cv identifier reference string, e.g. UO for unit obo 
				void setCVIdentifierRef(const String& cv_identifier_ref);

				/// returns the cv identifier reference string
				const String& getCVIdentifierRef() const;
				//@}
			
				/** @name Predicates
				*/
				//@{
				/// equality operator 
				bool operator == (const CVMappingTerm& rhs) const;
				
				/// inequality operator 
				bool operator != (const CVMappingTerm& rhs) const;
				//}

				protected:
				
				String accession_;

				bool use_term_name_;

				bool use_term_;

				String term_name_;

				bool is_repeatable_;

				bool allow_children_;

				String cv_identifier_ref_;
			};
						
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H
