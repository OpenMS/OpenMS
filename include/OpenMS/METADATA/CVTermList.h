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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_CVTERMLIST_H
#define OPENMS_METADATA_CVTERMLIST_H

#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

namespace OpenMS
{
	/**
		@brief Representation of controlled vocabulary term list
		
		This class should be used to inherit from, to allow to add
		an arbitrary number of CV terms to the inherited class
		
		@ingroup Metadata
	*/
			///Represenation of a CV term used by CVMappings
			class OPENMS_DLLAPI CVTermList
				:	public MetaInfoInterface
			{
				public:
			
				/// Defaults constructor
				CVTermList();

				/// Copy constructor
				CVTermList(const CVTermList& rhs);

				/// Destructor
				virtual ~CVTermList();

				/// Assignment operator
				CVTermList& operator = (const CVTermList& rhs);

				/** @name Accessors
				*/
				//@{
				/// sets the CV terms 
				void setCVTerms(const std::vector<CVTerm>& terms);

				/// returns the accession string of the term
				const Map<String, std::vector<CVTerm> >& getCVTerms() const;

				/// adds a CV term
				void addCVTerm(const CVTerm& term);

				/// checks whether the spellings of the CV terms stored are correct
				bool checkCVTerms(const ControlledVocabulary& cv) const;

				/// corrects the CVTerm names, according to the loaded CV
				void correctCVTermNames();
				//@}
			
				/** @name Predicates
				*/
				//@{
				/// equality operator 
				bool operator == (const CVTermList& cv_term_list) const;

				/// inequality operator 
				bool operator != (const CVTermList& cv_term_list) const;
				
				/// checks whether the term has a value
				bool hasCVTerm(const String& accession) const;
				
				/// checks whether the stored terms fullfil a given CVMappingRule
				bool checkCVTerms(const CVMappingRule& rule, const ControlledVocabulary& cv) const;

				/// return true if no terms are available
				bool empty() const;
				//}

				protected:
				
				Map<String, std::vector<CVTerm> > cv_terms_;

			};
						
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVTERMLIST_H
