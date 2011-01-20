// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/CVMappingTerm.h>
#include <vector>

namespace OpenMS
{
	/**
		@brief Representation of a CV Mapping rule used by CVMappings
		
		Representation of a controlled vocabulary mapping rule.

		@ingroup Datastructures
	*/
			class OPENMS_DLLAPI CVMappingRule
			{
				public:

				/// enum to specify the requirement level
				enum RequirementLevel
				{
					MUST = 0,
					SHOULD = 1,
					MAY = 2
				};

				/// enum to specify the combination operator
				enum CombinationsLogic
				{
					OR = 0,
					AND = 1,
					XOR = 2
				};
								
				/// Default constructor
				CVMappingRule();

				/// Copy constructor
				CVMappingRule(const CVMappingRule& rhs);

				/// Destructor
				virtual ~CVMappingRule();

				/// Assignment operator
				CVMappingRule& operator = (const CVMappingRule& rhs);

				/** @name Accessors
				*/
				//@{
				/// sets the identifier of the rule
				void setIdentifier(const String& identifier);

				/// returns the identifier of the rule
				const String& getIdentifier() const;

				/// sets the path of the element, where this rule is allowed
				void setElementPath(const String& element_path);

				/// returns the path of the element, where this rule is allowed
				const String& getElementPath() const;
	
				/// sets the requirement level of this rule
				void setRequirementLevel(RequirementLevel level);
				
				/// returns the requirement level of this rule
				RequirementLevel getRequirementLevel() const; 
			
				/// sets the combination operator of the rule
				void setCombinationsLogic(CombinationsLogic combinations_logic);

				/// returns the combinations operator of the rule
				CombinationsLogic getCombinationsLogic() const;
			
				/// sets the scope path of the rule
				void setScopePath(const String& path);

				/// returns the scope path of the rule
				const String& getScopePath() const;
				
				/// sets the terms which are allowed
				void setCVTerms(const std::vector<CVMappingTerm>& cv_terms);

				/// returns the allowed terms
				const std::vector<CVMappingTerm>& getCVTerms() const;

				/// adds a term to the allowed terms
				void addCVTerm(const CVMappingTerm& cv_terms);
				//@}

				/** @name Predicates
				*/
				//@{
				/// equality operator
				bool operator == (const CVMappingRule& rhs) const;

				/// inequality operator
				bool operator != (const CVMappingRule& rhs) const;
				//@}
				
				protected:

				String identifier_;

				String element_path_;

				RequirementLevel requirement_level_;

				String scope_path_;

				CombinationsLogic combinations_logic_;

				std::vector<CVMappingTerm> cv_terms_;
			};

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGRULE_H
