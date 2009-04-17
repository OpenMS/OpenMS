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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CVMAPPINGS_H
#define OPENMS_FORMAT_CVMAPPINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief Representation of controlled vocabulary mapping rules (for PSI formats)
		
		This file serves as object for the controlled vocabulary term usage definitions 
		used in CV-Mapping files. All the supported attributes supported in the 
		mapping file are supported by this class.
		
		@ingroup Format
	*/
	class OPENMS_DLLAPI CVMappings
	{
		public:

			///Represenation of a CV term used by CVMappings
			class OPENMS_DLLAPI CVTerm
			{
				public:
				
				/// Defaults constructor
				CVTerm();

				/// Copy constructor
				CVTerm(const CVTerm& rhs);

				/// Destructor
				virtual ~CVTerm();

				/// Assignment operator
				CVTerm& operator = (const CVTerm& rhs);

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
				
				protected:
				
				String accession_;

				bool use_term_name_;

				bool use_term_;

				String term_name_;

				bool is_repeatable_;

				bool allow_children_;

				String cv_identifier_ref_;
			};
						
			/** @brief Representation of a CV Mapping rule used by CVMappings
				
					Representation of a controlled vocabulary mapping rule.
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
				void setCVTerms(const std::vector<CVTerm>& cv_terms);

				/// returns the allowed terms
				const std::vector<CVTerm>& getCVTerms() const;

				/// adds a term to the allowed terms
				void addCVTerm(const CVTerm& cv_terms);
				
				protected:

				String identifier_;

				String element_path_;

				RequirementLevel requirement_level_;

				String scope_path_;

				CombinationsLogic combinations_logic_;

				std::vector<CVTerm> cv_terms_;
			};

			/* @brief Controlled Vocabulary Reference
			
				Reference to a controlled vocabulary, defined in the first section of a mapping file.

			*/
			class OPENMS_DLLAPI CVReference
			{
				public:

				/// Default constructor
				CVReference();

				/// Copy constructor
				CVReference(const CVReference& rhs);

				/// Destructor
				virtual ~CVReference();

				/// Assignment operator
				CVReference& operator = (const CVReference& rhs);

				/// sets the name of the CV reference
				void setName(const String& name);

				/// returns the name of the CV reference
				const String& getName() const;

				/// sets the CV identifier which is referenced
				void setIdentifier(const String& identifier);

				/// returns the CV identifier which is referenced
				const String& getIdentifier() const;
				
				protected:

				String name_;

				String identifier_;
			};
			
			
			/*
					@brief Representation of a complete mapping file

					This class represents the complete mapping file, consisting of CV-references and
					mapping rules.
			*/
			/// Default constructor
			CVMappings();
	
			/// Copy constructor
			CVMappings(const CVMappings& rhs);
	
			/// Destructor
			virtual ~CVMappings();
	
			/// Assignment operator
			CVMappings& operator = (const CVMappings& rhs);
	
			/// sets the mapping rules of the mapping file
			void setMappingRules(const std::vector<CVMappingRule>& cv_mapping_rules);
	
			/// returns the mapping rules
			const std::vector<CVMappingRule>& getMappingRules() const;
	
			/// adds a mapping rule
			void addMappingRule(const CVMappingRule& cv_mapping_rule);
	
			/// sets the CV references
			void setCVReferences(const std::vector<CVReference>& cv_references);
	
			/// adds a CV reference
			void addCVReference(const CVReference& cv_reference);
	
			/// returns true if a CV reference is given
			bool hasCVReference(const String& identifier);

		protected:

			std::vector<CVMappingRule> mapping_rules_;

			Map<String, CVReference> cv_references_;
	};
} // namespace OpenMS

#endif // OPENMS_FORMAT_CVMAPPINGS_H
