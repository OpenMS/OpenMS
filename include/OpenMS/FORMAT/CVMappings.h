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
// $Maintainer: Andreas Bertsch $
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
		
		@todo Docu (Andreas)
		
		@ingroup Format
	*/
	class CVMappings
	{
		public:

			///Represenation of a CV term used by CVMappings
			class CVTerm
			{
				public:
								
				CVTerm();

				CVTerm(const CVTerm& rhs);

				virtual ~CVTerm();

				CVTerm& operator = (const CVTerm& rhs);

				void setAccession(const String& accession);

				const String& getAccession() const;

				void setUseTermName(bool use_term_name);

				bool getUseTermName() const;

				void setUseTerm(bool use_term);

				bool getUseTerm() const;

				void setTermName(const String& term_name);

				const String& getTermName() const;

				void setIsRepeatable(bool is_repeatable);

				bool getIsRepeatable() const;

				void setAllowChildren(bool allow_children);

				bool getAllowChildren() const;

				void setCVIdentifierRef(const String& cv_identifier_ref);

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
						
			/// Representation of a CV Mapping rule used by CVMappings
			class CVMappingRule
			{
				public:

				enum RequirementLevel
				{
					MUST = 0,
					SHOULD = 1,
					MAY = 2
				};

				enum CombinationsLogic
				{
					OR = 0,
					AND = 1,
					XOR = 2
				};
								
				CVMappingRule();

				CVMappingRule(const CVMappingRule& rhs);

				virtual ~CVMappingRule();

				CVMappingRule& operator = (const CVMappingRule& rhs);

				void setIdentifier(const String& identifier);

				const String& getIdentifier() const;

				void setElementPath(const String& element_path);

				const String& getElementPath() const;

				void setRequirementLevel(RequirementLevel level);
				
				RequirementLevel getRequirementLevel() const; 
			
				void setCombinationsLogic(CombinationsLogic combinations_logic);

				CombinationsLogic getCombinationsLogic() const;
			
				void setScopePath(const String& path);

				const String& getScopePath() const;
				
				void setCVTerms(const std::vector<CVTerm>& cv_terms);

				const std::vector<CVTerm>& getCVTerms() const;

				void addCVTerm(const CVTerm& cv_terms);
				
				protected:

				String identifier_;

				String element_path_;

				RequirementLevel requirement_level_;

				String scope_path_;

				CombinationsLogic combinations_logic_;

				std::vector<CVTerm> cv_terms_;
			};

			class CVReference
			{
				public:

				CVReference();

				CVReference(const CVReference& rhs);

				virtual ~CVReference();

				CVReference& operator = (const CVReference& rhs);

				void setName(const String& name);

				const String& getName() const;

				void setIdentifier(const String& identifier);

				const String& getIdentifier() const;
				
				protected:

				String name_;

				String identifier_;
			};
			
		
			CVMappings();
	
			CVMappings(const CVMappings& rhs);
	
			virtual ~CVMappings();
	
			CVMappings& operator = (const CVMappings& rhs);
	
			void setMappingRules(const std::vector<CVMappingRule>& cv_mapping_rules);
	
			const std::vector<CVMappingRule>& getMappingRules() const;
	
			void addMappingRule(const CVMappingRule& cv_mapping_rule);
	
			void setCVReferences(const std::vector<CVReference>& cv_references);
	
			void addCVReference(const CVReference& cv_reference);
	
			bool hasCVReference(const String& identifier);

		protected:

			std::vector<CVMappingRule> mapping_rules_;

			Map<String, CVReference> cv_references_;
	};
} // namespace OpenMS

#endif // OPENMS_FORMAT_CVMAPPINGS_H
