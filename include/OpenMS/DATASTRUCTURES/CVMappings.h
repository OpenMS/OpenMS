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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGS_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGS_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/DATASTRUCTURES/CVReference.h>

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

			/// Default constructor
			CVMappings();
	
			/// Copy constructor
			CVMappings(const CVMappings& rhs);
	
			/// Destructor
			virtual ~CVMappings();
	
			/// Assignment operator
			CVMappings& operator = (const CVMappings& rhs);
	
			/** @name Accessors
			*/
			//@{
			/// sets the mapping rules of the mapping file
			void setMappingRules(const std::vector<CVMappingRule>& cv_mapping_rules);
	
			/// returns the mapping rules
			const std::vector<CVMappingRule>& getMappingRules() const;
	
			/// adds a mapping rule
			void addMappingRule(const CVMappingRule& cv_mapping_rule);
	
			/// sets the CV references
			void setCVReferences(const std::vector<CVReference>& cv_references);

			/// returns the CV references
			const std::vector<CVReference>& getCVReferences() const;
	
			/// adds a CV reference
			void addCVReference(const CVReference& cv_reference);
			//@}
	
			/** @name Predicates
			*/
			//@{
			/// returns true if a CV reference is given
			bool hasCVReference(const String& identifier);

			/// equality operator 
			bool operator == (const CVMappings& rhs) const;
			
			/// inequality operator 
			bool operator != (const CVMappings& rhs) const;
			//@}

		protected:

			std::vector<CVMappingRule> mapping_rules_;

			Map<String, CVReference> cv_references_;

			std::vector<CVReference> cv_references_vector_;
	};
} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGS_H
