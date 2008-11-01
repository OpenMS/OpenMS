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

#ifndef OPENMS_FORMAT_SEMANTICVALIDATOR_H
#define OPENMS_FORMAT_SEMANTICVALIDATOR_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/FORMAT/CVMappings.h>


namespace OpenMS 
{
	class ControlledVocabulary;
	
  /**
    @brief Semantically validates XML files using CVMappings and a ControlledVocabulary
		
  	@ingroup Format
  */
  class SemanticValidator
		: protected Internal::XMLHandler,
			public Internal::XMLFile
  {
    public:
    	/// Description of a CV term and its path in the XML instance file
			struct ValiationLocation
			{
				///Path in the XML instance file
				String path;
				/// CV term accession
				String accession;
				/// CV term name
				String name;
				/// CV term value
				String value;
			};
			
			/// Output container for validation results
			struct ValidationOutput
			{
				///Terms used that are not defined in the CV
				std::vector<ValiationLocation> unknown_terms;
					
				///Obolete terms used
				std::vector<ValiationLocation> obsolete_terms;
				
				///Terms used in the wrong schema location
				std::vector<ValiationLocation> invalid_location;

				///Terms used in locations for which no mapping rule exists
				std::vector<ValiationLocation> no_mapping;

				///Identifiers of violated rules (requirement level or combination logic)
				std::vector<String> violated;
				
				///Identifiers of violated rules (number of repeats)
				std::vector<String> violated_repeats;
			};
			
      /**
      	@brief Constructor
      
				@param mapping The mapping rules
				@param cv @em All controlled vocabularies required for the mapping 
			*/
      SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv);
			
			/// Destructor
			virtual ~SemanticValidator();
		  
			/**
				@brief semantically validates an XML file.
				
				@param filename The file to validate.s
				@param output If the validation failed, the errors are listed in this output parameter.
				
				@return true if the validation was successfull, false otherwise.
				
				@exception Exception::FileNotFound is thrown if the file could not be opened
			*/
	    bool validate(const String& filename, ValidationOutput& output);
			
			/// Sets the CV parameter tag name (default: 'cvParam')
			void setTag(const String& tag);
			
			/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'accession')
			void setAccessionAttribute(const String& accession);

			/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'name')
			void setNameAttribute(const String& name);

			/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'value')
			void setValueAttribute(const String& value);
			
		protected:

			// Docu in base class
			void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
			void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

      // Docu in base class
      void characters(const XMLCh* const chars, const unsigned int /*length*/);

			/// Reference to the mappings
			const CVMappings& mapping_;
			
			/// Reference to the CVs
			const ControlledVocabulary& cv_;
			
			/// Validation result
			ValidationOutput output_;
			
			/// List of open tags
			StringList open_tags_;
			
			/// Flag that indicates if the instance document is valid
			bool valid_;
			
			/// Rules (location => rule)
			Map<String, std::vector<CVMappings::CVMappingRule> > rules_;
			
			/// Fulfilled rules (location => rule ID => term ID => term count )
			/// When a tag is closed, the fulfilled rules of the current location are checked against the required rules
			/// The fulfilled rules for that location are then deleted.
			Map<String, Map< String, Map< String, UInt > > > fulfilled_;


			///@name Tag and attribute names
			//@{
			String cv_tag_;
			String accession_att_;
			String name_att_;
			String value_att_;
			//@}

		private:
			
			/// Not implemented
			SemanticValidator();
			
			/// Not implemented
			SemanticValidator(const SemanticValidator& rhs);

			/// Not implemented
			SemanticValidator& operator = (const SemanticValidator& rhs);

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_SEMANTICVALIDATOR_H
