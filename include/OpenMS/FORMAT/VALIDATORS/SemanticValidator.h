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
	namespace Internal
	{
		
	  /**
	    @brief Semantically validates XML files using CVMappings and a ControlledVocabulary.
	    
	    This is the general validator.
	    @n Specialized validators for specific file formats are derived from this class.
	  */
	  class SemanticValidator
			: protected Internal::XMLHandler,
				public Internal::XMLFile
	  {
	    public:
	      /**
	      	@brief Constructor
	      
					@param mapping The mapping rules
					@param cv @em All controlled vocabularies required for the mapping 
				*/
	      SemanticValidator(const CVMappings& mapping, const ControlledVocabulary& cv);
				
				/// Destructor
				virtual ~SemanticValidator();
			  
				/**
					@brief Semantically validates an XML file.
					
					@param filename The file to validate
					@param errors Errors during the validation are returned in this output parameter.
					@param warnings Warnings during the validation are returned in this output parameter.
					
					@return @em true if the validation was successfull, @em false otherwise.
					
					@exception Exception::FileNotFound is thrown if the file could not be opened
				*/
		    bool validate(const String& filename, StringList& errors, StringList& warnings);
				
				/// Sets the CV parameter tag name (default: 'cvParam')
				void setTag(const String& tag);
				
				/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'accession')
				void setAccessionAttribute(const String& accession);
	
				/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'name')
				void setNameAttribute(const String& name);
	
				/// Sets the name of the attribute for accessions in the CV parameter tag name (default: 'value')
				void setValueAttribute(const String& value);
				
			protected:
				
				///Representation of a parsed CV term
				struct CVTerm
				{
					String accession;
					String name;
					String value;
				};
					
				// Docu in base class
				void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
				
				// Docu in base class
				void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
	
	      // Docu in base class
	      void characters(const XMLCh* const chars, const unsigned int /*length*/);
				
				/// Returns the current element path
				virtual String getPath_(UInt remove_from_end = 0) const;
				
				/// Parses the CV term accession (required), name (required) and value (optional) from the XML attributes
				virtual void getCVTerm_(const xercesc::Attributes& attributes, CVTerm& parsed_term);
				
				/// Handling of the term
				virtual void handleTerm_(const String& path, const CVTerm& parsed_term); 

				/// Reference to the mappings
				const CVMappings& mapping_;
				
				/// Reference to the CVs
				const ControlledVocabulary& cv_;
				
				/// Validation erros
				StringList errors_;
				
				/// Validation warnings
				StringList warnings_;
				
				/// List of open tags
				StringList open_tags_;
				
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
 
	} // namespace Internal
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_SEMANTICVALIDATOR_H
