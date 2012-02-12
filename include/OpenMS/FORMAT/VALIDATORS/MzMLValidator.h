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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_VALIDATORS_MZMLVALIDATOR_H
#define OPENMS_FORMAT_VALIDATORS_MZMLVALIDATOR_H


#include <OpenMS/FORMAT/VALIDATORS/SemanticValidator.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>


namespace OpenMS 
{
	class ControlledVocabulary;
	namespace Internal
	{
		
	  /**
	    @brief Semantically validates MzXML files.
	  */
	  class OPENMS_DLLAPI MzMLValidator
			: public SemanticValidator
	  {
	    public:
	      /**
	      	@brief Constructor
	      
					@param mapping The mapping rules
					@param cv @em All controlled vocabularies required for the mapping 
				*/
	      MzMLValidator(const CVMappings& mapping, const ControlledVocabulary& cv);
				
				/// Destructor
				virtual ~MzMLValidator();
				
			protected:
				
				// Docu in base class
				void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

				// Docu in base class
				virtual String getPath_(UInt remove_from_end = 0) const;
				
				// Docu in base class
				virtual void handleTerm_(const String& path, const CVTerm& parsed_term); 
				
				///CV terms which can have a value (term => value type)
				Map<String,std::vector<CVTerm> > param_groups_;
				
				///Current referenceableParamGroup identifier
				String current_id_;
				
				///Binary data array name
				String binary_data_array_;
				///Binary data array type
				String binary_data_type_;
				
			private:
				
				/// Not implemented
				MzMLValidator();
				
				/// Not implemented
				MzMLValidator(const MzMLValidator& rhs);
	
				/// Not implemented
				MzMLValidator& operator = (const MzMLValidator& rhs);
	
	  };
 
	} // namespace Internal
 
} // namespace OpenMS

#endif

