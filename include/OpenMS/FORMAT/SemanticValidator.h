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

namespace OpenMS 
{
	class CVMappings;
	class ControlledVocabulary;
	
  /**
    @brief Semantically validates XML files using a CVMapping 

  	@ingroup Format
  */
  class SemanticValidator
		: protected Internal::XMLHandler,
			public Internal::XMLFile
  {
    public:
						
      /// Default constructor
      SemanticValidator();
			
			/// Destructor
			virtual ~SemanticValidator();
		  
			/**
				@brief semantically validates an XML file
				
				@param filename The file to validate
				@param mapping The mapping rules
				@param cv @em All controlled vocabularies required for the mapping 
			*/
	    bool validate(const String& filename, const CVMappings& mapping, const ControlledVocabulary& cv);
			
		protected:

			// Docu in base class
			void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
			void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

      // Docu in base class
      void characters(const XMLCh* const chars, const unsigned int /*length*/);

		private:
			
			/// Not implemented
			SemanticValidator(const SemanticValidator& rhs);

			/// Not implemented
			SemanticValidator& operator = (const SemanticValidator& rhs);

  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMXMLFILE_H
