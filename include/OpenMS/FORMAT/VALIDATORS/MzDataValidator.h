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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_VALIDATORS_MZDATAVALIDATOR_H
#define OPENMS_FORMAT_VALIDATORS_MZDATAVALIDATOR_H


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
	  class OPENMS_DLLAPI MzDataValidator
			: public SemanticValidator
	  {
	    public:
	      /**
	      	@brief Constructor
	      
					@param mapping The mapping rules
					@param cv @em All controlled vocabularies required for the mapping 
				*/
	      MzDataValidator(const CVMappings& mapping, const ControlledVocabulary& cv);
				
				/// Destructor
				virtual ~MzDataValidator();
				
			protected:
				
				//Docu in base class
				void handleTerm_(const String& path, const CVTerm& parsed_term);
				
			private:
				
				/// Not implemented
				MzDataValidator();
				
				/// Not implemented
				MzDataValidator(const MzDataValidator& rhs);
	
				/// Not implemented
				MzDataValidator& operator = (const MzDataValidator& rhs);
	
	  };
 
	} // namespace Internal
 
} // namespace OpenMS

#endif

