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

#ifndef OPENMS_FORMAT_VALIDATORS_MZIDENTMLVALIDATOR_H
#define OPENMS_FORMAT_VALIDATORS_MZIDENTMLVALIDATOR_H


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
	  class OPENMS_DLLAPI MzIdentMLValidator
			: public SemanticValidator
	  {
	    public:
	      /**
	      	@brief Constructor
	      
					@param mapping The mapping rules
					@param cv @em All controlled vocabularies required for the mapping 
				*/
	      MzIdentMLValidator(const CVMappings& mapping, const ControlledVocabulary& cv);
				
				/// Destructor
				virtual ~MzIdentMLValidator();
				
			private:
				
				/// Not implemented
				MzIdentMLValidator();
				
				/// Not implemented
				MzIdentMLValidator(const MzIdentMLValidator& rhs);
	
				/// Not implemented
				MzIdentMLValidator& operator = (const MzIdentMLValidator& rhs);
	
	  };
 
	} // namespace Internal
 
} // namespace OpenMS

#endif

