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

#ifndef OPENMS_FORMAT_CVMAPPINGFILE_H
#define OPENMS_FORMAT_CVMAPPINGFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/CVMappings.h>

#include <vector>

namespace OpenMS 
{
	class String;
	
  /**
    @brief Used to load CvMapping files
 
		This file contains the mapping of CV terms to the schema, which
		is used by PSI standard formats to semantically validate files.
  
  	@ingroup FileIO
  */
  class CVMappingFile
		: protected Internal::XMLHandler,
			public Internal::XMLFile
  {
    public:
						
      /// Default constructor
      CVMappingFile();
			
			/// Destructor
			virtual ~CVMappingFile();
		  
			/**
				@brief loads CvMappings from the given file
				
				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
				@param strip_namespaces if enable, namespace definitions of the paths are eliminated, e.g. 'pf:cvParam' -> 'cvParam'
			*/
	    void load(const String& filename, CVMappings& cv_mappings, bool strip_namespaces = false);
			
		protected:

			// Docu in base class
			void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
			void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

      // Docu in base class
      void characters(const XMLCh* const chars, const unsigned int /*length*/);

		private:
			
			///Not implemented
			CVMappingFile(const CVMappingFile& rhs);

			///Not implemented
			CVMappingFile& operator = (const CVMappingFile& rhs);

			String tag_;

			bool strip_namespaces_;
			
			CVMappings::CVMappingRule actual_rule_;

			std::vector<CVMappings::CVMappingRule> rules_;

			std::vector<CVMappings::CVReference> cv_references_;
			
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_CVMAPPINGFILE_H
