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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_UNIMODXMLFILE_H
#define OPENMS_FORMAT_UNIMODXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS 
{
	class ResidueModification;

  /**
    @brief Used to load XML files from unimod.org files
    
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI UnimodXMLFile
  	: public Internal::XMLFile
  {
    public:
						
      /// Default constructor
      UnimodXMLFile();
			
			/// Destructor
			virtual ~UnimodXMLFile();
		  /**
		    @brief loads data from unimod.xml file
		   
				@param filename the filename were the unimod xml file should be read from
				@param modifications the modifications which are read from the file
				@throw FileNotFound is thrown if the file could not be found
				@throw ParseError is thrown if the given file could not be parsed
				
		  	@ingroup FileIO
		  */
	    void load(const String& filename, std::vector<ResidueModification*>& modifications);

		private:
			
			///Not implemented
			UnimodXMLFile(const UnimodXMLFile& rhs);
			///Not implemented
			UnimodXMLFile& operator = (const UnimodXMLFile& rhs);
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_UNIMODXMLFILE_H
