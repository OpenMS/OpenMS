// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTXMLFILE_H
#define OPENMS_FORMAT_MASCOTXMLFILE_H

#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
	class String;
	
  /**
    @brief Used to load MascotXML files
    
    This class is used to load documents that implement 
    the schema of MascotXML files.
  
  	@ingroup FileIO
  */
  class MascotXMLFile
  {
    public:
      /// Constructor
      MascotXMLFile();
      
		  /**
		    @brief loads data from a MascotXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT

		    This class serves to read in a MascotXML file. The information can be 
		    retrieved via the load function.      
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<IdentificationData>& id_data ) const throw (Exception::FileNotFound, Exception::ParseError);
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTXMLFILE_H
