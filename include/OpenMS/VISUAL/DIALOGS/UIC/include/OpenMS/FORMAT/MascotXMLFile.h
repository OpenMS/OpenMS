// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTXMLFILE_H
#define OPENMS_FORMAT_MASCOTXMLFILE_H

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS 
{
	class ProteinIdentification;

  /**
    @brief Used to load MascotXML files
    
    This class is used to load documents that implement 
    the schema of MascotXML files.
  
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI MascotXMLFile
  	: public Internal::XMLFile
  {
    public:
      /// Constructor
      MascotXMLFile();
      
		  /**
		    @brief loads data from a MascotXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT
			 	@exception Exception::FileNotFound is thrown if the file does not exists.
				@exception Exception::ParseError is thrown if the file does not suit to the standard.
				
		    This method serves to read in a MascotXML file. The information can be 
		    retrieved via the load function.      
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data);
      					 
		  /**
		    @brief loads data from a MascotXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT
		    @param peptides a map of modified peptides identified by the String title
			 	@exception Exception::FileNotFound is thrown if the file does not exists.
				@exception Exception::ParseError is thrown if the file does not suit to the standard.

		    This method serves to read in a MascotXML file. The information can be 
		    retrieved via the load function.      
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, std::map< String, std::vector<AASequence> >& peptides);
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTXMLFILE_H
