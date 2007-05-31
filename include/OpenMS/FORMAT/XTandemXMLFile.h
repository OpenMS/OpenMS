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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMXMLFILE_H
#define OPENMS_FORMAT_XTANDEMXMLFILE_H

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS 
{
	class String;
	
  /**
    @brief Used to load XTandemXML files
    
    This class is used to load documents that implement 
    the schema of XTandemXML files.
  
  	@todo Handle Modifications (Andreas)
  	
  	@ingroup FileIO
  */
  class XTandemXMLFile
  {
    public:
						
      /// Default constructor
      XTandemXMLFile();
			
			/// Destructor
			virtual ~XTandemXMLFile();
		  /**
		    @brief loads data from a XTandemXML file
		    
		    @param filename the file to be loaded
		    @param protein_identification protein identifications belonging to the whole experiment
		    @param id_data the identifications with m/z and RT

		    This class serves to read in a XTandemXML file. The information can be 
		    retrieved via the load function.      
		  	
		  	@ingroup FileIO
		  */
	    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data) const throw (Exception::FileNotFound, Exception::ParseError);

		protected:

			XTandemXMLFile(const XTandemXMLFile& rhs);

			XTandemXMLFile& operator = (const XTandemXMLFile& rhs);
      					 
  };
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMXMLFILE_H
