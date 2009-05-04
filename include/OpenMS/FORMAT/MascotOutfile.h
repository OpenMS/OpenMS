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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MASCOTOUTFILE_H
#define OPENMS_FORMAT_MASCOTOUTFILE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
	class ProteinIdentification;

  /**
    @brief Representation of a Mascot output file
    
		!WARNING!

		This class is only provided for convenience. The MascotAdapter reads the 
		output from Mascot via MascotXML and pepXML. E.g. This class does not 
		provide support for modifications. MascotXML and pepXML can be exported
		from Mascot using export*.pl scripts.
		
    This class serves to read in a Mascot outfile. The information can be 
    retrieved via the load function.
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI MascotOutfile
  {
    public:
      
      /// Constructor
      MascotOutfile();

		  /**
		    @brief loads data from a Mascot outfile
		    
		    @param filename the file to be loaded
		    @param protein_identification the protein identification
		    @param peptide_identifications the peptide identifications
				@param p the significance level (for the protein hit scores)
				@throw ParseError is thrown if the file could not be parsed

		    This class serves to read in a Mascot outfile. The information can be 
		    retrieved via the load function. This class is only contained to be compatible with previous versions.
		    You should use the MascotXMLFile instead. 
		  */
			void load(String filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& peptide_identifications, Real p = 0.05);

    protected:

   };

} //namespace OpenMS

#endif // OPENMS_FORMAT_MASCOTOUTFILE_H
