// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_CONSENSUSXMLFILE_H
#define OPENMS_FORMAT_CONSENSUSXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>

namespace OpenMS
{
  /**
    @brief This class provides Input functionality for ConsensusMaps and Output functionality for 
    StarAlignments.
     
		This class can be used to load the content of a consensusXML file into a ConsensusMap 
		or to save the content of a StarAlignment object into an XML file.
		
		A documented schema for this format can be found at http://open-ms.sourceforge.net/schemas/.
		
		@todo Implement and test all PeakFileOptions (Hiwi, Marc, Clemens)
		
    @ingroup FileIO
  */
  class ConsensusXMLFile 
  	: public Internal::XMLFile
  {
    public:
      ///Default constructor
      ConsensusXMLFile();
      ///Destructor
      ~ConsensusXMLFile();


      /// Loads a consenus map from a ConsensusXML file.
      void load(const String& filename, ConsensusMap& map) throw (Exception::FileNotFound, Exception::ParseError)
      {
        map.clear(); // clear map
        Internal::ConsensusXMLHandler handler(map,filename,schema_version_);
        handler.setOptions(options_);
        parse_(filename, &handler);
      }

      /**
      	@brief Stores a staralignment object into consensusXML format.
      
      	@exception Exception::UnableToCreateFile is thrown if the file name is not writable
      	@exception Exception::IllegalArgument is thrown if the consensus map is not valid
      */
      void store(const String& filename, const ConsensusMap& map)
      {
      	if (!map.isValid())
      	{
      		throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Invalid consensus map cannot be stored!");
      	}
        Internal::ConsensusXMLHandler handler(const_cast<ConsensusMap&>(map),filename,schema_version_);
        save_(filename, &handler);
      }

      /// Mutable access to the options for loading/storing 
      PeakFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing 
      const PeakFileOptions& getOptions() const;
  
		protected:
		
			/// options for reading / writing
			PeakFileOptions options_;
		
  };
} // namespace OpenMS

#endif // OPENMS_FOMAT_CONSENSUSXMLFILE_H
