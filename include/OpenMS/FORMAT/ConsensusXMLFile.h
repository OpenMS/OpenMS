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
// $Maintainer: Eva Lange$
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
      template <typename ElementT>
      void load(const String& filename, ConsensusMap<ElementT>& map, bool load_element_maps = true) throw (Exception::FileNotFound, Exception::ParseError)
      {
        map.clear(); // clear map
        Internal::ConsensusXMLHandler< StarAlignment<ElementT> > handler(map,filename,load_element_maps);
        parse_(filename, &handler);
      }

      /// Stores a staralignment object into consensusXML format.
      template <typename AlignmentT>
      void store(const String& filename, const AlignmentT& alignment)
      const throw (Exception::UnableToCreateFile)
      {
        Internal::ConsensusXMLHandler<AlignmentT> handler(alignment,filename);
        save_(filename, &handler);
      }

  };
} // namespace OpenMS

#endif // OPENMS_FOMAT_CONSENSUSXMLFILE_H
