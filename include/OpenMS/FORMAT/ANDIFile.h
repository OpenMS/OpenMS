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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_ANDIFILE_H
#define OPENMS_FORMAT_ANDIFILE_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/HANDLERS/ANDIHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <fstream>

namespace OpenMS
{
	class String;

  /**
  	@brief File adapter for ANDI/MS files

  	@improvement Fix bugs on 64 bit systems (Hiwi)

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI ANDIFile
  	:	public ProgressLogger
  {
    public:
      ///Default constructor
      ANDIFile();
      ///Destructor
      ~ANDIFile();

			/**
      	@brief Loads a map from a ANDI/MS file.

      	@p map has to be a MSExperiment or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      template <typename MapType>
      void load(const String& filename, MapType& map)
      {
      	//try to open file
				std::ifstream is(filename.c_str());
		    if (!is)
		    {
		      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		    }

				map.reset();

				//set DocumentIdentifier
				map.setLoadedFileType(filename);
				map.setLoadedFilePath(filename);

				Internal::ANDIHandler<MapType> handler(map,*this);
				handler.parse(filename);
      }
  };

} // namespace OpenMS

#endif // OPENMS_FOMAT_ANDIFILE_H
