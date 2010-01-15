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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PTMXMLFILE_H
#define OPENMS_FORMAT_PTMXMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>

#include <map>
#include <vector>

namespace OpenMS 
{
	/**
		@brief Used to load and store PTMXML files
		
		This class is used to load and store documents that implement
		the schema of PTMXML files.
		
		@ingroup FileIO
	*/
	class OPENMS_DLLAPI PTMXMLFile
		: public Internal::XMLFile
	{
		public:
			/// Constructor
			PTMXMLFile();
			
			/**
				@brief Loads the informations of a PTMXML file
			
				@param filename The name of the file
				@param ptm_informations the PTM information from the file are stored herein
				@throw FileNotFound is thrown if the given file could not be found
				@throw ParseError is thrown if the given file could not be parsed
				The information is read in and stored in the corresponding variables
			*/
			void load(const String& filename, std::map< String, std::pair< String, String > >& ptm_informations);
			
			/**
				@brief Stores the data in an PTMXML file
				
				@throw UnableToCreateFile is thrown if the given filename could not be created

				The data is read in and stored in the file 'filename'.
			*/
			void store(String filename, std::map< String, std::pair< String, String > >& ptm_informations) const;
	};
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_PTMXMLFILE_H
