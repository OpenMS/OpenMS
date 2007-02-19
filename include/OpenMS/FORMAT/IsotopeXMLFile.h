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
// $Maintainer: Martin Langwisch $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_IsotopeXMLFILE_H
#define OPENMS_FORMAT_IsotopeXMLFILE_H

#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS 
{
	/**
		@brief Used to load and store IsotopeXML files
		
		This class is used to load and store documents that implement
		the schema of IsotopeXML files.
		
		@ingroup FileIO
	*/
	class IsotopeXMLFile
	{
		public:
			/// Constructor
			IsotopeXMLFile();
			
			/**
				@brief Loads the informations of a IsotopeXML file
				
				The information is read in and stored in the corresponding variables
			*/
			void load(const String& filename, std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& isotope_informations) const throw (Exception::FileNotFound, Exception::ParseError);
			
			/**
				@brief Stores the data in an IsotopeXML file
				
				The data is read in and stored in the file 'filename'.
			*/
			void store(String filename, const std::map< String, std::vector< std::pair< DoubleReal, DoubleReal > > >& isotope_informations) const throw (Exception::UnableToCreateFile);
	};
 
} // namespace OpenMS

#endif // OPENMS_FORMAT_IsotopeXMLFILE_H
