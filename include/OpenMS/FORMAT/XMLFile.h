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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XMLFILE_H
#define OPENMS_FORMAT_XMLFILE_H

// OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{	
	namespace Internal
	{
		class XMLHandler;
		
		///Base class for loading/storing XML files that have a handler derived from XMLHandler.
		class OPENMS_DLLAPI XMLFile
		{
			
			public:
				
				///Default constructor
				XMLFile();
				/// Constructor that sets the schema location
				XMLFile(const String& schema_location, const String& version);
				///Destructor
				~XMLFile();
				
				/**
					@brief Checks if a file validates against the XML schema
					
					Error messages are printed to the error stream, unless redirected with the attribute @p os .
				
			  	@exception Exception::FileNotFound is thrown if the file cannot be found
					@exception Exception::NotImplemented is thrown if there is no schema available for the file type
				*/
				bool isValid(const String& filename, std::ostream& os = std::cerr);
				
				///return the version of the schema
				const String& getVersion() const;
			
			protected:
				/**
				  @brief Parses the XML file given by @p filename using the handler given by @p handler.

				  @exception Exception::FileNotFound is thrown if the file is not found
				  @exception Exception::ParseError is thrown if an error occurred during the parsing
				*/
				void parse_(const String& filename, XMLHandler* handler);
	
				/**
				  @brief Stores the contents of the XML handler given by @p handler in the file given by @p filename.

				  @exception Exception::UnableToCreateFile is thrown if the file cannot be created
				*/
				void save_(const String& filename, XMLHandler* handler) const;
				
				/// XML schema file location
				String schema_location_;
				
				/// Version string
				String schema_version_;
		};
	}
} // namespace OpenMS

#endif // OPENMS_FOMAT_XMLFILE_H
