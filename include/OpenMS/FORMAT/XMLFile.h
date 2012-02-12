// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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

#include <xercesc/framework/XMLFormatter.hpp>

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
				virtual ~XMLFile();
				
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

  	// implementation of an XMLFormatTarget
  	class OPENMS_DLLAPI OpenMSXMLFormatTarget : public xercesc::XMLFormatTarget
  	{

			public:

    	OpenMSXMLFormatTarget(std::string &str)
      	: XMLFormatTarget(),
    			str_(str)
			{
    	}

    	virtual void writeChars(const XMLByte* const toWrite, const XMLSize_t count, xercesc::XMLFormatter* const /*formatter*/)
    	{
      	str_.append((const char*const)toWrite,count);
    	}

    	std::string &str_;
  	};
	
 	 /**
	    @brief Escapes a string to be storable into an XML File
    
			Some characters must be escaped which are allowed in user params. E.g. > and & are not in XML and 
  	  need to be escaped. Parsing those escaped strings is automatically done by xerces
 	 */
	 	void OPENMS_DLLAPI writeXMLEscape(const String& to_escape, std::ostream& os);

		/** 
	  	@brief Escapes a string and returns the escaped string

			Some characters must be escaped which are allowed in user params. E.g. > and & are not in XML and 
    	need to be escaped. Parsing those escaped strings is automatically done by xerces
		*/
		String OPENMS_DLLAPI writeXMLEscape(const String& to_escape);

	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FOMAT_XMLFILE_H
