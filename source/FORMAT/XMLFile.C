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

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <fstream>
#include <iomanip> // setprecision etc.

namespace OpenMS
{
	namespace Internal
	{
		XMLFile::XMLFile()
		{
		}
		
		XMLFile::XMLFile(const String& schema_location, const String& version)
			: schema_location_(schema_location),
				schema_version_(version)
		{
		}
		
		XMLFile::~XMLFile()
		{
		}
	
		void XMLFile::parse_(const String& filename, XMLHandler* handler) 
		{
			//try to open file
			if (!File::exists(filename))
			{
				throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}
			
			// initialize parser
			try 
			{
				xercesc::XMLPlatformUtils::Initialize();
			}
			catch (const xercesc::XMLException& toCatch) 
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + StringManager().convert(toCatch.getMessage()) );
			}
	
			xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
			parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
			parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
	
			parser->setContentHandler(handler);
			parser->setErrorHandler(handler);
			
			// try to parse file
			xercesc::LocalFileInputSource source(StringManager().convert(filename.c_str()));
				
			try 
			{
				parser->parse(source);
				delete(parser);
			}
			catch (const xercesc::XMLException& toCatch) 
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + StringManager().convert(toCatch.getMessage()) );
			}
			catch (const xercesc::SAXException& toCatch) 
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + StringManager().convert(toCatch.getMessage()) );
			}
			catch (const XMLHandler::EndParsingSoftly& /*toCatch*/)
			{
				//nothing to do here, as this exception is used to softly abort the parsing for whatever reason.
			}
		}
	
		void XMLFile::save_(const String& filename, XMLHandler* handler) const 
		{
			std::ofstream os(filename.c_str());
			
			//set high precision for writing of floating point numbers
			os.precision(writtenDigits(DoubleReal()));
			
			if (!os)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}
	
			// write data and close stream
			handler->writeTo(os);
			os.close();
		}

		bool XMLFile::isValid(const String& filename, std::ostream& os) 
		{
			if (schema_location_.empty())
			{
				throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			String current_location = File::find(schema_location_);
			return XMLValidator().isValid(filename,current_location,os);
		}

		const String& XMLFile::getVersion() const
		{
			return schema_version_;
		}

	} // namespace Internal
} // namespace OpenMS
