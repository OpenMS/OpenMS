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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/XMLValidator.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <fstream>

namespace OpenMS
{
	namespace Internal
	{
		XMLFile::XMLFile()
		{
		}
		
		XMLFile::XMLFile(const String& schema_location)
			: schema_location_(schema_location)
		{
		}
		
		XMLFile::~XMLFile()
		{
		}
	
		void XMLFile::parse_(const String& filename, XMLHandler* handler) throw (Exception::FileNotFound, Exception::ParseError)
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
	
		void XMLFile::save_(const String& filename, XMLHandler* handler) const throw (Exception::UnableToCreateFile)
		{
			std::ofstream os(filename.c_str());
			if (!os)
			{
				throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
			}
	
			// write data and close stream
			handler->writeTo(os);
			os.close();
		}

		bool XMLFile::isValid(const String& filename) throw (Exception::NotImplemented)
		{
			if (schema_location_.empty())
			{
				throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			return XMLValidator().isValid(filename,schema_location_);
		}

	} // namespace Internal
} // namespace OpenMS
