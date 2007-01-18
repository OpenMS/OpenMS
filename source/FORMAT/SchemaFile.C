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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SchemaFile.h>
#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>

#include <fstream>

namespace OpenMS
{
	namespace Internal
	{
		SchemaFile::SchemaFile() {}
		SchemaFile::~SchemaFile() {}
	
		void SchemaFile::parse_(const String& filename, SchemaHandler* handler) throw (Exception::FileNotFound, Exception::ParseError)
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
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Error during initialization: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
			}
	
			xercesc::SAX2XMLReader* parser = xercesc::XMLReaderFactory::createXMLReader();
			parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpaces,false);
			parser->setFeature(xercesc::XMLUni::fgSAX2CoreNameSpacePrefixes,false);
	
			parser->setContentHandler(handler);
			parser->setErrorHandler(handler);
			
			// try to parse file
			xercesc::LocalFileInputSource source( xercesc::XMLString::transcode(filename.c_str()) );
			try 
			{
				parser->parse(source);
				delete(parser);
			}
			catch (const xercesc::XMLException& toCatch) 
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("XMLException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
			}
			catch (const xercesc::SAXException& toCatch) 
			{
				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("SAXException: ") + xercesc::XMLString::transcode(toCatch.getMessage()) );
			}
		}
	
		void SchemaFile::save_(const String& filename, SchemaHandler* handler) const throw (Exception::UnableToCreateFile)
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
	} // namespace Internal
} // namespace OpenMS
