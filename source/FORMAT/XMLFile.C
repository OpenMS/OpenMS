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

#include <OpenMS/FORMAT/CompressedInputSource.h>

#include <xercesc/sax2/SAX2XMLReader.hpp>
#include <xercesc/framework/LocalFileInputSource.hpp>
#include <xercesc/sax2/XMLReaderFactory.hpp>
#include <fstream>
#include <iomanip> // setprecision etc.

using namespace std;

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
			
			//is it bzip2 or gzip compressed?
			std::ifstream file(filename.c_str());
			char bz[2];
			file.read(bz,2);
			xercesc::InputSource *source;
			char g1 = 0x1f;
			char g2 = 0;
			g2 |= 1 << 7;
			g2 |= 1 <<3;
			g2  |=1  <<1;
			g2 |=1 <<0;
			//g2 = static_cast<char>(0x8b); // can make troubles if it is casted to 0x7F which is the biggest number signed char can save
			if((bz[0] == 'B' && bz[1] =='Z' ) || 	(bz[0] == g1 && bz[1] == g2))
			{
				source = new CompressedInputSource(StringManager().convert(filename.c_str()), bz);
			}
			else
			{
				source = new xercesc::LocalFileInputSource(StringManager().convert(filename.c_str()));
			}	
			// try to parse file	
			try 
			{
				parser->parse(*source);
				delete(parser);
				delete source;
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

		void writeXMLEscape(const String& to_escape, ostream& os)
		{
			XMLCh* xmlch = xercesc::XMLString::transcode(to_escape.c_str());

			std::string out;
    	OpenMSXMLFormatTarget ft(out);
    	xercesc::XMLFormatter f("UTF-8", /* XMLUni::fgVersion1_1 */ "1.1", &ft);
  	  f << xercesc::XMLFormatter::StdEscapes << xmlch;
			os << out;
			xercesc::XMLString::release(&xmlch);
			return;

		}

		String writeXMLEscape(const String& to_escape)
		{
			stringstream ss;
			writeXMLEscape(to_escape, ss);
			return String(ss.str());
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
