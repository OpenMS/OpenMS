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

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <vector>
#include <string>
#include <netinet/in.h> //network format

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{

		XMLHandler::XMLHandler(const String& filename)
		: error_message_(""),
			file_(filename)
		{
		}
		
		XMLHandler::~XMLHandler()
		{
		}
	
		void XMLHandler::fatalError(const SAXParseException& exception)
		{
			fatalError(String(XMLString::transcode(exception.getMessage())));
		}
	
		void XMLHandler::error(const SAXParseException& exception)
		{
			error(String(XMLString::transcode(exception.getMessage())));
		}
		
		void XMLHandler::warning(const SAXParseException& exception)
		{
			warning(String(XMLString::transcode(exception.getMessage())));
		}
		
		void XMLHandler::fatalError(const String& msg)
		{
			error_message_ = "Fatal Error: " + msg;
			
			const xercesc::Locator* loc = 0;
			setDocumentLocator(loc);
			
			appendLocation_(loc, error_message_);
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, file_, error_message_);
		}
	
		void XMLHandler::error(const String& msg)
		{
			error_message_ = "Error: " + msg;
			
			const xercesc::Locator* loc = 0;
			setDocumentLocator(loc);
			
			appendLocation_(loc, error_message_);
			
			cerr << error_message_ << endl;
		}
		
		void XMLHandler::warning(const String& msg)
		{
			error_message_ = "Warning: " + msg;
			
			const xercesc::Locator* loc = 0;
			setDocumentLocator(loc);
			
			appendLocation_(loc, error_message_);
			
			cerr << error_message_ << endl;
		}
		
		void XMLHandler::characters(const XMLCh* const /*chars*/, unsigned int /*length*/)
		{
		}
		
		void XMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*localname*/, const XMLCh* const /*qname*/, const Attributes& /*attrs*/)
		{
		}
	
		void XMLHandler::endElement( const XMLCh* const /*uri*/, const XMLCh* const /*localname*/, const XMLCh* const /*qname*/)
		{
		}
	
		String XMLHandler::errorString()
		{
			return error_message_;
		}

		StringManager::StringManager()
			: xml_strings_(100),
				c_strings_(100)
		{
		}
		
		StringManager::~StringManager()
		{
			clear();
		}
		
		void StringManager::clear()
		{
			for(UInt i=0; i< xml_strings_.size(); ++i)
			{
				xercesc::XMLString::release(&xml_strings_[i]);
			}
			xml_strings_.clear();

			for(UInt i=0; i< c_strings_.size(); ++i)
			{
				xercesc::XMLString::release(&c_strings_[i]);
			}
			c_strings_.clear();
		}
		
		XMLCh* StringManager::convert(const char* str) const
		{
			XMLCh* result = xercesc::XMLString::transcode(str);
			xml_strings_.push_back(result);
			return result;
		}
		
		XMLCh* StringManager::convert(const std::string& str) const
		{
			XMLCh* result = xercesc::XMLString::transcode(str.c_str());
			xml_strings_.push_back(result) ;
			return result;
		}

		XMLCh* StringManager::convert(const String& str) const
		{
			XMLCh* result = xercesc::XMLString::transcode(str.c_str());
			xml_strings_.push_back(result) ;
			return result;
		}

		char* StringManager::convert(const XMLCh* str) const
		{
			char* result = xercesc::XMLString::transcode(str);
			c_strings_.push_back(result) ;
			return result;
		}

	} // namespace Internal
} // namespace OpenMS
