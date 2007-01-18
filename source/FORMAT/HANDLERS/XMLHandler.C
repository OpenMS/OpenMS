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
			// Thow error here as well?????
		}
		
		void XMLHandler::warning(const String& msg)
		{
			error_message_ = "Warning: " + msg;
			
			const xercesc::Locator* loc = 0;
			setDocumentLocator(loc);
			
			appendLocation_(loc, error_message_);
		}
		
		void XMLHandler::characters(const XMLCh* const /*chars*/, const unsigned int /*length*/)
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

	} // namespace Internal
} // namespace OpenMS
