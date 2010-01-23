// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{

		XMLHandler::XMLHandler(const String& filename, const String& version)
		: error_message_(""),
			file_(filename),
			version_(version)
		{
		}
		
		XMLHandler::~XMLHandler()
		{
		}
	
		void XMLHandler::fatalError(const SAXParseException& exception)
		{
			fatalError(LOAD, sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
	
		void XMLHandler::error(const SAXParseException& exception)
		{
			error(LOAD, sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
		
		void XMLHandler::warning(const SAXParseException& exception)
		{
			warning(LOAD, sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
		
		void XMLHandler::fatalError(ActionMode mode, const String& msg, UInt line, UInt column) const
		{
			if (mode==LOAD)  error_message_ =  String("While loading '") + file_ + "': " + msg;
			if (mode==STORE) error_message_ =  String("While storing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
			LOG_FATAL_ERROR << error_message_ << "\n";
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, file_, error_message_);
		}
	
		void XMLHandler::error(ActionMode mode, const String& msg, UInt line, UInt column) const
		{
			if (mode==LOAD)  error_message_ =  String("Non-fatal error while loading '") + file_ + "': " + msg;
			if (mode==STORE) error_message_ =  String("Non-fatal error while storing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
			LOG_ERROR << error_message_ << "\n";
		}
		
		void XMLHandler::warning(ActionMode mode, const String& msg, UInt line, UInt column) const
		{
			if (mode==LOAD)  error_message_ =  String("While loading '") + file_ + "': " + msg;
			if (mode==STORE) error_message_ =  String("While storing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
			LOG_WARN << error_message_ << "\n";
		}
		
		void XMLHandler::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
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

		void XMLHandler::writeUserParam_(const String& tag_name, std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
		{
			std::vector<String> keys;
			meta.getKeys(keys);
			
			for (Size i = 0; i!=keys.size();++i)
			{
				os << String(indent,'\t') << "<" << writeXMLEscape(tag_name) << " type=\"";
				
				DataValue d = meta.getMetaValue(keys[i]);
				//determine type
				if (d.valueType()==DataValue::INT_VALUE)
				{
					os << "int";
				}
				else if (d.valueType()==DataValue::DOUBLE_VALUE)
				{
					os << "float";
				}
				else //string or lists are converted to string
				{
					os << "string"; 
				}
				os << "\" name=\"" << keys[i] << "\" value=\"" << writeXMLEscape((String)(d)) << "\"/>" << "\n";
			}
		}
		
		
		//*******************************************************************************************************************

		StringManager::StringManager()
			: xml_strings_(0),
				c_strings_(0)
		{
		}
		
		StringManager::~StringManager()
		{
			clear();
		}
		
		void StringManager::clear()
		{
			for (Size i=0; i< xml_strings_.size(); ++i)
			{
				XMLString::release(&xml_strings_[i]);
			}
			xml_strings_.clear();

			for (Size i=0; i< c_strings_.size(); ++i)
			{
				XMLString::release(&c_strings_[i]);
			}
			c_strings_.clear();
		}
		
		XMLCh* StringManager::convert(const char* str) const
		{
			XMLCh* result = XMLString::transcode(str);
			xml_strings_.push_back(result);
			return result;
		}
		
		XMLCh* StringManager::convert(const std::string& str) const
		{
			XMLCh* result = XMLString::transcode(str.c_str());
			xml_strings_.push_back(result) ;
			return result;
		}

		XMLCh* StringManager::convert(const String& str) const
		{
			XMLCh* result = XMLString::transcode(str.c_str());
			xml_strings_.push_back(result) ;
			return result;
		}

		char* StringManager::convert(const XMLCh* str) const
		{
			char* result = XMLString::transcode(str);
			c_strings_.push_back(result) ;
			return result;
		}

	} // namespace Internal
} // namespace OpenMS
