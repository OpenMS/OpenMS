// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
			fatalError(sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
	
		void XMLHandler::error(const SAXParseException& exception)
		{
			error(sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
		
		void XMLHandler::warning(const SAXParseException& exception)
		{
			warning(sm_.convert(exception.getMessage()),exception.getLineNumber(),exception.getColumnNumber());
		}
		
		void XMLHandler::fatalError(const String& msg, UInt line, UInt column) const
		{
			error_message_ = String("Fatal error while parsing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
			cerr << error_message_ << endl;
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, file_, error_message_);
		}
	
		void XMLHandler::error(const String& msg, UInt line, UInt column) const
		{
			error_message_ = String("Non-fatal error while parsing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
			cerr << error_message_ << endl;
		}
		
		void XMLHandler::warning(const String& msg, UInt line, UInt column) const
		{
			error_message_ = String("Warning while parsing '") + file_ + "': " + msg;
			if (line!=0 || column!=0) error_message_ += String("( in line ") + line + " column " + column + ")";
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

		void XMLHandler::writeCVS_(std::ostream& os, int value, int map, const String& acc, const String& name, int indent)
		{
			//abort when receiving a wrong map index
			if (map>=cv_terms_.size())
			{
				warning(String("Cannot find map '") + map + "' needed to write CV term '" + name + "' with accession '" + acc + "'.");
				return;
			}
			//abort when receiving a wrong term index
			if (value>=cv_terms_[map].size())
			{
				warning(String("Cannot find value '") + value + "' needed to write CV term '" + name + "' with accession '" + acc + "'.");
				return;
			}
			XMLHandler::writeCVS_(os, cv_terms_[map][value], acc, name, indent);
		}

		void XMLHandler::writeUserParam_(const String& tag_name, std::ostream& os, const MetaInfoInterface& meta, UInt indent) const
		{
			std::vector<String> keys;
			meta.getKeys(keys);
			
			for (UInt i = 0; i!=keys.size();++i)
			{
				os << String(indent,'\t') << "<" << tag_name << " type=\"";
				
				DataValue d = meta.getMetaValue(keys[i]);
				//determine type
				if (d.valueType()==DataValue::STRING_VALUE)
				{
					os << "string\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
				if (d.valueType()==DataValue::INT_VALUE)
				{
					os << "int\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
				if (d.valueType()==DataValue::DOUBLE_VALUE)
				{
					os << "float\" name=\"" << keys[i] << "\" value=\"" << (String)(d) << "\"/>" << endl;
				}
			}
		}
		
		
		//*******************************************************************************************************************

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
				XMLString::release(&xml_strings_[i]);
			}
			xml_strings_.clear();

			for(UInt i=0; i< c_strings_.size(); ++i)
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
