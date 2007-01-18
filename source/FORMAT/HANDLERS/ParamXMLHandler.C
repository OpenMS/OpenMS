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

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <iostream>

#include <xercesc/sax2/Attributes.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
	namespace Internal
	{

	ParamXMLHandler::ParamXMLHandler(map<std::string,DataValue>& values, const String& filename)
		: XMLHandler(filename),
			values_(values)
	{

	}
	
	ParamXMLHandler::~ParamXMLHandler()
	{
		
	}

	void ParamXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{
		if (String("ITEM") == XMLString::transcode(qname))
		{
			SignedInt type_index = attributes.getIndex(XMLString::transcode("type"));
			SignedInt name_index = attributes.getIndex(XMLString::transcode("name"));
			SignedInt value_index = attributes.getIndex(XMLString::transcode("value"));
			
			//check if attributes are present
			if (type_index==-1 || name_index==-1)
			{
				const Locator* loc = 0;
				setDocumentLocator(loc);
				String message = String("Missing attribure type or name in ITEM");
				error(SAXParseException(XMLString::transcode(message.c_str()), *loc ));
			}		
			
			//parse value/type
			String type = XMLString::transcode(attributes.getValue(type_index));
			String value;
			if (value_index!=-1)
			{
				value = XMLString::transcode(attributes.getValue(value_index));
			}
			
			//cout << "  Type: '" << type << "' Name: '" << XMLString::transcode(attributes.getValue(name_index)) <<"' Value: '" << value << "'"<< endl;
			
			if (type == "int")
			{
				values_[path_+XMLString::transcode(attributes.getValue(name_index))]=DataValue(asSignedInt_(value));
			}
			else if (type == "string")
			{
				values_[path_+XMLString::transcode(attributes.getValue(name_index))]=DataValue(value);
			}
			else if (type == "float")
			{
				values_[path_+XMLString::transcode(attributes.getValue(name_index))]=DataValue(asFloat_(value));
			}
			else if (type == "double")
			{
				values_[path_+XMLString::transcode(attributes.getValue(name_index))]=DataValue(asDouble_(value));
			}
			else
			{
				cout << "Warning: Ignoring entry '" << path_+XMLString::transcode(attributes.getValue(name_index)) << "' because of unknown type '"<< type << "'" << endl;
			}
		}
		
		if (String("NODE") == XMLString::transcode(qname))
		{
			String tmp;
			for(UnsignedInt n = 0 ; n < attributes.getLength() ; n++ )
			{
			  String attributesValue = XMLString::transcode(attributes.getValue(n));
			  String attributesName = XMLString::transcode(attributes.getQName(n));
		
			  if( attributesName == "name" )
			  {
			    nodes_.push_back(attributesValue);
					tmp += attributesValue + ":";
				}
			}
			path_ += tmp;
		}
	}

	void ParamXMLHandler::endElement( const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		if (String("NODE") == XMLString::transcode(qname))
		{
			nodes_.pop_back();
			//renew path
			path_ = "";
			for (vector<string>::iterator it = nodes_.begin(); it != nodes_.end();++it)
			{
				path_ += *it+":";
			}
			
		}		
	}


	} // namespace Internal
} // namespace OpenMS
