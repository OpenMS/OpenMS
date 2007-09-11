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

	ParamXMLHandler::ParamXMLHandler(Param& param, const String& filename)
		: XMLHandler(filename),
			param_(param)
	{
	}
	
	ParamXMLHandler::~ParamXMLHandler()
	{
	}

	void ParamXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{
		if (String("ITEM") == sm_.convert(qname))
		{
			Int type_index = attributes.getIndex(sm_.convert("type"));
			Int name_index = attributes.getIndex(sm_.convert("name"));
			Int value_index = attributes.getIndex(sm_.convert("value"));
				
			//check if attributes are present
			if (type_index==-1 || name_index==-1)
			{
				const Locator* loc = 0;
				setDocumentLocator(loc);
				String message = String("Missing attribute type or name in ITEM");
				error(SAXParseException(sm_.convert(message.c_str()), *loc ));
			}		
			
			//parse value/type
			String type = sm_.convert(attributes.getValue(type_index));
			String value;
			if (value_index!=-1)
			{
				value = sm_.convert(attributes.getValue(value_index));
			}
			
			//parse description, if present
			String description = "";
			Int description_index = attributes.getIndex(sm_.convert("description"));
			if(description_index!=-1)
			{
				description = sm_.convert(attributes.getValue(description_index));
				description.substitute("#br#","\n");
			}

			//advanced parameters
			bool advanced = false;
			Int advanced_index = attributes.getIndex(sm_.convert("advanced"));
			if(advanced_index!=-1)
			{
				String value = sm_.convert(attributes.getValue(advanced_index));
				if (value=="true") 
				{
					advanced = true;
				}
			}

			if (type == "int")
			{
				param_.setValue(path_+sm_.convert(attributes.getValue(name_index)), asInt_(value), description, advanced);
			}
			else if (type == "string")
			{
				param_.setValue(path_+sm_.convert(attributes.getValue(name_index)), value, description, advanced);
			}
			else if (type == "float")
			{
				param_.setValue(path_+sm_.convert(attributes.getValue(name_index)), asFloat_(value), description, advanced);
			}
			else if (type == "double")
			{
				param_.setValue(path_+sm_.convert(attributes.getValue(name_index)), asDouble_(value), description, advanced);
			}
			else
			{
				cout << "Warning: Ignoring entry '" << path_+sm_.convert(attributes.getValue(name_index)) << "' because of unknown type '"<< type << "'" << endl;
			}
		}
		else if (String("NODE") == sm_.convert(qname))
		{
			//parse name
			Int name_index = attributes.getIndex(sm_.convert("name"));
			String name = sm_.convert(attributes.getValue(name_index));
	    nodes_.push_back(name);
	    path_ += name + ":";

			//parse description, if present
			Int description_index = attributes.getIndex(sm_.convert("description"));
			if(description_index!=-1)
			{
				String description = sm_.convert(attributes.getValue(description_index));
				description.substitute("#br#","\n");
				descriptions_[path_.substr(0,-1)] = description;
			}
		}
	}

	void ParamXMLHandler::endElement( const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		if (String("NODE") == sm_.convert(qname))
		{
			nodes_.pop_back();
			//renew path
			path_ = "";
			for (vector<String>::iterator it = nodes_.begin(); it != nodes_.end();++it)
			{
				path_ += *it+":";
			}
		}
		else if (String("PARAMETERS") == sm_.convert(qname))
		{
			//set all descriptions (now the nodes exist...)
			for(map<String,String>::const_iterator it = descriptions_.begin(); it !=descriptions_.end(); ++it)
			{
				param_.setSectionDescription(it->first,it->second);
			}
			descriptions_.clear();
		}
	}


	} // namespace Internal
} // namespace OpenMS
