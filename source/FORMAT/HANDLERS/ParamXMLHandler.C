// -*- Mode: C++; tab-width: 2; -*-
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
				error("Missing attribute type or name in ITEM", 0, 0);
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
			
			String name = path_+sm_.convert(attributes.getValue(name_index));
			
			//type
			if (type == "int")
			{
				param_.setValue(name, asInt_(value), description, advanced);
			}
			else if (type == "string")
			{
				param_.setValue(name, value, description, advanced);
			}
			else if (type == "float" || type == "double")
			{
				param_.setValue(name, asDouble_(value), description, advanced);
			}
			else
			{
				cout << "Warning: Ignoring entry '" << name << "' because of unknown type '"<< type << "'" << endl;
			}
			
			//restrictions
			Int restrictions_index = attributes.getIndex(sm_.convert("restrictions"));
			if(restrictions_index!=-1)
			{
				String value = sm_.convert(attributes.getValue(restrictions_index));
				std::vector<String> parts;
				if (type == "int")
				{
					value.split('-', parts);
					if (parts[0]!="")
					{
						param_.setMinInt(name,parts[0].toInt());
					}
					if (parts[1]!="")
					{
						param_.setMaxInt(name,parts[1].toInt());
					}
				}
				else if (type == "string")
				{
					value.split(',', parts);
					param_.setValidStrings(name,parts);
				}
				else if (type == "float" || type == "double")
				{
					value.split('-', parts);
					if (parts[0]!="")
					{
						param_.setMinFloat(name,parts[0].toDouble());
					}
					if (parts[1]!="")
					{
						param_.setMaxFloat(name,parts[1].toDouble());
					}
				}
			}
		}
		else if (String("NODE") == sm_.convert(qname))
		{
			//parse name
			Int name_index = attributes.getIndex(sm_.convert("name"));
			String name = sm_.convert(attributes.getValue(name_index));
	    open_tags_.push_back(name);
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
			open_tags_.pop_back();
			//renew path
			path_ = "";
			for (vector<String>::iterator it = open_tags_.begin(); it != open_tags_.end();++it)
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
