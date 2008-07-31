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

#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

#include <iostream>

#include <xercesc/sax2/Attributes.hpp>

using namespace xercesc;
using namespace std;

namespace OpenMS
{
	namespace Internal
	{

	ParamXMLHandler::ParamXMLHandler(Param& param, const String& filename, const String& version)
		: XMLHandler(filename,version),
			param_(param)
	{
	}
	
	ParamXMLHandler::~ParamXMLHandler()
	{
	}

	void ParamXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const Attributes& attributes)
	{
		String element = sm_.convert(qname);
		if (element == "ITEM")
		{
			//parse value/type
			String type = attributeAsString_(attributes,"type");
			String name = path_ + attributeAsString_(attributes,"name");
			String value = attributeAsString_(attributes,"value");
			//parse description, if present
			String description;
			optionalAttributeAsString_(description,attributes,"description");
			description.substitute("#br#","\n");
			//advanced parameters
			bool advanced = false;
			String advanced_string;
			optionalAttributeAsString_(advanced_string,attributes,"advanced");
			if (advanced_string=="true") 
			{
				advanced = true;
			}			
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
				warning(String("Ignoring entry '") + name + "' because of unknown type '" + type + "'");
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
		else if (element == "NODE")
		{
			//parse name
			String name = attributeAsString_(attributes,"name");
	    open_tags_.push_back(name);
	    path_ += name + ":";
			//parse description
			String description;
			optionalAttributeAsString_(description,attributes,"description");
			if (description!="")
			{
				description.substitute("#br#","\n");
				descriptions_[path_.substr(0,-1)] = description;
			}
		}
		else if (element == "ITEMLIST")
		{
			//parse name/type
		  list_.type = attributeAsString_(attributes,"type");
			list_.name = path_ + attributeAsString_(attributes,"name");
			//parse description, if present
			list_.description = "";
			optionalAttributeAsString_(list_.description,attributes,"description");
			list_.description.substitute("#br#","\n");
			//advanced parameters
			list_.advanced = false;
			String advanced_string;
			optionalAttributeAsString_(advanced_string,attributes,"advanced");
			if (advanced_string=="true") 
			{
				list_.advanced = true;
			}			
		}
		else if (element == "LISTITEM")
		{
			list_.list.push_back(attributeAsString_(attributes,"value"));
		}
		else if (element == "PARAMETERS")
		{
			//check file version against schema version
			String file_version="1.0";
			optionalAttributeAsString_(file_version,attributes,"version");
			if (file_version.toDouble()>version_.toDouble())
			{
				warning("The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
			}
		}
	}

	void ParamXMLHandler::endElement( const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String element = sm_.convert(qname);
		if (element == "NODE")
		{
			open_tags_.pop_back();
			//renew path
			path_ = "";
			for (vector<String>::iterator it = open_tags_.begin(); it != open_tags_.end();++it)
			{
				path_ += *it+":";
			}
		}
		else if (element == "ITEMLIST")
		{
			if (list_.type=="string")
			{
				param_.setValue(list_.name, list_.list, list_.description, list_.advanced);
			}
			else
			{
				warning(String("Ignoring list entry '") + list_.name + "' because of unknown type '" + list_.type + "'");				
			}
			list_.list.clear();
		}
		else if (element == "PARAMETERS")
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
