// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
			//tags
			String tags_string;
			optionalAttributeAsString_(tags_string,attributes,"tags");
			StringList tags = StringList::create(tags_string);
			//advanced (for downward compatibility with old Param files)
			String advanced_string;
			optionalAttributeAsString_(advanced_string,attributes,"advanced");
			if (advanced_string=="true") 
			{
				tags.push_back("advanced");
			}
			//type
			if (type == "int")
			{
				param_.setValue(name, asInt_(value), description, tags);
			}
			else if (type == "string")
			{
				param_.setValue(name, value, description, tags);
			}
			else if (type == "float" || type == "double")
			{
				param_.setValue(name, asDouble_(value), description, tags);
			}
			else
			{
				warning(LOAD, String("Ignoring entry '") + name + "' because of unknown type '" + type + "'");
			}
			
			//restrictions
			Int restrictions_index = attributes.getIndex(sm_.convert("restrictions"));
			if(restrictions_index!=-1)
			{
				String value = sm_.convert(attributes.getValue(restrictions_index));
				std::vector<String> parts;
				if (type == "int")
				{
					value.split(':', parts);
					if (parts.size()!=2) value.split('-', parts); //for downward compatibility
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
					value.split(':', parts);
					if (parts.size()!=2) value.split('-', parts); //for downward compatibility
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
				descriptions_[path_.chop(1)] = description;
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
			//tags
			String tags_string;
			optionalAttributeAsString_(tags_string,attributes,"tags");
			list_.tags = StringList::create(tags_string);
			//advanced (for downward compatibility with old Param files)
			String advanced_string;
			optionalAttributeAsString_(advanced_string,attributes,"advanced");
			if (advanced_string=="true") 
			{
				list_.tags.push_back("advanced");
			}
			list_.restrictions_index = attributes.getIndex(sm_.convert("restrictions"));
			if(list_.restrictions_index!=-1)
			{
				list_.restrictions = sm_.convert(attributes.getValue(list_.restrictions_index));
			}
		}
		else if (element == "LISTITEM")
		{
			if(list_.type == "string")
			{
				list_.stringlist.push_back(attributeAsString_(attributes,"value"));
			}
			else if(list_.type == "int")
			{
				list_.intlist.push_back(asInt_(attributeAsString_(attributes,"value")));
			}
			else if(list_.type =="float" || list_.type == "double")
			{
				list_.doublelist.push_back(asDouble_(attributeAsString_(attributes,"value")));
			}
			
		}
		else if (element == "PARAMETERS")
		{
			//check file version against schema version
			String file_version="";
			optionalAttributeAsString_(file_version,attributes,"version");
			if (file_version=="") file_version="1.0"; //default version is 1.0
			if (file_version.toDouble()>version_.toDouble())
			{
				warning(LOAD, "The XML file (" + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
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
			std::vector<String> parts;

			if (list_.type=="string")
			{
				param_.setValue(list_.name, list_.stringlist, list_.description, list_.tags);
				if(list_.restrictions_index!=-1)
				{
					list_.restrictions.split(',', parts);
					param_.setValidStrings(list_.name,parts);
				}
			}
			else if(list_.type=="int")
			{
				param_.setValue(list_.name,list_.intlist,list_.description,list_.tags);
				if(list_.restrictions_index!=-1)
				{
					list_.restrictions.split(':', parts);
					if (parts.size()!=2) list_.restrictions.split('-', parts); //for downward compatibility
					if (parts[0]!="")
					{
						param_.setMinInt(list_.name,parts[0].toInt());
					}
					if (parts[1]!="")
					{
						param_.setMaxInt(list_.name,parts[1].toInt());
					}
				}
			}
			else if(list_.type == "float" || list_.type == "double")
			{
				param_.setValue(list_.name,list_.doublelist,list_.description,list_.tags);
				if(list_.restrictions_index!=-1)
				{
					list_.restrictions.split(':', parts);
					if (parts.size()!=2) list_.restrictions.split('-', parts); //for downward compatibility
					if (parts[0]!="")
					{
						param_.setMinFloat(list_.name,parts[0].toDouble());
					}
					if (parts[1]!="")
					{
						param_.setMaxFloat(list_.name,parts[1].toDouble());
					}
				}
			}
			else
			{
				warning(LOAD, String("Ignoring list entry '") + list_.name + "' because of unknown type '" + list_.type + "'");				
			}
			list_.stringlist.clear();
			list_.intlist.clear();
			list_.doublelist.clear();
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
