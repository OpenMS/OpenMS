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

#include <OpenMS/FORMAT/HANDLERS/SchemaHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLSchemes.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
	namespace Internal
	{
	SchemaHandler::SchemaHandler(const String& filename)
  : XMLHandler(filename),
  	is_parser_in_tag_(),
		str2enum_array_(), enum2str_array_(),
		schema_(0)
	{

	}

	SchemaHandler::SchemaHandler(UInt tag_num, UInt map_num, const String& filename)
  : XMLHandler(filename),
  	is_parser_in_tag_(tag_num,false),
		str2enum_array_(map_num), 
		enum2str_array_(map_num),
		schema_(0)
	{

	}

	SchemaHandler::~SchemaHandler()
	{
		
	}

	void SchemaHandler::skipTag_()
	{
		skip_tag_.pop();
		skip_tag_.push(true);
	}
	
	UInt SchemaHandler::leaveTag(const XMLCh* const qname)
	{
		int tag = str2enum_(tag_map_, sm_.convert(qname),"closing tag"); // index of current tag
		is_parser_in_tag_[tag] = false;
		atts_ = 0;
		skip_tag_.pop();
		return tag;
	}
	
	UInt SchemaHandler::enterTag(const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		int tag = str2enum_(tag_map_, sm_.convert(qname),"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;
		atts_ = &attributes;
		if (!skip_tag_.empty())
		{
			skip_tag_.push(skip_tag_.top());
		}
		else
		{
			skip_tag_.push(false);
		}
		return tag;
	}
	
	UInt SchemaHandler::str2enum_(UInt index, const String& value, const char* message)
	{
		//std::cout << "looking up key \"" << value << "\" in map nr. " << index << "..." << std::endl;

		String2EnumMap::const_iterator it =  str2enum_array_[index].find(value);
		if (it == str2enum_array_[index].end()) // no enum-value for string defined
		{
			std::cout << "Warning: Unhandled object \"" << message << "\"=\"" << value << "\" parsed in " << file_ << std::endl;
		}
		else
		{
			return it->second;
		}
		
		return 0;
	}

	const String& SchemaHandler::enum2str_(UInt index, UInt value)
	{
		return enum2str_array_[index][value];
	}

	void SchemaHandler::fillMaps_(const String* schema)
	{
		for (UInt i= 0; i<str2enum_array_.size(); i++)
		{
			//i=0 contains scheme name -> i+1
			schema[i+1].split(';',enum2str_array_[i]);
			fillMap_(str2enum_array_[i], enum2str_array_[i]);
		}
	}

	void SchemaHandler::fillMap_(String2EnumMap& str2enum, const Enum2StringMap& enum2str)
	{
		for (UInt i=0; i<enum2str.size(); i++)
		{
			str2enum[ enum2str[i] ] = i;
		}
	}

	void SchemaHandler::writeCVS2_(std::ostream& os, int value, int map, const String& acc, const String& name, int indent)
	{
		XMLHandler::writeCVS_(os, enum2str_(map,value), acc, name, indent);
	}

	void SchemaHandler::checkAttribute_(UInt attribute, const String& required, const String& required_alt)
	{
		XMLCh* tmp = sm_.convert(enum2str_(att_map_, attribute));
		if (tmp==0) // no value
		{
			return;
		}
		if (atts_->getIndex(tmp)==-1) //not present
		{
			return;
		}
		//convert to String
		String value = sm_.convert(atts_->getValue(tmp));
		if (value!=required && value!=required_alt)
		{
			error("Invalid value \"" + value + "\" for attribute \"" + enum2str_(att_map_, attribute) + "\"");
		}
	}

	String SchemaHandler::getAttributeAsString_(UInt attribute, bool is_required, const XMLCh* tag)
	{
		XMLCh* tmp = sm_.convert(enum2str_(att_map_, attribute));
		if (atts_->getIndex(tmp)==-1) 
		{
			if (is_required)
			{
				error(String("Required attribute '") + enum2str_(att_map_, attribute) + "' missing in tag '" + sm_.convert(tag) + "'!");
			}
			return "";
		}
		return sm_.convert(atts_->getValue(tmp));
	}
	
	void SchemaHandler::setMaps_(UInt tagmap, UInt attmap)
	{
		tag_map_ = tagmap;
		att_map_ = attmap;
	}
		
	} // namespace Internal
} // namespace OpenMS

