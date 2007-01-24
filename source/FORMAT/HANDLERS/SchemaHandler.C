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

	SchemaHandler::SchemaHandler(Size tag_num, Size map_num, const String& filename)
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
	
	UnsignedInt SchemaHandler::leaveTag(const XMLCh* const qname)
	{
		int tag = str2enum_(tag_map_, xercesc::XMLString::transcode(qname),"closing tag"); // index of current tag
		is_parser_in_tag_[tag] = false;
		atts_ = 0;
		
		skip_tag_.pop();
		
		return tag;
	}
	
	UnsignedInt SchemaHandler::enterTag(const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		String tmp_str;
		int tag = str2enum_(tag_map_, xercesc::XMLString::transcode(qname),"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;
		atts_ = &attributes;
		
		if (!skip_tag_.empty())
			skip_tag_.push(skip_tag_.top());
		else
			skip_tag_.push(false);
		
		return tag;
	}
	
	UnsignedInt SchemaHandler::str2enum_(UnsignedInt index, const String& value, const char* message)
	{
		String2EnumMap::const_iterator it =  str2enum_array_[index].find(value);
		if (it == str2enum_array_[index].end()) // no enum-value for string defined
		{  
			const xercesc::Locator* loc = 0;
			setDocumentLocator(loc);
			std::cout << "Warning: Unhandled " << message << " \"" << value << "\" parsed in " << file_ << std::endl;
		}	
		else
		{
			return it->second;
		}
		return 0;
	}

	const String& SchemaHandler::enum2str_(UnsignedInt index, UnsignedInt value)
	{
		return enum2str_array_[index][value];
	}

	void SchemaHandler::fillMaps_(const String* schema)
	{
		for (Size i= 0; i<str2enum_array_.size(); i++)
		{
			//i=0 contains scheme name -> i+1
			schema[i+1].split(';',enum2str_array_[i]);
			fillMap_(str2enum_array_[i], enum2str_array_[i]);
		}
	}

	void SchemaHandler::fillMap_(String2EnumMap& str2enum, const Enum2StringMap& enum2str)
	{
		for (Size i=0; i<enum2str.size(); i++)
		{
			str2enum[ enum2str[i] ] = i;
		}
	}

	void SchemaHandler::setAddInfo_(	MetaInfoInterface& info, const String& name, const String& value, const String& description)
	{
		info.setMetaValue(info.metaRegistry().registerName(name, description),value);
	}

	void SchemaHandler::writeCVS_(std::ostream& os, float value, const String& acc, const String& name, int indent)
	{
		if (value)
			os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:"
					<< acc << "\" name=\""
					<< name << "\" value=\""
					<< value << "\"/>\n";
	}

	void SchemaHandler::writeCVS_(std::ostream& os, const String& value, const String& acc, const String& name, int indent)
	{
		if (value!="")
			os << String(indent,'\t') << "<cvParam cvLabel=\"psi\" accession=\"PSI:"
					<< acc << "\" name=\""
					<< name << "\" value=\""
					<< value << "\"/>\n";
	}

	void SchemaHandler::writeCVS_(std::ostream& os, int value, int map, const String& acc, const String& name, int indent)
	{
		writeCVS_(os, enum2str_(map,value), acc, name, indent);
	}

	void SchemaHandler::writeUserParam_(std::ostream& os, const MetaInfoInterface& meta, int indent)
	{
		std::vector<std::string> keys;  // Vector to hold keys to meta info
		meta.getKeys(keys);

		for (std::vector<std::string>::const_iterator it = keys.begin(); it!=keys.end(); ++it)
			if ( (*it)[0] != '#')  // internally used meta info start with '#'
				os << String(indent,'\t') << "<userParam name=\""
						<< *it << "\" value=\""
						<< meta.getMetaValue(*it) << "\"/>\n";
	}

	void SchemaHandler::checkAttribute_(UnsignedInt attribute, const String& required, const String& required_alt)
	{
		//TODO improve performace
		const XMLCh* tmp = xercesc::XMLString::transcode(enum2str_(att_map_, attribute).c_str());
		if (tmp==0) return;
		if (atts_->getIndex(tmp)==-1) return;
		String value = xercesc::XMLString::transcode(atts_->getValue(tmp));
		if (value!=required && value!=required_alt)
		{
			error("Invalid value \"" + value + "\" for attribute \"" + enum2str_(att_map_, attribute) + "\"");
		}
	}

	String SchemaHandler::getAttributeAsString_(UnsignedInt attribute)
	{
		//TODO improve performace
		const XMLCh* tmp = xercesc::XMLString::transcode(enum2str_(att_map_, attribute).c_str());
		if (atts_->getIndex(tmp)==-1) 
		{
			return "";
		}
		return xercesc::XMLString::transcode(atts_->getValue(tmp));
	}
	
	void SchemaHandler::setMaps_(UnsignedInt tagmap, UnsignedInt attmap)
	{
		tag_map_ = tagmap;
		att_map_ = attmap;
	}
		
	} // namespace Internal
} // namespace OpenMS

