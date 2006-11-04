// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

	UnsignedInt SchemaHandler::leaveTag(UnsignedInt tagmap, const XMLCh* const qname)
	{
		int tag = str2enum_(tagmap, xercesc::XMLString::transcode(qname),"closing tag"); // index of current tag
		is_parser_in_tag_[tag] = false;
		atts_ = 0;
		return tag;
	}
	
	UnsignedInt SchemaHandler::enterTag(UnsignedInt tagmap, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		String tmp_str;
		int tag = str2enum_(tagmap, xercesc::XMLString::transcode(qname),"opening tag");	// index of current tag
		is_parser_in_tag_[tag] = true;
		atts_ = &attributes;
		return tag;
	}
	
	} // namespace Internal
} // namespace OpenMS






