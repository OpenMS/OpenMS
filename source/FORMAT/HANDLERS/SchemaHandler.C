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
	SchemaHandler::SchemaHandler()
  : is_parser_in_tag_(),
		str2enum_array_(), enum2str_array_(),
		schema_(0)
	{
		file_ = __FILE__;
	}

	SchemaHandler::SchemaHandler(Size tag_num, Size map_num)
  : is_parser_in_tag_(tag_num,false),
		str2enum_array_(map_num), enum2str_array_(map_num),
		schema_(0)
	{
		file_ = __FILE__;
	}

	SchemaHandler::~SchemaHandler()
	{}

	} // namespace Internal
} // namespace OpenMS






