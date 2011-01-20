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
// $Maintainer: $
// $Authors: Chris Bielow, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ToolDescriptionFile.h>
//#include <OpenMS/FORMAT/VALIDATORS/ToolDescriptorValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/HANDLERS/ToolDescriptionHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

	ToolDescriptionFile::ToolDescriptionFile()
		: XMLFile("/SCHEMAS/ToolDescriptor_1_0.xsd","1.0.0")
	{
	}

	ToolDescriptionFile::~ToolDescriptionFile()
	{
	}

  void ToolDescriptionFile::load(const String& filename, std::vector<Internal::ToolDescription>& tds)
  {
  	Internal::ToolDescriptionHandler handler(filename, schema_version_);
    parse_(filename, &handler);
    tds = handler.getToolDescriptions();
  }

  void ToolDescriptionFile::store(const String& filename, const std::vector<Internal::ToolDescription>& tds) const
  {
  	Internal::ToolDescriptionHandler handler(filename, schema_version_);
    handler.setToolDescriptions(tds);
    save_(filename, &handler);
  }


}// namespace OpenMS

