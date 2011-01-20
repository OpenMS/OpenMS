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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/TraMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/HANDLERS/TraMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>

namespace OpenMS
{

	TraMLFile::TraMLFile()
		: XMLFile("/SCHEMAS/TraML0.9.3.xsd","0.9.3")
	{
	}

	TraMLFile::~TraMLFile()
	{
	}

  void TraMLFile::load(const String& filename, TargetedExperiment& exp)
  {
  	Internal::TraMLHandler handler(exp, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void TraMLFile::store(const String& filename, const TargetedExperiment& exp) const
  {
  	Internal::TraMLHandler handler(exp, filename, schema_version_, *this);
    save_(filename, &handler);
  }

	bool TraMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	{
		//load mapping
		CVMappings mapping;
		CVMappingFile().load(File::find("/MAPPING/TraML-mapping.xml"),mapping);
		
		//load cvs
		ControlledVocabulary cv;
		cv.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		
		//validate
		Internal::TraMLValidator v(mapping, cv);
		bool result = v.validate(filename, errors, warnings);
		
		return result;
	}


}// namespace OpenMS

