// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Andreas Bertsch, Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzIdentMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/HANDLERS/MzIdentMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

	MzIdentMLFile::MzIdentMLFile()
		: XMLFile("/SCHEMAS/mzIdentML1.1.0.xsd","1.1.0")
	{
	}

	MzIdentMLFile::~MzIdentMLFile()
	{
	}

  void MzIdentMLFile::load(const String& filename, Identification& id)
  {
  	Internal::MzIdentMLHandler handler(id, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void MzIdentMLFile::store(const String& filename, const Identification& id) const
  {
  	Internal::MzIdentMLHandler handler(id, filename, schema_version_, *this);
    save_(filename, &handler);
  }

  void MzIdentMLFile::store(const String& filename, const std::vector<ProteinIdentification>& poid, const std::vector<PeptideIdentification>& peid) const
  {
  	Internal::MzIdentMLHandler handler(poid, peid, filename, schema_version_, *this);
    save_(filename, &handler);
  }

	bool MzIdentMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	{
		//load mapping
		CVMappings mapping;
		CVMappingFile().load(File::find("/MAPPING/mzIdentML-mapping.xml"),mapping);
		
		//load cvs
		ControlledVocabulary cv;
		cv.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("BTO",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		
		//validate
		Internal::MzIdentMLValidator v(mapping, cv);
		bool result = v.validate(filename, errors, warnings);
		
		return result;
	}

}// namespace OpenMS

