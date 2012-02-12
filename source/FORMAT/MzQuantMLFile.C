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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/HANDLERS/MzQuantMLHandler.h>
#include <OpenMS/SYSTEM/File.h>
//~ #include <OpenMS/FORMAT/VALIDATORS/MzQuantMLValidator.h>

namespace OpenMS
{

	MzQuantMLFile::MzQuantMLFile()
		: XMLFile("/SCHEMAS/mzQuantML0.1.7.xsd","0.1.7")
	{
	}

	MzQuantMLFile::~MzQuantMLFile()
	{
	}

	void MzQuantMLFile::load(const String& filename, ConsensusMap& cm)
	{
		Internal::MzQuantMLHandler handler(cm, filename, schema_version_, *this);
		parse_(filename, &handler);
	}

	void MzQuantMLFile::store(const String& filename, const ConsensusMap& cm) const
	{
		Internal::MzQuantMLHandler handler(cm, filename, schema_version_, *this);
		save_(filename, &handler);
	}

	//~ TODO
	//~ bool MzQuantMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	//~ {
		//~ //load mapping
		//~ CVMappings mapping;
		//~ CVMappingFile().load(File::find("/MAPPING/mzQuantML-mapping.xml"),mapping);

		//~ //load cvs
		//~ ControlledVocabulary cv;
		//~ cv.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		//~ cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		//~ cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		//~ cv.loadFromOBO("BTO",File::find("/CV/brenda.obo"));
		//~ cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));

		//~ //validate
		//~ Internal::MzQuantMLValidator v(mapping, cv);
		//~ bool result = v.validate(filename, errors, warnings);

		//~ return result;
	//~ }

}// namespace OpenMS

