// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/AnalysisXMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/AnalysisXMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{

	AnalysisXMLFile::AnalysisXMLFile()
		: XMLFile("/SCHEMAS/analysisXML_1_0.xsd","1.10")
	{
	}

	AnalysisXMLFile::~AnalysisXMLFile()
	{
	}

  void AnalysisXMLFile::load(const String& filename, Identification& id)
  {
  	Internal::AnalysisXMLHandler handler(id, filename, schema_version_, *this);
    parse_(filename, &handler);
  }

  void AnalysisXMLFile::store(const String& filename, const Identification& id) const
  {
  	Internal::AnalysisXMLHandler handler(id, filename, schema_version_, *this);
    save_(filename, &handler);
  }

	bool AnalysisXMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	{
		//load mapping
		CVMappings mapping;
		CVMappingFile().load(File::find("/MAPPING/pi-mapping.xml"),mapping);
		
		//load cvs
		ControlledVocabulary cv;
		cv.loadFromOBO("PI",File::find("/CV/psi-pi.obo"));
		cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("BTO",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		
		//validate
		Internal::AnalysisXMLValidator v(mapping, cv);
		bool result = v.validate(filename, errors, warnings);
		
		return result;
	}

}// namespace OpenMS

