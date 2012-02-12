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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>
#include <OpenMS/FORMAT/VALIDATORS/XMLValidator.h>
#include <OpenMS/FORMAT/TextFile.h>


namespace OpenMS
{

	MzMLFile::MzMLFile()
		: XMLFile("/SCHEMAS/mzML_1_10.xsd","1.1.0"),
			indexed_schema_location_("/SCHEMAS/mzML_idx_1_10.xsd")
	{
	}

	MzMLFile::~MzMLFile()
	{
	}

	PeakFileOptions& MzMLFile::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& MzMLFile::getOptions() const
  {
  	return options_;
  }
	
	//reimplemented in order to handle index MzML
	bool MzMLFile::isValid(const String& filename, std::ostream& os) 
	{
		//determine if this is indexed mzML or not
		bool indexed = false;
		TextFile file(filename,true,4);
		if (file.concatenate().hasSubstring("<indexedmzML"))
		{
			indexed = true;
		}
		// find the corresponding schema
		String current_location;
		if (indexed)
		{
			current_location = File::find(indexed_schema_location_);
		}
		else
		{
			current_location = File::find(schema_location_);
		}
		
		return XMLValidator().isValid(filename,current_location,os);
	}
	
	bool MzMLFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	{
		//load mapping
		CVMappings mapping;
		CVMappingFile().load(File::find("/MAPPING/ms-mapping.xml"),mapping);
		
		//load cvs
		ControlledVocabulary cv;
		cv.loadFromOBO("MS",File::find("/CV/psi-ms.obo"));
		cv.loadFromOBO("PATO",File::find("/CV/quality.obo"));
		cv.loadFromOBO("UO",File::find("/CV/unit.obo"));
		cv.loadFromOBO("BTO",File::find("/CV/brenda.obo"));
		cv.loadFromOBO("GO",File::find("/CV/goslim_goa.obo"));
		
		//validate
		Internal::MzMLValidator v(mapping, cv);
		bool result = v.validate(filename, errors, warnings);
		
		return result;
	}

}// namespace OpenMS

