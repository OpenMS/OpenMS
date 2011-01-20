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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
#include <OpenMS/FORMAT/CVMappingFile.h>

namespace OpenMS
{
	MzDataFile::MzDataFile()
		: XMLFile("/SCHEMAS/mzData_1_05.xsd","1.05"),
			options_()
	{
	}

	MzDataFile::~MzDataFile()
	{
	}

	PeakFileOptions& MzDataFile::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& MzDataFile::getOptions() const
  {
  	return options_;
  }

	bool MzDataFile::isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings)
	{
		//load mapping
		CVMappings mapping;
		CVMappingFile().load(File::find("/MAPPING/mzdata-mapping.xml"),mapping);
		
		//load cvs
		ControlledVocabulary cv;
		cv.loadFromOBO("PSI",File::find("/CV/psi-mzdata.obo"));
		
		//validate
		Internal::MzDataValidator v(mapping, cv);
		bool result = v.validate(filename, errors, warnings);
		
		return result;
	}
	
}// namespace OpenMS

