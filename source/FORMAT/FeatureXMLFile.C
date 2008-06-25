// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include<OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

namespace OpenMS 
{
	FeatureXMLFile::FeatureXMLFile()
		: Internal::XMLFile("/SCHEMAS/FeatureXML_1_2.xsd","1.2"),
			options_()
	{
	}
	
	FeatureXMLFile::~FeatureXMLFile()
	{
	}

	void FeatureXMLFile::load(String filename, FeatureMap<>& feature_map) throw (Exception::FileNotFound, Exception::ParseError)
	{
		feature_map.clear();
		Internal::FeatureXMLHandler handler(feature_map,filename,schema_version_);
		handler.setOptions(options_);
		parse_(filename, &handler);
	}

	void FeatureXMLFile::store(String filename, const FeatureMap<>& feature_map) const throw (Exception::UnableToCreateFile)
	{
		Internal::FeatureXMLHandler handler(feature_map,filename,schema_version_);
		save_(filename, &handler);
	}

	PeakFileOptions& FeatureXMLFile::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& FeatureXMLFile::getOptions() const
  {
  	return options_;
  }

}
