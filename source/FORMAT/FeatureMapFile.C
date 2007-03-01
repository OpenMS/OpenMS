// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include<OpenMS/FORMAT/FeatureMapFile.h>
#include <OpenMS/FORMAT/HANDLERS/FeatureMapHandler.h>

namespace OpenMS 
{
	FeatureMapFile::FeatureMapFile()
		: SchemaFile()
	{
	}
	FeatureMapFile::~FeatureMapFile()
	{
	}

	void FeatureMapFile::load(String filename, FeatureMap<>& feature_map) throw (Exception::FileNotFound, Exception::ParseError)
	{
		feature_map.clear();
		Internal::FeatureMapHandler<2> handler(feature_map,filename);
		handler.setOptions(options_);
		parse_(filename, &handler);
	}

	void FeatureMapFile::store(String filename, const FeatureMap<>& feature_map) const throw (Exception::UnableToCreateFile)
	{
		if (feature_map.empty()) return;
		Internal::FeatureMapHandler<2> handler(feature_map,filename);
		save_(filename, &handler);
	}
}
