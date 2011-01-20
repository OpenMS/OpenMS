// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

namespace OpenMS
{
  FeatureFinder::FeatureFinder()
		: flags_()
	{			
	}

  FeatureFinder::~FeatureFinder()
	{
	}

	Param FeatureFinder::getParameters(const String& algorithm_name) const
	{
		Param tmp;
		if (algorithm_name!="none")
		{
			FeatureFinderAlgorithm<Peak1D, Feature>* a = Factory<FeatureFinderAlgorithm<Peak1D, Feature> >::create(algorithm_name);
			tmp.insert("", a->getDefaultParameters());
			delete(a);
		}	
		return tmp;
	}
}

