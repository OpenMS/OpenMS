// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

//Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

namespace OpenMS
{
	//register products here
	void MapAlignmentAlgorithm::registerChildren()
	{
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmSpectrumAlignment::getProductName(), &MapAlignmentAlgorithmSpectrumAlignment::create );
		Factory<MapAlignmentAlgorithm>::registerProduct ( MapAlignmentAlgorithmPoseClustering::getProductName(), &MapAlignmentAlgorithmPoseClustering::create );
	}

	MapAlignmentAlgorithm::MapAlignmentAlgorithm()
		: FactoryProduct("MapAlignmentAlgorithm")
	{
	}

	MapAlignmentAlgorithm::~MapAlignmentAlgorithm()
	{
	}

	void MapAlignmentAlgorithm::alignPeakMaps(std::vector< MSExperiment<> >&, std::vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}

	void MapAlignmentAlgorithm::alignFeatureMaps(std::vector< FeatureMap<> >&, std::vector<TransformationDescription>&)
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);				
	}

} 
