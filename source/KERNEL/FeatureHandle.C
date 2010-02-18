// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

namespace OpenMS
{
	
	FeatureHandle::FeatureHandle(UInt64 map_index, const ConsensusFeature& point)
		: Peak2D(point),
		  UniqueIdInterface(point),
			map_index_(map_index),
			charge_(point.getCharge())
	{
	}

  std::ostream& operator << (std::ostream& os, const FeatureHandle& cons)
  {
    os  << "---------- FeatureHandle -----------------\n"
		    << "RT: " << cons.getRT()<< std::endl
		    << "m/z: " << cons.getMZ()<< std::endl
		    << "Intensity: " << cons.getIntensity() << std::endl
		    << "Map Index: " << cons.getMapIndex() << std::endl
		    << "Element Id: " << cons.getUniqueId() << std::endl;
    return os;
  }
} 
