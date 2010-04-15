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
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LabeledPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>

#include <set>

namespace OpenMS
{

	BaseGroupFinder::BaseGroupFinder()
		: DefaultParamHandler("BaseGroupFinder")
	{
	}

	BaseGroupFinder::~BaseGroupFinder()
	{
	}

	void BaseGroupFinder::registerChildren()
  {
    Factory< BaseGroupFinder>::registerProduct(SimplePairFinder::   getProductName(), &SimplePairFinder::   create);
    Factory< BaseGroupFinder>::registerProduct(LabeledPairFinder::  getProductName(), &LabeledPairFinder::  create);
    Factory< BaseGroupFinder>::registerProduct(StablePairFinder::   getProductName(), &StablePairFinder::   create);
  }

	void BaseGroupFinder::checkIds_(const std::vector<ConsensusMap>& maps) const
	{
		std::set<Size> used_ids;
		for (Size i=0; i< maps.size(); ++i)
		{
			const ConsensusMap& map = maps[i];
			for (ConsensusMap::FileDescriptions::const_iterator it = map.getFileDescriptions().begin(); it!=map.getFileDescriptions().end(); ++it)
			{
				if (used_ids.find(it->first)!=used_ids.end())
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"file ids have to be unique");
				}
				else
				{
					used_ids.insert(it->first);
				}
			}
		}
	}

}
