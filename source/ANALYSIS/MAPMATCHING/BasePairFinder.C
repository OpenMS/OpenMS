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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>

namespace OpenMS
{

	BasePairFinder::BasePairFinder()
		: FactoryProduct("BasePairFinder")
	{
	}

	BasePairFinder::~BasePairFinder()
	{
	}

	void BasePairFinder::registerChildren()
  {
    Factory< BasePairFinder>::registerProduct(SimplePairFinder::getProductName(), &SimplePairFinder::create);
    Factory< BasePairFinder>::registerProduct(DelaunayPairFinder::getProductName(), &DelaunayPairFinder::create);
    Factory< BasePairFinder>::registerProduct(PairMatcher::getProductName(), &PairMatcher::create);
  }

} 
