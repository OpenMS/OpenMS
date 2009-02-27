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
// $Maintainer: Katharina Albers, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/CaapEvalAlgorithmRecall.h>

namespace OpenMS
{
	//register products here
	void CaapEvalAlgorithm::registerChildren()
	{
		Factory<CaapEvalAlgorithm>::registerProduct ( CaapEvalAlgorithmPrecision::getProductName(), &CaapEvalAlgorithmPrecision::create );
		Factory<CaapEvalAlgorithm>::registerProduct ( CaapEvalAlgorithmRecall::getProductName(), &CaapEvalAlgorithmRecall::create );
	}

	// TODO consider using (RT,MZ,IT) as a unique identifier ?
	bool CaapEvalAlgorithm::isSameHandle(const FeatureHandle & lhs, const FeatureHandle & rhs) //const		//geht nicht: fehler:/home/bude/albers/RAID/cmakeOpenMS/source/ANALYSIS/MAPMATCHING/CaapEvalAlgorithm.C:46: error: non-member function 'bool OpenMS::isSameHandle(const OpenMS::FeatureHandle&, const OpenMS::FeatureHandle&)' cannot have cv-qualifier
	{
#if 1
		// use (RT,MZ,IT) as "unique" identifier 
		if ( fabs ( lhs.getRT() - rhs.getRT() ) > 0.1 ) return false;  // TODO MAGIC_ALERT
		if ( fabs ( lhs.getMZ() - rhs.getMZ() ) > 0.1 ) return false;  // TODO MAGIC_ALERT
		if ( fabs ( lhs.getIntensity() - rhs.getIntensity() ) > 100 ) return false;  // TODO MAGIC_ALERT
		return true;
#else
		// use (map index, element index) as unique identifier
		return lhs.getMapIndex() == rhs.getMapIndex() && lhs.getElementIndex() == rhs.getElementIndex();
#endif
}

	CaapEvalAlgorithm::CaapEvalAlgorithm()
		: FactoryProduct("CaapEvalAlgorithm")
	{
	}

	CaapEvalAlgorithm::~CaapEvalAlgorithm()
	{
	}

} 
