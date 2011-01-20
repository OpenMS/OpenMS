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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>

// Derived classes are included here
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

namespace OpenMS
{
	//register products here
	void MapAlignmentEvaluationAlgorithm::registerChildren()
	{
		Factory<MapAlignmentEvaluationAlgorithm>::registerProduct ( MapAlignmentEvaluationAlgorithmPrecision::getProductName(), &MapAlignmentEvaluationAlgorithmPrecision::create );
		Factory<MapAlignmentEvaluationAlgorithm>::registerProduct ( MapAlignmentEvaluationAlgorithmRecall::getProductName(), &MapAlignmentEvaluationAlgorithmRecall::create );
	}

	// TODO consider using (RT,MZ,IT) as a unique identifier ?
	bool MapAlignmentEvaluationAlgorithm::isSameHandle(const FeatureHandle & lhs, const FeatureHandle & rhs, const DoubleReal& rt_dev, const DoubleReal& mz_dev, const Peak2D::IntensityType& int_dev, const bool use_charge)
	{
#if 1
		// use (RT,MZ,IT) as "unique" identifier 
		if ( fabs ( lhs.getRT() - rhs.getRT() ) > rt_dev ) return false;  // TODO MAGIC_ALERT
		if ( fabs ( lhs.getMZ() - rhs.getMZ() ) > mz_dev ) return false;  // TODO MAGIC_ALERT
		if ( fabs ( lhs.getIntensity() - rhs.getIntensity() ) > int_dev ) return false;  // TODO MAGIC_ALERT
		if ( use_charge && (lhs.getCharge() != rhs.getCharge())) return false;
		return true;
#else
		// use (map index, element index) as unique identifier
		return lhs.getMapIndex() == rhs.getMapIndex() && lhs.getElementIndex() == rhs.getElementIndex();
#endif
}

	MapAlignmentEvaluationAlgorithm::MapAlignmentEvaluationAlgorithm()
	{
	}

	MapAlignmentEvaluationAlgorithm::~MapAlignmentEvaluationAlgorithm()
	{
	}

} 
