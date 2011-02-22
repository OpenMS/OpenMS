// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERINGMETHOD_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERINGMETHOD_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DataSubset.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>

namespace OpenMS
{

/**
 * @brief Base class for all distance computing hierarchical clustering methods, which are used in HashClustering.
 * @see HashClustering
 * @ingroup SpectraClustering
 */

class OPENMS_DLLAPI ClusteringMethod {

public:
	/**
	 * @brief scale the distances in RT direction to get more oval clusters
	 */
	DoubleReal rt_scaling;
	/**
	 * @brief default constructor
	 */
	ClusteringMethod();
	/**
	 * @brief detailed constructor
	 * @param rt_scaling value for scaling the distances in RT direction
	 */
	ClusteringMethod(DoubleReal rt_scaling_);
	/**
	 * @brief gets the distance between two DataSubsets (calculation may be different than distance calculation of two DataPoints)
	 * @param subset1 first subset
	 * @param subset2 second subset
	 */
	virtual DoubleReal getDistance(DataSubset& subset1,DataSubset& subset2) = 0;
	/**
	 * @brief gets the distance between two DataPoints
	 * @param DataPoint1 first data point
	 * @param DataPoint1 second data point
	 */
	virtual DoubleReal getDistance(DataPoint& point1,DataPoint& point2) = 0;
};
}
#endif /* CLUSTERINGMETHOD_H_ */
