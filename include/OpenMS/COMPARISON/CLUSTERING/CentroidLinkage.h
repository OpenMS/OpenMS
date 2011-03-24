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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass$
// --------------------------------------------------------------------------

#ifndef OPENMS_COMPARISON_CLUSTERING_CENTROIDLINKAGE_H
#define OPENMS_COMPARISON_CLUSTERING_CENTROIDLINKAGE_H

#include <OpenMS/COMPARISON/CLUSTERING/ClusteringMethod.h>

namespace OpenMS
{

/**
 * @brief Clustering method with computes the distances of two clusters as the distance of their cluster centroids.
 *
 * Also known as <i>Unweighted Pair-Group Method using Centroids</i>, or UPGMC.
 * @ingroup SpectraClustering
 */

  class OPENMS_DLLAPI CentroidLinkage : public ClusteringMethod {
  public:

  /**
   * @brief default constructor
   */
   CentroidLinkage();

  /**
   * @brief destructor
   */
   ~CentroidLinkage();

  /**
	 * @brief detailed constructor
	 * @param rt_scaling_ value fpr sclaling distances in RT direction
	 */
   CentroidLinkage(DoubleReal rt_scaling_);

  /**
	 * @brief gets the distance between two DataSubsets
	 * @param subset1 first subset
	 * @param subset2 second subset
	 */
   DoubleReal getDistance(DataSubset& subset1,DataSubset& subset2);

  /**
	 * @brief gets the distance between two DataPoints
	 * @param DataPoint1 first data point
	 * @param DataPoint1 second data point
	 */
   DoubleReal getDistance(DataPoint& point1,DataPoint& point2);

  };
}

#endif /* CENTROIDLINKAGE_H_ */
