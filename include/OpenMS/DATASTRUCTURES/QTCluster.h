// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_QTSUBSET_H
#define OPENMS_DATASTRUCTURES_QTSUBSET_H

#include<OpenMS/DATASTRUCTURES/GridFeature.h>

namespace OpenMS {

/**
	@brief A representation of a QT cluster.

	A QT cluster holds all elements, which are clustered together during feature linking.
	To adjust the quality of a cluster, edit the getQuality() method.
	There exists the possibilty to use cluster members for the quality assignment as well as neighbor points whose distances are below a specified threshold.

	@ingroup Datastructures
*/

class OPENMS_DLLAPI QTCluster {
private:
	/**
	 * @brief pointer to the cluster center
	 */
	GridFeature* center_point;
	/**
	 * @brief cluster members. Each member is stored at the map index from which it comes from.
	 */
	std::map<Size,GridFeature*> cluster_members;
	/**
	 * @brief minimal distances from the cluster center to the element.
	 * The index corresponds to the map index the element comes from.
	 */
	std::vector<DoubleReal> min_distances;
	/**
	 * @brief cluster neighbors whose distances are below the specified maximal distance
	 */
	std::map<Size,std::set<GridFeature*> > cluster_neighbors;
	/**
	 * @biref maximal distance of two points
	 */
	DoubleReal max_distance_;
	/**
	 * @brief base constructor. Not accessible
	 */
	QTCluster();
public:
	/**
	 * @brief detailed constructor
	 * @param center_point_ pointer to the center point
	 * @param maps_size number of maps used in feature linking
	 * @param max_distance maximal distance of two points
	 */
	QTCluster(GridFeature* center_point_,Size maps_size,DoubleReal max_distance);
	/**
	 * @brief destructor
	 */
	virtual ~QTCluster();
	/**
	 * @brief returns the RT value of the cluster center
	 */
	DoubleReal getCenterRT();
	/**
	 * @brief returns the m/z valus of the cluster center
	 */
	DoubleReal getCenterMZ();
	/**
	 * @brief returns the size of the cluster
	 */
	Size size() const;
	/**
	 * @brief comparator
	 */
	bool operator<(const QTCluster &cluster) const;
	/**
	 * @brief adds a new element to the cluster
	 * @param element the element to be added
	 * @param distance of the element to the center point
	 */
	void add(GridFeature* element,DoubleReal distance);
	/**
	 * @brief non-mutable access to the cluster members
	 */
	std::map<Size,GridFeature*> getClusterMembers() const;
	/**
	 * @brief returns the cluster quality
	 */
	DoubleReal getQuality() const;
};
}
#endif /* QTSUBSET_H_ */
