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


#ifndef OPENMS_DATASTRUCTURES_QTCLUSTER_H
#define OPENMS_DATASTRUCTURES_QTCLUSTER_H

#include<OpenMS/DATASTRUCTURES/DataPoint.h>

namespace OpenMS {
/**
 * @brief QT clusters are used in QTClustering. They consists of a set of data points, while data point defines the center of the cluster.
 * QT clusters compute a two-dimensional diameter. The RT diameter corresponds to the maximal gap in RT direction of cluster.
 * The m/z diameter corresponds to the maximal cluster extent in m/z direction.
 */
class OPENMS_DLLAPI QTSILACCluster {
private:
	/**
	 * @brief the cluster center
	 */
	DataPoint* center_point;
	/**
	 * @brief members of the cluster
	 */
	std::set<DataPoint*> cluster_members;
	/**
	 * default constructor
	 */
	QTSILACCluster();
public:
	/**
	 * @brief detailed constructor
	 * @param center_point cluster center
	 */
	QTSILACCluster(DataPoint* center_point_);
	/**
	 * @brief destructor
	 */
	virtual ~QTSILACCluster();
	/**
	 * @brief gets the center RT position of the cluster
	 */
	DoubleReal getCenterRT();
	/**
	 * @brief gets the center m/z position of the cluster
	 */
	DoubleReal getCenterMZ();
	/**
	 * @brief gets the size of the cluster
	 */
	Size size() const;
	/**
	 * @brief cluster comparator
	 */
	bool operator<(const QTSILACCluster &cluster) const;
	/**
	* @brief adds an element to the cluster
	* @param element the element to be added
	 */
	void add(DataPoint* element);
	/**
	 * @brief non-mutable access to the cluster members
	 */
	std::set<DataPoint*> getClusterMembers();
	/**
	 * @brief gets the diameter pair of the cluster
	 * first: rt_diameter, second: mz_diameter
	 * The diameters are computed by considering a further data point, which is a candidate to be added to the cluster.
	 * @param point data point, which should be added to the cluster
	 */
	std::pair<DoubleReal,DoubleReal> getDiameters(DataPoint* point);
	/**
	 * @brief checks if an element is contained in the cluster
	 */
	bool contains(DataPoint* element);
};
}
#endif /* QTSUBSET_H_ */
