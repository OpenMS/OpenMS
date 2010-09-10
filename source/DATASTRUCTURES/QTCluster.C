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


#include <OpenMS/DATASTRUCTURES/QTCluster.h>

namespace OpenMS {

QTCluster::QTCluster()
{

}
QTCluster::QTCluster(GridFeature* center_point_,Size maps_size,DoubleReal max_distance) : center_point(center_point_), min_distances(maps_size,max_distance),max_distance_(max_distance)
{
	Size map_index=center_point_->getMapIndex();
	min_distances[map_index]=0.0;
	cluster_members[map_index]=center_point;
}

QTCluster::~QTCluster() {
	// TODO Auto-generated destructor stub
}

DoubleReal QTCluster::getCenterRT()
{
	return center_point->rt;
}

DoubleReal QTCluster::getCenterMZ()
{
	return center_point->mz;
}

Size QTCluster::size() const
{
	return cluster_members.size();
}

bool QTCluster::operator<(const QTCluster &cluster) const
{
	return (this->getQuality() < cluster.getQuality());
}

void QTCluster::add(GridFeature* element,DoubleReal distance)
{
	Size map_index=element->getMapIndex();
	if (distance < min_distances[map_index])
	{
		min_distances[map_index]=distance;
		std::map<Size,GridFeature*>::iterator feature_position=cluster_members.find(map_index);
		if (feature_position!=cluster_members.end())
		{
			cluster_neighbors[map_index].insert(feature_position->second);
		}
		cluster_members[map_index]=element;
	}
}

std::map<Size,GridFeature*> QTCluster::getClusterMembers() const
{
	return cluster_members;
}


DoubleReal QTCluster::getQuality() const
{
	DoubleReal internal_distance=0.0;
	for (std::vector<DoubleReal>::const_iterator it=min_distances.begin();it!=min_distances.end();++it)
	{
		internal_distance+=*it;
	}
	internal_distance/=min_distances.size()-1;
	DoubleReal quality=(max_distance_-internal_distance)/max_distance_;
	return quality;
}

}
