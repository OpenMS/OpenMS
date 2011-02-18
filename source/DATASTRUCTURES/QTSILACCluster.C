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


#include <OpenMS/DATASTRUCTURES/QTSILACCluster.h>

namespace OpenMS {

QTSILACCluster::QTSILACCluster()
{

}
QTSILACCluster::QTSILACCluster(DataPoint* center_point_) : center_point(center_point_)
{
	cluster_members.insert(center_point);
}

QTSILACCluster::~QTSILACCluster() {
}

DoubleReal QTSILACCluster::getCenterRT()
{
	return center_point->rt;
}

DoubleReal QTSILACCluster::getCenterMZ()
{
	return center_point->mz;
}

Size QTSILACCluster::size() const
{
	return cluster_members.size();
}

bool QTSILACCluster::operator<(const QTSILACCluster &cluster) const
{
	return (this->size() < cluster.size());
}

void QTSILACCluster::add(DataPoint* element)
{
	cluster_members.insert(element);
}

bool QTSILACCluster::contains(DataPoint* element)
{
	std::set<DataPoint*>::iterator pos = cluster_members.find(element);
	return pos != cluster_members.end();
}

std::set<DataPoint*> QTSILACCluster::getClusterMembers()
{
	return cluster_members;
}

std::pair<DoubleReal,DoubleReal> QTSILACCluster::getDiameters(DataPoint* point)
{
	DoubleReal point_mz=point->mz;
	DoubleReal point_rt=point->rt;
	DoubleReal rt_diameter=std::numeric_limits<Real>::max();
	DoubleReal mz_diameter=0.0;
	for (std::set<DataPoint*>::iterator it=cluster_members.begin();it!=cluster_members.end();++it)
	{
		DoubleReal rt_dist=std::abs((*it)->rt-point_rt);
		if (rt_dist < rt_diameter)
			rt_diameter=rt_dist;
		DoubleReal mz_dist=std::abs((*it)->mz-point_mz);
		if (mz_dist > mz_diameter)
			mz_diameter=mz_dist;
	}
	return std::make_pair(rt_diameter,mz_diameter);
}


}
