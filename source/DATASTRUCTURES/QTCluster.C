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
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/QTCluster.h>

using namespace std;

namespace OpenMS 
{

	QTCluster::QTCluster()
	{
	}

	QTCluster::QTCluster(GridFeature* center_point, Size num_maps, 
											 DoubleReal max_distance) :
		center_point_(center_point), neighbors_(), max_distance_(max_distance), 
		num_maps_(num_maps), quality_(0.0), changed_(false)
	{
	}

	QTCluster::~QTCluster()
	{
	}

	QTCluster& QTCluster::operator=(const QTCluster& rhs)
	{
		center_point_ = rhs.center_point_;
		neighbors_ = rhs.neighbors_;
		max_distance_ = rhs.max_distance_;
		num_maps_ = rhs.num_maps_;
		quality_ = rhs.quality_;
		changed_ = rhs.changed_;
		return *this;
	}

	DoubleReal QTCluster::getCenterRT()
	{
		return center_point_->rt;
	}

	DoubleReal QTCluster::getCenterMZ()
	{
		return center_point_->mz;
	}

	Size QTCluster::size() const
	{
		return neighbors_.size() + 1; // + 1 for the center
	}

	bool QTCluster::operator<(QTCluster &cluster)
	{
		return (this->getQuality() < cluster.getQuality());
	}

	void QTCluster::add(GridFeature* element, DoubleReal distance)
	{
		Size map_index = element->getMapIndex();
		if (map_index != center_point_->getMapIndex())
		{
			neighbors_[map_index].insert(make_pair(distance, element));
			changed_ = true;
		}
	}

	void QTCluster::getElements(map<Size, GridFeature*>& elements) const
	{
		elements.clear();
		elements[center_point_->getMapIndex()] = center_point_;
		for (NeighborMap::const_iterator it = neighbors_.begin(); 
				 it != neighbors_.end(); ++it)
		{
			elements[it->first] = it->second.begin()->second;
		}
	}

	bool QTCluster::update(const map<Size, GridFeature*>& removed)
	{
		// check if the cluster center was removed:
		for (map<Size, GridFeature*>::const_iterator rm_it = removed.begin(); 
				 rm_it != removed.end(); ++rm_it)
		{
			if (rm_it->second == center_point_) return false;
		}
		// update the cluster contents:
		for (map<Size, GridFeature*>::const_iterator rm_it = removed.begin(); 
				 rm_it != removed.end(); ++rm_it)
		{
			NeighborMap::iterator pos = neighbors_.find(rm_it->first);
			if (pos == neighbors_.end()) continue; // no points from this map
			for (multimap<DoubleReal, GridFeature*>::iterator feat_it = 
						 pos->second.begin(); feat_it != pos->second.end(); ++feat_it)
			{
				if (feat_it->second == rm_it->second) // remove this neighbor
				{
					pos->second.erase(feat_it);
					changed_ = true;
					break;
				}
			}
			if (pos->second.empty()) // only neighbor from this map was just removed
			{
				neighbors_.erase(pos);
			}
		}
		return true;
	}

	DoubleReal QTCluster::getQuality()
	{
		if (changed_)
		{
			computeQuality_();
			changed_ = false;
		}
		return quality_;
	}

	void QTCluster::computeQuality_()
	{
		DoubleReal internal_distance = 0.0;
		Size counter = 0;
		for (NeighborMap::iterator it = neighbors_.begin(); it != neighbors_.end();
				 ++it)
		{
			internal_distance += it->second.begin()->first;
			counter++;
		}
		Size num_other = num_maps_ - 1;
		// add max. distance for missing cluster elements:
		internal_distance += (num_other - counter) * max_distance_;
		// normalize:
		internal_distance /= num_other;
		quality_ = (max_distance_ - internal_distance) / max_distance_;
	}

}
