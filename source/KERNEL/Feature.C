// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Feature.h>

using namespace std;

namespace OpenMS
{
	Feature& Feature::operator = (const Feature& rhs)
	{
		if (this==&rhs) return *this;
		
		RichPeak2D::operator=(rhs);
		overall_quality_  			= rhs.overall_quality_;
		copy(rhs.qualities_,rhs.qualities_+2,qualities_);
		model_desc_							= rhs.model_desc_;
		convex_hulls_						= rhs.convex_hulls_;
		convex_hulls_modified_	= rhs.convex_hulls_modified_;
		convex_hull_      			= rhs.convex_hull_;
		charge_       					= rhs.charge_;
		identifications_				= rhs.identifications_;
		subordinates_						= rhs.subordinates_;
		
		return *this;
	}

	bool Feature::operator == (const Feature& rhs) const
	{
		return (RichPeak2D::operator == (rhs) 
						&& (overall_quality_   == rhs.overall_quality_)
						&& (charge_ == rhs.charge_)
						&& equal(qualities_, qualities_+2, rhs.qualities_)
						&& (model_desc_ == rhs.model_desc_)
						&& (convex_hulls_ == rhs.convex_hulls_)
						&& (subordinates_  == rhs.subordinates_))
            ;
	}
	
	bool Feature::encloses(DoubleReal rt, DoubleReal mz) const
	{
		ConvexHull2D::PointType tmp(rt,mz);
		for (vector<ConvexHull2D>::const_iterator	it=convex_hulls_.begin(); it!=convex_hulls_.end(); ++it)
		{
			if (it->encloses(tmp)) return true;
		}
		return false;
	}

  ConvexHull2D& Feature::getConvexHull() const
  {
  	//recalculate convex hull if necessary
  	if (convex_hulls_modified_)
  	{
  		//only one mass trace convex hull => use it as overall convex hull
  		if (convex_hulls_.size()==1)
  		{
  			convex_hull_ = convex_hulls_[0];
  		}
  		else
  		{
				ConvexHull2D::PointArrayType all_points;
				for (Size hull=0; hull<convex_hulls_.size(); ++hull)
				{
					all_points.insert(all_points.end(), convex_hulls_[hull].getPoints().begin(), convex_hulls_[hull].getPoints().end());
				}
				convex_hull_ = all_points;
  		}
  		
  		convex_hulls_modified_ = false;
  	}

  	return convex_hull_;
  }

}
