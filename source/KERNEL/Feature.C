// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Feature.h>

namespace OpenMS
{
	Feature& Feature::operator = (const Feature& rhs)
	{
		if (this==&rhs) return *this;
		
		Peak2D::operator = (rhs);
		overall_quality_  = rhs.overall_quality_;
		std::copy(rhs.qualities_,rhs.qualities_+2,qualities_);
		model_desc_       = rhs.model_desc_;
		convex_hulls_     = rhs.convex_hulls_;
		charge_           = rhs.charge_;
		
		return *this;
	}

	bool Feature::operator == (const Feature& rhs) const
	{
		return (Peak2D::operator == (rhs) 
						&& (overall_quality_   == rhs.overall_quality_)
						&& (charge_ == rhs.charge_)
						&& std::equal(qualities_, qualities_+2, rhs.qualities_)
						&& (model_desc_ == rhs.model_desc_)
						&& (convex_hulls_ == rhs.convex_hulls_));
	}

	DBoundingBox<2> Feature::getBoundingBox() const
	{
		DBoundingBox<2> bb, tmp;
		
		for (ConvexHullVector::const_iterator	it=convex_hulls_.begin(); it!=convex_hulls_.end(); ++it)
		{
			tmp = it->getBoundingBox();
			bb.enlarge(tmp.min());
			bb.enlarge(tmp.max());
		}
		
		return bb;
	}
    
}
