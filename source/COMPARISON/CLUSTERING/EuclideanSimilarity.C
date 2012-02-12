// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/EuclideanSimilarity.h>

namespace OpenMS
{
	EuclideanSimilarity::EuclideanSimilarity() : scale_(1)
	{
	}

	EuclideanSimilarity::EuclideanSimilarity(const EuclideanSimilarity& source) : scale_(source.scale_)
	{
	}

	EuclideanSimilarity::~EuclideanSimilarity()
	{
	}

	EuclideanSimilarity& EuclideanSimilarity::operator = (const EuclideanSimilarity& source)
	{
		if (this != &source)
		{
			scale_ = source.scale_;
		}
		return *this;
	}

	Real EuclideanSimilarity::operator () (const std::pair<Real,Real>& c) const
	{
		return operator () (c, c);
	}

	// calculates euclidean distance between two points
	Real EuclideanSimilarity::operator () (const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) const
	{
		if(scale_==0)
		{
			//unapplicable scaling
			throw Exception::DivisionByZero(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		return 1-(sqrt( (a.first-b.first) * (a.first-b.first) + (a.second-b.second) * (a.second-b.second) )/scale_);
	}

	void EuclideanSimilarity::setScale(Real x)
	{
		scale_ = x;
	}

}
