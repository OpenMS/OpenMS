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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H
#define OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H

#include <cmath>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

	/**
		@brief CompareFunctor for 2Dpoints

		each 2D point as a pair of Real holds a Real coordinate for each Dimension

		@ingroup DummyComparison
	*/

	class OPENMS_DLLAPI EuclideanSimilarity
	{
	private:

		Real scale_;

	public:

		/// default constructor
		EuclideanSimilarity();

		/// copy constructor
		EuclideanSimilarity(const EuclideanSimilarity& source);

		/// destructor
		virtual ~EuclideanSimilarity();

		/// assignment operator
		EuclideanSimilarity& operator = (const EuclideanSimilarity& source);


	/**
		@brief calculates similarity between two points in euclidean space

		@param a a pair of Real, giving the x and the y coordinates of the first point
		@param b a pair of Real, giving the x and the y coordinates of the second point

		calculates similarity from the euclidean distance between given 2D points, scaled in [0,1] @see setScale
	*/
	Real operator () (const std::pair<Real,Real>& a, const std::pair<Real,Real>& b) const;

	/**
		@brief calculates self similarity, will yield 0

		@param c a pair of Real, giving the x and the y coordinates

	*/
	Real operator () (const std::pair<Real,Real>& c) const;

	/**
		@brief clusters the indices according to their respective element distances

		@param x Real value to scale the result
		@throw Exception::DivisionByZero if scaling is unapplicable because it is 0

		sets the scale so that similarities can be correctly calculated from distances. Should be set so that the greatest distance in a chosen set will be scales to 1 (i.e. @p x = greatest possible distance in the set)
*/
	void setScale (Real x);

	};

}
#endif //OPENMS_COMPARISON_CLUSTERING_EUCLIDEANSIMILARITY_H
