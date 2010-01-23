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

#ifndef OPENMS_MATH_STATISTICS_AVERAGEPOSITION_H
#define OPENMS_MATH_STATISTICS_AVERAGEPOSITION_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

namespace OpenMS
{

	namespace Math
	{
		/**
			 @brief Maintain an average position by summing up positions with
			 weights.
			 
			 @improvement A lot of convenience methods could be added.  For example,
			 we could overload operator + and -, consequently +=, -=, even < and ==,
			 inherit the whole thing from an extended DPosition with similar
			 methods, initialize from iterator ranges, ... but there are not
			 concrete use cases at the moment.  Please contact the maintainer if you
			 need something for your application field. (Clemens)
		*/
		template < UInt D >
		class AveragePosition
		{
		 public:

			/// Dimensionality
			enum
				{
					DIMENSION = D
				};

			/// Position type (a D-dimensional position)
			typedef DPosition<DIMENSION> PositionType;

			/// Weight type (for weighted average - a scalar type)
			typedef typename PositionType::CoordinateType CoordinateType;

		 public:

			/// Default constructor
			AveragePosition ()
				: position_(),
					position_weighted_sum_(),
					weight_sum_()
			{
			}
			
			/// Copy constructor
			AveragePosition ( AveragePosition const & rhs )
				: position_               ( rhs.position_),
					position_weighted_sum_  ( rhs.position_weighted_sum_),
					weight_sum_             ( rhs.weight_sum_)
			{
			}

			/// Returns the current average position.
			PositionType const & getPosition() const
			{
				return position_;
			}

			/// Returns the total weight.
			CoordinateType const & getWeight() const
			{
				return weight_sum_;
			}

			/// Reset everything.  (Note that \c update() will cause a division by zero after that.)
			void clear ()
			{
				position_.clear();
				position_weighted_sum_.clear();
				weight_sum_ = 0;
				return;
			}

			/// Add a position.
			void add ( PositionType position, CoordinateType const weight = 1 )
			{
				
				weight_sum_ += weight;
				position *= weight;
				position_weighted_sum_ += position;
				
				// if the sum of weights is 0, set all coordinates to zero as well
				if (weight_sum_ == 0)
				{
					position_.clear();
				}
				else 
				{
					position_ = position_weighted_sum_;
					position_ /= weight_sum_;
				}
				return;
			}

		 protected:

			PositionType position_;
			PositionType position_weighted_sum_;
			CoordinateType weight_sum_;

		};

	} // namespace Math

} // namespace OpenMS

#endif //  OPENMS_MATH_STATISTICS_AVERAGEPOSITION_H
