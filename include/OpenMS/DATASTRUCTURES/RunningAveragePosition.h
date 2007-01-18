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

#ifndef OPENMS_DATASTRUCTURES_RUNNINGAVERAGEPOSITION_H
#define OPENMS_DATASTRUCTURES_RUNNINGAVERAGEPOSITION_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

namespace OpenMS
{

  /**
  	@brief A running average position.  Positions can be added and removed.
           Positions have weights.

		Template parameter Position_ will usually be an instance of
		DPosition, but any type that has begin(), end(), clear(),
		iterators, a typedef CoordinateType and an enum DIMENSION should work (maybe
		I forgot a few concept requirements here).

		@todo A lot of convenience methods could be added.  For
		example, we could overload operator + and -, consequently +=, -=, even <
		and ==, inherit the whole thing from an extended DPosition with similar
		methods, initialize from iterator ranges, and and and. (Clemens)
	*/
  template < typename Position_ >
  class RunningAveragePosition
  {
   public:

    typedef Position_ PositionType;
    typedef typename PositionType::CoordinateType CoordinateType;
		enum { DIMENSION = PositionType::DIMENSION };

  public:

		RunningAveragePosition ()
			: position_(), position_weight_sum_(), weight_sum_()
		{}

		RunningAveragePosition ( RunningAveragePosition const & _arg )
			: position_             (_arg.position_),
				position_weight_sum_  (_arg.position_weight_sum_),
				weight_sum_           (_arg.weight_sum_)
		{}

		/// Returns the current running average position. */
    inline PositionType const & getPosition() const throw() { return position_; }

		/// Returns the total weight. */
		inline CoordinateType const & getWeight() const throw() { return weight_sum_; }

		/** Dimensionality of the underlying position type, might be useful in a polymorphic context.*/
		static Size const getDimension() { return DIMENSION; }

		/** Reset everything.  (Note that \c update() will cause a division by zero after that.) */
		void clear ()
		{
			position_.clear();
			position_weight_sum_.clear();
			weight_sum_ = 0;
		}

		/** Add a position. */
		void add ( PositionType const & position, CoordinateType const weight = 1 )
		{
			for ( Size i = 0; i < DIMENSION; ++i )
			{
				position_weight_sum_[i] += position[i] * weight;
			}
			weight_sum_ += weight;
			update();
		}

		/**@brief Subtract a position.  <code>subtract ( pos, weight )</code> is
			 equivalent to <code>add ( pos, -weight )</code>, but might be faster.
		*/
		void subtract ( PositionType const & position, CoordinateType const weight = 1 )
		{
			for ( Size i = 0; i < DIMENSION; ++i )
			{
				position_weight_sum_[i] -= position[i] * weight;
			}
			weight_sum_ -= weight;
			update();
		}

   protected:

		/**
			@brief updates the current average
		   
		  If the sum of weights is zero, the average will be set to zero as well.
		*/
		void update ()
		{
						
			// if the sum of weights is 0, set all
			// coordinates to zero as well
			if (weight_sum_ == 0)
			{
				for ( Size i = 0; i < DIMENSION; ++i )
				{
					position_[i] = 0;
				}
			}
			else 
			{
				// if not, proceed as usual.
				for ( Size i = 0; i < DIMENSION; ++i )
				{
					position_[i] = position_weight_sum_[i] / weight_sum_;
				}
			}
		}

    PositionType position_;
    PositionType position_weight_sum_;
		CoordinateType weight_sum_;
  };

} // namespace OpenMS

#endif //  OPENMS_DATASTRUCTURES_RUNNINGAVERAGEPOSITION_H
