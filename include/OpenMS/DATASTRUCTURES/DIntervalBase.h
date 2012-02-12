// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DINTERVALBASE_H
#define OPENMS_DATASTRUCTURES_DINTERVALBASE_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <utility> // std::pair

namespace OpenMS
{
	namespace Internal
	{
		/**
			@brief A base class for D-dimensional interval.

			See DIntervalBase for a closed interval and DRange for a half-open interval class.

			@invariant All methods maintain the invariant that minPosition() is geometrically less or equal maxPosition() 
			           i.e. minPosition()[x] <= maxPosition()[x].
		*/
		template <UInt D>
		class DIntervalBase
		{
			public:
			
				/**	
					@name Type definitions 
				*/	
				//@{
				/// Dimensions
				enum { DIMENSION = D };
				/// Position type
				typedef DPosition<D>							PositionType;
				/// Coordinate type of the positions
				typedef typename PositionType::CoordinateType CoordinateType;
				//@}
		
				/**	@name Constructors and Destructor */
				//@{
				
				/** 
					@brief Default constructor.  
					
					Creates an empty interval with corners at infinity.
				*/
				DIntervalBase()
					: min_(PositionType::maxPositive()),
						max_(PositionType::minNegative())
				{
				}
				
				/// Copy constructor
				DIntervalBase(const DIntervalBase& rhs)
					: min_(rhs.min_),
						max_(rhs.max_)
				{
				}
				
				/// Assignment operator
				DIntervalBase & operator=(const DIntervalBase & rhs)
				{
					min_ = rhs.min_;
					max_ = rhs.max_;
					return *this;
				}
				
				/// Destructor
				~DIntervalBase() 
				{
				}
		
				/**
					@brief This constructor sets min_ and max_ directly.
				*/
				DIntervalBase( PositionType const & minimum, PositionType const & maximum )
					: min_(minimum),
						max_(maximum)
				{
					normalize_();
				}
		
				//@}

				/**	@name Accessors */
				//@{
				
				/// Accessor to minimum position
				inline PositionType const & minPosition() const
				{
					return min_;
				}
		
				/// Accessor to maximum position
				inline PositionType const & maxPosition() const
				{
					return max_;
				}
				
				/**
					@brief Mutator for minimum position
				
					@note The minimum position given here will be returned my minPosition() after the method.
					      If necessary the value returned by maxPosition() will be adjusted.					
				*/
				void setMin(PositionType const & position)
				{
					min_ = position;
					for ( UInt i = 0; i < DIMENSION; ++i )
					{
						if (min_[i]>max_[i]) max_[i] = min_[i];
					}
				}
	
				/**
					@brief Mutator for maximum position
				
					@note The maximum position given here will be returned my maxPosition() after the method.
					      If necessary the value returned by minPosition() will be adjusted.
				*/
				void setMax(PositionType const & position)
				{
					max_ = position;
					for ( UInt i = 0; i < DIMENSION; ++i )
					{
						if (min_[i]>max_[i]) min_[i] = max_[i];
					}
				}

				/**
					@brief Mutator for minimum and maximum position
				*/
				void setMinMax(PositionType const & min, PositionType const & max)
				{
					min_ = min;
					max_ = max;
					normalize_();
				}
				
				/**
					@brief Assignment from a DIntervalBase of different dimensions.
					
					Only the dimensions 0 upto min(D,D2)-1 are copied.
				*/
				template <UInt D2>
				void assign(const DIntervalBase<D2> rhs)
				{
					for (UInt i=0; i<std::min(D,D2); ++i)
					{
						min_[i] = rhs.minPosition()[i];
						max_[i] = rhs.maxPosition()[i];
					}
				}
				
				//}@
		
				/**	@name Predicates */
				//@{
				/// Equality operator
				bool operator == (const DIntervalBase& rhs) const 
				{
					return (min_ == rhs.min_) && (max_ == rhs.max_);
				}
				
				/// Equality operator
				bool operator != (const DIntervalBase& rhs) const 
				{
					return !(operator==(rhs));
				}
		    
				/// Make the interval empty
				inline void clear() 
				{ 
					*this = empty; 
				}
				
				//@}
	
				/**	@name Misc */
				//@{
	
				///Returns the center of the interval
				PositionType center() const
				{
					PositionType center(min_);
					center += max_;
					center /= 2;
					return center;
				}

				/// Returns the diagonal of the area, i.e. max_ - min_.
				PositionType diagonal() const
				{
					return max_ - min_;
				}
	
				/// empty instance
				static DIntervalBase const empty;
				/// instance with all positions zero
				static DIntervalBase const zero;
		
				//}@

				/**	@name Accessors for 2D-intervals (for convenience) */
				//@{
		
				/// Accessor for min_ coordinate minimum
		    inline CoordinateType minX() const
		    {
		      return min_[0];  
		    }
				
				/// Accessor for max_ coordinate minimum
		    inline CoordinateType minY() const
		    {
					return min_[1];
		    }
		
				/// Accessor for min_ coordinate maximum
				inline CoordinateType maxX() const
		    {
		      return max_[0];
		    }
		
				/// Accessor for max_ coordinate maximum
				inline CoordinateType maxY() const
		    {
					return max_[1];
		    }	
		
				/// Mutator for min_ coordinate of the smaller point 
				void setMinX(CoordinateType const c)
		    {
		      min_[0] = c;
					if (min_[0]>max_[0]) max_[0] = min_[0];
		    }
				
				/// Mutator for max_ coordinate of the smaller point 
				void setMinY(CoordinateType const c)
		    {
		      min_[1] = c;
					if (min_[1]>max_[1]) max_[1] = min_[1];
				}
		
				/// Mutator for min_ coordinate of the larger point.
				void setMaxX(CoordinateType const c)
		    {
		      max_[0] = c;
					if (min_[0]>max_[0]) min_[0] = max_[0];
		    }
		
				/// Mutator for max_ coordinate of the larger point.
				void setMaxY(CoordinateType const c)
		    {
		      max_[1] = c;
					if (min_[1]>max_[1]) min_[1] = max_[1];
		    }
		
				/// Returns the width of the area i.e. the difference of dimension zero (X).
				inline CoordinateType width() const
				{
					return max_[0] - min_[0];
				}
				
				/// Returns the height of the area i.e. the difference of dimension one (Y).
				inline CoordinateType height() const
				{
					return max_[1] - min_[1];
				}
		
				//@}
	
			protected:
				
				/// lower left point
				PositionType min_; 
				
				/// upper right point
				PositionType max_; 
				
				/// normalization to keep all dimensions in the right geometrical order (min_[X] < max_[X])
				void normalize_()
				{
					for ( UInt i = 0; i < DIMENSION; ++i )
					{
						if (min_[i]>max_[i]) 
						{
							std::swap(min_[i],max_[i]);
						}
					}
				}

				///Protected constructor for the construction of static instances
				DIntervalBase(const std::pair< PositionType , PositionType >& pair )
					: min_(pair.first),
						max_(pair.second)
				{
					
				}
		
			};
		
			template <UInt D>
			DIntervalBase<D> const DIntervalBase<D>::zero
			= DIntervalBase<D>( DIntervalBase<D>::PositionType::zero(), DIntervalBase<D>::PositionType::zero());
			
			template <UInt D>
			DIntervalBase<D> const DIntervalBase<D>::empty 
			= DIntervalBase<D>(std::make_pair(DIntervalBase<D>::PositionType::maxPositive(), DIntervalBase<D>::PositionType::minNegative()));
			
			///Print the contents to a stream.
			template <UInt D>
			std::ostream& operator << (std::ostream& os, const DIntervalBase<D>& rhs)
			{
				os << "--DIntervalBase BEGIN--"<<std::endl;
				os << "MIN --> " << rhs.minPosition() << std::endl;
				os << "MAX --> " << rhs.maxPosition() << std::endl;
				os << "--DIntervalBase END--"<<std::endl;
				return os;
			}
	
	} // namespace Internal
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DINTERVALBASE_H
