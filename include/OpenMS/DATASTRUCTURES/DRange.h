// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DRANGE_H
#define OPENMS_DATASTRUCTURES_DRANGE_H

#include <OpenMS/DATASTRUCTURES/DIntervalBase.h>

namespace OpenMS
{
	/**	
		@brief A D-dimensional half-open interval.
	
		This class describes a range in D-dimensional space delimited
		by two points (i.e. a D-dimensional hyper-rectangle). The 
		two points define the lower left and the upper right corner
		in 2D and analogous points in higer dimensions.
		
		A range is a pair of positions in D-space represented by DPosition.
		The two limiting points are accessed as min() and max().
		
		A range denotes a semi-open interval. A lower coordinate of each
		dimension is part the range, the higher coordinate is not.
		
		@ingroup Datastructures
	*/
	template <UInt D>
	class DRange
		:	public Internal::DIntervalBase<D>
	{
	 public:

		/**	
			@name Type definitions 
		*/	
		//@{
		/// Dimensions
		enum { DIMENSION = D };
		/// Base class type
		typedef Internal::DIntervalBase<D> Base; 
		/// Position type
		typedef typename Base::PositionType PositionType;
		/// Coordinate type of the positions
		typedef typename Base::CoordinateType CoordinateType;
		///Types that describe the kind of intersection between two ranges
		enum DRangeIntersection 
		{
			Disjoint,   ///< No intersection
			Intersects, ///< Intersection
			Inside      ///< One contains the other
		};

		//@}
		
		using Base::min_;
		using Base::max_;

		/**	@name Constructors and Destructor */
		//@{
		/** 
			@brief Default constructor.
			
			Creates a range with all coordinates zero.
		*/
		DRange()
			: Base()
		{
		}
		
		/// Constructor that takes two Points and constructs a range.
		DRange(const PositionType& lower, const PositionType& upper)
			:Base(lower, upper)
		{
		}
		
		/// Copy constructor.
		DRange(const DRange& range)
			: Base(range)
		{
		}

		/// Copy constructor for the base class
		DRange(const Base& range)
			: Base(range)
		{
		}

		///Convenient constructor for DRange<2>
		DRange(CoordinateType minx,CoordinateType miny, CoordinateType maxx,CoordinateType maxy)
		{
			OPENMS_PRECONDITION(D == 2, "DRange<D>:DRange(minx, miny, maxx, maxy): index overflow!");
			min_[0]=minx;
			min_[1]=miny;
			max_[0]=maxx;
			max_[1]=maxy;
		}

		/// Assignement operator
		DRange & operator=(const DRange& rhs)
		{
			Base::operator=(rhs);
			return *this;
		}

		/// Assignement operator for the base class
		DRange & operator=(const Base& rhs)
		{
			Base::operator=(rhs);
			return *this;
		}

		/// Destuctor 
		~DRange()
		{
		}
		//@}

		/**	@name Predicates */
		//@{
		///Equality operator
		bool operator == (const DRange& rhs) const 
		{
			return Base::operator==(rhs);
		}
		/// Equality operator
		bool operator == (const Base& rhs) const 
		{
			return Base::operator==(rhs);
		}

		/**
			 @brief Checks whether this range contains a certain point.
			 
			 @param position The point's position.
			 @returns true if point lies inside this area.
		*/
		bool encloses(const PositionType& position) const
		{
      for (UInt i = 0; i != D; i++)
      {
        if (position[i]<min_[i]) return false;
        if (position[i]>=max_[i]) return false;
      }
			return true;
		}

		///@brief 2D-version of encloses for convenience only
		bool encloses(CoordinateType x, CoordinateType y) const
		{
      if (x<min_[0]) return false;
      if (x>=max_[0]) return false;
      if (y<min_[1]) return false;
      if (y>=max_[1]) return false;
			return true;
		}
		
		/// Returns the smallest range containing this range and @p other_range
		DRange united(const DRange<D>& other_range) const
		{
			PositionType united_min;
			PositionType united_max;
			DRange<D> united_range = DRange<D>::empty;
			
			PositionType other_min = other_range.min();
			PositionType other_max = other_range.max();
			
			for (Size i = 0; i != D; ++i)
			{
				united_min[i] = min_[i] < other_min[i] ? min_[i] : other_min[i];
				united_max[i] = max_[i] > other_max[i] ? max_[i] : other_max[i];
			}
			united_range.setMinMax(united_min, united_max);
			
			return united_range;
		}
		
		/**
			 @brief Checks how this range intersects with another @p range.
			 
			 @param range The max_ range.
		*/
		DRangeIntersection intersects(const DRange& range) const
		{
		  //check if r.min_ is in this area
		  if (encloses(range.min_))
		  {
		  	//check if r.max_ in this area => Inside / Intersects
	      for (Size i = 0; i != D; i++)
	      {
	        if (range.max_[i]>max_[i])
	        {
	        	return Intersects;
	      	}
	      }
		  	return Inside;
		  }
		  // => r.min_ is not inside this area
	  	//check if any r.min_ >= max_ => Disjoint
      for (Size i = 0; i != D; i++)
      {
        if (range.min_[i]>=max_[i])
        {
        	return Disjoint;
      	}
      }
      // => some coordinate of r.min_ has to be smaller than the one of min_
  		//check if all coords of r are smaller than the those of the range
      for (Size i = 0; i != D; i++)
      {
        if (range.max_[i]<=min_[i])
        {
        	return Disjoint;
      	}
      }
  		return Intersects;
	  }         

		/**
			 @brief Checks whether this range intersects with another @p range.
			 
			 @param range The max_ range.
			 @returns True if the areas intersect (i.e. they intersect or one contains the other).
		*/
		bool isIntersected(const DRange& range) const
		{
		  //check if r.min_ is in this area
		  if (encloses(range.min_))
		  {
		  	return true;
		  }
		  
		  // => r.min_ is not inside this area
	  	//check if any r.min_ >= max_ => Disjoint
      for (Size i = 0; i != D; i++)
      {
        if (range.min_[i]>=max_[i])
        {
        	return false;
      	}
      }
      // => some coordinate of r.min_ has to be smaller than the one of min_
  		//check if all coords of r are smaller than the those of the range
      for (Size i = 0; i != D; i++)
      {
        if (range.max_[i]<=min_[i])
        {
        	return false;
      	}
      }
  		return true;
    }
    
		/// Checks if the range is empty
		bool isEmpty() const 
		{ 
      for (UInt i = 0; i != D; i++)
      {
        if (max_[i]<=min_[i])
        {
        	return true;
      	}
      }
			return false; 
		}
		
		//@}
	};

	///Print the contents to a stream.
	template <UInt D>
	std::ostream& operator << (std::ostream& os, const DRange<D>& area)
	{
		os << "--DRANGE BEGIN--"<<std::endl;
		os << "MIN --> " << area.min_ << std::endl;
		os << "MAX --> " << area.max_ << std::endl;
		os << "--DRANGE END--"<<std::endl;
		return os;
	}

	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DRANGE_H
