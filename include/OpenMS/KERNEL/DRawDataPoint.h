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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_DRAWDATAPOINT_H
#define OPENMS_KERNEL_DRAWDATAPOINT_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <sstream>

namespace OpenMS
{

	/**
		@brief	D-dimensional raw data point.
	 
		This datastructure is intended for continous data.
		If you want to handle picked data use DPeak or DPickedPeak.

		@ingroup Kernel
	*/
	template <Size D>
	class DRawDataPoint 	
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
		/// Number of dimenstions
		enum { DIMENSION = D };
		/// Intensity type
		typedef DoubleReal IntensityType;
		/// Coordinate type (of the position)
		typedef DoubleReal CoordinateType;
		/// Position type
		typedef DPosition<D> PositionType;
		//@}

		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor
		DRawDataPoint() : position_(), intensity_(0) {}
		/// Copy constructor
		DRawDataPoint(const DRawDataPoint& p) 
			: position_(p.position_), intensity_(p.intensity_)
		{}
		/**@brief Destructor

		\note The destructor is non-virtual although many classes are derived from
		  DRawDataPoint.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance - but DRawDataPoints are
		  used in great amounts for storing raw data. Normally you should not derive other classes from
		  DRawDataPoint (unless you know what you are doing, of course).
		*/
		~DRawDataPoint() {}
		//@}
		
		/**	
			@name Accessors
		 */
		//@{
		/// Non-mutable access to the data point intensity (height)
		const IntensityType& getIntensity() const { return intensity_; }
		/// Mutable access to the data point intensity (height)
		IntensityType& getIntensity() { return intensity_; }
		/// Non-mutable access to the data point intensity (height)
		void setIntensity(const IntensityType& intensity) { intensity_ = intensity; }

		/// Non-mutable access to the data point position (multidimensional)
		const PositionType& getPos() const { return position_; }
		/// Mutable access to the data point position (multidimensional)
		PositionType& getPos() { return position_; }
		/// Mutable access to the data point position (multidimensional)
		void setPos(PositionType const& position) { position_ = position; }

		//@}

		/// Assignment operator
		DRawDataPoint& operator = (const DRawDataPoint& rhs)
		{
			if (this==&rhs) return *this;
			
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
		
			return *this;
		}
		
		/// Equality operator
		bool operator == (const DRawDataPoint& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		bool operator != (const DRawDataPoint& rhs) const
		{
			return !( operator==(rhs) );
		}

									
 		/**	@name	Comparator classes.
				These classes implement binary predicates that can be used 
				to compare two peaks with respect to their intensities, positions.
				They are employed by the sort methods in container classes such as DPeakArray.
		*/
		//@{

		/// Compare by getIntensity()
		struct IntensityLess
			: std::binary_function < DRawDataPoint, DRawDataPoint, bool >
		{
			inline bool operator () ( DRawDataPoint const & left, DRawDataPoint const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( DRawDataPoint const & left, IntensityType const & right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType const & left, DRawDataPoint const & right ) const
			{
				return ( left< right.getIntensity() );
			}
			inline bool operator () ( IntensityType const & left, IntensityType const & right ) const
			{
				return ( left < right );
			}
		};
		
		/**
			@brief Comparator for the i-th coordinate of the position.
		*/
		template <UnsignedInt i>
		struct NthPositionLess
			: std::binary_function <DRawDataPoint, DRawDataPoint, bool>
		{
			enum { DIMENSION = i };
			
			/// comparison of two DRawDataPoints
			inline bool operator () ( DRawDataPoint const & left, DRawDataPoint const & right ) const throw()
			{
				return (left.getPos()[i] < right.getPos()[i]);
			}
			
			/// comparison of a DRawDataPoint with a CoordinateType
			inline bool operator () ( DRawDataPoint const & left, CoordinateType const & right ) const throw()
			{
				return (left.getPos()[i] < right );
			}
			
			/// comparison of a CoordinateType with a DRawDataPoint
			inline bool operator () ( CoordinateType const & left, DRawDataPoint const & right ) const throw()
			{
				return (left < right.getPos()[i] );
			}

			/**
				@brief Operator to check if comparison is done increasing or decreasing.
				
				Sometimes we need a way to find out which way the CoordinateType is
				sorted and adding this overload seems to be the best way to achieve that goal.
			*/
			inline bool operator () ( CoordinateType const & left, CoordinateType const & right ) const throw()
			{
				return (left < right );
			}

		};

		/**
			@brief Comparator for the position.
			
			Lexicographical comparison from dimension 0 to dimension D-1 is done.
		*/
		struct PositionLess
			: public std::binary_function <DRawDataPoint, DRawDataPoint, bool>
		{
			inline bool operator () (const DRawDataPoint& a, const DRawDataPoint& b) const
			{
				return (a.getPos() < b.getPos());
			}
		};
		
		//@}
		
		protected:
		/// The data point position
		PositionType	position_;
		/// The data point intensity
		IntensityType intensity_;
	};

	///Print the contents to a stream.
	template <Size D>
	std::ostream& operator << (std::ostream& os, const DRawDataPoint<D>& point)
	{
		os << "POS: "<< point.getPos() << " INT: "<<point.getIntensity();
		
		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_DRAWDATAPOINT_H
