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

#ifndef OPENMS_KERNEL_RAWDATAPOINT1D_H
#define OPENMS_KERNEL_RAWDATAPOINT1D_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

	/**
		@brief	1-dimensional raw data point.
	 
		This datastructure is intended for continuous data.
		If you want to handle picked data use Peak1D or PickedPeak1D.

		@ingroup Kernel
	*/
	class RawDataPoint1D 	
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
    /// Dimension
    enum { DIMENSION = 1 };
		/// Intensity type
		typedef DoubleReal IntensityType;
		/// Position type
		typedef DPosition<1> PositionType;
		/// Coordinate type
		typedef DoubleReal CoordinateType;
		//@}

		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor
		inline RawDataPoint1D() 
      : position_(), 
        intensity_(0) 
    {
    }
		/// Copy constructor
		inline RawDataPoint1D(const RawDataPoint1D& p) 
			: position_(p.position_),
        intensity_(p.intensity_)
		{
    }
		/**@brief Destructor

		  @note The destructor is non-virtual although many classes are derived from
		  RawDataPoint1D.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance - but RawDataPoint1Ds are
		  used in great amounts for storing raw data. Normally you should not derive other classes from
		  RawDataPoint1D (unless you know what you are doing, of course).
		*/
		~RawDataPoint1D()
    {
    }
		//@}
		
		/**	
			@name Accessors
		*/
		//@{
		/// Non-mutable access to the data point intensity (height)
		inline IntensityType getIntensity() const { return intensity_; }
		/// Mutable access to the data point intensity (height)
		inline void setIntensity(IntensityType intensity) { intensity_ = intensity; }

		/// Non-mutable access to m/z
		inline CoordinateType getMZ() const
		{
			return position_[0];
		}
		/// Mutable access to m/z
		inline void setMZ(CoordinateType mz)
		{
			position_[0] = mz;
		}

		/// Alias for getMZ()
		inline CoordinateType getPos() const
		{
			return position_[0];
		}
		/// Alias for setMZ()
		inline void setPos(CoordinateType pos)
		{
			position_[0] = pos;
		}

    /// Non-mutable access to the position
    inline PositionType const & getPosition() const
    {
      return position_;
    }
    /// Mutable access to the position
    inline PositionType & getPosition()
    {
      return position_;
    }
    /// Mutable access to the position
    inline void setPosition(PositionType const& position)
    {
      position_ = position;
    }
    //@}

		/// Assignment operator
		inline RawDataPoint1D& operator = (const RawDataPoint1D& rhs)
		{
			if (this==&rhs) return *this;
		
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
	
			return *this;
		}
		
		/// Equality operator
		inline bool operator == (const RawDataPoint1D& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		inline bool operator != (const RawDataPoint1D& rhs) const
		{
			return !( operator==(rhs) );
		}

									
 		/**	@name	Comparator classes.
				These classes implement binary predicates that can be used 
				to compare two peaks with respect to their intensities, positions.
				They are employed by the sort methods in container classes such as PeakArray.
		*/
		//@{

		/// Compare by getIntensity()
		struct IntensityLess
			: std::binary_function < RawDataPoint1D, RawDataPoint1D, bool >
		{
			inline bool operator () ( RawDataPoint1D const & left, RawDataPoint1D const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( RawDataPoint1D const & left, IntensityType const & right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType const & left, RawDataPoint1D const & right ) const
			{
				return ( left< right.getIntensity() );
			}
			inline bool operator () ( IntensityType const & left, IntensityType const & right ) const
			{
				return ( left < right );
			}
		};

		///Comparator for the position.
		struct PositionLess
			: public std::binary_function <RawDataPoint1D, RawDataPoint1D, bool>
		{
			inline bool operator () (const RawDataPoint1D& a, const RawDataPoint1D& b) const
			{
				return (a.getMZ() < b.getMZ());
			}
		
			/// comparison of a RawDataPoint2D with a CoordinateType
			inline bool operator () ( RawDataPoint1D const & left, CoordinateType right ) const throw()
			{
				return (left.getMZ() < right );
			}
			
			/// comparison of a CoordinateType with a RawDataPoint2D
			inline bool operator () ( CoordinateType left, RawDataPoint1D const & right ) const throw()
			{
				return (left < right.getMZ() );
			}

			/**
				@brief Operator to check if comparison is done increasing or decreasing.
				
				Sometimes we need a way to find out which way the CoordinateType is
				sorted and adding this overload seems to be the best way to achieve that goal.
			*/
			inline bool operator () ( CoordinateType left, CoordinateType right ) const throw()
			{
				return (left < right );
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
	std::ostream& operator << (std::ostream& os, const RawDataPoint1D& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_RAWDATAPOINT1D_H
