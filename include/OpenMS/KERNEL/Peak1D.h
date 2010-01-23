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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_PEAK1D_H
#define OPENMS_KERNEL_PEAK1D_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

	/**
		@brief A 1-dimensional raw data point or peak.
	 
		This datastructure is intended for continuous data or peak data.
		If wou want to annotated single peaks with meta data, use RichPeak1D instead.

		@ingroup Kernel
	*/
	class OPENMS_DLLAPI Peak1D 	
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
    /// Dimension
    enum { DIMENSION = 1 };
		/// Intensity type
		typedef Real IntensityType;
		/// Position type
		typedef DPosition<1> PositionType;
		/// Coordinate type
		typedef DoubleReal CoordinateType;
		//@}

		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor
		inline Peak1D() 
      : position_(), 
        intensity_(0) 
    {
    }
		/// Copy constructor
		inline Peak1D(const Peak1D& p) 
			: position_(p.position_),
        intensity_(p.intensity_)
		{
    }
		/**@brief Destructor

		  @note The destructor is non-virtual although many classes are derived from
		  Peak1D.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance. Normally you should not derive other classes from
		  Peak1D (unless you know what you are doing, of course).
		*/
		~Peak1D()
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
		inline Peak1D& operator = (const Peak1D& rhs)
		{
			if (this==&rhs) return *this;
		
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
	
			return *this;
		}
		
		/// Equality operator
		inline bool operator == (const Peak1D& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		inline bool operator != (const Peak1D& rhs) const
		{
			return !( operator==(rhs) );
		}

									
 		/**	@name	Comparator classes.
				These classes implement binary predicates that can be used 
				to compare two peaks with respect to their intensities, positions.
		*/
		//@{
		/// Comparator by intensity
		struct IntensityLess
			: std::binary_function < Peak1D, Peak1D, bool >
		{
			inline bool operator () ( Peak1D const & left, Peak1D const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( Peak1D const & left, IntensityType right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType left, Peak1D const & right ) const
			{
				return ( left< right.getIntensity() );
			}
			inline bool operator () ( IntensityType left, IntensityType right ) const
			{
				return ( left < right );
			}
		};
		///Comparator by m/z position.
		struct MZLess
			: public std::binary_function <Peak1D, Peak1D, bool>
		{
			inline bool operator () (const Peak1D& left, const Peak1D& right) const
			{
				return (left.getMZ() < right.getPos());
			}
			inline bool operator () ( Peak1D const & left, CoordinateType right ) const 
			{
				return (left.getMZ() < right );
			}
			inline bool operator () ( CoordinateType left, Peak1D const & right ) const 
			{
				return (left < right.getMZ() );
			}
			inline bool operator () ( CoordinateType left, CoordinateType right ) const 
			{
				return (left < right );
			}
		};
		/// Comparator by position. As this class has dimension 1, this is basically an alias for MZLess.
		struct PositionLess
			: public std::binary_function <Peak1D, Peak1D, bool>
		{
			inline bool operator () ( const Peak1D& left, const Peak1D& right) const
			{
				return (left.getPosition() < right.getPosition());
			}
			inline bool operator () ( const Peak1D& left, const PositionType& right ) const 
			{
				return (left.getPosition() < right );
			}
			inline bool operator () ( const PositionType& left, const Peak1D& right ) const 
			{
				return (left < right.getPosition() );
			}
			inline bool operator () ( const PositionType& left, const PositionType& right ) const 
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
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const Peak1D& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_PEAK1D_H
