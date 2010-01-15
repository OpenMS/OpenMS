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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_CHROMATOGRAMPEAK_H
#define OPENMS_KERNEL_CHROMATOGRAMPEAK_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

	/**
		@brief A 1-dimensional raw data point or peak for chromatograms.
	 
		This datastructure is intended for chromatograms.
		If wou want to annotated single peaks with meta data, use RichChromatogramPeak instead.

		@ingroup Kernel
	*/
	class OPENMS_DLLAPI ChromatogramPeak 	
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
		inline ChromatogramPeak() 
      : position_(), 
        intensity_(0) 
    {
    }
		/// Copy constructor
		inline ChromatogramPeak(const ChromatogramPeak& p) 
			: position_(p.position_),
        intensity_(p.intensity_)
		{
    }
		/**@brief Destructor

		  @note The destructor is non-virtual although many classes are derived from
		  ChromatogramPeak.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance. Normally you should not derive other classes from
		  ChromatogramPeak (unless you know what you are doing, of course).
		*/
		~ChromatogramPeak()
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
		inline CoordinateType getRT() const
		{
			return position_[0];
		}
		/// Mutable access to m/z
		inline void setRT(CoordinateType rt)
		{
			position_[0] = rt;
		}

		/// Alias for getRT()
		inline CoordinateType getPos() const
		{
			return position_[0];
		}
		/// Alias for setRT()
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
		inline ChromatogramPeak& operator = (const ChromatogramPeak& rhs)
		{
			if (this==&rhs) return *this;
		
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
	
			return *this;
		}
		
		/// Equality operator
		inline bool operator == (const ChromatogramPeak& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		inline bool operator != (const ChromatogramPeak& rhs) const
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
			: std::binary_function < ChromatogramPeak, ChromatogramPeak, bool >
		{
			inline bool operator () ( ChromatogramPeak const & left, ChromatogramPeak const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( ChromatogramPeak const & left, IntensityType right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType left, ChromatogramPeak const & right ) const
			{
				return ( left< right.getIntensity() );
			}
			inline bool operator () ( IntensityType left, IntensityType right ) const
			{
				return ( left < right );
			}
		};
		///Comparator by RT position.
		struct RTLess
			: public std::binary_function <ChromatogramPeak, ChromatogramPeak, bool>
		{
			inline bool operator () (const ChromatogramPeak& left, const ChromatogramPeak& right) const
			{
				return (left.getRT() < right.getPos());
			}
			inline bool operator () ( ChromatogramPeak const & left, CoordinateType right ) const 
			{
				return (left.getRT() < right );
			}
			inline bool operator () ( CoordinateType left, ChromatogramPeak const & right ) const 
			{
				return (left < right.getRT() );
			}
			inline bool operator () ( CoordinateType left, CoordinateType right ) const 
			{
				return (left < right );
			}
		};
		/// Comparator by position. As this class has dimension 1, this is basically an alias for RTLess.
		struct PositionLess
			: public std::binary_function <ChromatogramPeak, ChromatogramPeak, bool>
		{
			inline bool operator () ( const ChromatogramPeak& left, const ChromatogramPeak& right) const
			{
				return (left.getPosition() < right.getPosition());
			}
			inline bool operator () ( const ChromatogramPeak& left, const PositionType& right ) const 
			{
				return (left.getPosition() < right );
			}
			inline bool operator () ( const PositionType& left, const ChromatogramPeak& right ) const 
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
	OPENMS_DLLAPI std::ostream& operator << (std::ostream& os, const ChromatogramPeak& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_CHROMATOGRAMPEAK_H
