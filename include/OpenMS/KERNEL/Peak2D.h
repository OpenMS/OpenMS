// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_KERNEL_PEAK2D_H
#define OPENMS_KERNEL_PEAK2D_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

	/**
		@brief A 2-dimensional raw data point or peak.
	 
		This datastructure is intended for continuous data or peak data.
		If wou want to annotated single peaks with meta data, use RichPeak2D instead.

		@ingroup Kernel
	*/
	class Peak2D
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
    
    /// Intensity type
    typedef Real IntensityType;
    /// Coordinate type (of the position)
		typedef DoubleReal CoordinateType;
		/// Position type
		typedef DPosition<2> PositionType;
    //@}

		/// Dimensions
		//@{
		
		/// This enum maps the symbolic names of the dimensions to numbers
		enum DimensionDescription
			{
				RT = 0, ///< Mass-to-charge dimension id (0 if used as a const int)
				MZ = 1, ///< Retention time dimension id (1 if used as a const int)
				DIMENSION = 2 ///< Number of dimensions
			};
		
    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionName(UInt const dim);
    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionNameRT();
    /// Short name of the dimension (abbreviated form)
    static char const * shortDimensionNameMZ();

    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionName(UInt const dim);
    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionNameRT();
    /// Full name of the dimension (self-explanatory form)
    static char const * fullDimensionNameMZ();

    /// Unit of measurement (abbreviated form)
		static char const * shortDimensionUnit(UInt const dim);
    /// Unit of measurement (abbreviated form)
		static char const * shortDimensionUnitRT();
    /// Unit of measurement (abbreviated form)
		static char const * shortDimensionUnitMZ();

    /// Unit of measurement (self-explanatory form)
		static char const * fullDimensionUnit(UInt const dim);
    /// Unit of measurement (self-explanatory form)
		static char const * fullDimensionUnitRT();
    /// Unit of measurement (self-explanatory form)
		static char const * fullDimensionUnitMZ();

	 protected:
    /// Short name of the dimension (abbreviated form)
		static char const * const dimension_name_short_[DIMENSION];
		
    /// Full name of the dimension (self-explanatory form)
    static char const * const dimension_name_full_[DIMENSION];
		
    /// Unit of measurement (abbreviated form)
    static char const * const dimension_unit_short_[DIMENSION];
		
    /// Unit of measurement (self-explanatory form)
    static char const * const dimension_unit_full_[DIMENSION];
    
		//@}

	 public:

		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor
		inline Peak2D() 
      : position_(), 
        intensity_(0) 
    {
    }
		/// Copy constructor
		inline Peak2D(const Peak2D& p) 
			: position_(p.position_), 
        intensity_(p.intensity_)
		{
    }
		
    /**@brief Destructor

		  @note The destructor is non-virtual although many classes are derived from
		  Peak2D.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance. Normally you should not derive other classes from
		  Peak2D (unless you know what you are doing, of course).
		*/
		~Peak2D()
    {
    }
		//@}
		
		/**	
			@name Accessors
		 */
		//@{
		/// Non-mutable access to the data point intensity (height)
		inline IntensityType getIntensity() const { return intensity_; }
		/// Non-mutable access to the data point intensity (height)
		inline void setIntensity(IntensityType intensity) { intensity_ = intensity; }

		/// Non-mutable access to the position
		inline PositionType const & getPosition() const 
    { 
      return position_; 
    }
		/// Mutable access to the position
		inline PositionType& getPosition() 
    {
      return position_; 
    }
		/// Mutable access to the position
		inline void setPosition(const PositionType& position) 
    { 
      position_ = position; 
    }

    /// Returns the m/z coordinate (index 1)
    inline CoordinateType getMZ() const 
    { 
      return position_[1]; 
    }
    /// Mutable access to the m/z coordinate (index 1)
    inline void setMZ(CoordinateType coordinate) 
    { 
      position_[1] = coordinate; 
    }

    /// Returns the RT coordinate (index 0)
    inline CoordinateType getRT() const 
    { 
      return position_[0]; 
    }
    /// Mutable access to the RT coordinate (index 0)
    inline void setRT(CoordinateType coordinate) 
    { 
      position_[0] = coordinate; 
    }
    
		//@}

		/// Assignment operator
		inline Peak2D& operator = (const Peak2D& rhs)
		{
			if (this==&rhs) return *this;
		
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
	
			return *this;
		}
				
		/// Equality operator
		inline bool operator == (const Peak2D& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		inline bool operator != (const Peak2D& rhs) const
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
			: std::binary_function < Peak2D, Peak2D, bool >
		{
			inline bool operator () ( Peak2D const & left, Peak2D const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( Peak2D const & left, IntensityType const & right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType const & left, Peak2D const & right ) const
			{
				return ( left< right.getIntensity() );
			}
			inline bool operator () ( IntensityType const & left, IntensityType const & right ) const
			{
				return ( left < right );
			}
		};
		
		/**
			@brief Comparator for the n-th coordinate of the position.
		*/
		template <UInt i>
		struct NthPositionLess
			: std::binary_function <Peak2D, Peak2D, bool>
		{
			enum { DIMENSION = i };
			
			/// comparison of two Peak2Ds
			inline bool operator () ( Peak2D const & left, Peak2D const & right ) const 
			{
				return (left.getPosition()[i] < right.getPosition()[i]);
			}
			
			/// comparison of a Peak2D with a CoordinateType
			inline bool operator () ( Peak2D const & left, CoordinateType right ) const 
			{
				return (left.getPosition()[i] < right );
			}
			
			/// comparison of a CoordinateType with a Peak2D
			inline bool operator () ( CoordinateType left, Peak2D const & right ) const 
			{
				return (left < right.getPosition()[i] );
			}

			/**
				@brief Operator to check if comparison is done increasing or decreasing.
				
				Sometimes we need a way to find out which way the CoordinateType is
				sorted and adding this overload seems to be the best way to achieve that goal.
			*/
			inline bool operator () ( CoordinateType left, CoordinateType right ) const 
			{
				return (left < right );
			}

		};

		/// Comparator with respect to retention time
		typedef NthPositionLess < RT > LessRT;
		
		/// Comparator with respect to mass-to-charge
		typedef NthPositionLess < MZ > LessMZ;
		

		/**
			@brief Comparator for the position.
			
			Lexicographical comparison from dimension 0 to dimension D-1 is done.
		*/
		struct PositionLess
			: public std::binary_function <Peak2D, Peak2D, bool>
		{
			inline bool operator () (const Peak2D& a, const Peak2D& b) const
			{
				return (a.getPosition() < b.getPosition());
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
	std::ostream& operator << (std::ostream& os, const Peak2D& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_PEAK2D_H
