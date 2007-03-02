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

#ifndef OPENMS_KERNEL_RAWDATAPOINT2D_H
#define OPENMS_KERNEL_RAWDATAPOINT2D_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <ostream>
#include <functional>

namespace OpenMS
{

	/**
		@brief	2-dimensional raw data point.
	 
		This datastructure is intended for continuous data.
		If you want to handle picked data use Peak2D or PickedPeak2D.

		@ingroup Kernel, Serialization
	*/
	class RawDataPoint2D
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
    
    /// Intensity type
    typedef DoubleReal IntensityType;
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
    static char const * const shortDimensionName(UnsignedInt const dim);
    /// Short name of the dimension (abbreviated form)
    static char const * const shortDimensionNameRT();
    /// Short name of the dimension (abbreviated form)
    static char const * const shortDimensionNameMZ();

    /// Full name of the dimension (self-explanatory form)
    static char const * const fullDimensionName(UnsignedInt const dim);
    /// Full name of the dimension (self-explanatory form)
    static char const * const fullDimensionNameRT();
    /// Full name of the dimension (self-explanatory form)
    static char const * const fullDimensionNameMZ();

    /// Unit of measurement (abbreviated form)
		static char const * const shortDimensionUnit(UnsignedInt const dim);
    /// Unit of measurement (abbreviated form)
		static char const * const shortDimensionUnitRT();
    /// Unit of measurement (abbreviated form)
		static char const * const shortDimensionUnitMZ();

    /// Unit of measurement (self-explanatory form)
		static char const * const fullDimensionUnit(UnsignedInt const dim);
    /// Unit of measurement (self-explanatory form)
		static char const * const fullDimensionUnitRT();
    /// Unit of measurement (self-explanatory form)
		static char const * const fullDimensionUnitMZ();

	 protected:
    /// Short name of the dimension (abbreviated form)
		static char const * const dimension_name_short [DIMENSION];
		
    /// Full name of the dimension (self-explanatory form)
    static char const * const dimension_name_full  [DIMENSION];
		
    /// Unit of measurement (abbreviated form)
    static char const * const dimension_unit_short [DIMENSION];
		
    /// Unit of measurement (self-explanatory form)
    static char const * const dimension_unit_full  [DIMENSION];
    
		//@}

	 public:

		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor
		RawDataPoint2D() 
      : position_(), 
        intensity_(0) 
    {
    }
		/// Copy constructor
		RawDataPoint2D(const RawDataPoint2D& p) 
			: position_(p.position_), 
        intensity_(p.intensity_)
		{
    }
		
    /**@brief Destructor

		  @note The destructor is non-virtual although many classes are derived from
		  RawDataPoint2D.  This is intentional, since otherwise we would "waste"
		  space for a vtable pointer in each instance - but RawDataPoint2Ds are
		  used in great amounts for storing raw data. Normally you should not derive other classes from
		  RawDataPoint2D (unless you know what you are doing, of course).
		*/
		~RawDataPoint2D()
    {
    }
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

		/// Non-mutable access to the position
		PositionType const & getPos() const 
    { 
      return position_; 
    }
		/// Mutable access to the position
		PositionType& getPos() 
    {
      return position_; 
    }
		/// Mutable access to the position
		void setPos(const PositionType& position) 
    { 
      position_ = position; 
    }

    /// Returns the m/z coordinate (index 1)
    CoordinateType const & getMZ() const 
    { 
      return position_[1]; 
    }
    /// Mutable access to the m/z coordinate (index 1)
    void setMZ(const CoordinateType& coordinate) 
    { 
      position_[1] = coordinate; 
    }

    /// Returns the RT coordinate (index 0)
    CoordinateType const & getRT() const 
    { 
      return position_[0]; 
    }
    /// Mutable access to the RT coordinate (index 0)
    void setRT(const CoordinateType& coordinate) 
    { 
      position_[0] = coordinate; 
    }
    
		//@}

		/// Assignment operator
		RawDataPoint2D& operator = (const RawDataPoint2D& rhs)
		{
			if (this==&rhs) return *this;
			
			intensity_ = rhs.intensity_;
			position_ = rhs.position_;
		
			return *this;
		}
		
		/// Equality operator
		bool operator == (const RawDataPoint2D& rhs) const
		{
			return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
		}

		/// Equality operator
		bool operator != (const RawDataPoint2D& rhs) const
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
			: std::binary_function < RawDataPoint2D, RawDataPoint2D, bool >
		{
			inline bool operator () ( RawDataPoint2D const & left, RawDataPoint2D const & right ) const
			{
				return ( left.getIntensity() < right.getIntensity() );
			}
			inline bool operator () ( RawDataPoint2D const & left, IntensityType const & right ) const
			{
				return ( left.getIntensity() < right );
			}
			inline bool operator () ( IntensityType const & left, RawDataPoint2D const & right ) const
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
		template <UnsignedInt i>
		struct NthPositionLess
			: std::binary_function <RawDataPoint2D, RawDataPoint2D, bool>
		{
			enum { DIMENSION = i };
			
			/// comparison of two RawDataPoint2Ds
			inline bool operator () ( RawDataPoint2D const & left, RawDataPoint2D const & right ) const throw()
			{
				return (left.getPos()[i] < right.getPos()[i]);
			}
			
			/// comparison of a RawDataPoint2D with a CoordinateType
			inline bool operator () ( RawDataPoint2D const & left, CoordinateType const & right ) const throw()
			{
				return (left.getPos()[i] < right );
			}
			
			/// comparison of a CoordinateType with a RawDataPoint2D
			inline bool operator () ( CoordinateType const & left, RawDataPoint2D const & right ) const throw()
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

		/// Comparator with respect to retention time
		typedef NthPositionLess < RT > LessRT;
		
		/// Comparator with respect to mass-to-charge
		typedef NthPositionLess < MZ > LessMZ;
		

		/**
			@brief Comparator for the position.
			
			Lexicographical comparison from dimension 0 to dimension D-1 is done.
		*/
		struct PositionLess
			: public std::binary_function <RawDataPoint2D, RawDataPoint2D, bool>
		{
			inline bool operator () (const RawDataPoint2D& a, const RawDataPoint2D& b) const
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


		///@name Serialization
		//@{
	 private:
		/// Serialization interface
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* version */ )
		{
			ar & boost::serialization::make_nvp("pos",this->position_);
			ar & boost::serialization::make_nvp("it",this->intensity_);
		}
		//@}

		/// Serialization
		friend class boost::serialization::access;

	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const RawDataPoint2D& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_DRAWDATAPOINT_H
