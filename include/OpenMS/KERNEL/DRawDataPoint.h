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

		@ingroup Kernel, Serialization
	*/
	template <Size D, typename Traits = KernelTraits>
	class DRawDataPoint 	
	{
	 public:
		
		/** @name Type definitions
		 */
		//@{
		/// Number of dimenstions
		enum { DIMENSION = D };
		/// Traits types
		typedef Traits TraitsType;
		/// Intensity type
		typedef typename Traits::IntensityType IntensityType;
		/// Coordinate type (of the position)
		typedef typename Traits::CoordinateType CoordinateType;
		/// Position type
		typedef DPosition<D, Traits> PositionType;
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
		const PositionType& getPosition() const { return position_; }
		/// Mutable access to the data point position (multidimensional)
		PositionType& getPosition() { return position_; }
		/// Mutable access to the data point position (multidimensional)
		void setPosition(PositionType const& position) { position_ = position; }

		/// Non-mutable access to the zero-th data point position (D==1 only)
		CoordinateType const & getPos() const
		{
			// static int dimension_must_be_one_ [2-DIMENSION];
			return position_[0];
		}
		/// Mutable access to the zero-th data point position (D==1 only)
		CoordinateType & getPos()
		{
			// static int dimension_must_be_one_  [2-DIMENSION];
			return position_[0];
		}
		/// Mutable access to the zero-th data point position (D==1 only)
		void setPos(CoordinateType const& coordinate)
		{
			// static int dimension_must_be_one_ [2-DIMENSION];
			position_[0] = coordinate;
		}

		/// Non-mutable access to the i-th data point dimension
		CoordinateType const & getPos(Size const i) const { return position_[i]; }
		/// Mutable access to the i-th data point dimension
		CoordinateType& getPos(Size const i) { return position_[i]; }
		/// Mutable access to the i-th data point dimension
		void setPos(Size const i, const CoordinateType& coordinate) { position_[i] = coordinate; }

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
				return (left.getPosition()[i] < right.getPosition()[i]);
			}
			
			/// comparison of a DRawDataPoint with a CoordinateType
			inline bool operator () ( DRawDataPoint const & left, CoordinateType const & right ) const throw()
			{
				return (left.getPosition()[i] < right );
			}
			
			/// comparison of a CoordinateType with a DRawDataPoint
			inline bool operator () ( CoordinateType const & left, DRawDataPoint const & right ) const throw()
			{
				return (left < right.getPosition()[i] );
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
				return (a.getPosition() < b.getPosition());
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
	template <Size D, typename Traits>
	std::ostream& operator << (std::ostream& os, const DRawDataPoint<D, Traits>& point)
	{
		os << "POS: "<< point.getPosition() << " INT: "<<point.getIntensity();
		
		return os;
	}

} // namespace OpenMS

#endif // OPENMS_KERNEL_DRAWDATAPOINT_H
