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

#ifndef OPENMS_KERNEL_DPEAK_H
#define OPENMS_KERNEL_DPEAK_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <sstream>

namespace OpenMS
{
	
	/**
		@brief	D-dimensional raw data point or peak.
		
		This datastructure is intended for continuous data or peak data.
		If wou want to annotated single peaks with meta data, use DRichPeak instead.
		
		@ingroup Kernel
	*/
	template <UInt D>
	class DPeak
	{
		public:
			
			///@name Type definitions
			//@{
			/// Dimensionality
			enum
			{
				DIMENSION = D
			};
			/// Intensity type
			typedef Real IntensityType;
			/// Coordinate type (of the position)
			typedef DoubleReal CoordinateType;
			/// Position type
			typedef DPosition<D> PositionType;
			//@}
			
			///@name Constructors and Destructor
			//@{
			/// Default constructor
			DPeak() : position_(), intensity_(0)
			{}
			/// Copy constructor
			DPeak(const DPeak& p)
			: position_(p.position_), intensity_(p.intensity_)
			{}
			/**@brief Destructor
			
			 @note The destructor is non-virtual although many classes are derived from
			 DPeak.  This is intentional, since otherwise we would "waste"
			 space for a vtable pointer in each instance. Normally you should not derive
			 other classes from DPeak (unless you know what you are doing, of course).
			 */
			~DPeak()
			{}
			//@}
			
			/**
			 @name Accessors
			 */
			//@{
			
			/// Non-mutable access to the data point intensity (height)
			IntensityType getIntensity() const
			{ return intensity_; }
			
			/// Non-mutable access to the data point intensity (height)
			void setIntensity(IntensityType intensity)
			{ intensity_ = intensity; }
			
			/// Non-mutable access to the data point position (multidimensional)
			const PositionType& getPosition() const
			{ return position_; }
			
			/// Mutable access to the data point position (multidimensional)
			PositionType& getPosition()
			{ return position_; }
			
			/// Mutable access to the data point position (multidimensional)
			void setPosition(PositionType const& position)
			{ position_ = position; }
			
			//@}
			
			/// Assignment operator
			DPeak& operator = (const DPeak& rhs)
			{
				if (this==&rhs) return *this;
				
				intensity_ = rhs.intensity_;
				position_ = rhs.position_;
				
				return *this;
			}
			
			/// Equality operator
			bool operator == (const DPeak& rhs) const
			{
				return  intensity_ == rhs.intensity_ && position_ == rhs.position_ ;
			}
			
			/// Equality operator
			bool operator != (const DPeak& rhs) const
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
			: std::binary_function < DPeak, DPeak, bool >
			{
				inline bool operator () ( DPeak const & left, DPeak const & right ) const
				{
					return ( left.getIntensity() < right.getIntensity() );
				}
				inline bool operator () ( DPeak const & left, IntensityType const & right ) const
				{
					return ( left.getIntensity() < right );
				}
				inline bool operator () ( IntensityType const & left, DPeak const & right ) const
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
			template <UInt i>
			struct NthPositionLess
			: std::binary_function <DPeak, DPeak, bool>
			{
				enum
				{ DIMENSION = i };
				
				/// comparison of two DPeaks
				inline bool operator () ( DPeak const & left, DPeak const & right ) const 
				{
					return (left.getPosition()[i] < right.getPosition()[i]);
				}
				
				/// comparison of a DPeak with a CoordinateType
				inline bool operator () ( DPeak const & left, CoordinateType right ) const 
				{
					return (left.getPosition()[i] < right );
				}
				
				/// comparison of a CoordinateType with a DPeak
				inline bool operator () ( CoordinateType left, DPeak const & right ) const 
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
			
			/**
			 @brief Comparator for the position.
			
			 Lexicographical comparison from dimension 0 to dimension D-1 is done.
			 */
			struct PositionLess
			: public std::binary_function <DPeak, DPeak, bool>
			{
				inline bool operator () (const DPeak& a, const DPeak& b) const
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
	template <UInt D>
	std::ostream& operator << (std::ostream& os, const DPeak<D>& point)
	{
		os << "POS: "<< point.getPosition() << " INT: "<<point.getIntensity();
		
		return os;
	}
	
} // namespace OpenMS

#endif // OPENMS_KERNEL_DPEAK_H
