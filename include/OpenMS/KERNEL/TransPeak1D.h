// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_TRANSPEAK1D_H
#define OPENMS_KERNEL_TRANSPEAK1D_H

#include <OpenMS/KERNEL/Peak1D.h>

namespace OpenMS
{

	/**
		@brief A 1-dimensional raw data point or peak that stores a transformed intensity additionally.
	 
		This data structure is intended for continuous data or peak data.
		In addition this peak holds a transformed (e.g. by a wavelet) intensity. 
		For performance reasons this has not been implemented via RichTransPeak1D. 
		If you want to annotated single peaks with meta data, use RichTransPeak1D instead.

		@ingroup Kernel
	*/
	class TransPeak1D : public Peak1D 	
	{
		public:
		
			/** @name Constructors and Destructor
			 */
			//@{
			/// Default constructor
			inline TransPeak1D() 
				: Peak1D(), trans_intensity_(0) 
			{
			}

			/// Copy constructor
			inline TransPeak1D(const TransPeak1D& p) 
				: Peak1D(p), trans_intensity_(p.trans_intensity_)
			{
			}

			/**@brief Destructor

				@note The destructor is non-virtual although many classes are derived from
				TransPeak1D.  This is intentional, since otherwise we would "waste"
				space for a vtable pointer in each instance. Normally you should not derive other classes from
				TransPeak1D (unless you know what you are doing, of course).
			*/
			virtual ~TransPeak1D()
			{
			}
			//@}
			
			/**	
				@name Accessors
			*/
			//@{
			/// Non-mutable access to the data point's transformed intensity (height)
			inline IntensityType getTransIntensity() const { return trans_intensity_; }
			/// Mutable access to the data point's transformed intensity (height)
			inline void setTransIntensity(IntensityType trans_intensity) { trans_intensity_ = trans_intensity; }

			/// Assignment operator
			inline TransPeak1D& operator = (const TransPeak1D& rhs)
			{
				if (this==&rhs) return *this;
			
				intensity_ = rhs.intensity_;
				trans_intensity_ = rhs.trans_intensity_;
				position_ = rhs.position_;
		
				return *this;
			}
			
			/// Equality operator
			inline bool operator == (const TransPeak1D& rhs) const
			{
				return  intensity_ == rhs.intensity_ && position_ == rhs.position_ && trans_intensity_ == rhs.trans_intensity_;
			}

			/// Equality operator
			inline bool operator != (const TransPeak1D& rhs) const
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
				: std::binary_function < TransPeak1D, TransPeak1D, bool >
			{
				inline bool operator () ( TransPeak1D const & left, TransPeak1D const & right ) const
				{
					return ( left.getIntensity() < right.getIntensity() );
				}
				inline bool operator () ( TransPeak1D const & left, IntensityType const & right ) const
				{
					return ( left.getIntensity() < right );
				}
				inline bool operator () ( IntensityType const & left, TransPeak1D const & right ) const
				{
					return ( left< right.getIntensity() );
				}
				inline bool operator () ( IntensityType const & left, IntensityType const & right ) const
				{
					return ( left < right );
				}
			};
		
			struct TransIntensityLess
				: std::binary_function < TransPeak1D, TransPeak1D, bool >
			{
				inline bool operator () ( TransPeak1D const & left, TransPeak1D const & right ) const
				{
					return ( left.getTransIntensity() < right.getTransIntensity() );
				}
				inline bool operator () ( TransPeak1D const & left, IntensityType const & right ) const
				{
					return ( left.getTransIntensity() < right );
				}
				inline bool operator () ( IntensityType const & left, TransPeak1D const & right ) const
				{
					return ( left< right.getTransIntensity() );
				}
				inline bool operator () ( IntensityType const & left, IntensityType const & right ) const
				{
					return ( left < right );
				}
			};

			///Comparator for the position.
			struct PositionLess
				: public std::binary_function <TransPeak1D, TransPeak1D, bool>
			{
				inline bool operator () (const TransPeak1D& a, const TransPeak1D& b) const
				{
					return (a.getMZ() < b.getMZ());
				}
			
				/// comparison of a Peak2D with a CoordinateType
				inline bool operator () ( TransPeak1D const & left, CoordinateType right ) const 
				{
					return (left.getMZ() < right );
				}
				
				/// comparison of a CoordinateType with a Peak2D
				inline bool operator () ( CoordinateType left, TransPeak1D const & right ) const 
				{
					return (left < right.getMZ() );
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
			
			//@}
		
		protected:
			
			/// The data point position
			PositionType	position_;
			/// The original data point intensity
			IntensityType intensity_;
			/// The transformed intensity
			IntensityType trans_intensity_;
	};

	///Print the contents to a stream.
	std::ostream& operator << (std::ostream& os, const TransPeak1D& point);

} // namespace OpenMS

#endif // OPENMS_KERNEL_PEAK1D_H
