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

#ifndef OPENMS_KERNEL_RANGEMANAGER_H
#define OPENMS_KERNEL_RANGEMANAGER_H

#include <OpenMS/DATASTRUCTURES/DRange.h>

namespace OpenMS
{
	/**	
		@brief Handles the managment of a position and intensity range.
		
		This is needed for all peak and feature container like Spectrum, MSExperiment and FeatureMap.
	*/
	template <UInt D>
	class RangeManager
	{
		public:
			/// Dimension of the position range
			enum { DIMENSION = D };
			/// Position range type			
			typedef DRange<D> PositionRangeType;
			/// Position Type
			typedef DPosition<D> PositionType;
			/// Intensity range type			
			typedef DRange<1> IntensityRangeType;
			
			/// Default constructor
			RangeManager()
				: int_range_(),
					pos_range_()
			{
			}

			/// Copy constructor
			RangeManager(const RangeManager& rhs)
				: int_range_(rhs.int_range_),
					pos_range_(rhs.pos_range_)
			{
			}
			
			/// Destructor
			virtual ~RangeManager()
			{
			}
	
			/// Assignment operator
			RangeManager& operator = (const RangeManager& rhs)
			{
				if (this==&rhs) return *this;
				
				int_range_ = rhs.int_range_;
				pos_range_ = rhs.pos_range_;
				
				return *this;
			}
	
			/// Equality operator
			bool operator == (const RangeManager& rhs) const
			{
				return
					int_range_ == rhs.int_range_ &&
					pos_range_ == rhs.pos_range_
					;				
			}
			
			/// Equality operator
			bool operator != (const RangeManager& rhs) const
			{
				return !(operator==(rhs));
			}

			
			/**	
				@name Range methods
			
				@note The range values are not updated automatically. Call updateRanges() to update the values!
			*/
			//@{
	
			/// Returns the minimum position
			const PositionType& getMin() const	
			{ 
				return pos_range_.minPosition(); 
			}
			/// Returns the maximum position
		  const PositionType& getMax() const 
		  { 
		  	return pos_range_.maxPosition(); 
		  }
	
			/// Returns the minimum intensity
			DoubleReal getMinInt() const	
			{ 
				return int_range_.minPosition()[0]; 
			}
			/// Returns the maximum intensity
		  DoubleReal getMaxInt() const 
		  { 
		  	return int_range_.maxPosition()[0]; 
		  }
			
			/**
				@brief Updates minimum and maximum position/intensity.
				
				This method is usually implemented by calling clearRanges() and
				updateRanges_().
			*/
			virtual void updateRanges() = 0;

			/// Resets the ranges
			void clearRanges()
			{
				int_range_ = IntensityRangeType::empty;
				pos_range_ = PositionRangeType::empty;
			}

			//@}

		protected:
			/// Intensity range (1-dimensional)
			IntensityRangeType int_range_;
			/// Position range (D-dimensional)
			PositionRangeType pos_range_;
			
			/// Updates the range using data points in the iterator range. 
			template <class PeakIteratorType>
			void updateRanges_(const PeakIteratorType& begin, const PeakIteratorType& end)
			{
				
				//prevent invalid range by empty container
				if (begin==end)
				{
					return;
				}
				
				PositionType min = pos_range_.minPosition();
				PositionType max = pos_range_.maxPosition();
				
				DoubleReal it_min = int_range_.minPosition()[0];
				DoubleReal it_max = int_range_.maxPosition()[0];
				
				for (PeakIteratorType it = begin; it != end; ++it)
				{
					//update position
					for (UInt i = 0; i < D; ++i)
					{
						DoubleReal tmp = it->getPosition()[i];
						if (tmp < min[i])
						{
							min[i] = tmp;
						}
						if (tmp > max[i])
						{
							max[i] = tmp;
						}
					}
				
					//update intensity
					DoubleReal tmp = it->getIntensity();
					if (tmp < it_min)
					{
						it_min = tmp;
					}
					if (tmp > it_max)
					{
					  it_max = tmp;
					}
				}
				
				pos_range_.setMin(min);
				pos_range_.setMax(max);
				
				int_range_.setMinX(it_min);
				int_range_.setMaxX(it_max);
      }
      
	};  // class
  
}  // namespace OpenMS

#endif  // OPENMS_KERNEL_DRANGE_H
