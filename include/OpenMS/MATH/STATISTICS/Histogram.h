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

#ifndef OPENMS_MATH_STATISTICS_HISTOGRAM_H
#define OPENMS_MATH_STATISTICS_HISTOGRAM_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

//STL
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>

namespace OpenMS
{
	namespace Math
	{
	
		/**
			@brief Representation of a histogram
			
			The first template argument gives the Type of the
			values that are stored in the bins. The second argument
			gives the type for the bin size and range.
			
			@ingroup Math
		*/
		template < typename ValueType = UnsignedInt, typename BinSizeType = float>
		class Histogram
		{
		 public:

			/// Non-mutable iterator of the bins
			typedef typename std::vector<ValueType>::const_iterator ConstIterator;
	
			/** @name Constructors and Destructors
			 */
			//@{
			///default constructor
			Histogram()
				:	min_(0),
					max_(0),
					bin_size_(0)
			{
			}
	
			///copy constructor
			Histogram(const Histogram& histogram)
				: min_(histogram.min_),
					max_(histogram.max_),
					bin_size_(histogram.bin_size_),
					bins_(histogram.bins_)
			{
			}
	
			///constructor with min, max and bin size
			Histogram(BinSizeType min, BinSizeType max, BinSizeType bin_size) throw(Exception::OutOfRange)
				: min_(min),
					max_(max),
					bin_size_(bin_size)
			{
				if (bin_size_ <= 0)
				{
					throw(Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__));
				}
				else
				{
					// if max_ == min_ there is only one bin
					if (max_ != min_)
					{
						bins_ = std::vector<ValueType>(UnsignedInt(ceil((double(max_)-double(min_))/double(bin_size_))),0);
					}
					else
					{
						bins_ = std::vector<ValueType>(1, 0);
					}
				}
			}
	
			///destructor
			~Histogram()
			{
			}
			//@}
	
			/** @name Accessors
			 */
			//@{
			///returns to the lower bound
			BinSizeType min() const
			{
				return min_;
			}
	
			///returns the upper bound
			BinSizeType max() const
			{
				return max_;
			}
	
			///returns the highest value of all bins
			ValueType maxValue() const
			{
				return *(std::max_element(bins_.begin(), bins_.end()));
			}
	
			///returns the lowest value of all bins
			ValueType minValue() const
			{
				return *(std::min_element(bins_.begin(), bins_.end()));
			}
	
			///returns the bin size
			BinSizeType binSize() const
			{
				return bin_size_;
			}
	
			///returns the number of bins
			Size size() const
			{
				return bins_.size();
			}
			
			///returns the value of bin @p index
			ValueType operator [] (UnsignedInt index) const throw(Exception::IndexOverflow)
			{
				if (index >= bins_.size())
				{
					throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				return bins_[index];
			}
	
			///returns the value of bin corresponding to the value @p val
			ValueType binValue(BinSizeType val) const throw(Exception::OutOfRange)
			{
					return bins_[valToBin_(val)];
			}
	
			///increases the bin corresponding to value @p val by @p increment
			void inc(BinSizeType val, ValueType increment=1) throw(Exception::OutOfRange)
			{
			  bins_[valToBin_(val)]+=increment;
			}
	
			///resets the histogram with the given range and bin size
			void set(BinSizeType min, BinSizeType max, BinSizeType bin_size) throw(Exception::OutOfRange)
			{
				if (bin_size <= 0)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					min_ = min;
					max_ = max;
					bin_size_ = bin_size;
					
					bins_.clear();
					bins_.resize(UnsignedInt(ceil((max_-min_)/bin_size_)),0);
				}
			}
			//@}
	
			/** @name Predicates
			 */
			//@{
			/// equality operator
			bool operator == (const Histogram& histogram) const
			{
				return (min_ == histogram.min_ &&
								max_ == histogram.max_ &&
								bin_size_ == histogram.bin_size_ &&
								bins_ == histogram.bins_);
			}
	
			/// inequality operator
			bool operator != (const Histogram& histogram) const
			{
				return !operator==(histogram);
			}
			//@}
	
			/** @name Assignment
			 */
			Histogram& operator = (const Histogram& histogram)
			{
				if (&histogram == this) return *this;
				
				min_ = histogram.min_;
				max_ = histogram.max_;
				bin_size_ = histogram.bin_size_;
				bins_ = histogram.bins_;
				
				return *this;
			}
			//@}
	
			/** @name Iterators
			 */
			//@{
			/// Non-mutable iterator pointing to the first bin
			inline ConstIterator begin() const { return bins_.begin(); }
	
			/// Non-mutable iterator pointing after the last bin
			inline ConstIterator end() const { return bins_.end(); }
			//@}

			/// Transforms the bin values with f(x)=multiplier*log(x+1) 	 
			void applyLogTransformation(float multiplier) 	 
			{ 	 
				for (typename std::vector<ValueType>::iterator it = bins_.begin(); it!=bins_.end(); ++it) 	 
				{ 	 
					*it = (ValueType)(multiplier*log((float)(*it+1.0f))); 	 
				} 	 
			}
			
		protected:
			/// Lower bound
			BinSizeType min_;
			/// Upper bound
			BinSizeType max_;
			/// Bin size
			BinSizeType bin_size_;
			/// Vector of bins
			std::vector<ValueType> bins_;
			///
			UnsignedInt valToBin_(BinSizeType val) const throw (Exception::OutOfRange)
			{
				if (val < min_ || val > max_)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				if (val == max_)
				{
					return (bins_.size()-1);
				}
				else
				{
					return (UnsignedInt) floor ( (double(val)-double(min_)) / (double(max_)-double(min_)) * bins_.size() );
				}				
			}
		};

} // namespace Math

} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_HISTOGRAM_H
