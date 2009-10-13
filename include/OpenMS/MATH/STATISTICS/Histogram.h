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
// $Maintainer: $
// $Authors: Marc Sturm $
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
#include <algorithm>

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
		template<typename ValueType=UInt, typename BinSizeType=DoubleReal>
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
	
			/**
			  @brief constructor with min, max and bin size

			  @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
			*/
			Histogram(BinSizeType min, BinSizeType max, BinSizeType bin_size)
				: min_(min),
					max_(max),
					bin_size_(bin_size)
			{
				if (bin_size_ <= 0)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					// if max_ == min_ there is only one bin
					if (max_ != min_)
					{
						bins_ = std::vector<ValueType>(Size(ceil((max_-min_)/bin_size_)),0);
					}
					else
					{
						bins_ = std::vector<ValueType>(1, 0);
					}
				}
			}
	
			///destructor
			virtual ~Histogram()
			{
			}
			//@}
	
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
			
			/**
			  @brief returns the value of bin @p index

			  @exception Exception::IndexOverflow is thrown for invalid indices
			*/
			ValueType operator [] (Size index) const
			{
				if (index >= bins_.size())
				{
					throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				return bins_[index];
			}
	
			/**
			  @brief returns the center position of the bin with the index @p bin_index

			  @exception Exception::IndexOverflow is thrown for invalid indices
			*/
			BinSizeType centerOfBin(Size bin_index) const
			{
				if (bin_index >= bins_.size())
				{
					throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				
				return (BinSizeType)(min_+((BinSizeType)bin_index+0.5)*bin_size_);
			}

			/**
			  @brief returns the value of bin corresponding to the value @p val

			  @exception Exception::OutOfRange is thrown if the value is out of valid range
			*/
			ValueType binValue(BinSizeType val) const
			{
				return bins_[valToBin_(val)];
			}
	
			/**
			  @brief increases the bin corresponding to value @p val by @p increment
				
				@return The index of the increased bin.
				
			  @exception Exception::OutOfRange is thrown if the value is out of valid range
			*/
			Size inc(BinSizeType val, ValueType increment=1)
			{
				Size bin_index = valToBin_(val);
			  bins_[bin_index]+=increment;
			  return bin_index;
			}
	
			/**
			  @brief resets the histogram with the given range and bin size

			  @exception Exception::OutOfRange is thrown if @p bin_size negative or zero
			*/
			void reset(BinSizeType min, BinSizeType max, BinSizeType bin_size)
			{
				min_ = min;
				max_ = max;
				bin_size_ = bin_size;
				bins_.clear();
					
				if (bin_size_ <= 0)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					// if max_ == min_ there is only one bin
					if (max_ != min_)
					{
						bins_.resize(Size(ceil((max_-min_)/bin_size_)),0);
					}
					else
					{
						bins_.resize(1, 0);
					}
				}
			}
	
			/** @name Assignment and equality operators
			 */
			//@{
			///Equality operator
			bool operator == (const Histogram& histogram) const
			{
				return (min_ == histogram.min_ &&
								max_ == histogram.max_ &&
								bin_size_ == histogram.bin_size_ &&
								bins_ == histogram.bins_);
			}
	
			///Inequality operator
			bool operator != (const Histogram& histogram) const
			{
				return !operator==(histogram);
			}
	
			///Assignment
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
			void applyLogTransformation(BinSizeType multiplier) 	 
			{ 	 
				for (typename std::vector<ValueType>::iterator it = bins_.begin(); it!=bins_.end(); ++it) 	 
				{ 	 
					*it = (ValueType)(multiplier*log((BinSizeType)(*it+1.0f))); 	 
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
			/**
			  @brief Returns the bin a given value belongs to

			  @exception Exception::OutOfRange is thrown if the value is out of valid range
			*/
			Size valToBin_(BinSizeType val) const
			{
				//std::cout << "val: " << val << " (min: " << min_ << " max: " << max_ << ")" << std::endl;
				if (val < min_ || val > max_)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				if (val == max_)
				{
					return Size(bins_.size()-1);
				}
				else
				{
					return (Size) floor ( (val-min_) / (max_-min_) * bins_.size() );
				}				
			}
		};

		///Print the contents to a stream.
		template<typename ValueType, typename BinSizeType>
		std::ostream& operator << (std::ostream& os, const Histogram<ValueType,BinSizeType>& hist)
		{
			for (Size i=0; i<hist.size(); ++i)
			{
				os << hist.centerOfBin(i) << "	" << hist[i] << std::endl;
			}
			return os;
		}

	} // namespace Math

} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_HISTOGRAM_H
