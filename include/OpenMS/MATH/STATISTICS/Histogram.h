// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Histogram.h,v 1.9 2006/06/08 15:51:32 marc_sturm Exp $
// $Author: marc_sturm $
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
			@brief Representation of a Histogram
			
			The first template argument gives the Type of the
			values that are stored in the bins. The second argument
			gives the type of the binsize. This type is also used for
			the range of the stored values.
			
			@ingroup Math
		*/
		template < typename T = UnsignedInt, typename BinSizeType = float>
		class Histogram
		{
		 public:
	
			/** @name Typedefs
			 */
			//@{
			typedef T ValueType;
			typedef BinSizeType BinType;
			typedef typename std::vector<T>::const_iterator ConstIterator;
			typedef typename std::vector<T>::const_iterator const_iterator;
			typedef typename std::vector<T>::iterator Iterator;
			typedef typename std::vector<T>::iterator iterator;
			//@}
	
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
						bins_ = std::vector<T>(UnsignedInt(ceil((double(max_)-double(min_))/double(bin_size_))),0);
					}
					else
					{
						bins_ = std::vector<T>(1, 0);
					}
				}
			}
	
			///destructor
			virtual ~Histogram()
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
			T maxValue() const
			{
				T max = std::numeric_limits<T>::min();
				for (Size i=0;i!=bins_.size();++i)
				{
					if (bins_[i] > max)
					{
						max = bins_[i];
					}
				}
				return max;
			}
	
			///returns the lowest value of all bins
			T minValue() const
			{
				T min = std::numeric_limits<T>::max();
				for (Size i=0;i!=bins_.size();++i)
				{
					if (bins_[i] < min)
					{
						min = bins_[i];
					}
				}
				return min;
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
	
			///returns the value of bin <i>index</i>
			T bin(UnsignedInt index) const throw(Exception::IndexOverflow)
			{
				if (index>=bins_.size())
				{
					throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
	
				return bins_[index];
			}
	
			///returns the value of bin <i>index</i>
			T bin(SignedInt index) const throw(Exception::IndexOverflow,Exception::IndexUnderflow)
			{
				if (index < 0)
				{
					throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					if ((UnsignedInt)index >= bins_.size())
					{
						throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
					}
				}
				return binValue((UnsignedInt)index);
			}
	
			T operator [] (UnsignedInt index) const throw(Exception::IndexOverflow)
			{
				if (index >= bins_.size())
				{
					throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				return bins_[index];
			}
	
			T operator [] (SignedInt index) const throw(Exception::IndexOverflow, Exception::IndexUnderflow)
			{
				if (index < 0)
				{
					throw Exception::IndexUnderflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					if ((UnsignedInt)index >= bins_.size())
					{
						throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__);
					}
				}
				return bins_[(UnsignedInt)index];
			}
	
			///returns the value of bin corresponding to the value <i>val</i>
			T binValue(BinSizeType val) const throw(Exception::OutOfRange)
			{
				if (val < min_ || val > max_)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				return bins_[UnsignedInt((val-min_)/bin_size_)];
			}
	
			///increases the bin corresponding to value <i>val</i> by <i>increment</i> and returns the value of the bin
			T inc(BinSizeType val, T increment=1) throw(Exception::OutOfRange)
			{
				if (val < min_ || val > max_)
				{
					throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				}
				else
				{
					return (bins_[UnsignedInt(ceil((double(val)-double(min_))/double(bin_size_)))]+=increment);
				}
			}
	
			///resets the histogram to the given values (can be used with the default constructur)
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
					// (cg) 2005-04-01 suggest the following:
					// bins_.clear();
					// bins_.resize(UnsignedInt(ceil((max_-min_)/bin_size_)),0);
					// (cg) 2005-04-01 instead of this:
					bins_ = std::vector<T>(UnsignedInt(ceil((max_-min_)/bin_size_)),0);
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
				return (min_ != histogram.min_ ||
								max_ != histogram.max_ ||
								bin_size_ != histogram.bin_size_ ||
								bins_ != histogram.bins_);
			}
			//@}
	
			/** @name Assignment
			 */
			Histogram& operator = (const Histogram& histogram)
			{
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
			/// constant iterator pointing to the min value of the histogram
			inline ConstIterator begin() const { return bins_.begin(); }
	
			/// constant iterator pointing to one position behind the max value of the histogram
			inline ConstIterator end() const { return bins_.end(); }
	
			/// non-constant iterator pointing to the min value of the histogram
			inline Iterator begin() { return bins_.begin(); }
	
			/// non-constant iterator pointing to one position behind the max value of the histogram
			inline Iterator end() { return bins_.end(); }
			//@}
	
		 protected:
	
			BinSizeType min_;
	
			BinSizeType max_;
	
			BinSizeType bin_size_;
	
			std::vector<T> bins_;
		};

} // namespace Math

} // namespace OpenMS

#endif // OPENMS_MATH_STATISTICS_HISTOGRAM_H
