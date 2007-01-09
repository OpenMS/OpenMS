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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_PEAKITERATOR_H
#define OPENMS_KERNEL_PEAKITERATOR_H

// OpenMS includes
#include <OpenMS/CONCEPT/Types.h>

// STL includes
#include <iterator>

namespace OpenMS
{
	/**
		@brief Adaptor class for bidirectional iterator on objects of DPeak<1>
		
		This iterator allows us to move through the data structure in a linear
		manner i.e. we don't need to jump to the next spectrum manually.
		
		The class has a member  DPeakArray<>::iterator pointing
		to the current peak. The class also remembers the retention time of the current
		scan.
	*/
	template<class ValueT, class ReferenceT, class PointerT, class ExperimentT>
	class PeakIterator : public std::iterator<std::bidirectional_iterator_tag,  ValueT>
	{
		typedef PeakIterator<ValueT, ReferenceT, PointerT, ExperimentT> Self;
	
	public:
		typedef double CoordinateType;
		typedef ValueT PeakType;
		typedef ExperimentT ExperimentType;
		
		/** @name Typedefs for STL compliance
		*/
		//@{
		/// The iterator's value type
		typedef ValueT value_type;
		/// The reference type as returned by operator*()
		typedef ReferenceT reference;
		/// The pointer type as returned by operator->()
		typedef PointerT pointer;
		/// The difference type
		typedef unsigned int difference_type;
		//@}

		/// Default constructor
		PeakIterator()
			: peak_index_(), rt_(), scan_index_(), exp_()
		{}

		/// Constructor
		PeakIterator(UnsignedInt pind, CoordinateType & co, UnsignedInt sind, ExperimentType& exp)
			: peak_index_(pind), rt_(co), scan_index_(sind), exp_(&exp)
		{}

		/// Destructor
		~PeakIterator()
		{}

		/// Copy constructor
		PeakIterator(const PeakIterator& rhs)
			: peak_index_(rhs.peak_index_), rt_(rhs.rt_),
				scan_index_(rhs.scan_index_), exp_(rhs.exp_)
		{}

		/// Assignment operator
		Self& operator=(const PeakIterator& rhs)
		{
			if (&rhs == this) return *this;

			peak_index_ = rhs.peak_index_;
			rt_         = rhs.rt_;
			scan_index_ = rhs.scan_index_;
			exp_        = rhs.exp_;

			return (*this);
		}

		/// Test for equality
		bool operator==(const PeakIterator& rhs)
		{
			return ( peak_index_ == rhs.peak_index_ &&
			         rt_ == rhs.rt_ &&
			         scan_index_ == rhs.scan_index_ );
		}

		/// Test for inequality
		bool operator!=(const PeakIterator& rhs)
		{
			return !(*this == rhs);
		}

		/// Step forward by one (prefix operator)
		Self& operator++()
		{
			++peak_index_;
			// test whether we arrived at the end of the current scan
			if ( peak_index_ >= (*exp_)[scan_index_].size() && scan_index_ != ( (*exp_).size() - 1) )
			{
				// we are at the end of a scan, but this scan is not the very last one
				// so we can jump into the next scan
				peak_index_ = 0;
				++scan_index_;
				rt_ = (*exp_)[scan_index_].getRetentionTime();
			}
			return (*this);
		}

		/// Step backward by one (prefix operator)
		Self& operator--()
		{
			// test whether we are at the start of a scan
			if (peak_index_  == 0)
			{
				// update scan index and move to end of previous scan
				if (scan_index_ == 0)
				{
					std::cout << "PeakIterator: In first scan and moving backwards ! " << std::endl;
					return (*this);
				}
				--scan_index_;
				peak_index_  = ( (*exp_)[scan_index_].size() -1) ;
				rt_          = (*exp_)[scan_index_].getRetentionTime();
			}
			else
			{
				// simply one step backwards
				--peak_index_;
			}
			return (*this);
		}

		/// Step forward by one (postfix operator)
		Self operator++(int)
		{
			PeakIterator tmp(*this);
			++(*this);
			return tmp;
		}

		/// Step backward by one (postfix operator)
		Self operator--(int)
		{
			PeakIterator tmp(*this);
			--(*this);
			return tmp;
		}

		/// Dereferencing of this pointer yields the underlying peak
		reference operator*()
		{
			return (*exp_)[scan_index_][peak_index_];
		}

		/// Dereferencing of this pointer yields the underlying peak
		pointer operator->()
		{
			return &((*exp_)[scan_index_][peak_index_]);
		}
		
		/** @name Accesssors
		*/
		//@{
		/// Returns the current retention time (mutable)
		CoordinateType& getRt() { return rt_; }
		/// Returns the current retention time (not mutable)
		const CoordinateType& getRt() const { return rt_; }
		/// Returns the index of the peak this iterator points to 
		/// NOTE: Call updateRanges() before using this function
		UnsignedInt getPeakNumber()
		{
			if (scan_index_ > 0)
				return (exp_->spectra_lengths_[ (scan_index_-1) ] + peak_index_);
			else
				return peak_index_;
		}
		//@}

	private:
		/// Points to the current peak
		UnsignedInt peak_index_;
		/// Retention time of the current spectrum
		CoordinateType rt_;
		/// Index of the current spectrum
		UnsignedInt scan_index_;
		/// Pointer to the experiment
		ExperimentType* exp_;
	};
}

#endif
