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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#ifndef OPENMS_KERNEL_AREAITERATOR_H
#define OPENMS_KERNEL_AREAITERATOR_H

// OpenMS includes
#include <OpenMS/CONCEPT/Types.h>

// STL includes
#include <iterator>

namespace OpenMS
{
namespace Internal
{
	/**
		@brief Forward iterator for an area of peaks in an experiment
		
		This iterator allows us to move through the data structure in a linear
		manner i.e. we don't need to jump to the next spectrum manually.
		
		@todo Add filter for MS-level and Int-Range (Marc)
	*/
	template<class ValueT, class ReferenceT, class PointerT, class SpectrumIteratorT, class PeakIteratorT>
	class AreaIterator : public std::iterator<std::bidirectional_iterator_tag, ValueT>
	{
		typedef AreaIterator<ValueT, ReferenceT, PointerT, SpectrumIteratorT, PeakIteratorT> Self;

	public:
		typedef double CoordinateType;
		typedef ValueT PeakType;
		typedef SpectrumIteratorT SpectrumIteratorType;
		typedef PeakIteratorT PeakIteratorType;
		
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
		AreaIterator()
		{}

		/// Constructor
		AreaIterator(SpectrumIteratorType begin, SpectrumIteratorType end, CoordinateType low_mz, CoordinateType high_mz)
			: current_scan_(begin), last_scan_(end), low_mz_(low_mz), high_mz_(high_mz)
		{
			nextScan_();
		}

		/// Destructor
		~AreaIterator()
		{}

		/// Copy constructor
		AreaIterator(const AreaIterator& rhs)
			: current_scan_(rhs.current_scan_),
				last_scan_(rhs.last_scan_),
				current_peak_(rhs.current_peak_),
				last_peak_(rhs.last_peak_),
				low_mz_(rhs.low_mz_),
				high_mz_(rhs.high_mz_)
		{}

		/// Assignment operator
		AreaIterator& operator=(const AreaIterator& rhs)
		{
			if (&rhs == this) return *this;

			current_scan_ = rhs.current_scan_;
			last_scan_ = rhs.last_scan_;
			current_peak_ = rhs.current_peak_;
			last_peak_ = rhs.last_peak_;
			low_mz_ = rhs.low_mz_;
			high_mz_ = rhs.high_mz_;

			return (*this);
		}

		/// Test for equality
		bool operator==(const AreaIterator& rhs)
		{
			// only test for equality to the end iterator
			return current_scan_ == last_scan_ && current_peak_ == last_peak_ &&
			       rhs.current_scan_ == rhs.last_scan_ && rhs.current_peak_ == rhs.last_peak_;
			
/*			return ( current_scan_ == rhs.current_scan_ &&
			         last_scan_ == rhs.last_scan_ &&
			         current_peak_ == rhs.current_peak_ &&
			         last_peak_ == rhs.last_peak_ &&
			         low_mz_ == rhs.low_mz_ &&
			         high_mz == rhs.high_mz_ );*/
		}

		/// Test for inequality
		bool operator!=(const AreaIterator& rhs)
		{
			return !(*this == rhs);
		}

		/// Step forward by one (prefix operator)
		AreaIterator& operator++()
		{
			++current_peak_;
			// test whether we arrived at the end of the current scan
			if (current_peak_ == current_scan_->end())
			{
				++current_scan_;
				nextScan_();
			}
			return (*this);
		}

		/// Step forward by one (postfix operator)
		AreaIterator operator++(int)
		{
			AreaIterator tmp(*this);
			++(*this);
			return tmp;
		}

		/// Dereferencing of this pointer yields the underlying peak
		reference operator*()
		{
			return current_peak_.operator*();
		}

		/// Dereferencing of this pointer yields the underlying peak
		pointer operator->()
		{
			return current_peak_.operator->();
		}
		
	private:
		void nextScan_()
		{
			if (current_scan_ != last_scan_)
			{
				current_peak_ = current_scan_->MZBegin(low_mz_);
				last_peak_ = current_scan_->MZBegin(high_mz_);
			}
			else
			{
				current_peak_ = last_peak_ = PeakIteratorType();
			}
		}
	
		/// Iterator to the current spectrum
		SpectrumIteratorType current_scan_;
		/// Iterator to the last spectrum
		SpectrumIteratorType last_scan_;
		/// Iterator to the current peak
		PeakIteratorType current_peak_;
		/// Iterator to the last peak in the current spectrum
		PeakIteratorType last_peak_;
		
		CoordinateType low_mz_;
		CoordinateType high_mz_;
	};
}
}

#endif
