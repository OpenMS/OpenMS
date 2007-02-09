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
			
			This iterator iterates over spectra with MS level 1 only!
		*/
		template<class ValueT, class ReferenceT, class PointerT, class SpectrumIteratorT, class PeakIteratorT>
		class AreaIterator : public std::iterator<std::forward_iterator_tag, ValueT>
		{
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
		
				/// Constructor for the begin iterator
				AreaIterator(SpectrumIteratorType begin, SpectrumIteratorType end, CoordinateType low_mz, CoordinateType high_mz)
					: current_scan_(begin), 
						end_scan_(end), 
						low_mz_(low_mz), 
						high_mz_(high_mz)
				{
					nextScan_();
				}
		
				/// Constructor for the end iterator
				AreaIterator( SpectrumIteratorType spectrum_end, PeakIteratorType peak_end )
					: current_scan_(spectrum_end), 
						end_scan_(spectrum_end), 
						current_peak_(peak_end), 
						end_peak_(peak_end)
				{
				}
		
				/// Destructor
				~AreaIterator()
				{
				}
		
				/// Copy constructor
				AreaIterator(const AreaIterator& rhs)
					: current_scan_(rhs.current_scan_),
						end_scan_(rhs.end_scan_),
						current_peak_(rhs.current_peak_),
						end_peak_(rhs.end_peak_),
						low_mz_(rhs.low_mz_),
						high_mz_(rhs.high_mz_)
				{
				}
		
				/// Assignment operator
				AreaIterator& operator=(const AreaIterator& rhs)
				{
					if (&rhs == this) return *this;
		
					current_scan_ = rhs.current_scan_;
					end_scan_ = rhs.end_scan_;
					current_peak_ = rhs.current_peak_;
					end_peak_ = rhs.end_peak_;
					low_mz_ = rhs.low_mz_;
					high_mz_ = rhs.high_mz_;
		
					return *this;
				}
		
				/// Test for equality
				bool operator==(const AreaIterator& rhs) const
				{
					return ( &(*current_peak_) == &(*(rhs.current_peak_))) //Equality of pointed to peak adresses
								 ||
								 ( current_peak_ == end_peak_ && //Equality to the end iterator
								   current_scan_ == end_scan_ &&
								   rhs.current_scan_ == rhs.end_scan_ &&  
								   rhs.current_peak_ == rhs.end_peak_ ) ;
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
					if (current_peak_ == end_peak_)
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
				reference operator*() const
				{
					return current_peak_.operator*();
				}
		
				/// Dereferencing of this pointer yields the underlying peak
				pointer operator->() const
				{
					return current_peak_.operator->();
				}
	
				/// returns the retention time of the current scan
				const CoordinateType& getRetentionTime() const
				{
					return current_scan_->getRetentionTime();
				}
				
			private:
				//Advances to the iterator to the next valid peak in the next valid spectrum
				void nextScan_()
				{
					while (true)
					{
						//if (current_scan_ != end_scan_) std::cout << "RT: " << current_scan_->getRetentionTime() << std::endl;
						while (current_scan_ != end_scan_ && current_scan_->getMSLevel()!=1)
						{
							++current_scan_;
						}
						if (current_scan_ == end_scan_)
						{
							current_peak_ = end_peak_ = PeakIteratorType();
							return;
						}
						current_peak_ = current_scan_->MZBegin(low_mz_);
						end_peak_ = current_scan_->MZEnd(high_mz_);
						if (current_peak_!=end_peak_)
						{
							return;
						}
						++current_scan_;
					}
				}
			
				/// Iterator to the current spectrum
				SpectrumIteratorType current_scan_;
				/// Past-the-end iterator of spectra
				SpectrumIteratorType end_scan_;
				/// Iterator to the current peak
				PeakIteratorType current_peak_;
				/// Past-the-end iterator of peaks in the current spectrum
				PeakIteratorType end_peak_;
				/// low m/z boundary
				CoordinateType low_mz_;
				/// high m/z boundary
				CoordinateType high_mz_;
			
			private:
				//Hidden default constructor
				AreaIterator();
		};

	}
}

#endif
