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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_SCANINDEXMSEXPERIMENT_H
#define OPENMS_DATASTRUCTURES_SCANINDEXMSEXPERIMENT_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

	/** @brief Adaptor class for peak container for rapid navigation between scans.
	
		This class supports class MSExperiment and all OpenMS datastructures with
		the same (or similar interfaces).

		This data structure allows us to move rapidly between different scans by storing begin and end iterator
		of a scan and performing a binary search for the m/z of the interator range.
	*/
	template < typename PeakContainer_, typename PeakContainerIterator_ = typename PeakContainer_::PIterator >
	class ScanIndexMSExperiment : public std::vector < PeakContainerIterator_ >
	{
	 public:

		typedef PeakContainer_ PeakContainerType;
		typedef typename PeakContainerType::PeakType PeakType;			
		typedef typename PeakType::CoordinateType CoordinateType;
		typedef PeakContainerIterator_ PeakIterator;
		
		typedef std::vector < CoordinateType > ScanPositionContainerType;
		typedef std::vector < PeakIterator > ScanBeginContainerType;
		typedef typename PeakType::template NthPositionLess<0> MZless;
    	    
		using ScanBeginContainerType::clear;
		using ScanBeginContainerType::begin;
		using ScanBeginContainerType::end;
		using ScanBeginContainerType::size;
		using ScanBeginContainerType::back;
		using ScanBeginContainerType::push_back;
		
		/// Constructor
		ScanIndexMSExperiment()
			:  std::vector < PeakContainerIterator_ >(),
				last_rt_(), last_rank_()
		{
		}
		
		/// Copy Constructor
		ScanIndexMSExperiment(const ScanIndexMSExperiment& rhs )
			:  std::vector < PeakContainerIterator_ >(rhs),
				last_rt_(rhs.last_rt_), last_rank_(rhs.last_rank_)
		{
		}
		
		/// Destructor
		~ScanIndexMSExperiment() 
		{
		}	
		
		/// Assignment 
		ScanIndexMSExperiment& operator = (const ScanIndexMSExperiment& rhs)
		{
			if (this==&rhs) return *this;
			
			std::vector < PeakContainerIterator_ >::operator = (rhs);
			last_rt_     = rhs.last_rt_;
			last_rank_ = rhs.last_rank_;
			
			return *this;
		}

		/// Test for equality 
		bool operator == (const ScanIndexMSExperiment& rhs) const
		{
			return (std::vector < PeakContainerIterator_ >::operator == (rhs) &&
					    		last_rt_     == rhs.last_rt_ &&
									last_rank_ == rhs.last_rank_);	
		}

		/// Test for inequality
		bool operator != (const ScanIndexMSExperiment& rhs) const
		{
			return !(operator == (rhs));
		}	
		
		/// We throw this exception if the next (previous) peak is requested for a peak in the last (first) scan.
		class NoSuccessor
     : public Exception::Base
     {
     public:
       NoSuccessor(const char* file, int line, const char* function, const UnsignedInt& index) throw()
       :	Base(file, line, function, "NoSuccessor", "no successor/predecessor"), index_(index)
			{
				what_ = "there is no successor/predecessor for the given Index: " + index_;
				Exception::globalHandler.setMessage(what_);
			}
       
       virtual ~NoSuccessor() throw(){}
       
     protected:
       UnsignedInt index_;  // index without successor/predecessor
       
     }; // end of class NoSuccessor
		
		/** @brief Initialize the container.
		  
		    Precondition: The range is sorted with respect to
				ComparatorType.  Otherwise the result will be garbage (like with std::unique_copy()).
				Note that <code>this->size() == this->scan_position_.size() + 1</code>
				because the first one deals with ranges (pairs).
				
				@NOTE: Use the peak iterators peakBegin() and peakEnd() to initialize the DS.
				
		*/
		void init ( PeakIterator _begin, PeakIterator const _end ) throw ()
		{
				clear();
				push_back ( _begin );
      
      	// iterate over range, save the iterator if
      	// a new scan starts.
      	while ( ++_begin != _end )
	  		{
					if (  back().getRt() < _begin.getRt() ) 
					{
						push_back ( _begin );
					}
			}		
			scan_position_.clear();
			scan_position_.reserve ( size() );
			
			for ( typename ScanBeginContainerType::const_iterator peak1 = begin();
						peak1 != end();
						++peak1
					)
			{
				scan_position_.push_back ( (*peak1) .getRt() );
			}
			
			push_back ( _end ); // we will need the end() of the last scan as well
		}
		
		/** @brief Moves to the next scan.
		    
		    Retrieves the peak in the next scan whose m/z is closest
		    to @p peak.
		*/
		PeakIterator getNextRt(const CoordinateType& current_rt, const CoordinateType& current_mz) throw (NoSuccessor)
		{
			int current_scan = 0;
			if (current_rt == last_rt_)
			{
				current_scan =   last_rank_;
			}
			else
			{
				current_scan	= getRank(current_rt);
				last_rank_     = current_scan;
			}
			if (current_scan >= int(size()-2)) throw NoSuccessor(__FILE__, __LINE__, "getNextRt()", current_scan);
		
			// determine start and end of the next scan
			PeakIterator scan_begin = (*this)[current_scan+1];
	 		PeakIterator scan_end   = (*this)[current_scan+2];	// Seems to be dangerous, but does work.
	 		
			return searchInScan_(scan_begin,scan_end,current_mz);			
		}
		
		/** @brief Moves to the previous scan.
		    
		    Retrieves the peak in the previous scan whose m/z is closest
		    to @p peak.
		*/
		PeakIterator getPrevRt(const CoordinateType& current_rt, const CoordinateType& current_mz) throw (NoSuccessor)
		{			
			
			int current_scan = 0;
			if (current_rt == last_rt_)
			{
				current_scan =   last_rank_;
			}
			else
			{
				current_scan	= getRank(current_rt);
				last_rank_     = current_scan;
			}
			// if we are already in the first scan, there will be no predeccessor....
			if (current_scan == 0) throw NoSuccessor(__FILE__, __LINE__, "getPrevRt()", current_scan);
		
			// determine start and end of the next scan
			PeakIterator scan_begin = (*this)[current_scan-1];
			PeakIterator scan_end   = (*this)[current_scan];
		
			// binary search
			return searchInScan_(scan_begin,scan_end,current_mz);
		}
		     
    	 /// Returns the rank of position \p coord, starting with 0, usually this will give us the scan number.	 
		typename ScanPositionContainerType::size_type getRank ( CoordinateType const & coord ) const throw()
		{
			// perform binary search to retreive rank of this retention time
			return std::lower_bound
				( scan_position_.begin(),
					scan_position_.end(),
					coord
				) - scan_position_.begin();
		}
		
		
		protected:	
		
		/// Performs binary search on an iterator range to find the
		/// peak with m/z that comes closest to peak at @p current_mz.
		PeakIterator searchInScan_(PeakIterator scan_begin, PeakIterator scan_end, CoordinateType current_mz) const
		{
		
			// perform binary search to find the neighbour in rt dimension
			PeakIterator insert_iter = std::lower_bound(scan_begin,scan_end,current_mz,MZless());	
					
			// In some cases (i.e. for picked peaks) the binary search does
			// not give us the correct position. Therefore we check the peaks
			// to our left and right hand side and verify which one has a m/z
			// closer to the original peak.
			if ( insert_iter == scan_end ) // only one choice
			{
		 			return --insert_iter;
			}
			else
			{
				// if the found peak is at the beginning of the spectrum,
				// there is not much we can do.
				if ( insert_iter == scan_begin ) 
				{
					return insert_iter;
				}
				else // see if the next smaller one fits better
				{
					CoordinateType delta_mz = (insert_iter->getPos() - current_mz);
					--insert_iter;
									
					if ( current_mz - insert_iter->getPos() > delta_mz )
					{
						return insert_iter; // peak to the right is closer (in m/z dimension)
					}
					else
					{
						return ++insert_iter;    // peak to the left is closer
					}
				}
			}
		} // end of searchInScan_ 
		
		/// Records the retention time for a given scan
		ScanPositionContainerType scan_position_;  	
		CoordinateType last_rt_;
		int last_rank_;
		
	};
	  
} // namespace OpenMS

#endif //  OPENMS_DATASTRUCTURES_SCANINDEXMSEXPERIMENT_H
