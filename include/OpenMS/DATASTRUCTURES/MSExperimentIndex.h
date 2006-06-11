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
// $Id: MSExperimentIndex.h,v 1.3 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_MSEXPERIMENTINDEX_H
#define OPENMS_DATASTRUCTURES_MSEXPERIMENTINDEX_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/KERNEL/DimensionDescription.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{

	/**
		@brief This is an adaptor class for MSExperiment.

	*/
	template < typename PeakT_  >
	class MSExperimentIndex : public std::vector < typename MSExperiment< PeakT_ >::PeakIterator >
	{
	 public:	 
	 	enum DimensionId
			{ 
				RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
				MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
			};
	 	
		typedef PeakT_ PeakType;
		typedef typename MSExperiment< PeakT_ >::PeakIterator PeakIterator;
		typedef typename PeakType::CoordinateType CoordinateType;
			
		typedef std::vector < CoordinateType > ScanPositionContainerType;
		typedef std::vector < PeakIterator > ScanBeginContainerType;
    
		using ScanBeginContainerType::clear;
		using ScanBeginContainerType::begin;
		using ScanBeginContainerType::end;
		using ScanBeginContainerType::size;
		using ScanBeginContainerType::back;
		using ScanBeginContainerType::push_back;
		
		/// We throw this exception if the next (previous) peak is requested for a peak in the last (first) scan.
		class NoSuccessor
     : public Exception::Base
     {
     public:
       NoSuccessor(const char* file, int line, const char* function, const UnsignedInt& index) throw()
       :	Base(file, line, function, "NoSuccessor", "no successor/predecessor"), index_(index)
			{
				what_ = "there is no successor/predecessor for the given Index: " + String(index_);
				OpenMS::Exception::globalHandler.setMessage(what_);
			}
       
       virtual ~NoSuccessor() throw(){}
       
     protected:
       UnsignedInt index_;  // index without successor/predecessor
       
     }; // end of class NoSuccessor
		
		/**
			@brief Initialize the container.
		  
		  Precondition: The range is sorted with respect to
			ComparatorType.  Otherwise the result will be garbage (like with std::unique_copy()).
			Note that <code>this->size() == this->scan_position_.size() + 1</code>
			because the first one deals with ranges (pairs).
		*/
		void init ( PeakIterator _begin, PeakIterator _end ) throw ()
		{
			clear();
			push_back ( _begin );
      
      		// iterate over range and save the iterator if a new scan starts.
			for ( ; _begin != _end; ++_begin )
			{
				if ( this->back().getRt() < _begin.getRt() ) push_back ( _begin );
			}		
			
			scan_position_.clear();
			scan_position_.reserve ( size() );
			
			for ( typename ScanBeginContainerType::const_iterator iter = begin();
						iter != end();
						++iter
					)
			{
				scan_position_.push_back ( iter->getRt() );
				//std::cout << "Inserted rt: " <<  iter->getRt() << std::endl;
			}
			
			push_back ( _end ); // we will need the end() of the last scan as well
		}
		
		/**
			@brief Moves to the next scan.
		    
	    Retrieves the peak in the next scan whose m/z is closest
	    to @p peak.
		*/
		PeakIterator getNextRt(DRawDataPoint<2> const & peak) const throw (NoSuccessor)
		{
			//std::cout << "In getNextRt() " << std::endl;
			CoordinateType current_mz = peak.getPosition()[MZ];
			CoordinateType current_rt   = peak.getPosition()[RT];
			
			int const current_scan    = getRank(current_rt);
			//std::cout << "Rank: " << current_scan << std::endl;
			if (current_scan >= int(size()-2)) throw NoSuccessor(__FILE__, __LINE__, "getNextRt()", current_scan);
		
			// determine start and end of the next scan
			PeakIterator scan_begin = (*this)[current_scan+1];
	 		PeakIterator scan_end   = (*this)[current_scan+2];	// Seems to be dangerous, but does work.
	 		
			return searchInScan_(scan_begin,scan_end,current_mz);
			
		}
		
		/**
			@brief Moves to the previous scan.
		    
	    Retrieves the peak in the previous scan whose m/z is closest
	    to @p peak.
		*/
		PeakIterator getPrevRt(DRawDataPoint<2> const & peak) const throw (NoSuccessor)
		{			
			CoordinateType current_mz = peak.getPosition()[MZ];
			CoordinateType current_rt   = peak.getPosition()[RT];
			
			int const current_scan    = getRank(current_rt);
			// if we are already in the first scan, there will be no predeccessor....
			if (current_scan == 0) throw NoSuccessor(__FILE__, __LINE__, "getPrevRt()", current_scan);
		
			// determine start and end of the next scan
			PeakIterator scan_begin = (*this)[current_scan-1];
			PeakIterator scan_end   = (*this)[current_scan];
		
			// binary search
			return searchInScan_(scan_begin,scan_end,current_mz);
		}
		     
    /// Returns the scan number of retention time @p rt (starting with 0). 
		typename ScanPositionContainerType::size_type getRank ( CoordinateType const & rt ) const throw()
		{
			return std::lower_bound
				( scan_position_.begin(),
					scan_position_.end(),
					rt
				) - scan_position_.begin();
		}
		
		template<typename Peak >
		class MZless
		{
				typedef Peak PeakType;

				public : 
					/// Check if comparison is done increasing or decreasing.
					bool operator () ( CoordinateType const & left, CoordinateType const &  right ) const throw()
					{
						return (left < right );
					}
					
					bool operator () (PeakType const &  left, PeakType const& right ) const throw()
					{
						return (left.getPosition()[0] < right.getPosition()[0]  );
					}
					
					bool operator () (PeakType const & left, CoordinateType const & right ) const throw()
					{
						return (left.getPosition()[0] < right );
					}
					
					bool operator () (CoordinateType const &   left, PeakType const &  right ) const throw()
					{
						return (left < right.getPosition()[0] );
					}
			};
		
		
		protected:	
		
		/// Performs binary search on an iterator range to find the
		/// peak with the m/z coordinate that comes closest to the starting peak.
		PeakIterator searchInScan_(PeakIterator scan_begin, 
															 PeakIterator scan_end , 
															 CoordinateType current_mz) const
		{
			// perform binary search to find the neighbour in rt dimension
			PeakIterator insert_iter = std::lower_bound(scan_begin,scan_end,current_mz,MZless<PeakType>());	
			
			// In some cases (i.e. for picked peaks) the binary search does
			// not give us the correct position. Therefore we check the peaks
			// to our left and right hand side and verify which one has a m/z
			// closer to the original peak.
			if ( insert_iter == scan_end ) // only one choice
			{
					//std::cout << "At end of scan" << std::endl;
					return --insert_iter;					
			}
			else
			{
				// if the found peak is at the beginning of the spectrum,
				// there is not much we can do.
				if ( insert_iter == scan_begin ) 
				{
					//std::cout << "At Begin of scan" << std::endl;
					return insert_iter;
				}
				else // see if the next smaller one fits better
				{
					CoordinateType delta_mz = (insert_iter->getPosition()[0] - current_mz);
					--insert_iter;
// 					std::cout << "Testing next one." << std::endl;	
// 					std::cout << "delta_mz : " << delta_mz << std::endl;	
// 					std::cout << "Diff to right " << (current_mz - insert_iter->getPosition()[0]) << std::endl;
					if ( (current_mz - insert_iter->getPosition()[0]) < delta_mz )
					{
// 						std::cout << "Returning right pea: " << delta_mz << std::endl;		
						return insert_iter; // peak to the right is closer (in m/z dimension)
					}
					else
					{
// 						std::cout << "Returning left peak." << delta_mz << std::endl;		
						return ++insert_iter;    // peak to the left is closer
					}
				}
			}
		} // end of searchInScan_ 
		
		/// @brief Records the retention time for a given scan
		ScanPositionContainerType scan_position_;  	

	};
	  
} // namespace OpenMS

#endif //  OPENMS_DATASTRUCTURES_SCANINDEX_H
