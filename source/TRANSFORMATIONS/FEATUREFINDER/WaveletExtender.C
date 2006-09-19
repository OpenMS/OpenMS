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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/WaveletExtender.h>

namespace OpenMS
{

WaveletExtender::WaveletExtender() 
	: BaseExtender(), is_initialized_(false)
{
    name_ = WaveletExtender::getName();
		
    param_ = defaults_;
}

WaveletExtender::~WaveletExtender()
{}

const IndexSet& WaveletExtender::extend(const UnsignedInt /*seed_index*/)
{
	
	if (!is_initialized_)
	{	
		std::cout << "Starting WaveletExtender..." << std::endl;
		
		/// Very ugly. Should be removed in final version (ost)
		DPeakArrayNonPolymorphic<2, DRawDataPoint<2> > peaks_ = traits_->getAllPeaks();	
		scan_index_ = traits_->getScanIndex();
		MSExperiment<DRawDataPoint<2>  > exp;
		exp.set2DData(peaks_);
		
		std::cout << "Copying of data finished..." << std::endl;
		
		IsotopeFinder<MSExperiment<DRawDataPoint<2>  > > finder;
		finder.setData(exp);
				
		// setting params
		finder.setWtCutOff (0);			// threshold for intensities in wavelet transform
   		finder.setScoreCutOff (0);		// scores are ignored
   		finder.setRTVotesCutOff (3); 	// we need isotopic patterns in at least six consecutive scans
		
		std::cout << "Starting detection: " << std::endl;
		
		hash_ = finder.findFeatures(7, (exp.size()-1), true);
		
		hash_iter = hash_.begin();		
		is_initialized_ = true;	
		
		 avMZSpacing_ = finder.getAvMZSpacing();	
		
		 exp.updateRanges();
		 min_mass_ = exp.getMin().Y();
		 
		 exp.clear();
	}
	
	region_.clear();	// empty last region 
	
	if (hash_iter == hash_.end() || hash_.size() == 0 ) throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);
	 
	 std::cout << "m/z range: ";
	 std::cout << (min_mass_ + (hash_iter->first-1)*avMZSpacing_) << " ";
	 std::cout << (min_mass_ + (hash_iter->first)*avMZSpacing_) << " " << std::endl;
			
	CoordinateType mass_to_find = min_mass_ + (hash_iter->first-1)*avMZSpacing_;
	std::cout << "I am searching for m/z : " << 	mass_to_find << std::endl;
	
	// check all scans that support this isotopic pattern
	for (std::list<double>::const_iterator iter_cl2 = hash_iter->second.first.begin(); 
		  iter_cl2 != hash_iter->second.first.end(); 
		   ++iter_cl2)
	{
				CoordinateType rt_to_find = *iter_cl2;
				
				std::cout << "Searching for rt: " << rt_to_find << std::endl;
				
				unsigned int current_scan = scan_index_.getRank(rt_to_find);
				
				if (current_scan >= (scan_index_.size()-1) )
				{
					std::cout << "Wrong scan number:" << current_scan  << std::endl;
					break;
				}
				
// 				std::cout << "Searching for " << current_scan << std::endl;
// 				std::cout << "Scan index " << scan_index_.size() << std::endl;
							
				PeakIterator scan_begin = scan_index_[current_scan];
				PeakIterator scan_end   = scan_index_[current_scan+1];				
					
				PeakIterator insert_iter = std::lower_bound(scan_begin,scan_end,mass_to_find,MZless());	
				UnsignedInt peak_index = (insert_iter - traits_->getAllPeaks().begin());		
	
				std::cout << "Adding peak at mass " << traits_->getPeakMz(peak_index) << std::endl;
				
				if ( (peak_index-1 ) >= 0) region_.add(peak_index-1);
				std::cout << "Adding peak at mass " << traits_->getPeakMz(peak_index-1) << std::endl;
				if ( (peak_index-2 ) >= 0) region_.add(peak_index-2);
				std::cout << "Adding peak at mass " << traits_->getPeakMz(peak_index-2) << std::endl;
						
				region_.add(peak_index);
				
				CoordinateType mass_distance = 0;
				CoordinateType miso_mass      = traits_->getPeakMz(peak_index);
				
				unsigned int nr_peaks = traits_->getNumberOfPeaks();
				//we added the first peak (hopefully the monoisotopic one), now
				// we want to walk for about 10 Thompson (or Dalton or whatever)
				// into positive  (increasing) m/z direction
				
				while (mass_distance < 10 && peak_index < nr_peaks)
				{
					++peak_index;
					region_.add(peak_index);
					std::cout << "Adding peak " << peak_index << std::endl;
					
					mass_distance = (traits_->getPeakMz(peak_index) - miso_mass);
					std::cout << "Current mass distance : " << mass_distance << std::endl;
				} 
				std::cout << "This scan is done." << std::endl;
								
	}	// for (std::list...)
			
	 ++hash_iter;
	
	 std::cout << "Extension done. Size of region: " << region_.size() << std::endl;
	return region_;
}

} // end of namespace OpenMS

