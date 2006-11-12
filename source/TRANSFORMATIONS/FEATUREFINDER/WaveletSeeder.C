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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/WaveletSeeder.h>


#include <iostream>

namespace OpenMS
{

	WaveletSeeder::WaveletSeeder(): 
		BaseSeeder(), is_initialized_(false)
	{
		name_ = WaveletSeeder::getName();
		
		defaults_.setValue("rtvotes_cutoff",6);
		defaults_.setValue("wt_cutoff",0.0);
		defaults_.setValue("score_cutoff",0.0);

		param_ = defaults_;
	}
	
	WaveletSeeder::~WaveletSeeder(){}
	

  Index WaveletSeeder::nextSeed() throw (NoSuccessor)
	{
		if (!is_initialized_) 
		{
			std::cout << "Initializing wavelet seeder..." << std::endl;
		
			/// Very ugly. Should be removed in final version (ost)
			DPeakArray<2, DRawDataPoint<2> > peaks_ = traits_->getAllPeaks();	
			scan_index_ = traits_->getScanIndex();
			MSExperiment<DRawDataPoint<2>  > exp;
			
			exp.set2DData(peaks_);
		
			std::cout << "Copying of data finished..." << std::endl;
		
			IsotopeFinder<MSExperiment<DRawDataPoint<2>  > > finder(exp);
		
			int votes_cutoff        = param_.getValue("rtvotes_cutoff");
			double wt_cutoff      = param_.getValue("wt_cutoff");
			double score_cutoff = param_.getValue("score_cutoff");
		
			std::cout << "RT votes cut off: " << votes_cutoff << std::endl;
				
			// setting params
			finder.setWtCutOff (wt_cutoff);			// threshold for intensities in wavelet transform
   		finder.setScoreCutOff (score_cutoff);		// scores are ignored
   		finder.setRTVotesCutOff (votes_cutoff); 	// we need isotopic patterns in at least six consecutive scans
				
			std::cout << "Starting detection: " << std::endl;
		
			hash_ = finder.findFeatures(0, (exp.size()-1), true);
		
			hash_iter = hash_.begin();		
			is_initialized_ = true;	
		
			avMZSpacing_ = finder.getAvMZSpacing();	
		
			exp.updateRanges();
			min_mass_ = exp.getMin().Y();
		 
			// delete temporary map
			exp.clear();		
			
		} // end if(....)
		
		if (hash_iter == hash_.end() || hash_.size() == 0 ) throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);
	 	
		
		
		// compute mass we are searching for
		CoordinateType mass_to_find = min_mass_ + (hash_iter->first-1)*avMZSpacing_;
		IndexSet region_;		
	
		// check all scans that support this isotopic pattern
		for (std::list<double>::const_iterator iter_cl2 = hash_iter->second.first.begin(); 
		  		iter_cl2 != hash_iter->second.first.end(); 
		   		++iter_cl2)
				{
				
				CoordinateType rt_to_find = *iter_cl2;
				
// 				std::cout << "Searching for rt: " << rt_to_find << std::endl;
				unsigned int current_scan = scan_index_.getRank(rt_to_find);
				
				if (current_scan >= (scan_index_.size()-1) )
				{
// 					std::cout << "Wrong scan number:" << current_scan  << std::endl;
					break;
				}
				
// 				std::cout << "Searching for " << current_scan << std::endl;
// 				std::cout << "Scan index " << scan_index_.size() << std::endl;
							
				PeakIterator scan_begin = scan_index_[current_scan];
				PeakIterator scan_end   = scan_index_[current_scan+1];				
					
				PeakIterator insert_iter = std::lower_bound(scan_begin,scan_end,mass_to_find,MZless());	
				int peak_index = (insert_iter - traits_->getAllPeaks().begin());		
	
// 				std::cout << "Adding peak at mass " << traits_->getPeakMz(peak_index) << std::endl;
				
				CoordinateType miso_mass = traits_->getPeakMz(peak_index);
				CoordinateType miso_rt       =	traits_->getPeakRt(peak_index);			
				
				// walk a bit to the left from the monoisotopic peak
				for (int p=0; p<=10; ++p)
				{
					if ( (peak_index - p) > 0 
					     && traits_->getPeakFlag( (peak_index - p) ) == FeaFiTraits::UNUSED
						 && traits_->getPeakRt( (peak_index - p) ) == miso_rt
						 && traits_->getPeakMz( peak_index) - traits_->getPeakMz( (peak_index - p) ) < 2  )
					{
					 	region_.add( (peak_index - p) );
					 	traits_->getPeakFlag(  (peak_index - p) ) = FeaFiTraits::INSIDE_FEATURE;
					}
				}
		
				CoordinateType mass_distance = 0;							
				int nr_peaks = traits_->getNumberOfPeaks();
				
				//we added the first peak (hopefully the monoisotopic one), now
				// we want to walk for about 10 Thompson (or Dalton or whatever)
				// into positive  (increasing) m/z direction
				while (mass_distance < 7
				           && peak_index < (nr_peaks-2)
						   && traits_->getPeakRt(peak_index) == miso_rt)
				{
					++peak_index;
					if (traits_->getPeakFlag(peak_index) == FeaFiTraits::UNUSED) region_.add(peak_index);
// 					std::cout << "Adding peak " << peak_index << std::endl;
					mass_distance = (traits_->getPeakMz(peak_index+1) - miso_mass);
// 					std::cout << "Current mass distance : " << mass_distance << std::endl;
				} 
// 				std::cout << "This scan is done." << std::endl;
								
	}	// for (std::list...)
	
// 	std::cout << "Done Size of region: " << region_.size() << std::endl;	
	
	region_.sort();
	 
	 ++hash_iter;
	
	 // find maximum intensity
	 UnsignedInt max_index           = 0;
	 CoordinateType max_intensity = 0.0;
	 
	// find maximum of region
	for (IndexSet::const_iterator citer = region_.begin();
	      citer != region_.end(); 
		  ++citer)
	{
		if (traits_->getPeakIntensity(*citer) > max_intensity) 
		{
			max_index = *citer;
			max_intensity = traits_->getPeakIntensity(*citer);
		}	
	
	} // end of for (Index in region)
	traits_->getPeakFlag(max_index) = FeaFiTraits::SEED;
	
	return max_index;
}


} // end of namespace OpenMS

