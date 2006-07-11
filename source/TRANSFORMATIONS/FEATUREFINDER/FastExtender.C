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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FastExtender.h>

#include<OpenMS/SYSTEM/StopWatch.h>

namespace OpenMS
{

	FastExtender::FastExtender() : BaseExtender(),
	first_seed_seen_(false), last_extracted_(0), 
	nr_peaks_seen_(0), intensity_sum_(0), moving_avg_(), 
	last_avg_(0)
	{
		name_ = FastExtender::getName();
		defaults_.setValue("tolerance_rt",2.0f);
		defaults_.setValue("tolerance_mz",0.5f);
		defaults_.setValue("dist_mz_up",6.0f);
		defaults_.setValue("dist_mz_down",2.0f);
		defaults_.setValue("dist_rt_up",5.0f);
		defaults_.setValue("dist_rt_down",5.0f);
		defaults_.setValue("intensity_factor",0.001f);
		defaults_.setValue("intensity_avg_tolerance",500.0f);
		defaults_.setValue("extension_baseline",300.0f);
		defaults_.setValue("moving_avg_size",10);
				
		param_ = defaults_;
		
		std::cout << "Starting FastExtender..." << std::endl;
		
		init_();
	}

	FastExtender::~FastExtender(){}
	
	void FastExtender::init_()
	{
		float tol_rt       = param_.getValue("tolerance_rt");
		float tol_mz       = param_.getValue("tolerance_mz");
		intensity_factor_  = param_.getValue("intensity_factor");
		intensity_avg_tol_ = param_.getValue("intensity_avg_tolerance");
			
		dist_mz_up_   = param_.getValue("dist_mz_up");
		dist_mz_down_ = param_.getValue("dist_mz_down");
		dist_rt_up_   = param_.getValue("dist_rt_up");
		dist_rt_down_ = param_.getValue("dist_rt_down");
		
		extension_baseline_ = param_.getValue("extension_baseline");			
		
		// initialise distributions to interpolate the priority
		score_distribution_rt_.getData().push_back(1.0);
		score_distribution_rt_.setScale(tol_rt);
		score_distribution_rt_.setOffset(0);
		
		score_distribution_mz_.getData().push_back(1.0);
		score_distribution_mz_.setScale(tol_mz);
		score_distribution_mz_.setOffset(0);
				
	}

  const IndexSet& FastExtender::extend(const UnsignedInt seed_index)
  {  	
		// empty region and boundary datastructures
		region_.clear();
		while (boundary_.size() > 0) boundary_.pop();
		priorities_.clear();
		running_avg_.clear();
		
		nr_peaks_seen_  = 0;
		intensity_sum_    = 0;
		last_avg_            = 0;
				
		// remember last peak to be extracted from the boundary
		// in this case, the seed
		last_extracted_ = seed_index;
		
		float prior = computePeakPriority_(seed_index);
		IndexWithPriority seed(seed_index,prior);
		
		// at the beginning, the boundary contains only the seed
		boundary_.push(seed);
		priorities_[seed_index] = prior;
			
		while (!boundary_.empty())
		{			
			// remove peak with highest priority
			IndexWithPriority const index_priority = boundary_.top();
			boundary_.pop();
			
			++nr_peaks_seen_;
									
			UnsignedInt const current_index       = index_priority.index;
			IntensityType const current_intensity = traits_->getPeakIntensity(current_index);
						
			// remember last peak to be extracted from the boundary
			last_extracted_ = current_index;
		
			double intensity_thr = intensity_factor_ * intensity_sum_;
			if (current_intensity < intensity_thr && current_intensity < extension_baseline_) 
			{
				std::cout << "Skipping peak : current intensity: " << current_intensity << " threshold: " << intensity_thr << std::endl;	
				continue;
			}				
			
			IntensityType intens_avg = 0;
			for (unsigned int i=0; i< moving_avg_.size(); ++i)
			{
				intens_avg += moving_avg_.at(i); 	
			}
			if ( moving_avg_.size() > 0) 	intens_avg /= moving_avg_.size();
			else intens_avg = current_intensity;
						
			if (intens_avg < (last_avg_ + intensity_avg_tol_) 
			    &&  intens_avg > (last_avg_ - intensity_avg_tol_) 
			    && current_intensity < extension_baseline_)
			{
				std::cout << "Current intensity average: " << intens_avg << " last avg: " << last_avg_ << std::endl;
				std::cout << "Skipping peak because of insufficient change in moving avg. " << std::endl;
				continue;
			}
			last_avg_ = intens_avg;
		
			// Now we explore the neighbourhood of the current peak. Peaks in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.
			
			// get position of current peak
			FeaFiTraits::PositionType const curr_pos = traits_->getPeak(current_index).getPosition();
			// and add it to the current average of positions weighted by intensity
			running_avg_.add(curr_pos);
			
			// explore neighbourhood of current peak
			moveMzUp_(current_index);
			moveMzDown_(current_index);
			moveRtUp_(current_index);
			moveRtDown_(current_index);
			
			// set flag
			if (traits_->getPeakFlag(current_index) != FeaFiTraits::SEED)
				traits_->getPeakFlag(current_index) = FeaFiTraits::INSIDE_FEATURE;
						
			region_.add(current_index);
			
			intensity_sum_ += current_intensity;
			
			// update moving average of collected intensities
			if (moving_avg_.size() == (unsigned int)param_.getValue("moving_avg_size") )
			{
				moving_avg_.erase(moving_avg_.begin());
			}
			moving_avg_.push_back(current_intensity);
		
		} // end of while ( !boundary_.empty() )		
		
		std::cout << "Feature region size: " << region_.size() << std::endl;
		
		region_.sort();		
		return region_;								
		
	} // end of extend(Unsigned int seed_index)

	/**
	 * \brief Checks whether the current peak is to far away from the seed
	 * 
	 * */
	bool FastExtender::isTooFarFromCentroid_(UnsignedInt current_peak) 
	{
				
		CoordinateType const current_mz   = traits_->getPeakMz(current_peak);
		CoordinateType const current_rt   = traits_->getPeakRt(current_peak);
		
		FeaFiTraits::PositionType const curr_mean = running_avg_.getPosition();
		
		if ( current_mz > curr_mean[MZ] + dist_mz_up_   ||
				 current_mz < curr_mean[MZ] - dist_mz_down_ ||
				 current_rt > curr_mean[RT] + dist_rt_up_   ||
				 current_rt < curr_mean[RT] - dist_rt_down_ )
			{
				return true;
				std::cout << "FastExtender: Too far from centroid !! " << std::endl;
			}
	
		return false;
	}	
	
	void FastExtender::moveMzUp_(UnsignedInt current_index)
	{
		
		UnsignedInt current_scan = traits_->getPeakScanNr(current_index);

		try // moving up in m/z direction
			{			
				while (true) 
				{			
					current_index = traits_->getNextMz(current_index); // take next peak
						
					// stop if we've left the current scan
					if (current_scan != traits_->getPeakScanNr(current_index)
					    || isTooFarFromCentroid_(current_index) ) break;
					
					// check this neighbour for insertion into the boundary
					checkNeighbour_(current_index);
														
				} // end of while (true)
				 
		} catch(NoSuccessor) {}
	}
	
	void FastExtender::moveMzDown_(UnsignedInt current_index)
	{
		
		UnsignedInt current_scan = traits_->getPeakScanNr(current_index);
		
		try // moving down in m/z direction
			{			
				while (true) 
				{			
					current_index    = traits_->getPrevMz(current_index);	// take next peak
									
					// stop if we've left the current scan
					if (current_scan != traits_->getPeakScanNr(current_index)
					    || isTooFarFromCentroid_(current_index)) break;
					
					// check this neighbour for insertion into the boundary
					checkNeighbour_(current_index);
					
				} // end of while (true)
				 
		} catch(NoSuccessor) {}
	}
	
	void FastExtender::moveRtUp_(UnsignedInt current_index)
	{
				
		try // moving up in retention time
			{			
				while (true) 
				{				
					current_index = traits_->getNextRt(current_index); // take next peak
					
					if (isTooFarFromCentroid_(current_index)) break;
										
					// check this neighbour for insertion into the boundary
					checkNeighbour_(current_index);
					
				} // end of while (true)
				 
		} catch(NoSuccessor) {}
	}
	
	void FastExtender::moveRtDown_(UnsignedInt current_index)
	{
				
		try // moving up in retention time
			{			
								
				while (true) 
				{				
					current_index = traits_->getPrevRt(current_index); // take next peak
					
					if (isTooFarFromCentroid_(current_index)) break;

					// check this neighbour for insertion into the boundary
					checkNeighbour_(current_index);
					
				} // end of while (true)
				 
		} catch(NoSuccessor) {}
	}
		
	double FastExtender::computePeakPriority_(UnsignedInt current_peak) 
	{
						
		IntensityType curr_intens = traits_->getPeakIntensity(current_peak);
		
		CoordinateType curr_mz = traits_->getPeakMz(current_peak);
		CoordinateType last_mz = traits_->getPeakMz(last_extracted_);
		
		CoordinateType curr_rt = traits_->getPeakRt(current_peak);
		CoordinateType last_rt = traits_->getPeakRt(last_extracted_);
				
		return curr_intens *
					 score_distribution_rt_.value(curr_rt-last_rt) *
					 score_distribution_mz_.value(curr_mz-last_mz);
				
	}
		
	void FastExtender::checkNeighbour_(UnsignedInt current_index)
	{
		// we don't care about points with intensity zero
		if (traits_->getPeakIntensity(current_index) == 0) return;
			
		if (traits_->getPeakFlag(current_index)== FeaFiTraits::UNUSED) 
		{
			// is this point allready contained in the boundary??
			std::map<UnsignedInt, double>::iterator piter = priorities_.find(current_index);
			double pr_new = computePeakPriority_(current_index); 
													
			if (piter == priorities_.end()) // not yet in boundary
			{
					priorities_[current_index] = pr_new;
					IndexWithPriority peak(current_index,pr_new);
					boundary_.push(peak);	// add to boundary
			} 
			else std::cout << "Already in boundary !! " << std::endl;
			
		} // end if traits_->...		
	}
	

} // end of class FastExtender
