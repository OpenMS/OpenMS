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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>

namespace OpenMS
{

	SimpleExtender::SimpleExtender() : BaseExtender(),
																		 first_seed_seen_(false), intensity_threshold_(0),
																		 last_pos_extracted_(), intensity_sum_(0)
	{
    name_ = SimpleExtender::getName();
    defaults_.setValue("tolerance_rt",2.0f);
    defaults_.setValue("tolerance_mz",0.5f);
    defaults_.setValue("dist_mz_up",6.0f);
    defaults_.setValue("dist_mz_down",2.0f);
    defaults_.setValue("dist_rt_up",5.0f);
    defaults_.setValue("dist_rt_down",5.0f);
    defaults_.setValue("priority_thr",0.0f);
    defaults_.setValue("intensity_factor",0.03f);
		defaults_.setValue("min_intensity_contribution",0.01f);

    param_ = defaults_;
	}

	SimpleExtender::~SimpleExtender()
	{}

	const IndexSet& SimpleExtender::extend(const IndexSet& seed_region)
	{
    if (!first_seed_seen_)
    {
			CoordinateType tol_rt      = param_.getValue("tolerance_rt");
			CoordinateType tol_mz   = param_.getValue("tolerance_mz");
			intensity_factor_            = param_.getValue("intensity_factor");

			dist_mz_up_     = param_.getValue("dist_mz_up");
			dist_mz_down_ = param_.getValue("dist_mz_down");
			dist_rt_up_        = param_.getValue("dist_rt_up");
			dist_rt_down_    = param_.getValue("dist_rt_down");

			priority_threshold_              = param_.getValue("priority_thr");
			min_intensity_contribution_ = param_.getValue("min_intensity_contribution");

			// initialise priority distributions
			score_distribution_rt_.getData().push_back(1.0);
			score_distribution_rt_.setScale(tol_rt);
			score_distribution_rt_.setOffset(0);

			score_distribution_mz_.getData().push_back(1.0);
			score_distribution_mz_.setScale(tol_mz);
			score_distribution_mz_.setOffset(0);
			
			first_seed_seen_ = true;			
    }

    // empty region and boundary datastructures
    region_.clear();
		priorities_.clear();
    running_avg_.clear();
		while (boundary_.size() > 0) boundary_.pop();
		
		CoordinateType max_intensity = 0.0;
		UnsignedInt seed_index = 0;
		
		// find maximum of region and set seed flag
    for (IndexSet::const_iterator citer = seed_region.begin();
            citer != seed_region.end();
            ++citer)
    {
				PeakType p	 = traits_->getPeak(*citer);
				
        if (p.getIntensity() > max_intensity)
        {
            seed_index     = *citer;
            max_intensity = p.getIntensity();						
				}
				
				ProbabilityType pi = computePeakPriority_(p);
    		IndexWithPriority iwp(*citer,pi);
				
				boundary_.push(iwp);

    } // end of for (Index in region)
    traits_->getPeakFlag(seed_index) = FeaFiTraits::SEED;
    		
		 // remember last peak to be extracted from the boundary
    // (in this case the seed !)
		PeakType seed        = traits_->getPeak(seed_index);
		last_pos_extracted_  = seed.getPosition();
		
		intensity_threshold_ = intensity_factor_ * seed.getIntensity();
		intensity_sum_        = seed.getIntensity();

    ProbabilityType prior = computePeakPriority_(seed);
    IndexWithPriority seed_p(seed_index,prior);

    // at the beginning, the boundary contains only the seed
    boundary_.push(seed_p);
    priorities_[seed_index] = prior;

    while (!boundary_.empty())
    {
// 			std::cout << "Size of boundary: " << boundary_.size() << std::endl;
// 			std::cout << "Size of region: " << region_.size() << std::endl;
			
			// remove peak with highest priority
			IndexWithPriority const index_priority = boundary_.top();
			boundary_.pop();
			
			const UnsignedInt  current_index = index_priority.index;
			const PeakType current_peak  		= traits_->getPeak(current_index);
			
			// skip this point if its intensity is too large
			if (current_peak.getIntensity() <  intensity_threshold_)	continue;
				
			if (current_peak.getIntensity() < (intensity_sum_ * min_intensity_contribution_) )
			{
// 				std::cout << "Skipping point because of low intensity contribution. " << std::endl;
// 				std::cout << current_peak.getIntensity() << " " << (intensity_sum_ * min_intensity_contribution_) << std::endl;
				continue;			 
			}
			
			// remember last peak to be extracted from the boundary
			last_pos_extracted_ = current_peak.getPosition();

			// Now we explore the neighbourhood of the current peak. Peaks in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.
			
			// Add position to the current average of positions weighted by intensity
			running_avg_.add(last_pos_extracted_,current_peak.getIntensity());

			// explore neighbourhood of current peak
			moveMzUp_(current_index);
			moveMzDown_(current_index);
			moveRtUp_(current_index);
			moveRtDown_(current_index);

			// check peak flags (if data point is already inside a feature region or used as seed, we discard it)
			if (traits_->getPeakFlag(current_index) == FeaFiTraits::SEED || traits_->getPeakFlag(current_index) == FeaFiTraits::UNUSED)
			{
				traits_->getPeakFlag(current_index) = FeaFiTraits::INSIDE_FEATURE;
				region_.add(current_index);
				intensity_sum_ += current_peak.getIntensity();
// 				std::cout << "Added point to region. Intensity sum is now: " << intensity_sum_ << std::endl;
// 				std::cout << "Intensity of the added point is : " << current_peak.getIntensity() << std::endl;
			}
    } // end of while ( !boundary_.empty() )

    std::cout << "SimpleExtender: Feature region size: " << region_.size() << std::endl;

    region_.sort();
    return region_;

	} // end of extend(Unsigned int seed_index)

	/**
	 * \brief Checks whether the current peak is to far away from the seed
	 * 
	 * */
	bool SimpleExtender::isTooFarFromCentroid_(UnsignedInt current_peak)
	{

		PeakType p = traits_->getPeak(current_peak);
    FeaFiTraits::PositionType const curr_mean = running_avg_.getPosition();

    if ( p.getPosition()[MZ] > curr_mean[MZ] + dist_mz_up_   ||
				 p.getPosition()[MZ] < curr_mean[MZ] - dist_mz_down_ ||
				 p.getPosition()[RT] > curr_mean[RT] + dist_rt_up_   ||
				 p.getPosition()[RT] < curr_mean[RT] - dist_rt_down_ )
    {
			// close enough
			return true;
    }
		
		// too far away from centroid of region
    return false;
	}

	void SimpleExtender::moveMzUp_(UnsignedInt current_index)
	{
    try // moving up in m/z direction
    {
			while (true)
			{
				current_index = traits_->getNextMz(current_index); // take next peak
				
				// stop if we've left the current scan
				if (isTooFarFromCentroid_(current_index) ) break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveMzDown_(UnsignedInt current_index)
	{
    try // moving down in m/z direction
    {
			while (true)
			{
				current_index    = traits_->getPrevMz(current_index);	// take next peak
				
				// stop if we've left the current scan
				if (isTooFarFromCentroid_(current_index))	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveRtUp_(UnsignedInt current_index)
	{
    try // moving up in retention time
    {
			while (true)
			{
				current_index = traits_->getNextRt(current_index); // take next peak

				if (isTooFarFromCentroid_(current_index)) 	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveRtDown_(UnsignedInt current_index)
	{
    try // moving up in retention time
    {

			while (true)
			{
				current_index = traits_->getPrevRt(current_index); // take next peak

				if (isTooFarFromCentroid_(current_index)) 	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	SimpleExtender::ProbabilityType SimpleExtender::computePeakPriority_(PeakType& p)
	{
		return p.getIntensity() *
 			score_distribution_rt_.value(p.getPosition()[RT]-last_pos_extracted_[RT]) *
			score_distribution_mz_.value(p.getPosition()[MZ]-last_pos_extracted_[MZ]);
	}

	void SimpleExtender::checkNeighbour_(UnsignedInt current_index)
	{
		PeakType p = traits_->getPeak(current_index);
    // we don't care about points with intensity zero (or < 0 which might occur if TopHatFilter was applied)
    if (p.getIntensity() <= 0) 	return;

    if (traits_->getPeakFlag(current_index)== FeaFiTraits::UNUSED)
    {
			std::map<UnsignedInt, double>::iterator piter = priorities_.find(current_index);
			double pr_new = computePeakPriority_(p);

			if (piter == priorities_.end()) // not yet in boundary
			{
				if (pr_new > priority_threshold_) // check if priority larger than threshold
				{
					priorities_[current_index] = pr_new;
					IndexWithPriority peak(current_index,pr_new);
					boundary_.push(peak);	// add to boundary
				}
			} // end if (piter ...)
			
		} // end if traits_->...
	}


} // end of class SimpleExtender
