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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>

namespace OpenMS
{

	SimpleExtender::SimpleExtender() 
	: BaseExtender(),
		first_seed_seen_(false), 
		intensity_threshold_(0),
		last_pos_extracted_(),
		intensity_sum_(0)
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
	{
	}

	const FeaFiModule::IndexSet& SimpleExtender::extend(const IndexSet& seed_region)
	{
    if (!first_seed_seen_)
    {
			CoordinateType tol_rt = param_.getValue("tolerance_rt");
			CoordinateType tol_mz = param_.getValue("tolerance_mz");
			intensity_factor_ = param_.getValue("intensity_factor");

			dist_mz_up_ = param_.getValue("dist_mz_up");
			dist_mz_down_ = param_.getValue("dist_mz_down");
			dist_rt_up_ = param_.getValue("dist_rt_up");
			dist_rt_down_ = param_.getValue("dist_rt_down");

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
		IDX seed_index;
		
		// find maximum of region and set seed flag
    for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
    {	
    	FeaFiTraits::MapType::SpectrumType& scan = traits_->getData()[citer->first];
      if (scan[citer->second].getIntensity() > max_intensity)
      {
        seed_index     = *citer;
        max_intensity = scan[citer->second].getIntensity();						
			}
			
  		IndexWithPriority iwp(*citer,computePeakPriority_(*citer));
			boundary_.push(iwp);
    }
    
    traits_->getPeakFlag(seed_index) = FeaFiTraits::SEED;
    		
		 // remember last peak to be extracted from the boundary
    // (in this case the seed !)
		IDX seed = seed_index;
		last_pos_extracted_[RT] = traits_->getData()[seed_index.first].getRetentionTime();
		last_pos_extracted_[MZ] = traits_->getData()[seed_index.first][seed_index.second].getPos();
		
		std::cout << "New seed at " << seed.first << " " << seed.second << std::endl;
		
		intensity_threshold_ = intensity_factor_ * traits_->getPeakIntensity(seed);
		intensity_sum_        = traits_->getPeakIntensity(seed);

    ProbabilityType prior = computePeakPriority_(seed);
    IndexWithPriority seed_p(seed_index,prior);

    // at the beginning, the boundary contains only the seed
    boundary_.push(seed_p);
    priorities_[seed_index] = prior;

    while (!boundary_.empty())
    {
			
			// remove peak with highest priority
			IndexWithPriority const index_priority = boundary_.top();
			boundary_.pop();
			
			const IDX  current_index = index_priority.index;
			if (current_index.first >= traits_->getData().size() || current_index.second >= traits_->getData()[current_index.first].size())
			{
				std::cout << "WARNING: Index too large " << current_index.first << " " << current_index.second << std::endl;
				continue;
			}
			
			FeaFiTraits::MapType::SpectrumType& scan = traits_->getData()[current_index.first];
			
			// skip this point if its intensity is too large
			if (scan[current_index.second].getIntensity() <  intensity_threshold_)	continue;
				
			if (scan[current_index.second].getIntensity() < (intensity_sum_ * min_intensity_contribution_) )
			{
// 				std::cout << "Skipping point because of low intensity contribution. " << std::endl;
// 				std::cout << current_peak.getIntensity() << " " << (intensity_sum_ * min_intensity_contribution_) << std::endl;
				continue;			 
			}
			
			// remember last peak to be extracted from the boundary
			last_pos_extracted_[RT] = scan.getRetentionTime();
			last_pos_extracted_[MZ] = scan[current_index.second].getPos();

			// Now we explore the neighbourhood of the current peak. Peaks in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.
			
			// Add position to the current average of positions weighted by intensity
			running_avg_.add(last_pos_extracted_,scan[current_index.second].getIntensity());

			// explore neighbourhood of current peak
			moveMzUp_(current_index);
			moveMzDown_(current_index);
			moveRtUp_(current_index);
			moveRtDown_(current_index);

			// check peak flags (if data point is already inside a feature region or used as seed, we discard it)
			if (traits_->getPeakFlag(current_index) == FeaFiTraits::SEED || traits_->getPeakFlag(current_index) == FeaFiTraits::UNUSED)
			{
				traits_->getPeakFlag(current_index) = FeaFiTraits::INSIDE_FEATURE;
				region_.insert(current_index);
				intensity_sum_ += scan[current_index.second].getIntensity();
// 				std::cout << "Added point to region. Intensity sum is now: " << intensity_sum_ << std::endl;
// 				std::cout << "Intensity of the added point is : " << current_peak.getIntensity() << std::endl;
			}
    } // end of while ( !boundary_.empty() )

    std::cout << "SimpleExtender: Feature region size: " << region_.size() << std::endl;
    
    return region_;

	} // end of extend(Unsigned int seed_index)

	/**
	 * \brief Checks whether the current peak is to far away from the seed
	 * 
	 * */
	bool SimpleExtender::isTooFarFromCentroid_(const IDX& index)
	{
		// error check, this should not happen
		if (index.first >= traits_->getData().size() || index.second >= traits_->getData()[index.first].size())
		{
			std::cout << "WARNING: Index too large " << index.first << " " << index.second << " (SimpleExtender::isTooFarFromCentroid_)" << std::endl;
			return false;
		}
	
		FeaFiTraits::MapType::SpectrumType& scan = traits_->getData()[index.first];
	
    FeaFiTraits::PositionType2D const curr_mean = running_avg_.getPosition();

    if ( scan[index.second].getPos() > curr_mean[MZ] + dist_mz_up_   ||
				 scan[index.second].getPos() < curr_mean[MZ] - dist_mz_down_ ||
				 scan.getRetentionTime() > curr_mean[RT] + dist_rt_up_   ||
				 scan.getRetentionTime() < curr_mean[RT] - dist_rt_down_ )
    {
			// close enough
			return true;
    }
		
		// too far away from centroid of region
    return false;
	}

	void SimpleExtender::moveMzUp_(const IDX& index)
	{
    try // moving up in m/z direction
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getNextMz(tmp); // take next peak
				
				// stop if we've left the current scan
				if (isTooFarFromCentroid_(tmp) ) break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(tmp);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveMzDown_(const IDX& index)
	{
    try // moving down in m/z direction
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getPrevMz(tmp);	// take next peak
				
				// stop if we've left the current scan
				if (isTooFarFromCentroid_(tmp))	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(tmp);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveRtUp_(const IDX& index)
	{
    try // moving up in retention time
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getNextRt(tmp); // take next peak

				if (isTooFarFromCentroid_(tmp)) 	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(tmp);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveRtDown_(const IDX& index)
	{
    try // moving up in retention time
    {
			IDX tmp = index;
			while (true)
			{
				traits_->getPrevRt(tmp); // take next peak

				if (isTooFarFromCentroid_(tmp)) 	break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(tmp);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {
    }
	}

	SimpleExtender::ProbabilityType SimpleExtender::computePeakPriority_(const IDX& index)
	{
		return traits_->getData()[index.first][index.second].getIntensity() *
 			score_distribution_rt_.value(traits_->getData()[index.first].getRetentionTime()-last_pos_extracted_[RT]) *
			score_distribution_mz_.value(traits_->getData()[index.first][index.second].getPos()-last_pos_extracted_[MZ]);
	}

	void SimpleExtender::checkNeighbour_(const IDX& index)
	{
		// error check, this should not happen
		if (index.first >= traits_->getData().size() || index.second >= traits_->getData()[index.first].size())
		{
			std::cout << "WARNING: Index too large " << index.first << " " << index.second << " (SimpleExtender::checkNeighbour_)" << std::endl;
			return;
		}
	
		FeaFiTraits::MapType::SpectrumType& scan = traits_->getData()[index.first];
		
    // we don't care about points with intensity zero (or < 0 which might occur if TopHatFilter was applied)
    if (scan[index.second].getIntensity() <= 0) 	return;

    if (traits_->getPeakFlag(index)== FeaFiTraits::UNUSED)
    {
			double pr_new = computePeakPriority_(index);

			if (pr_new > priority_threshold_) // check if priority larger than threshold
			{
				std::map<IDX, double>::iterator piter = priorities_.find(index);
				if (piter == priorities_.end()) // not yet in boundary
				{
					priorities_[index] = pr_new;
					IndexWithPriority peak(index,pr_new);
					boundary_.push(peak);	// add to boundary
				}
			}
		}
	}


} // end of class SimpleExtender
