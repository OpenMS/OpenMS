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

using namespace std;

namespace OpenMS
{

	SimpleExtender::SimpleExtender() 
	: BaseExtender(),
		last_pos_extracted_()
	{
    setName(getProductName());
    
    defaults_.setValue("tolerance_rt",2.0f);
    defaults_.setValue("tolerance_mz",0.5f);
    defaults_.setValue("dist_mz_up",6.0f);
    defaults_.setValue("dist_mz_down",2.0f);
    defaults_.setValue("dist_rt_up",5.0f);
    defaults_.setValue("dist_rt_down",5.0f);
    defaults_.setValue("priority_thr",0.0f);
    defaults_.setValue("intensity_factor",0.03f);
		defaults_.setValue("min_intensity_contribution",0.01f);

    defaultsToParam_();
	}

	SimpleExtender::~SimpleExtender()
	{
	}

  SimpleExtender::SimpleExtender(const SimpleExtender& rhs)
    : BaseExtender(rhs)
  {
    updateMembers_();
  }
  
  SimpleExtender& SimpleExtender::operator= (const SimpleExtender& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseExtender::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

  void SimpleExtender::updateMembers_()
  {
		dist_mz_up_     = param_.getValue("dist_mz_up");
		dist_mz_down_ = param_.getValue("dist_mz_down");
		dist_rt_up_       = param_.getValue("dist_rt_up");
		dist_rt_down_   = param_.getValue("dist_rt_down");

		priority_threshold_ = param_.getValue("priority_thr");

		// initialise priority distributions
		if (score_distribution_rt_.getData().size()!=1)
		{
			score_distribution_rt_.getData().push_back(1.0);
			score_distribution_rt_.setOffset(0);
		}
		score_distribution_rt_.setScale(param_.getValue("tolerance_rt"));
		
		if (score_distribution_mz_.getData().size()!=1)
		{
			score_distribution_mz_.getData().push_back(1.0);
			score_distribution_mz_.setOffset(0);	
 		}
 		score_distribution_mz_.setScale(param_.getValue("tolerance_mz"));
  }


	const FeaFiModule::IndexSet& SimpleExtender::extend(const IndexSet& seed_region)
	{
    // empty region and boundary datastructures
    region_.clear();
		priorities_.clear();
    running_avg_.clear();
		boundary_ = priority_queue< IndexWithPriority, vector < IndexWithPriority > , IndexWithPriority::PriorityLess >();
		
		// find maximum of region (seed)
		CoordinateType max_intensity = 0.0;
		IDX seed;
    for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
    {	
      if (traits_->getPeakIntensity(*citer) > max_intensity)
      {
        seed = *citer;
        max_intensity = traits_->getPeakIntensity(seed);						
			}
    }
    traits_->getPeakFlag(seed) = FeaFiTraits::SEED;
    		
		// remember last extracted point (in this case the seed !)
		last_pos_extracted_[RT] = traits_->getPeakRt(seed);
		last_pos_extracted_[MZ] = traits_->getPeakMz(seed);
		
		// Add peaks in the region to boundary
    for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
    {	
    	ProbabilityType priority = computePeakPriority_(*citer);
    	priorities_[*citer] = priority;
			boundary_.push(IndexWithPriority(*citer,priority));
    }
		
		cout << "Extending from " << traits_->getPeakRt(seed) << "/" << traits_->getPeakMz(seed); 
		cout << " (" << seed.first << "/" << seed.second << ")" << endl;
		
		//compute intensity threshold and sum
		intensity_threshold_ = (double)param_.getValue("intensity_factor") * traits_->getPeakIntensity(seed);
		IntensityType intensity_sum = 0.0;
		IntensityType min_intensity_contribution = param_.getValue("min_intensity_contribution");
				
    while (!boundary_.empty())
    {
			// remove peak with highest priority
			const IDX  current_index = boundary_.top().index;
			boundary_.pop();
			
    	//Corrupt index
    	OPENMS_PRECONDITION(current_index.first<traits_->getData().size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(current_index.second<traits_->getData()[current_index.first].size(), "Peak index outside of scan!");
			
			// Here one could also devide the intensity_sum by max((IntensityType)region_.size(),1.0))
			if (traits_->getPeakIntensity(current_index) < intensity_sum  *  min_intensity_contribution )
			{
 				//cout << "Skipping point because of low intensity contribution. " << endl;
 				//cout << current_peak.getIntensity() << " " << (intensity_sum * min_intensity_contribution) << endl;
				continue;			 
			}
			
			// remember last extracted peak
			last_pos_extracted_[RT] = traits_->getPeakRt(current_index);
			last_pos_extracted_[MZ] = traits_->getPeakMz(current_index);

			// Now we explore the neighbourhood of the current peak. Peaks in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.
			
			// Add position to the current average of positions weighted by intensity
			running_avg_.add(last_pos_extracted_,traits_->getPeakIntensity(current_index));

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
				intensity_sum += traits_->getPeakIntensity(current_index);
	 			//cout << "Added point to region. Intensity sum is now: " << intensity_sum << endl;
	 			//cout << "Intensity of the added point is : " << current_peak.getIntensity() << endl;
			}
    } // end of while ( !boundary_.empty() )

    cout << "Feature region size: " << region_.size() << endl;
    
    return region_;

	} // end of extend


	bool SimpleExtender::isTooFarFromCentroid_(const IDX& index)
	{
	
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<traits_->getData().size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<traits_->getData()[index.first].size() , "Peak index outside of scan!");

     const FeaFiTraits::PositionType2D& curr_mean = running_avg_.getPosition();

    if ( traits_->getPeakMz(index) > curr_mean[MZ] + dist_mz_up_   ||
				 traits_->getPeakMz(index) < curr_mean[MZ] - dist_mz_down_ ||
				 traits_->getPeakRt(index) > curr_mean[RT] + dist_rt_up_   ||
				 traits_->getPeakRt(index) < curr_mean[RT] - dist_rt_down_ )
    {
    	//too far
			return true;
    }
		
		//close enough
		return false;
	}

	void SimpleExtender::moveMzUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getNextMz(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveMzDown_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getPrevMz(tmp);
				if (isTooFarFromCentroid_(tmp))	break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveRtUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;

			while (true)
			{
				traits_->getNextRt(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void SimpleExtender::moveRtDown_(const IDX& index)
	{
    try
    {
			IDX tmp = index;
			while (true)
			{
				traits_->getPrevRt(tmp);
				if (isTooFarFromCentroid_(tmp)) break;
				checkNeighbour_(tmp);
			}
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
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<traits_->getData().size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<traits_->getData()[index.first].size(), "Peak index outside of scan!");
		
    // skip this point if its intensity is too low
    if (traits_->getPeakIntensity(index) <= intensity_threshold_) 	return;

    if (traits_->getPeakFlag(index) == FeaFiTraits::UNUSED)
    {
			double pr_new = computePeakPriority_(index);
			
			if (pr_new > priority_threshold_)
			{
				map<IDX, double>::iterator piter = priorities_.find(index);
				if (piter == priorities_.end()) // not yet in boundary
				{
					priorities_[index] = pr_new;
					boundary_.push(IndexWithPriority(index,pr_new));
				}
			}
			
			// Note that Clemens used to update priorities in his FF algo version...
			// I don't think that this is necessary. 
			
		}
	}


} // end of class SimpleExtender
