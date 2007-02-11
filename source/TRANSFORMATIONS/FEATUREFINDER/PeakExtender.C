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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakExtender.h>

using namespace std;

namespace OpenMS
{

	PeakExtender::PeakExtender() 
	: BaseExtender(),
		seed_pos_()
	{
    setName(getProductName());
        
    defaults_.setValue("dist_mz_up",6.0f);
    defaults_.setValue("dist_mz_down",2.0f);
    
		defaults_.setValue("dist_rt_up",5.0f);
    defaults_.setValue("dist_rt_down",5.0f);
    
    defaults_.setValue("intensity_factor",0.03f);
		
		defaults_.setValue("min_intensity_contribution",0.01f);

    defaultsToParam_();
	}

	PeakExtender::~PeakExtender()
	{
	}

  PeakExtender::PeakExtender(const PeakExtender& rhs)
    : BaseExtender(rhs)
  {
    updateMembers_();
  }
  
  PeakExtender& PeakExtender::operator= (const PeakExtender& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseExtender::operator=(rhs);
    
    updateMembers_();
    
    return *this;
  }

  void PeakExtender::updateMembers_()
  {
		// initialize members
		dist_mz_up_     = param_.getValue("dist_mz_up");
		dist_mz_down_ = param_.getValue("dist_mz_down");
		dist_rt_up_       = param_.getValue("dist_rt_up");
		dist_rt_down_   = param_.getValue("dist_rt_down");
  }

	const FeaFiModule::IndexSet& PeakExtender::extend(const IndexSet& seed_region)
	{
    // empty region and boundary datastructures
    region_.clear();
		boundary_ = priority_queue< IDX >();
				
		// find maximum of region (seed) and add data points in seeding region to feature
		CoordinateType max_intensity = 0.0;
		IDX seed;
    for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
    {	
      if (traits_->getPeakIntensity(*citer) > max_intensity)
      {
        seed = *citer;
        max_intensity = traits_->getPeakIntensity(seed);						
			}
			region_.insert(*citer);
    }
  
		// remember last extracted point (in this case the seed !)
		seed_pos_[RT] = traits_->getPeakRt(seed);
		seed_pos_[MZ] = traits_->getPeakMz(seed);

		cout << "Extending from " << traits_->getPeakRt(seed) << "/" << traits_->getPeakMz(seed); 
		cout << " (" << seed.first << "/" << seed.second << ")" << endl;
		
		// compute intensity threshold and sum
		IntensityType intensity_sum = 0.0;
		IntensityType min_intensity_contribution = param_.getValue("min_intensity_contribution");
				
    while (!boundary_.empty())
    {
			// remove peak with highest priority
			const IDX  current_index = boundary_.top();
			boundary_.pop();
			
    	//Corrupt index
    	OPENMS_PRECONDITION(current_index.first<traits_->getData().size(), "Scan index outside of map!");
      OPENMS_PRECONDITION(current_index.second<traits_->getData()[current_index.first].size(), "Peak index outside of scan!");
		
			if (traits_->getPeakIntensity(current_index) < (intensity_sum/std::sqrt(std::max((IntensityType)region_.size(),1.0)) )  *  min_intensity_contribution )
			{
 				cout << "Skipping point because of low intensity contribution. " << endl;
 				cout << traits_->getPeakIntensity(current_index) << " " << (intensity_sum/std::sqrt(std::max((IntensityType) region_.size(),1.0) )  *  min_intensity_contribution ) << endl;
				continue;			 
			}
						
			// explore neighbourhood !
			moveMzUp_(current_index);
			moveMzDown_(current_index);
			moveRtUp_(current_index);
			moveRtDown_(current_index);

			// check flags (if data point is already inside a feature region or used as seed, we discard it)
			if (traits_->getPeakFlag(current_index) == FeaFiTraits::SEED || traits_->getPeakFlag(current_index) == FeaFiTraits::UNUSED)
			{
				traits_->getPeakFlag(current_index) = FeaFiTraits::INSIDE_FEATURE;
				region_.insert(current_index);
				intensity_sum += traits_->getPeakIntensity(current_index);
	 			cout << "Added point to region. Intensity sum is now: " << intensity_sum << endl;
	 			cout << "Intensity of included point is : " << traits_->getPeakIntensity(current_index) << endl;
			}
			
			cout << "Size of boundary: " << boundary_.size() << endl;
			
    } // end of while ( !boundary_.empty() )

    cout << "Feature region size: " << region_.size() << endl;
    
    return region_;

	} // end of extend


	bool PeakExtender::isTooFarFromSeed_(const IDX& index)
	{	
  	//Corrupt index
  	OPENMS_PRECONDITION(index.first<traits_->getData().size(), "Scan index outside of map!");
    OPENMS_PRECONDITION(index.second<traits_->getData()[index.first].size() , "Peak index outside of scan!");

    if ( traits_->getPeakMz(index) > seed_pos_[MZ] + dist_mz_up_   ||
				 traits_->getPeakMz(index) < seed_pos_[MZ] - dist_mz_down_ ||
				 traits_->getPeakRt(index) >  seed_pos_[RT] + dist_rt_up_   ||
				 traits_->getPeakRt(index) <  seed_pos_[RT] - dist_rt_down_ )
    {
    	//too far
			return true;
    }
		
		//close enough
		return false;
	}

	void PeakExtender::moveMzUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getNextMz(tmp);
				if (isTooFarFromSeed_(tmp)) break;
				
				if (traits_->getPeakFlag(index) != FeaFiTraits::INSIDE_FEATURE)
				{
					boundary_.push(tmp);
				}
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void PeakExtender::moveMzDown_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;
			while (true)
			{
				traits_->getPrevMz(tmp);
				if (isTooFarFromSeed_(tmp))	break;
				
				if (traits_->getPeakFlag(index) != FeaFiTraits::INSIDE_FEATURE)
				{
					boundary_.push(tmp);
				}
				
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void PeakExtender::moveRtUp_(const IDX& index)
	{
    try
    {
    	IDX tmp = index;

			while (true)
			{
				traits_->getNextRt(tmp);
				if (isTooFarFromSeed_(tmp)) break;
				
				if (traits_->getPeakFlag(index) != FeaFiTraits::INSIDE_FEATURE)
				{
					boundary_.push(tmp);
				}
				
			}
    }
    catch(NoSuccessor)
    {
    }
	}

	void PeakExtender::moveRtDown_(const IDX& index)
	{
    try
    {
			IDX tmp = index;
			while (true)
			{
				traits_->getPrevRt(tmp);
				if (isTooFarFromSeed_(tmp)) break;
				
				if (traits_->getPeakFlag(index) != FeaFiTraits::INSIDE_FEATURE)
				{
					boundary_.push(tmp);
				}
			
			}
    }
    catch(NoSuccessor)
    {
    }
	}
	

} // end of class PeakExtender
