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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SweepExtender.h>

namespace OpenMS
{

SweepExtender::SweepExtender() 
	: BaseExtender(), is_initialized_(false)
{
    name_ = SweepExtender::getName();
	
	defaults_.setValue("charge1_ub",1.3f);
	defaults_.setValue("charge1_lb",0.7f);
	defaults_.setValue("charge2_ub",0.8f);
	defaults_.setValue("charge2_lb",0.3f);
	defaults_.setValue("tolerance_mz",0.7f);
    param_ = defaults_;
}

SweepExtender::~SweepExtender()
{}

const IndexSet& SweepExtender::extend(const UnsignedInt /*seed_index*/)
{
	if (!is_initialized_) sweep_();
	
	region_.clear();
	
	is_initialized_ = true;
	
	if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 ) throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);
	
	// check if this cluster consists of one scan only
	if ((*curr_region_).second.scans_.size() == 1)
	{
		//std::cout << "Only one scan " << std::endl;
		++curr_region_;
		return region_;	// this isotopic pattern occurs in one scan only and is therefore unreliable => omit it	
	}	
	std::vector<UnsignedInt> next_region = (*curr_region_).second.peaks_;
		
	for (std::vector<UnsignedInt>::const_iterator cit = next_region.begin();
	       cit != next_region.end();
		   ++cit)
	{
		region_.add(*cit);
	}
	
	++curr_region_;
	
	region_.sort();
	
	std::cout << "Retrieved region with " << region_.size() << " peaks. " << std::endl;
	
    return region_;

} // end of extend(Unsigned int seed_index)

void SweepExtender::sweep_()
{
    // stores the monoisotopic peaks of isotopic clusters
    std::vector<double> iso_last_scan;
    std::vector<double> iso_curr_scan;

    unsigned int nr_peaks    = traits_->getNumberOfPeaks();
    CoordinateType last_rt  = traits_->getPeakRt(0);
	
	// retrieve values for accepted peaks distances
	charge1_ub_    = param_.getValue("charge1_ub");
	charge1_lb_     = param_.getValue("charge1_lb");
	
	charge2_ub_    = param_.getValue("charge2_ub");
	charge2_lb_     = param_.getValue("charge2_lb");
	
	CoordinateType tolerance_mz = param_.getValue("tolerance_mz");
	UnsignedInt current_charge     = 0;
	
    // sweep through scans
    for (unsigned int curr_peak = 0; curr_peak < nr_peaks; ++curr_peak)
    {
        CoordinateType current_rt =	traits_->getPeakRt(curr_peak);
        // check if new scan has begun
        if (current_rt != last_rt)
        {
            // copy cluster information of least scan
            iso_last_scan = iso_curr_scan;
            iso_curr_scan.clear();
			last_rt = current_rt;
        }
        // store the m/z of the current peak
        CoordinateType curr_mz         = traits_->getPeakMz(curr_peak);
		if ( (curr_peak+1) >= nr_peaks ) break;
        CoordinateType dist2nextpeak = ( traits_->getPeakMz(curr_peak+1) - curr_mz);
		
		// test for different charge states
		current_charge = testDistance2NextPeak_(dist2nextpeak);
			
		//std::cout << "Charge set to " << current_charge << std::endl;
		
        // check for charge 1 isotope
        if (current_charge > 0)
        {
            IsotopeCluster iso_clust; // stores scan and peaks

            if (iso_last_scan.size() > 0)
            {
                // there were some isotopic clustures in the last scan...
                std::vector<double>::iterator it = searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
                double delta_mz = (*it - curr_mz);

                if ( fabs(delta_mz) > tolerance_mz)
                {
                    // create new isotopic cluster
                    iso_clust.left_mz_ = curr_mz;
                    iso_clust.charge_  = current_charge;
                    iso_clust.scans_.push_back( traits_->getPeakRt(curr_peak) );
                }
                else
                {
                    std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
                    curr_mz = *it;
                    iso_clust = iso_map_[curr_mz];
                    iso_clust.scans_.push_back( traits_->getPeakRt(curr_peak) );
                    //std::cout << "Cluster with " << iso_clust.peaks_.size() << " peaks retrieved." << std::endl;
                }

            }
            else
            {
                //std::cout << "Last scan empty. New cluster." << std::endl;
			    iso_clust.left_mz_ = curr_mz;
                iso_clust.charge_  = current_charge;
                iso_clust.scans_.push_back( traits_->getPeakRt(curr_peak) );
            }

            //std::cout << "store found peak in current isotopic cluster" << std::endl;
            iso_clust.peaks_.push_back(curr_peak);
			iso_curr_scan.push_back(traits_->getPeakMz(curr_peak));
            ++curr_peak;
		    if (curr_peak == nr_peaks ) break;
			
            // std::cout << "and store next one as well" << std::endl;
            iso_clust.peaks_.push_back(curr_peak);
           // iso_curr_scan.push_back(traits_->getPeakMz(curr_peak));

		   //std::cout << "computing next distance.." << std::endl;
		   if ( (curr_peak+1) >= nr_peaks ) break;
		   
            dist2nextpeak = ( traits_->getPeakMz(curr_peak+1) -  traits_->getPeakMz(curr_peak));
			
			if (testDistance2NextPeak_(dist2nextpeak) != current_charge) continue;	// charge state should remain the same 
						
            while (current_charge > 0)
            {
                ++curr_peak;
                if (curr_peak == nr_peaks) break;
           		iso_clust.peaks_.push_back(curr_peak);		// save peak in cluster

               dist2nextpeak = ( traits_->getPeakMz(curr_peak+1) -  traits_->getPeakMz(curr_peak)); // get distance to next peak
			   current_charge = testDistance2NextPeak_(dist2nextpeak);
			   
                // store current cluster
                iso_map_[curr_mz] = iso_clust;
            
			} // end while(...)
			
			iso_map_[curr_mz] =   iso_clust; 	// update cluster for this mass
			
        } // end of if (charge > 0)
    
		current_charge = 0; // reset charge
	} // end for (...)

	curr_region_ = iso_map_.begin();
	
	std::cout << iso_map_.size() << " clusters were found ! " << std::endl;
	
} // end of void sweep_()


UnsignedInt SweepExtender::testDistance2NextPeak_(CoordinateType dist2nextpeak)
{

	if (dist2nextpeak < charge1_ub_ && dist2nextpeak > charge1_lb_) 
	{
		return 1;		
	} 
	else if (dist2nextpeak < charge2_ub_ && dist2nextpeak > charge2_lb_) 
	{
		return 2;
	}
	else
	{
	 	return 0;
	}
}


} // end of class SweepExtender
