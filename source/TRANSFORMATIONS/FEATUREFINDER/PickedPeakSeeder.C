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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PickedPeakSeeder.h>

using namespace std;

namespace OpenMS
{

	PickedPeakSeeder::PickedPeakSeeder()
		: BaseSeeder(), 
			is_initialized_(false)
	{
    setName(getProductName());

    // lower and upper bounds for distances between isotopic peaks (defaults)
    // charge 1
    defaults_.setValue("charge1_ub",1.1f);
    defaults_.setValue("charge1_lb",0.9f);
    // charge 2
    defaults_.setValue("charge2_ub",0.50f);
    defaults_.setValue("charge2_lb",0.4f);
    // charge 3
    defaults_.setValue("charge3_ub",0.4f);
    defaults_.setValue("charge3_lb",0.2f);
    // charge 4
    defaults_.setValue("charge4_ub",0.24f);
    defaults_.setValue("charge4_lb",0.26f);
    // charge 5
    defaults_.setValue("charge5_ub",0.21f);
    defaults_.setValue("charge5_lb",0.19f);

    // tolerance in m/z for an monoisotopic peak in the previous scan
    defaults_.setValue("tolerance_mz",1.0f);
		
		 // minimum number of scan per isotopic cluster
    defaults_.setValue("min_number_scans",3);
    // minimum number of peaks per cluster
    defaults_.setValue("min_number_peaks",6);

    defaultsToParam_();
	}
	
	PickedPeakSeeder::~PickedPeakSeeder()
	{
	}

  PickedPeakSeeder::PickedPeakSeeder(const PickedPeakSeeder& rhs)
    : BaseSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  PickedPeakSeeder& PickedPeakSeeder::operator= (const PickedPeakSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
	
	FeaFiModule::IndexSet PickedPeakSeeder::nextSeed() throw (NoSuccessor)
	{
    if (!is_initialized_)
    {
    	sweep_();
    	is_initialized_ = true;
    }
    
		if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 ) 
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(0u,0u));
		}
		
    IndexSet next_region = (*curr_region_).second.peaks_;

    ++curr_region_;

    return next_region;	
	}
	
	void PickedPeakSeeder::updateMembers_()
	{
		// retrieve values for accepted peaks distances
		charge1_ub_	= param_.getValue("charge1_ub");
		charge1_lb_	 = param_.getValue("charge1_lb");
	
		charge2_ub_	= param_.getValue("charge2_ub");
		charge2_lb_	 = param_.getValue("charge2_lb");
	
		charge3_ub_	= param_.getValue("charge3_ub");
		charge3_lb_	 = param_.getValue("charge3_lb");
	
		charge4_ub_	= param_.getValue("charge4_ub");
		charge4_lb_	 = param_.getValue("charge4_lb");
	
		charge5_ub_	= param_.getValue("charge5_ub");
		charge5_lb_	 = param_.getValue("charge5_lb");		
	}
	
	void PickedPeakSeeder::sweep_()
	{
		// stores the monoisotopic peaks of isotopic clusters
		vector<double> iso_last_scan;
		vector<double> iso_curr_scan;
	
		CoordinateType tolerance_mz = param_.getValue("tolerance_mz");
	
		UnsignedInt current_charge	 = 0;			// charge state of the current isotopic cluster
		CoordinateType mz_in_hash   = 0;			// used as reference to the current isotopic peak
	
		// sweep through scans
		for (UnsignedInt i=0; i < traits_->getData().size(); ++i)
		{
			const FeaFiTraits::MapType::SpectrumType& scan = traits_->getData()[i];
			
#ifdef DEBUG_FEATUREFINDER
			cout << "Next scan with rt: " << scan.getRetentionTime() << endl;
			cout << "---------------------------------------------------------------------------" << endl;
#endif
			// copy cluster information of last scan
			iso_last_scan = iso_curr_scan;
			iso_curr_scan.clear();
			
			for (UnsignedInt j=0; j < traits_->getData()[i].size()-1; ++j)
			{
				// test for different charge states
				current_charge = distanceToCharge_(scan[j+1].getPos() - scan[j].getPos());
				
				//remove false positives by looking at the intensity ratio of the peaks
				if ( fabs( scan[j].getIntensity()/scan[j+1].getIntensity() - 1.0) < 0.1)
				{
					current_charge = 0;
				}
				
				// charger = 0 => no isotope
				if (current_charge > 0) 
				{
#ifdef DEBUG_FEATUREFINDER
					cout << "Isotopic pattern found ! " << endl;
					cout << "We are at: " << scan.getRetentionTime() << " " << scan[j].getPos() << endl;
#endif
					// hash entry to write in
					TableType::iterator entry_to_insert;			
		
					if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
					{
						// there were some isotopic clustures in the last scan...
						std::vector<double>::iterator it = searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),scan[j].getPos());
						double delta_mz = fabs(*it - scan[j].getPos());
						
						// check if first peak of last cluster is close enough -> create new isotopic cluster
						if ( delta_mz > tolerance_mz)
						{
	#ifdef DEBUG_FEATUREFINDER
							cout << "Last peak cluster too far, creating new cluster" << endl;
	#endif
	
							mz_in_hash = scan[j].getPos(); // update current hash key
			
							IsotopeCluster isoclust;
							isoclust.charge_ = current_charge;
							isoclust.scans_.push_back(i);
	
							entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );		   
						}
						else
						{
	#ifdef DEBUG_FEATUREFINDER
							cout << "Found neighbouring peak with distance (m/z) " << delta_mz << endl;
	#endif
							mz_in_hash = *it;	// retrieve hash key
	   										
							pair<TableType::iterator, TableType::iterator> range = iso_map_.equal_range(mz_in_hash);
																	
							if (range.first != range.second)		// several peak cluster at this m/z found
							{
								// we want to find the previous scan
								UnsignedInt scan_wanted = (i - 1);
																									
								for (TableType::iterator iter = range.first; iter != range.second; ++iter)
								{
									// check if last scan of this cluster is the previous scan
									if (	*( iter->second.scans_.end() - 1) == scan_wanted )
									{
										entry_to_insert	= iter;
										continue;
									}
								}	
							}
							else	// only one cluster with this m/z
							{
								entry_to_insert	 = range.first;										
							}
							// save current rt and m/z
							entry_to_insert->second.scans_.push_back( i );
		
	#ifdef DEBUG_FEATUREFINDER
							cout << "Cluster with " << entry_to_insert->second.peaks_.size() << " peaks retrieved." << endl;
	#endif
						}
					}
					else // last scan did not contain any isotopic cluster
					{
	#ifdef DEBUG_FEATUREFINDER
						cout << "Last scan was empty => creating new cluster." << endl;
						cout << "Creating new cluster at m/z: " << scan[j].getPos() << endl;
	#endif
						mz_in_hash = scan[j].getPos(); // update current hash key
										
						IsotopeCluster isoclust;
						isoclust.charge_ = current_charge;
						isoclust.scans_.push_back(i);
	
						entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );		  
					}
		
	#ifdef DEBUG_FEATUREFINDER
					cout << "Storing found peak in current isotopic cluster" << endl;
	#endif
		
					// add next two peaks to current cluster
					iso_curr_scan.push_back(mz_in_hash);
					entry_to_insert->second.peaks_.insert(std::make_pair(i,j));
					++j;
					entry_to_insert->second.peaks_.insert(std::make_pair(i,j));
								
					// if distance to next peak does not correspond to a isotopic spacing, just add one more peak 
					// and continue
					if (distanceToCharge_(scan[j+1].getPos() - scan[j].getPos()) != current_charge)
					{
						entry_to_insert->second.peaks_.insert(std::make_pair(i,j+1));
						continue;
					}
		
					// loop until end of isotopic pattern in this scan
					while (j+1!=scan.size() && distanceToCharge_(scan[j+1].getPos() - scan[j].getPos()) == current_charge )
					{
						++j;
						entry_to_insert->second.peaks_.insert(std::make_pair(i,j));				// save peak in cluster
					}
					
				} // end of if (charge > 0)
						
			} // end keep loop
		
		} //end scan loop
		
		vector< TableType::iterator > toDelete;
	
		cout << iso_map_.size() << " isotopic clusters were found." << endl;
			
		UnsignedInt min_number_scans = param_.getValue("min_number_scans");
		UnsignedInt min_number_peaks = param_.getValue("min_number_peaks");
	
		// remove cluster having less than 6 peaks or less than 3 scans
		for (TableType::iterator iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			if (iter->second.scans_.size() < min_number_scans ||  iter->second.peaks_.size() < min_number_peaks)
			{
				toDelete.push_back(iter);
			}
		}
	
		for (unsigned int i=0; i<toDelete.size();++i)
		{
			iso_map_.erase(toDelete[i]);
		}
	
		curr_region_ = iso_map_.begin();
		cout << iso_map_.size() << " clusters remained after filtering." << endl;
	
	} // end of void sweep_()
	
	
	UnsignedInt PickedPeakSeeder::distanceToCharge_(const CoordinateType& dist)
	{
    if (dist < charge1_ub_ && dist > charge1_lb_)
    {
    	return 1;
    }
    else if (dist < charge2_ub_ && dist > charge2_lb_)
    {
    	return 2;
    }
    else if (dist < charge3_ub_ && dist > charge3_lb_)
    {
    	return 3;
    }
    else if (dist < charge4_ub_ && dist > charge4_lb_)
    {
    	return 4;
    }
    else if (dist < charge5_ub_ && dist > charge5_lb_)
    {
    	return 5;
    }
    else
    {
    	return 0;
    }
	}

} // end of namespace OpenMS
