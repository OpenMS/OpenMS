// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//				   OpenMS Mass Spectrometry Framework
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakSeeder.h>

using namespace std;

namespace OpenMS
{
	PeakSeeder::PeakSeeder()
			: BaseSeeder(), 
				is_initialized_(false)
	{
		setName(getProductName());
		
		// minimum sn ratio for "interesting" peaks
		defaults_.setValue("min_snratio",1.1);
		// tolerance in m/z for a peak in the previous scan
		defaults_.setValue("tolerance_mz",1.1f);
			// minimum number of scan per isotopic cluster
		defaults_.setValue("min_number_scans",5);
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_peaks",20);
		// maximum distance to neighbouring peaks in the same scan
		defaults_.setValue("max_peak_distance",1.2);
		
		defaultsToParam_();
	}
	
	PeakSeeder::~PeakSeeder()
	{
	}

  PeakSeeder::PeakSeeder(const PeakSeeder& rhs)
    : BaseSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  PeakSeeder& PeakSeeder::operator= (const PeakSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
	
	FeaFiModule::IndexSet PeakSeeder::nextSeed() throw (NoSuccessor)
	{
		if (!is_initialized_)
		{
			sweep_();
			is_initialized_ = true;
		}
		
		if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(0,0));
		}
		
		cout << "Retrieving next region... " << endl; 
			
		return (curr_region_++)->second.peaks_;
	}
	
	void PeakSeeder::sweep_()
	{
		// stores the monoisotopic peaks of isotopic clusters
		vector<CoordinateType> iso_last_scan;			// in the previous scan
		vector<CoordinateType> iso_curr_scan;			// in the last scan

		CoordinateType mz_in_hash = 0; // used as reference to the current isotopic peak
		CoordinateType current_rt = 0; // retention time of current scan
		
		// max. distance to nighbouring peaks
		CoordinateType max_dist_mz = param_.getValue("max_peak_distance");
			
		for (UnsignedInt currscan_index = 0; currscan_index < traits_->getData().size(); ++currscan_index)
		{
			SpectrumType current_scan = traits_->getData()[currscan_index];
			current_rt = current_scan.getRetentionTime();

			cout << "---------------------------------------------------------------------------" << endl;
			cout << "Processing scan " << currscan_index << " of " << traits_->getData().size() << endl;
			cout << "Retention time: " << current_rt << endl;

			// copy cluster coordinates
			iso_last_scan = iso_curr_scan;
			iso_curr_scan.clear();
			
			#ifdef DEBUG_FEATUREFINDER
			String fname = String("scan_") + current_rt;
			ofstream peakfile( fname.c_str() );
			for (SpectrumType::const_iterator it = current_scan.begin(); it != current_scan.end(); ++it)			
			{
				peakfile << it->getPos() << " " << it->getIntensity() << endl;
			}
			peakfile.close();
			#endif
						
			// search for maximal positions and extract potential peaks
			vector<int> local_maxima;
			filterAndComputeLocalMax_(current_scan, local_maxima
																								#ifdef DEBUG_FEATUREFINDER
																								, currscan_index
																								#endif
																								);
	
			int nr_maxima = local_maxima.size();
			cout << "# local maxima after filtering : " << nr_maxima << endl;
			
			// test for groups of local maxima resembling isotopic pattern
			for (int z = 0; z< nr_maxima; ++z)
			{
				// store the m/z of the current peak
				CoordinateType curr_mz = current_scan[ local_maxima[z] ].getPos();
				
				#ifdef DEBUG_FEATUREFINDER
				cout << "Peak found ! " << endl;
				cout << "We are at: " << current_rt << " " << curr_mz << endl;
				cout << "Storing found peak in current isotopic cluster" << endl;
				#endif
				
				// retrieve iterator for hash entry
				TableType::iterator entry_to_insert = retrieveHashIter_(curr_mz, mz_in_hash, iso_last_scan, currscan_index);				
				
				// store peak index, scan index and set flag  
				entry_to_insert->second.peaks_.insert(make_pair(currscan_index, local_maxima[z] ));	
				traits_->getPeakFlag( make_pair(currscan_index, local_maxima[z] ) ) = FeaFiTraits::SEED;
				entry_to_insert->second.scans_.push_back( currscan_index );		
			
				iso_curr_scan.push_back( mz_in_hash );			// store peak position for sweep line
								
				if (z < (nr_maxima - 1) ) 
				{
					// there are more local max, check how far away they are
					CoordinateType dist2nextpeak = (current_scan[ local_maxima[ (z+1) ]  ].getPos() - curr_mz);
				
					// collect all local maxima which are not too far away
					// and insert them into the same peak cluster (hash entry)
					while (dist2nextpeak < max_dist_mz && z < (nr_maxima - 1) )	 
					{						
						// skip peaks that have already been used
						if ( traits_->getPeakFlag( make_pair(currscan_index, local_maxima[z] ) ) != FeaFiTraits::UNUSED) 
						{
							++z;
							continue;
						}
						
						// store next local max in hash and set flag accordingly						
						curr_mz = current_scan[ local_maxima[z] ].getPos();
						entry_to_insert->second.peaks_.insert(make_pair(currscan_index, local_maxima[z] ));	
						traits_->getPeakFlag( make_pair(currscan_index, local_maxima[z] ) ) = FeaFiTraits::SEED;
											
						// add data points between local maxima
						if (z > 0)
						{
							int begin = local_maxima[z-1];
							int end   = local_maxima[z];
							
							for (int k = begin; k < end; ++k)
							{
								entry_to_insert->second.peaks_.insert( make_pair(currscan_index, k ) ); 
								traits_->getPeakFlag( make_pair(currscan_index, k ) ) = FeaFiTraits::SEED;
							}
						}
						// get next distance
						dist2nextpeak = (current_scan[ local_maxima[ (z+1) ]  ].getPos() - curr_mz);			
						++z;				
						
					} // end of if (distance2nextpek < max_dist)
				
			} // end for  (# maximum ))
	
			
		}
		
		cout << "Numer of potential seeds in this scan: " << iso_curr_scan.size() << endl;
		
		} // end scan loop
	
		typedef map<CoordinateType,IsotopeCluster>::iterator HashIterator;
		vector<HashIterator> toDelete;
	
		cout << iso_map_.size() << " isotopic clusters were found." << endl;
	
		UnsignedInt min_number_scans = param_.getValue("min_number_scans");
		UnsignedInt min_number_peaks = param_.getValue("min_number_peaks");
	
		// Remove cluster containing too few scans or peaks
		for (HashIterator iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
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
	
	
	PeakSeeder::TableType::iterator PeakSeeder::retrieveHashIter_(const CoordinateType& curr_mz, 
																																															CoordinateType& mz_in_hash, 
																																															const std::vector<CoordinateType>& iso_last_scan,
																																								 							const UnsignedInt& currscan_index  )
	{
		// hash entry to write in
		TableType::iterator entry_to_insert;		
		
		// mass tolerance for peak cluster in previous scans
		CoordinateType tolerance_mz = param_.getValue("tolerance_mz");
		
		if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
		{
				// there were some isotopic cluster in the last scan...
				vector<double>::const_iterator it = searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
				CoordinateType delta_mz = fabs(*it - curr_mz);
	
				if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
				{
					mz_in_hash = curr_mz; // update current hash key
	
					#ifdef DEBUG_FEATUREFINDER
					cout << "Last peak cluster too far, creating new cluster" << endl;
					#endif
	
					// create new isotopic cluster												
					IsotopeCluster isoclust;
					entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );   
				}
				else
				{						
					#ifdef DEBUG_FEATUREFINDER
					cout << "Found neighbouring peak with distance (m/z) " << delta_mz << endl;
					#endif
			
					mz_in_hash = *it;	// retrieve hash key
					bool entry_found = false;				
			
					pair<TableType::iterator, TableType::iterator> range = iso_map_.equal_range(mz_in_hash);
					if (range.first != range.second)		// several peak cluster at this m/z found
					{
						// we want to find the previous scan
						UnsignedInt scan_wanted = (currscan_index - 1);
						for (TableType::iterator iter = range.first; iter != range.second; ++iter)
						{
							// check if last scan of this cluster is the previous scan
							if (	*( iter->second.scans_.end() - 1) == scan_wanted )
							{
								entry_to_insert	= iter;
								entry_found = true;
								continue;
							}
						}									
				}
				else	// only one cluster with this m/z
				{
					entry_to_insert	 = range.first;		
					entry_found = true;								
				}
				
				if (!entry_found)
				{
					// This should not happen
					IsotopeCluster isoclust;
					entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );		  
				}
					
				#ifdef DEBUG_FEATUREFINDER
				cout << "Cluster with " << entry_to_insert->second.peaks_.size() << " peaks retrieved." << endl;
				#endif
			}
		}
		else // last scan did not contain any isotopic cluster
		{
			#ifdef DEBUG_FEATUREFINDER
			cout << "Last scan was empty => creating new cluster." << endl;
			cout << "Creating new cluster at m/z: " << curr_mz << endl;
			#endif	
			mz_in_hash = curr_mz; // update current hash key
					
			IsotopeCluster isoclust;
			entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );		  
	
		} // end if (iso_last_scan.size() > 0)
	
		return entry_to_insert;		
	}
	
	void PeakSeeder::filterAndComputeLocalMax_(const SpectrumType & vec, 
														 						 															std::vector<int>& localmax
																																			#ifdef DEBUG_FEATUREFINDER
																																			, const UnsignedInt& currscan_index
																																			#endif
																																			)
	{
	
		// minimum s/n ratio for an interesting peak
		IntensityType min_sn = param_.getValue("min_snratio");			
		
		DSignalToNoiseEstimatorMeanIterative<1, SpectrumType::const_iterator > estimator;
		estimator.init(vec.begin(), vec.end());
		
		int i = 0;
		for(SpectrumType::const_iterator it = vec.begin(); it != vec.end(); ++it)
		{
			// Check for maximum at position i 
			if(   (i > 0) && 
					(vec[i-1].getIntensity() - vec[i].getIntensity()   < 0) &&
					(vec[i].getIntensity() - vec[i+1].getIntensity() > 0) )
			{
				#ifdef DEBUG_FEATUREFINDER
				IDX id;
				id.first      =  currscan_index;
				id.second = i;
				String locmaxfile_name = String("sn_localmax_") + traits_->getPeakRt(id);
				ofstream lmfile( locmaxfile_name.c_str(), ios_base::app);
				lmfile << traits_->getPeakMz(id)  << "  " <<  traits_->getPeakIntensity(id) << endl;
				lmfile.close();
				#endif
				
				if (estimator.getSignalToNoise(it) > min_sn) localmax.push_back(i);
				
				}
				++i;
			}
		
	} // end of getMaxPositions_(...)

} // end of namespace OpenMS
