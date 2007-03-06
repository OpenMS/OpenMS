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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MarrWaveletSeeder.h>

using namespace std;

namespace OpenMS
{
	MarrWaveletSeeder::MarrWaveletSeeder()
			: BaseSeeder(), 
				is_initialized_(false)
	{
		setName(getProductName());
	
		// lower and upper bounds for distances between isotopic peaks (defaults)
		// charge 1
		defaults_.setValue("charge1_ub",1.5f);
		defaults_.setValue("charge1_lb",0.70f);
		// charge 2
		defaults_.setValue("charge2_ub",0.69f);
		defaults_.setValue("charge2_lb",0.41f);
		// charge 3
		defaults_.setValue("charge3_ub",0.4f);
		defaults_.setValue("charge3_lb",0.1f);
	
		// tolerance in m/z for an monoisotopic peak in the previous scan
		defaults_.setValue("tolerance_mz",1.0f);
	
		// params for the cwt
		defaults_.setValue("cwt_scale",0.1f);
		defaults_.setValue("noise_level_signal",1000);
		defaults_.setValue("noise_level_cwt",6000);
	
		// number of scans used for alignment
		defaults_.setValue("scans_to_sumup",5);
		// mass tolerance during scan alignment
		defaults_.setValue("mass_tolerance_alignment", 0.1);
		
		// minimum number of scan per isotopic cluster
		defaults_.setValue("min_number_scans",5);
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_peaks",20);
		
	
		defaultsToParam_();
	}
	
	MarrWaveletSeeder::~MarrWaveletSeeder()
	{
	}

  MarrWaveletSeeder::MarrWaveletSeeder(const MarrWaveletSeeder& rhs)
    : BaseSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  MarrWaveletSeeder& MarrWaveletSeeder::operator= (const MarrWaveletSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
	
	FeaFiModule::IndexSet MarrWaveletSeeder::nextSeed() throw (NoSuccessor)
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
		
		std::cout << "Retrieving next region with charge " << (*curr_region_).second.charge_; 
		std::cout << " and size " << (*curr_region_).second.peaks_.size() << std::endl;
			
		return (curr_region_++)->second.peaks_;
	}
	
	void MarrWaveletSeeder::updateMembers_()
	{
		// retrieve values for accepted peaks distances
		charge1_ub_	= param_.getValue("charge1_ub");
		charge1_lb_	 = param_.getValue("charge1_lb");
	
		charge2_ub_	= param_.getValue("charge2_ub");
		charge2_lb_	 = param_.getValue("charge2_lb");
	
		charge3_ub_	= param_.getValue("charge3_ub");
		charge3_lb_	 = param_.getValue("charge3_lb");
		
		mass_tolerance_alignment_ = param_.getValue("mass_tolerance_alignment");
	}
	
	void MarrWaveletSeeder::sweep_()
	{
		// stores the monoisotopic peaks of isotopic clusters
		vector<CoordinateType> iso_last_scan;			// in the previous scan
		vector<CoordinateType> iso_curr_scan;			// in the last scan

		CoordinateType tolerance_mz = param_.getValue("tolerance_mz");
		double cwt_scale                   = param_.getValue("cwt_scale");

		noise_level_signal_ 							 = param_.getValue("noise_level_signal");
		noise_level_cwt_                   = param_.getValue("noise_level_cwt");
	
		UInt current_charge = 0;															 // charge state of the current isotopic cluster
		CoordinateType mz_in_hash = 0; 								// used as reference to the current isotopic peak
		CoordinateType current_rt = 0; 										// retention time of current scan
	
		for (UInt currscan_index = 0; currscan_index < traits_->getData().size(); ++currscan_index)
		{
			//make a copy of the scan. This is necessary as the peak intensities have to be modified in the 
			//sumUp_ method
			SpectrumType current_scan = traits_->getData()[currscan_index];
			
			current_rt = current_scan.getRetentionTime();

			cout << "---------------------------------------------------------------------------" << endl;
			cout << "Processing scan " << currscan_index << " of " << traits_->getData().size() << endl;
			cout << "Retention time: " << current_rt << endl;

			// copy cluster information of last scan
			iso_last_scan = iso_curr_scan;
			iso_curr_scan.clear();
			
#ifdef DEBUG_FEATUREFINDER
			String fname = String("scan_") + current_rt;
			ofstream peakfile( fname.c_str() );
			for(unsigned k = 0; k<current_scan.size();++k)
			{
				peakfile << current_scan[k].getMZ() << " " << current_scan[k].getIntensity() << endl;
			}
			peakfile.close();
#endif
	
			// align and sum
			sumUp_(current_scan,currscan_index);
					
			double spacing_cwt = 0.0001; 		// spacing between sampling points of the wavelet
			double resolution_cwt = 1.0;			 // compute convolution for each data point
								
			cwt_.init(cwt_scale, spacing_cwt);
			cwt_.transform(current_scan.begin(), current_scan.end(),resolution_cwt);
	
#ifdef DEBUG_FEATUREFINDER
			fname = String("cwt_") + current_rt;
			ofstream gpfile( fname.c_str() );
			for (int i=0;i<cwt_.getSize(); ++i)
			{
				gpfile << (current_scan.begin() + i)->getMZ() << "  " << cwt_[i] << endl;
			}
			gpfile.close();
#endif
	
			if  (cwt_.getSize() == 0)
			{
				#ifdef DEBUG_FEATUREFINDER
				cout << "Empty cwt for this scan => continue." << endl;
				#endif
				continue;
			}
	
			// compute average intensity of cwt
			double cwt_avg = 0.0;
			for (int k = 0; k<cwt_.getSize();++k)
			{
				cwt_avg += cwt_[k];
			}
			cwt_avg /= cwt_.getSize();
			cout << "Average strength in cwt: " << cwt_avg<< endl;
	
			// search for maximal positions in the cwt and extract potential peaks
			vector<int> local_maxima;
			getMaxPositions_(current_scan.begin(), current_scan.end(), cwt_, local_maxima 
																	#ifdef DEBUG_FEATUREFINDER
											 						,current_rt 
																	#endif
											 						);
	
			int nr_maxima = local_maxima.size();
			cout << "# local maxima in cwt : " << nr_maxima << endl;
			
			// test for groups local maxima resembling isotopic pattern
			for (int z = 0; z< ( nr_maxima - 2); ++z)
			{
				// store the m/z of the current peak
				CoordinateType curr_mz = current_scan[ local_maxima[z] ].getMZ();
				CoordinateType dist2nextpeak = (current_scan[ local_maxima[ (z + 1) ] ].getMZ() - curr_mz);
	
				// test for different charge states
				current_charge = distanceToCharge_(dist2nextpeak);
	
				if (current_charge > 0) // charge = 0 <=> no isotope
				{
	
#ifdef DEBUG_FEATUREFINDER
					cout << "Isotopic pattern found ! " << endl;
					cout << "We are at: " << current_rt << " " << curr_mz << endl;
#endif
									
					// hash entry to write in
					TableType::iterator entry_to_insert;		
	
					if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
					{
						// there were some isotopic clustures in the last scan...
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
							isoclust.charge_ = current_charge;
							isoclust.scans_.push_back( currscan_index );

							std::cout << "Creating cluster at m/z " << mz_in_hash << std::endl;
							entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );   
						}
						else
						{
						#ifdef DEBUG_FEATUREFINDER
						cout << "Found neighbouring peak with distance (m/z) " << delta_mz << endl;
						#endif
						mz_in_hash = *it;	// retrieve hash key
   										
						pair<TableType::iterator, TableType::iterator> range = iso_map_.equal_range(mz_in_hash);
						
						std::cout << iso_map_.count(mz_in_hash) << " entries found for m/z " << mz_in_hash << std::endl;
					
						bool scan_found = false;								
					
						if (range.first != range.second)		// several peak cluster at this m/z found
						{
							// we want to find the previous scan
							UInt scan_wanted = (currscan_index - 1);
																								
							for (TableType::iterator iter = range.first; iter != range.second; ++iter)
							{

								// enumerate all scans
								// the scan number we are searching for is not necessarily the last one
								// in this cluster if there were other very close local maxima in the same scan.
								for (std::vector<UInt>::const_iterator it = iter->second.scans_.begin();
											 it != iter->second.scans_.end(); ++it)
								{
								
									if (	*it == scan_wanted )
									{
										scan_found = true;
										entry_to_insert	= iter;
										continue;
									}														
								}
							}								
						}
						else	// only one cluster with this m/z
						{
							entry_to_insert	 = range.first;										
						}
						
						if (!scan_found)
						{
							// corrupt hash map / isotope counter
							throw Exception::InvalidIterator(	__FILE__, __LINE__, "MarrWaveletSeeder::sweep_()");
						}
						
						// save current rt and m/z
						entry_to_insert->second.scans_.push_back( currscan_index );
	
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
						isoclust.charge_ = current_charge;
						isoclust.scans_.push_back( currscan_index );

						entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );		  
	
					} // end if (iso_last_scan.size() > 0)
	
#ifdef DEBUG_FEATUREFINDER
					cout << "Storing found peak in current isotopic cluster" << endl;
#endif
					// store position of the the potential isotopic pattern
					iso_curr_scan.push_back( mz_in_hash );
					// walk a bit to the left: we include the five peaks to the left of the
					// monisotopic one into the feature region
					int ind = local_maxima[z];
					for (int p=1; p<=5; ++p)
					{
						if ( (ind-p) > 0 && traits_->getPeakFlag( make_pair(currscan_index,ind-p) ) == FeaFiTraits::UNUSED)
						{
							entry_to_insert->second.peaks_.insert(make_pair(currscan_index,ind-p));
							traits_->getPeakFlag( make_pair(currscan_index,ind-p) ) = FeaFiTraits::SEED;
						}
					}
	
					++z;	// next maximum in cwt
					if ((local_maxima[z] - local_maxima[z-1]) > 1)
					{
						for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
						{
							entry_to_insert->second.peaks_.insert(make_pair(currscan_index,v));
							traits_->getPeakFlag( make_pair(currscan_index,v)) = FeaFiTraits::SEED;
						}
					}
	
					entry_to_insert->second.peaks_.insert( make_pair(currscan_index,local_maxima[z]));
					traits_->getPeakFlag( make_pair(currscan_index,local_maxima[z])) = FeaFiTraits::SEED;
	
					// check distance to next peak
					if ( (z+1) >= nr_maxima) break;
	
					dist2nextpeak = ( current_scan[ local_maxima[z+1] ].getMZ()  -   current_scan[ local_maxima[z] ].getMZ() );
	
					if (distanceToCharge_(dist2nextpeak) != current_charge)
					{
						// charge state has changed. Insert m/z of last maximum and continue.
						++z;
						if ((local_maxima[z] - local_maxima[z-1]) > 1)
						{
							for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
							{
								entry_to_insert->second.peaks_.insert( make_pair(currscan_index,v));
								traits_->getPeakFlag(  make_pair(currscan_index,v) ) = FeaFiTraits::SEED;
							}
						}
						entry_to_insert->second.peaks_.insert( make_pair(currscan_index, local_maxima[z]) );
						traits_->getPeakFlag(  make_pair(currscan_index,local_maxima[z])  ) = FeaFiTraits::SEED;
						continue;
					}
	
					if ( (z+1) >= nr_maxima) break;
	
					CoordinateType monoiso_mass = current_scan[ local_maxima[z+1] ].getMZ();
					CoordinateType mass_diff		 = 0.0;
	
					// loop until end of isotopic pattern in this scan
					while (distanceToCharge_(dist2nextpeak) == current_charge &&
							z < (nr_maxima-2) &&
							mass_diff < 6 &&
							traits_->getPeakFlag(make_pair(currscan_index, local_maxima[z])) == FeaFiTraits::UNUSED)
					{
						entry_to_insert->second.peaks_.insert( make_pair(currscan_index, local_maxima[z]) );				// save peak in cluster
						traits_->getPeakFlag(  make_pair(currscan_index, local_maxima[z])  ) = FeaFiTraits::SEED;
						++z;
	
						// between two local maxima in the cwt there might be several raw data points
						if ((local_maxima[z] - local_maxima[z-1]) > 1)
						{
							for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
							{
								entry_to_insert->second.peaks_.insert( make_pair(currscan_index,v) );
								traits_->getPeakFlag(  make_pair(currscan_index,v) ) = FeaFiTraits::SEED;
							}
						}
	
						dist2nextpeak = (current_scan[ local_maxima[z+1] ].getMZ() -  current_scan[ local_maxima[z] ].getMZ() ); // get distance to next peak
						mass_diff	   = ( current_scan[ local_maxima[z+1] ].getMZ() - monoiso_mass);
					} // end while(...)
	
				} // end of if (charge > 0)
				
				current_charge = 0; // reset charge
	
			} // end for each (local maximum in cwt)
	
			cout << "Numer of potential seeds in this scan: " << iso_curr_scan.size() << endl;
	
		} // end scan loop
	
		typedef map<CoordinateType,IsotopeCluster>::iterator HashIterator;
		vector<HashIterator> toDelete;
	
		cout << iso_map_.size() << " isotopic clusters were found." << endl;
	
		UInt min_number_scans = param_.getValue("min_number_scans");
		UInt min_number_peaks = param_.getValue("min_number_peaks");
	
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
		
		std::cout << iso_map_.size() << " clusters remained after filtering." << std::endl;
		
		#ifdef DEBUG_FEATUREFINDER 
		for (HashIterator iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			std::cout << "m/z " << iter->first << " charge: " << iter->second.charge_ << std::endl;
			
			for (std::vector<UInt>::const_iterator citer = 	iter->second.scans_.begin(); 
						citer != iter->second.scans_.end();
						++citer)
			{
				std::cout << "# scan : " << *citer << " (";
				IDX tmp;
				tmp.first = *citer;
				std::cout << traits_->getPeakRt(tmp) << ")" << std::endl;				
			}
		
		}
		#endif
	
	} // end of void sweep_()
	
	void MarrWaveletSeeder::getMaxPositions_( const SpectrumType::const_iterator&  first, 
																						const SpectrumType::const_iterator&  last, 
																						const ContinuousWaveletTransform& wt, 
																						vector<int>& localmax
																						#ifdef DEBUG_FEATUREFINDER 
																							, CoordinateType current_rt
																						#endif
																						)
	{
		if (wt.getSize() == 0) return;
	
		int zeros_left_index  = wt.getLeftPaddingIndex();
		int zeros_right_index = wt.getRightPaddingIndex();
			
		// Points to most intensive data point in the signal
		SpectrumType::const_iterator it_max_pos;
			
		IntensityType max_value = 0.0;	// intensity of the signal peak
		UInt max_index   = 0;		// index of the signal peak (counted from the very first data point)

		int start = zeros_left_index + 2;
		int end  = zeros_right_index - 1;
		
		#ifdef DEBUG_FEATUREFINDER 
		// remove debug file (below we append entries only,
		// so if the file exists, we will get some quite confusing results)
		String fname = String("cwt_localmax_") + current_rt;
		if (File::exists(fname) )
		{
			File::remove(fname);
		}		
		#endif
	
		int i=0, j=0;
		for(i=start; i<end; ++i)
		{
			// Check for maximum in cwt at position i with cwt intensity > noise
			if( ((wt[i-1] - wt[i]  ) < 0) &&
					((wt[i] - wt[i+1]) > 0)  &&
					( wt[i]  > noise_level_cwt_ ) )
			{
				#ifdef DEBUG_FEATUREFINDER
				String fname = String("cwt_localmax_") + current_rt;
				ofstream gpfile( fname.c_str(), ios_base::app);
				gpfile << (first + i)->getMZ()  << "  " << cwt_[i] << endl;
				gpfile.close();
				#endif
				max_value=(first +  i)->getIntensity();
	
				// search for the corresponding maximum in the signal (consider the "radius" left and right adjacent points)
				int radius = 3;  // search radius for peaks in data
				int start_intervall = (( i - (int)radius) < 0 ) ? 0 : ( i - (int)radius);
				int end_intervall  = (( i + (int)radius) >= distance(first,last)) ? 0 : (i + (int)radius);
												
				for(j = start_intervall; j <= end_intervall; ++j)
				{
					if((first + j)->getIntensity() > max_value)
					{
						max_value = (first + j)->getIntensity();
						max_index = j;
					}
				}
	
				// check if this intensity is higher than the signal intensity threshold
				if (max_value > noise_level_signal_)
				{
					localmax.push_back(j);
				}
			}
		}
	}
	
	void MarrWaveletSeeder::sumUp_(SpectrumType& scan, UInt current_scan_index)
	{
		//Sum up those following scans that exist
		for ( UInt i=current_scan_index + 1
					; i < current_scan_index + 1 + (UInt)(param_.getValue("scans_to_sumup")) && i < traits_->getData().size()
				 	; ++i
				)
		{
			AlignAndSum_(scan,traits_->getData()[i]);
		}
	}
	
	void MarrWaveletSeeder::AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour)
	{
		if (scan.size() == 0 || neighbour.size() == 0) return;
	
		
	
		UInt index_newscan = 0;
		for (UInt k=0; k<neighbour.size(); ++k)
		{
			PeakType p			   = neighbour[k];
			CoordinateType mass = p.getMZ();
	
			while (scan[index_newscan].getMZ() < mass && index_newscan < scan.size())
				++index_newscan;
	
			// This seems to happen more frequently than expected -> quit the loop
			if (index_newscan >= scan.size() ) break;
	
			if (index_newscan > 0)
			{
				double left_diff   = fabs(scan[index_newscan-1].getMZ() - mass);
				double right_diff = fabs(scan[index_newscan].getMZ() - mass);

				// check which neighbour is closer
				if (left_diff < right_diff && (left_diff < mass_tolerance_alignment_) )
				{
					scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() + p.getIntensity() );
				}
				else if (right_diff < mass_tolerance_alignment_)
				{
					scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() + p.getIntensity() );
				}
			}
			else // no left neighbour available
			{
				double right_diff = fabs(scan[index_newscan].getMZ() - mass);
				if (right_diff < mass_tolerance_alignment_)
				{
					scan[index_newscan].setIntensity( scan[index_newscan].getIntensity() + p.getIntensity() );
				}
			}
		} // end for (all peaks in neighbouring scan)
	}
	
	UInt MarrWaveletSeeder::distanceToCharge_(CoordinateType dist)
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
// 		else if (dist < charge4_ub_ && dist > charge4_lb_)
// 		{
// 			return 4;
// 		}
// 		else if (dist < charge5_ub_ && dist > charge5_lb_)
// 		{
// 			return 5;
// 		}
		return 0;
	}

} // end of namespace OpenMS
