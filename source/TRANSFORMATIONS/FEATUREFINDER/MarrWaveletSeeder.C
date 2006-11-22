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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/MarrWaveletSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

MarrWaveletSeeder::MarrWaveletSeeder()
        : BaseSeeder(), is_initialized_(false)
{
    name_ = MarrWaveletSeeder::getName();

    // lower and upper bounds for distances between isotopic peaks (defaults)
    // charge 1
    defaults_.setValue("charge1_ub",1.2f);
    defaults_.setValue("charge1_lb",0.8f);
    // charge 2
    defaults_.setValue("charge2_ub",0.65f);
    defaults_.setValue("charge2_lb",0.41f);
    // charge 3
    defaults_.setValue("charge3_ub",0.4f);
    defaults_.setValue("charge3_lb",0.2f);
		// charge 4
		defaults_.setValue("charge4_ub",0.27f);
    defaults_.setValue("charge4_lb",0.23f);
		// charge 5
		defaults_.setValue("charge5_ub",0.19f);
    defaults_.setValue("charge5_lb",0.21f);

    // tolerance in m/z for an monoisotopic peak in the previous scan
    defaults_.setValue("tolerance_mz",0.5f);

		// params for the cwt
    defaults_.setValue("cwt_scale",0.3f);
    defaults_.setValue("noise_level_signal",750);
    defaults_.setValue("noise_level_cwt",6000);
		defaults_.setValue("high_peak_intensity_factor",4);
		
		// number of scans used for alignment
		defaults_.setValue("scans_to_sumup",5);
		
		// minimum number of scan per isotopic cluster
		defaults_.setValue("min_number_scans",5);
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_peaks",20);

    param_ = defaults_;
}

MarrWaveletSeeder::~MarrWaveletSeeder()
{}

Index MarrWaveletSeeder::nextSeed() throw (NoSuccessor)
{ 
		if (!is_initialized_)  sweep_();		
		
    is_initialized_ = true;

    if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
        throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);

    std::vector<UnsignedInt> next_region = (*curr_region_).second.peaks_;
    ++curr_region_;

    // find maximum intensity in the decteted region
    UnsignedInt max_index           = 0;
    CoordinateType max_intensity = 0.0;

    // find maximum of region
    for (std::vector<UnsignedInt>::const_iterator citer = next_region.begin();
            citer != next_region.end();
            ++citer)
    {
        if (traits_->getPeakIntensity(*citer) > max_intensity)
        {
            max_index = *citer;
            max_intensity = traits_->getPeakIntensity(*citer);
        }

    } // end of for (Index in region)
    traits_->getPeakFlag(max_index) = FeaFiTraits::SEED;

    std::cout << "Found seed with intensity " << traits_->getPeakIntensity(max_index) ;
    std::cout << " at (" << traits_->getPeakRt(max_index) << "," << traits_->getPeakMz(max_index) << ")" << std::endl;

    return max_index;

}

void MarrWaveletSeeder::sweep_()
{
	  // stores the monoisotopic peaks of isotopic clusters
    std::vector<CoordinateType> iso_last_scan;			// in the previous scan
    std::vector<CoordinateType> iso_curr_scan;			// in the last scan
    
    // retrieve values for accepted peaks distances
    charge1_ub_    = param_.getValue("charge1_ub");
    charge1_lb_     = param_.getValue("charge1_lb");

    charge2_ub_    = param_.getValue("charge2_ub");
    charge2_lb_     = param_.getValue("charge2_lb");

    charge3_ub_    = param_.getValue("charge3_ub");
    charge3_lb_     = param_.getValue("charge3_lb");
		
		charge4_ub_    = param_.getValue("charge4_ub");
    charge4_lb_     = param_.getValue("charge4_lb");
		
		charge5_ub_    = param_.getValue("charge5_ub");
    charge5_lb_     = param_.getValue("charge5_lb");
		
		// thresholds for cwt
		CoordinateType tolerance_mz = param_.getValue("tolerance_mz");
    double cwt_scale										= param_.getValue("cwt_scale");
    noise_level_signal_                = param_.getValue("noise_level_signal");
    noise_level_cwt_                   = param_.getValue("noise_level_cwt");
		high_peak_intensity_factor_   = param_.getValue("high_peak_intensity_factor");

    UnsignedInt current_charge    = 0;			// charge state of the current isotopic cluster
    CoordinateType mz_in_hash  = 0;			// used as reference to the current isotopic peak
		CoordinateType current_rt      = 0;			// retention time of current scan
		
    ContainerType current_scan;						// stores the current scan of points
		UnsignedInt current_offset = 0; 				// size offset of the current scan
		
    // sweep through scans
		UnsignedInt nr_scans = traits_->getData().size();
		
    for (UnsignedInt currscan_index = 0; currscan_index < nr_scans; ++currscan_index)
    {
					current_scan  = traits_->getData()[currscan_index].getContainer();
					current_rt       = traits_->getData()[currscan_index].getRetentionTime();
					current_offset = traits_->getData().getSpectraLengths()[currscan_index];

          #ifdef DEBUG_FEATUREFINDER
          std::cout << "---------------------------------------------------------------------------" << std::endl;
					std::cout << "I am in " << currscan_index << " of " << nr_scans << " scans. " << std::endl;
          std::cout << "Retention time: " << current_rt << std::endl;
					std::cout << "Size offset: " << current_offset << std::endl;
					#endif

            String fname = "scan_" + String(current_rt);
            std::ofstream peakfile( fname.c_str() );
            for(unsigned k = 0; k<current_scan.size();++k)
            {
                peakfile << current_scan[k].getPosition()[0] << " " << current_scan[k].getIntensity() << std::endl;
            }
            peakfile.close();

            // align and sum 
					 sumUp_(current_scan,currscan_index,current_rt);

            // compute cwt for this scan, use spacing of 0.0001
						// TODO: estimate spacing for this scan.
            cwt_.init(cwt_scale, 0.0001);
            cwt_.transform(current_scan.begin(), current_scan.end(),1.0);

            fname = "cwt_" + String(current_rt);
            std::ofstream gpfile( fname.c_str() );
            for (int i=0;i<cwt_.getSize(); ++i)
            {
            	gpfile << (current_scan.begin() + i)->getPosition()[0] << "  " << cwt_[i] << std::endl;
            }
            gpfile.close();

            if  ( cwt_.getSize() == 0)
            {
                std::cout << "Empty cwt for this scan => continue." << std::endl;
                continue;
            }

           // compute average intensity of cwt
            double cwt_avg = 0.0;
						double cwt_var = 0.0;
            for (int k = 0; k<cwt_.getSize();++k)
            {
              cwt_avg += cwt_[k];
            }
            cwt_avg /= cwt_.getSize();
						for (int p=0; p<cwt_.getSize();++p)
						{
							cwt_var += (cwt_[p] - cwt_avg) * (cwt_[p] - cwt_avg) ;	
						}
						cwt_var /= (cwt_.getSize() - 1);
            std::cout << "Average strength in cwt: " << cwt_avg<< std::endl;
						std::cout << "Variance in cwt: " << cwt_var << std::endl;

						 // search for maximal positions in the cwt and extract potential peaks
            std::vector<int> local_maxima;
            getMaxPositions_(current_scan.begin(), current_scan.end(), cwt_, local_maxima,current_rt);

            int nr_maxima = local_maxima.size();
            std::cout << "# local maxima in cwt : " << nr_maxima << std::endl;

            for (int z = 0; z< ( nr_maxima - 1); ++z)
            {
                // store the m/z of the current peak
                CoordinateType curr_mz         = current_scan[ local_maxima[z] ].getPos();
                CoordinateType dist2nextpeak = (current_scan[ local_maxima[ (z + 1) ] ].getPos() - curr_mz);

                std::cout << "Dist2nextPeak = " << dist2nextpeak  << std::endl;

                // test for different charge states
                current_charge = testDistance2NextPeak_(dist2nextpeak);
								
								std::cout << "Predicted charge: " << current_charge << std::endl;

                if (current_charge > 0) // charger = 0 <=> no isotope
                {

										#ifdef DEBUG_FEATUREFINDER
                		std::cout << "Isotopic pattern found ! " << std::endl;
                    std::cout << "We are at: " << current_rt << " " << curr_mz << std::endl;
										#endif

                    if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
                    {
                        // there were some isotopic clustures in the last scan...
                        std::vector<double>::iterator it = searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
                        CoordinateType delta_mz = fabs(*it - curr_mz);

                        if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
                        {
                            mz_in_hash = curr_mz; // update current hash key

                            // create new isotopic cluster
														#ifdef DEBUG_FEATUREFINDER
                            std::cout << "Last peak cluster too far, creating new cluster" << std::endl;
														#endif

                            iso_map_[mz_in_hash] = IsotopeCluster();

                            iso_map_[mz_in_hash].charge_  = current_charge;
                            iso_map_[mz_in_hash].scans_.push_back( current_rt );
                        }
                        else
                        {
														#ifdef DEBUG_FEATUREFINDER
                            std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
														#endif

                            mz_in_hash = *it;	// retrieve new hash key
                            // save current rt and m/z
                            iso_map_[mz_in_hash].scans_.push_back(current_rt);

														#ifdef DEBUG_FEATUREFINDER
                            std::cout << "Cluster with " << iso_map_[mz_in_hash].peaks_.size() << " peaks retrieved." << std::endl;
														#endif

                        }

                    }
                    else // last scan did not contain any isotopic cluster
                    {
												#ifdef DEBUG_FEATUREFINDER
                        std::cout << "Last scan was empty => creating new cluster." << std::endl;
                        std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
												#endif

                        mz_in_hash = curr_mz; // update current hash key

                        // create new isotopic cluster
                        iso_map_[mz_in_hash] = IsotopeCluster();
                        iso_map_[mz_in_hash].charge_  = current_charge;
                        iso_map_[mz_in_hash].scans_.push_back( current_rt );

                    } // end if (iso_last_scan.size() > 0)

										#ifdef DEBUG_FEATUREFINDER
                    std::cout << "Storing found peak in current isotopic cluster" << std::endl;
										#endif

                    // walk a bit to the left: we include the five peaks to the left of the
                    // monisotopic one into the feature region
                    int ind = local_maxima[z];
                    for (int p=0; p<=5; ++p)
                    {
                        if ( (ind-p) > 0 && traits_->getPeakFlag( current_offset + (ind-p) ) == FeaFiTraits::UNUSED)
                        {
                            iso_map_[mz_in_hash].peaks_.push_back( current_offset + ind-p );
//                             traits_->getPeakFlag(ind-p) = FeaFiTraits::INSIDE_FEATURE;
                        }
                    }
										std::cout << "Storing m/z " << mz_in_hash << std::endl;
                    // store position of detected  cluster
                    iso_curr_scan.push_back(  mz_in_hash );
                    ++z;
                    if ((local_maxima[z] - local_maxima[z-1]) > 1)
                    {
                        for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
                        {
                            iso_map_[mz_in_hash].peaks_.push_back( current_offset + v );
//                             traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
                        }
                    }

                    iso_map_[mz_in_hash].peaks_.push_back( current_offset + local_maxima[z] );
//                     traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;

//                     iso_curr_scan.push_back(traits_->getPeakMz( local_maxima[z] ));
										iso_curr_scan.push_back( current_scan[ local_maxima[z] ].getPos() );

                    // check distance to next peak
                    if ( (z+1) >= nr_maxima) break;

                    dist2nextpeak = ( current_scan[ local_maxima[z+1] ].getPos()  -   current_scan[ local_maxima[z] ].getPos() );

                    if (testDistance2NextPeak_(dist2nextpeak) != current_charge)
                    {
                        // charge state has changed. Insert m/z of last maximum and continue.
                        ++z;
                        if ((local_maxima[z] - local_maxima[z-1]) > 1)
                        {
                            for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
                            {
                                iso_map_[mz_in_hash].peaks_.push_back( current_offset + v );
//                                 traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
                            }
                        }
                        iso_map_[mz_in_hash].peaks_.push_back( current_offset + local_maxima[z] );
//                         traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;
                        continue;
                    }

                    if ( (z+1) >= nr_maxima) break;
										
                    CoordinateType monoiso_mass = current_scan[ local_maxima[z+1] ].getPos();
                    CoordinateType mass_diff         = 0.0;

                    // loop until end of isotopic pattern in this scan
                    while (testDistance2NextPeak_(dist2nextpeak) == current_charge &&
                            z < (nr_maxima-2) &&
                            mass_diff < 6 &&
                            traits_->getPeakFlag(local_maxima[z]) == FeaFiTraits::UNUSED)
                    {
                        iso_map_[mz_in_hash].peaks_.push_back( current_offset + local_maxima[z] );				// save peak in cluster
//                         traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;
                        ++z;

                        // between two local maxima in the cwt there might be several raw data points
                        if ((local_maxima[z] - local_maxima[z-1]) > 1)
                        {
                            for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
                            {
                                iso_map_[mz_in_hash].peaks_.push_back( current_offset + v );
//                                 traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
                            }
                        }

                        dist2nextpeak = (current_scan[ local_maxima[z+1] ].getPos() -  current_scan[ local_maxima[z] ].getPos() ); // get distance to next peak
                        mass_diff       = ( current_scan[ local_maxima[z+1] ].getPos() - monoiso_mass);
                    } // end while(...)

                } // end of if (charge > 0)
                else if ( cwt_[ local_maxima[z] ] > (high_peak_intensity_factor_ *cwt_avg) )
                {
                    // We have reached a data point with high intensity in the wavelet
                    // transform. This might be an interesting signal even without isotopic pattern,
                    // so we collect some surrounding data points and see if we meet this point
                    // in the next scan again....
                    std::cout << "HIGH PEAK IN CWT !!! " << std::endl;

                    UnsignedInt this_peak = /*current_offset +*/ local_maxima[z];

                    // extend this data point
                    CoordinateType this_mass     = current_scan[ this_peak ].getPos();
                    CoordinateType this_intensity = current_scan[ this_peak ].getIntensity();
										
                    CoordinateType next_mass     = this_mass;
                    CoordinateType next_intensity = this_intensity;

                    iso_map_[this_mass].charge_  = 1;
                    iso_map_[this_mass].scans_.push_back( current_rt );
                    iso_map_[this_mass].peaks_.push_back(current_offset + this_peak);

                    // walk to the left
                    while ( fabs(next_mass - this_mass) < 4 && next_intensity > (this_intensity * 0.003)  && this_peak > 0)
                    {
                        --this_peak;
                        iso_map_[this_mass].peaks_.push_back(current_offset + this_peak);

                        next_mass     = current_scan[ this_peak ].getPos();  /*traits_->getPeakMz(this_peak);*/
                        next_intensity = current_scan[ this_peak ].getIntensity();
                    }

                    UnsignedInt nr_peaks = traits_->getNumberOfPeaks();
										 this_peak = local_maxima[z];
										 
                    // walk to the right
                    while ( fabs(next_mass - this_mass) < 4 && next_intensity > (this_intensity * 0.003)  && this_peak < nr_peaks)
                    {
                        ++this_peak;
                        iso_map_[this_mass].peaks_.push_back(this_peak);

                       next_mass     = current_scan[ this_peak ].getPos();  /*traits_->getPeakMz(this_peak);*/
                        next_intensity = current_scan[ this_peak ].getIntensity();
                    }

                }
                current_charge = 0; // reset charge

        } // end for each (local maximum in cwt)

				 // copy cluster information of last scan
				iso_last_scan = iso_curr_scan;
        iso_curr_scan.clear();
				
				std::cout << "Numer of potential seeds in this scan: " << iso_last_scan.size() << std::endl;

		 } // end for ( all scans ) => sweepline finished

    typedef std::map<CoordinateType,IsotopeCluster>::iterator HashIterator;
    std::vector<HashIterator> toDelete;

    std::cout << iso_map_.size() << " isotopic clusters were found." << std::endl;
		
		UnsignedInt min_number_scans = param_.getValue("min_number_scans");		
		UnsignedInt min_number_peaks = param_.getValue("min_number_peaks");

    // remove cluster having less than 6 peaks or less than 3 scans
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

} // end of void sweep_()

void MarrWaveletSeeder::getMaxPositions_( RawDataPointIterator first, RawDataPointIterator last, const ContinuousWaveletTransform& wt, std::vector<int>& localmax, CoordinateType current_rt)
{
		if (wt.getSize() == 0) return;
		
		int zeros_left_index  = wt.getLeftPaddingIndex();
    int zeros_right_index = wt.getRightPaddingIndex();

    	// Points to most intensive data point in the signal
    	RawDataPointIterator it_max_pos;
    	double max_value = 0.0;

    	// Given direction, start the search from left or right
    	int start = zeros_left_index + 2;
    	int end  = zeros_right_index - 1;

    	int i=0, j=0;
    	for(i=start; i<end; i+=1)
    	{
      		// Check for maximum in cwt at position i with cwt intensity > noise 
     		if( ((wt[i-1] - wt[i]  ) < 0) &&
			    	((wt[i] - wt[i+1]) > 0)  && 
						( wt[i]  > noise_level_cwt_ ) ) 
      			{					
					 		String fname = "cwt_localmax_" + String(current_rt);
					 		std::ofstream gpfile( fname.c_str(), std::ios_base::app); 
					 		gpfile << (first + i)->getPos()  << "  " << cwt_[i] << std::endl;
					 		gpfile.close();
							
// 							std::cout << "Inserting : " << (curr_peak + i) << std::endl;
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
          				}
        			}
					
					if (max_value > noise_level_signal_)
					{
						localmax.push_back(/*curr_peak +*/ i);			
					}

      		}
    	}
	
    
 }

void MarrWaveletSeeder::sumUp_(ContainerType& scan, UnsignedInt& current_scan_index, CoordinateType& current_rt)
{
	UnsignedInt scans_to_collect = param_.getValue("scans_to_sumup");
	ContainerType ascan;	
	
	std::cout << "Summing up " << scans_to_collect << " scans." << std::endl;
	
	if (current_scan_index + scans_to_collect > traits_->getData().size() )
	{
		UnsignedInt diff = (current_scan_index + scans_to_collect) - traits_->getData().size();
		
		if (diff > scans_to_collect)
			scans_to_collect = 0;
		else
			scans_to_collect -= diff;
		
		std::cout << "End of map. Summing up " << scans_to_collect << " scans. "<< std::endl;	
	}
		
	for (UnsignedInt i=current_scan_index; i < (current_scan_index + scans_to_collect); ++i)
	{
		ascan = traits_->getData()[i].getContainer();
		AlignAndSum_(scan,ascan);		
	} 
	
	#ifdef DEBUG_FEATUREFINDER
	// write scan to file
// 	String fname = "scan_summed_" + String(current_rt);
// 	std::ofstream outfile( fname.c_str() );
//   for (unsigned int i=0;i<scan.size(); ++i)
//   {
//   	outfile << scan[i].getPosition()[0] << "  " << scan[i].getIntensity() << std::endl;
//   }
//   outfile.close();
	#endif
	
}




void MarrWaveletSeeder::AlignAndSum_(ContainerType& scan, ContainerType& neighbour)
{
			if (scan.size() == 0 || neighbour.size() == 0) return;

			double mass_tolerance = 0.1;
			
			UnsignedInt index_newscan = 0;
			for (UnsignedInt k=0; k<neighbour.size(); ++k)
			{
				PeakType p               = neighbour[k];
				CoordinateType mass = p.getPosition()[0];
			
				while (scan[index_newscan].getPos() < mass && index_newscan < scan.size()) ++index_newscan;
	
				// This seems to happen more frequently than expected
				if (index_newscan >= scan.size() )  break;	// and quit the loop
					
				if (index_newscan > 0)
				{
					double left_diff   = fabs(scan[index_newscan-1].getPosition()[0] - mass);
					double right_diff = fabs(scan[index_newscan].getPosition()[0] - mass);			
// 					std::cout << "Checking neighbours: " << left_diff << " " << right_diff << std::endl;

					// check which neighbour is closer
					if (left_diff < right_diff && (left_diff < mass_tolerance) )
					{
// 						std::cout << "Left. Old intensity: " << scan[ (index_newscan-1) ].getIntensity() << std::endl;
						scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() + p.getIntensity() );
// 						std::cout << "Left. New intensity: " << scan[ (index_newscan-1) ].getIntensity() << std::endl;
					}
					else if (right_diff < mass_tolerance)
					{
// 						std::cout << "Right. Old intensity: " << scan[ (index_newscan) ].getIntensity() << std::endl;
						scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() + p.getIntensity() );
// 						std::cout << "Right. New intensity: " << scan[ (index_newscan) ].getIntensity() << std::endl;
					}
				} 
				else // no left neighbour available
				{
					double right_diff = fabs(scan[index_newscan].getPosition()[0] - mass);		
					if (right_diff < mass_tolerance)
					{
						scan[index_newscan].getIntensity() += p.getIntensity();
					}					
				} 
								
			} // end for (all peaks in neighbouring scan)
}

UnsignedInt MarrWaveletSeeder::testDistance2NextPeak_(CoordinateType dist2nextpeak)
{

    if (dist2nextpeak < charge1_ub_ && dist2nextpeak > charge1_lb_)
    {
        return 1;
    }
    else if (dist2nextpeak < charge2_ub_ && dist2nextpeak > charge2_lb_)
    {
        return 2;
    }
    else if (dist2nextpeak < charge3_ub_ && dist2nextpeak > charge3_lb_)
    {
        return 3;
    }
		else if (dist2nextpeak < charge4_ub_ && dist2nextpeak > charge4_lb_)
    {
        return 4;
    }
		else if (dist2nextpeak < charge5_ub_ && dist2nextpeak > charge5_lb_)
    {
        return 5;
    }
    
		return 0;
    
}

}
