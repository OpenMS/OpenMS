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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PPExtender.h>

namespace OpenMS
{

PPExtender::PPExtender()
        : BaseExtender(), is_initialized_(false)
{
    name_ = PPExtender::getName();

    // lower and upper bounds for distances between isotopic peaks (defaults)
    // charge 1
    defaults_.setValue("charge1_ub",1.3f);
    defaults_.setValue("charge1_lb",0.8f);
    // charge 2
    defaults_.setValue("charge2_ub",0.65f);
    defaults_.setValue("charge2_lb",0.41f);
    // charge 3
    defaults_.setValue("charge3_ub",0.4f);
    defaults_.setValue("charge3_lb",0.2f);

    // tolerance in m/z for an monoisotopic peak in the previous scan
    defaults_.setValue("tolerance_mz",1.5f);

    defaults_.setValue("cwt_scale",0.6f);

    defaults_.setValue("noise_level_signal",750);

    defaults_.setValue("noise_level_cwt",6000);
	
	defaults_.setValue("scans_to_sumup",5);

    param_ = defaults_;
}

PPExtender::~PPExtender()
{}

const IndexSet& PPExtender::extend(const UnsignedInt /*seed_index*/)
{
    if (!is_initialized_)
        sweep_();

    region_.clear();

    is_initialized_ = true;

    // check if this cluster consists of one scan only and if so, just move to the next one
    while ((*curr_region_).second.scans_.size() == 1 &&
            curr_region_ != iso_map_.end() )
    {
        ++curr_region_;
    }

    if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
        throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);

    std::vector<UnsignedInt> next_region = (*curr_region_).second.peaks_;

    for (std::vector<UnsignedInt>::const_iterator cit = next_region.begin();
            cit != next_region.end();
            ++cit)
    {
        region_.add(*cit);
    }

    ++curr_region_;

    region_.sort();

    return region_;

} // end of extend(Unsigned int seed_index)

void PPExtender::sweep_()
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

    charge3_ub_    = param_.getValue("charge3_ub");
    charge3_lb_     = param_.getValue("charge3_lb");

    CoordinateType tolerance_mz = param_.getValue("tolerance_mz");

    cwt_scale_     = param_.getValue("cwt_scale");

    noise_level_signal_  = param_.getValue("noise_level_signal");

    noise_level_cwt_ = param_.getValue("noise_level_cwt");

    UnsignedInt current_charge     = 0;			// charge state of the current isotopic cluster
    CoordinateType mz_in_hash   = 0;			// used as reference to the current isotopic peak

    RawDataArrayType current_scan;	// stores the current scan of points
	
    unsigned int first_peak_in_scan = 0;

    // sweep through scans
    for (unsigned int curr_peak = 0; curr_peak < (nr_peaks-1); ++curr_peak)
    {
        CoordinateType current_rt =	traits_->getPeakRt(curr_peak);

        // check if new scan has begun
        if (current_rt != last_rt)
        {
            ////#ifdef DEBUG_FEATUREFINDER
            std::cout << "Next scan with rt: " << current_rt << std::endl;
            std::cout << "---------------------------------------------------------------------------" << std::endl;
            ////#endif

			// sum scan
			sumUp_(current_scan,curr_peak);
			
            // compute cwt for this scan
            cwt_.init(cwt_scale_, 0.0001);
            cwt_.transform(current_scan.begin(), current_scan.end(),1.);
		
			String fname = "cwt_" + String(last_rt);
            std::ofstream gpfile( fname.c_str() );
             for (int i=0;i<cwt_.getSize(); ++i)
             {
            	gpfile << (current_scan.begin() + i)->getPosition()[0] << "  " << cwt_[i] << std::endl;
             }
             gpfile.close();

			if  ( cwt_.getSize() == 0) 
			{
				std::cout << "Empty cwt => continue." << std::endl;
				continue;
			} 

            // search for maximal positions in the cwt and extract potential peaks
            std::vector<int> local_maxima;
			double cwt_sum = 0.0;
			for (int k = 0; k<cwt_.getSize();++k)
			{
				cwt_sum += cwt_[k];
			}
			std::cout << "Average strength in cwt: " << ( cwt_sum / cwt_.getSize() ) << std::endl;

			if (cwt_.size() == 0) continue;

            getMaxPositions_(current_scan.begin(), current_scan.end(), cwt_, local_maxima,curr_peak);

			fname = "scan_" + String(last_rt);
          	std::ofstream peakfile( fname.c_str() );
            for(unsigned k = 0; k<current_scan.size();++k)
            {
            	peakfile << current_scan[k].getPosition()[0] << " " << current_scan[k].getIntensity() << std::endl;
            }
            peakfile.close();

            current_scan.clear();	// prepare for next scan                        
			int nr_maxima = local_maxima.size();
			std::cout << "Local maxima in cwt : " << nr_maxima << std::endl;
			
            for (int z = 0; z< ( nr_maxima - 1); ++z)
			{
                // store the m/z of the current peak
                CoordinateType curr_mz         = traits_->getPeakMz(local_maxima[z]);
                CoordinateType dist2nextpeak = ( traits_->getPeakMz( local_maxima[z + 1]) - curr_mz);
				
				//std::cout << "Dist2nextPeak = " << dist2nextpeak  << std::endl;

                // test for different charge states
                current_charge = testDistance2NextPeak_(dist2nextpeak);

                if (current_charge > 0) // charger = 0 <=> no isotope
                {

					#ifdef DEBUG_FEATUREFINDER
//                     std::cout << "Isotopic pattern found ! " << std::endl;
//                     std::cout << "We are at: " << traits_->getPeakRt( local_maxima[z] ) << " " << curr_mz << std::endl;
					#endif

                    if (iso_last_scan.size() > 0)  // Did we find any isotopic cluster in the last scan?
                    {
                        // there were some isotopic clustures in the last scan...
                        std::vector<double>::iterator it = searchInScan_(iso_last_scan.begin(),iso_last_scan.end(),curr_mz);
                        double delta_mz = fabs(*it - curr_mz);

                        if ( delta_mz > tolerance_mz) // check if first peak of last cluster is close enough
                        {
                            mz_in_hash = curr_mz; // update current hash key

                            // create new isotopic cluster
							#ifdef DEBUG_FEATUREFINDER
//                             std::cout << "Last peak cluster too far, creating new cluster" << std::endl;
							#endif

                            iso_map_[mz_in_hash] = IsotopeCluster();

                            iso_map_[mz_in_hash].charge_  = current_charge;
                            iso_map_[mz_in_hash].scans_.push_back( traits_->getPeakRt( local_maxima[z] ) );
                        }
                        else
                        {
							#ifdef DEBUG_FEATUREFINDER
                          	//std::cout << "Found neighbouring peak with distance (m/z) " << delta_mz << std::endl;
							#endif

                            mz_in_hash = *it;	// retrieve new hash key
                            // save current rt and m/z
                            iso_map_[mz_in_hash].scans_.push_back( traits_->getPeakRt( local_maxima[z] ) );

							#ifdef DEBUG_FEATUREFINDER
                            //std::cout << "Cluster with " << iso_map_[mz_in_hash].peaks_.size() << " peaks retrieved." << std::endl;
							#endif

                        }

                    }
                    else // last scan did not contain any isotopic cluster
                    {
						#ifdef DEBUG_FEATUREFINDER
//                         std::cout << "Last scan was empty => creating new cluster." << std::endl;
//                         std::cout << "Creating new cluster at m/z: " << curr_mz << std::endl;
						#endif

                        mz_in_hash = curr_mz; // update current hash key

                        // create new isotopic cluster
                        iso_map_[mz_in_hash] = IsotopeCluster();
                        iso_map_[mz_in_hash].charge_  = current_charge;
                        iso_map_[mz_in_hash].scans_.push_back( traits_->getPeakRt( local_maxima[z] ) );
											
                    } // end if (iso_last_scan.size() > 0)

					#ifdef DEBUG_FEATUREFINDER
//                     std::cout << "Storing found peak in current isotopic cluster" << std::endl;
					#endif
					
					// walk a bit to the left: we include the five peaks to the left of the
					// monisotopic one into the feature region
					int ind = local_maxima[z];
					for (int p=0; p<=5; ++p)
					{
						if ( (ind-p) > 0 && traits_->getPeakFlag( (ind-p) ) == FeaFiTraits::UNUSED)
						{
						 	iso_map_[mz_in_hash].peaks_.push_back( ind-p );
						 	traits_->getPeakFlag(ind-p) = FeaFiTraits::INSIDE_FEATURE;
						}
					}
									
                   // store position of detected  cluster
                    iso_curr_scan.push_back(  mz_in_hash );
                    ++z;
					if ((local_maxima[z] - local_maxima[z-1]) > 1)
					{
						for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
						{
							iso_map_[mz_in_hash].peaks_.push_back( v );
							traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
						}
					}

                    iso_map_[mz_in_hash].peaks_.push_back( local_maxima[z] );
					traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;
					
                    //iso_curr_scan.push_back(traits_->getPeakMz( local_maxima[z] ));

                    // check distance to next peak
                    if ( (z+1) >= nr_maxima) break;
					
                    dist2nextpeak = ( traits_->getPeakMz( local_maxima[z+1] ) -  traits_->getPeakMz( local_maxima[z] ));
								
                    if (testDistance2NextPeak_(dist2nextpeak) != current_charge)
                    {
						// charge state has changed. Insert m/z of last maximum and continue.
						++z; 
						if ((local_maxima[z] - local_maxima[z-1]) > 1)
						{
							for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
							{
								iso_map_[mz_in_hash].peaks_.push_back( v );
								traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
							}
						}
						iso_map_[mz_in_hash].peaks_.push_back( local_maxima[z] );
						traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;
                        continue;
                    }
					
					if ( (z+1) >= nr_maxima) break;
					CoordinateType monoiso_mass = traits_->getPeakMz(local_maxima[z+1]);
					CoordinateType mass_diff         = 0.0;

                    // loop until end of isotopic pattern in this scan
                    while (testDistance2NextPeak_(dist2nextpeak) == current_charge && 
					           z < (nr_maxima-2) &&
							   mass_diff < 6 &&
							   traits_->getPeakFlag(local_maxima[z]) != FeaFiTraits::INSIDE_FEATURE)
                    {
                        iso_map_[mz_in_hash].peaks_.push_back( local_maxima[z] );				// save peak in cluster
						traits_->getPeakFlag(local_maxima[z]) = FeaFiTraits::INSIDE_FEATURE;
                        ++z;
						
						// between two local maxima in the cwt there might be several raw data points
						if ((local_maxima[z] - local_maxima[z-1]) > 1)
						{
							for (int v = local_maxima[z-1];v <  local_maxima[z]; ++v)
							{
								iso_map_[mz_in_hash].peaks_.push_back( v );
								traits_->getPeakFlag(v) = FeaFiTraits::INSIDE_FEATURE;
							}
						}
						
                        dist2nextpeak = ( traits_->getPeakMz( local_maxima[z+1] ) -  traits_->getPeakMz(local_maxima[z] )); // get distance to next peak
						mass_diff       = ( traits_->getPeakMz( local_maxima[z+1] ) - monoiso_mass);
                    } // end while(...)

                } // end of if (charge > 0)

                current_charge = 0; // reset charge
            } // end for (local maxima in cwt)
			
			// copy cluster information of least scan
            iso_last_scan = iso_curr_scan;
            iso_curr_scan.clear();
            last_rt = current_rt;		
           
        } // end if (current_rt != last_rt) => current scan is finished
		
		 // collect data points until new scan starts
         DRawDataPoint<1> p;
         p.setIntensity( traits_->getPeakIntensity(curr_peak) );
         p.getPosition()[0] = traits_->getPeakMz(curr_peak);
         current_scan.push_back(p);

         first_peak_in_scan = curr_peak;
			
    } // end for ( all data points ) => sweepline finished

    typedef std::map<CoordinateType,IsotopeCluster>::iterator HashIterator;
    std::vector<HashIterator> toDelete;

    std::cout << iso_map_.size() << " isotopic clusters were found." << std::endl;

    // remove cluster having less than 6 peaks or less than 3 scans
    for (HashIterator iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
    {
        if (iter->second.scans_.size() < 3 ||  iter->second.peaks_.size() < 6)
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

void PPExtender::getMaxPositions_( RawDataPointIterator first, RawDataPointIterator last, const ContinuousWaveletTransform& wt, std::vector<int>& localmax, unsigned int curr_peak)
{
		if (wt.getSize() == 0) return;
		
    	/// empirical values => check if this is a good idea
		//double noise_level_signal        = peak_bound_;
		//double noise_level_cwt       = 16000;   

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
					 String fname = "cwt_localmax_" + String(traits_->getPeakRt(curr_peak));
					 std::ofstream gpfile( fname.c_str(), std::ios_base::app); 
					 gpfile << (first + i)->getPos()  << "  " << cwt_[i] << std::endl;
					 gpfile.close();
//					std::cout << "Inserting : " << (curr_peak + i) << std::endl;
					max_value=(first +  i)->getIntensity();
					
					int radius = 3;  // search radius for peaks in data
					
        			// search for the corresponding maximum in the signal (consider the radius left and right adjacent points)
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
						localmax.push_back(curr_peak + i);			
					}

      		}
    	}
	
    
 }


void PPExtender::sumUp_(RawDataArrayType& scan, unsigned int current_index)
{
	unsigned int scans_to_collect = param_.getValue("scans_to_sumup");
	
	std::cout << "Summing up " << scans_to_collect << " scans." << std::endl;
	
	unsigned int scans_collected  = 0;
	
	CoordinateType current_rt = traits_->getPeakRt(current_index);
	CoordinateType last_rt      = current_rt;
	
	unsigned int nr_peaks      = traits_->getNumberOfPeaks();
	
	RawDataArrayType ascan;
	
	for (unsigned int j = current_index; j<nr_peaks; ++j)
	{
		current_rt = traits_->getPeakRt(j);
		if (current_rt != last_rt)
		{
			// new scan; store collected peaks			
			if (scans_collected == scans_to_collect)
			{
				// we have enough scans
				std::cout << "Collected " << scans_collected << " scans. Stop. " << std::endl;
				String fname = "scan_summed_" + String(traits_->getPeakRt(current_index));
				std::ofstream outfile( fname.c_str() );
             	for (unsigned int i=0;i<scan.size(); ++i)
            	 {
            		outfile << scan[i].getPosition()[0] << "  " << scan[i].getIntensity() << std::endl;
             	}
             	outfile.close();
				
				break;
			} 
			// sum up the intensities of the neighbouring scan
			AlignAndSum_(scan,ascan);		
			
			// prepare for next scan
			ascan.clear();			
			last_rt = current_rt;		
			++scans_collected;	
		} // end of if (new scan)
		
		// collect data points until new scan starts
        DRawDataPoint<1> p;
        p.setIntensity( traits_->getPeakIntensity(j) );
        p.getPosition()[0] = traits_->getPeakMz(j);
        ascan.push_back(p);
		
			
	}	// end of for (all peaks)
		
	//std::cout << "End of data. Collected " << scans_collected << " scans." << std::endl;
}




void PPExtender::AlignAndSum_(RawDataArrayType& scan, RawDataArrayType& neighbour)
{
			double mass_tolerance = 0.2;
			
			unsigned int index_newscan = 0;
			for (unsigned int k=0; k<neighbour.size(); ++k)
			{
				DRawDataPoint<1> p = neighbour[k];
				double mass              = p.getPosition()[0];
			
				while (scan[index_newscan].getPosition()[0] < mass && index_newscan < scan.size()) ++index_newscan;
	
				if (index_newscan >= scan.size() ) 
				{
					//std::cout << "End of scan. Adding " << ( neighbour1.size() -k ) << " peaks." << std::endl;
					// if the end of the scan is reached, add the remaining peaks
					//// Check if this is really needed . ///////
// 					while (k < neighbour1.size())
// 					{
// 						new_scan.push_back(neighbour1[k++]);
// 					}
					break;	// and quit the loop
				}
			
				if (index_newscan > 0)
				{
					double left_diff   = fabs(scan[index_newscan-1].getPosition()[0] - mass);
					double right_diff = fabs(scan[index_newscan].getPosition()[0] - mass);			

					// check which neighbour is closer
					if (left_diff < right_diff && (left_diff < mass_tolerance) )
					{
						scan[ (index_newscan-1) ].getIntensity() += p.getIntensity();
					}
					else if (right_diff < mass_tolerance)
					{
						scan[index_newscan].getIntensity() += p.getIntensity();;
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

UnsignedInt PPExtender::testDistance2NextPeak_(CoordinateType dist2nextpeak)
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
    else
    {
        return 0;
    }
}

} // end of namespace OpenMS
