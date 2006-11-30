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
// $Maintainer: Ole Schulz-Trieglaff , Rene Hussong$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>


#include <iostream>

namespace OpenMS
{

IsotopeWaveletSeeder::IsotopeWaveletSeeder():
        BaseSeeder(), is_initialized_(false),
        peak_cut_off_(5), waveletLength_(0),
        avMZSpacing_(0), min_spacing_(0),
        rt_votes_cutoff_(3), max_charge_(5),
				intensity_factor_(0), avg_intensity_factor_(0)
{
    name_ = IsotopeWaveletSeeder::getName();

//     defaults_.setValue("rtvotes_cutoff",6);
		defaults_.setValue("max_charge",5);
		defaults_.setValue("intensity_factor",1.5);
		defaults_.setValue("avg_intensity_factor",3.0);
		
    param_ = defaults_;
}

IsotopeWaveletSeeder::~IsotopeWaveletSeeder()
{}


IndexSet IsotopeWaveletSeeder::nextSeed() throw (NoSuccessor)
{
    if (!is_initialized_)
    {
        // reading params
//         rt_votes_cutoff_         = param_.getValue("rtvotes_cutoff");
				intensity_factor_        = param_.getValue("intensity_factor");
				avg_intensity_factor_ = param_.getValue("avg_intensity_factor");
				max_charge_            = param_.getValue("max_charge");
				
        // store charge states
        for (UnsignedInt i=1;i<=max_charge_;++i) charges_.push_back(i);           

        // compute spacings
        computeSpacings_();
				
				#ifdef DEBUG_FEATUREFINDER
        std::cout << "Average m/z spacing: " << avMZSpacing_ << std::endl;
        std::cout << "Minimal m/z spacing: " << min_spacing_ << std::endl;
				#endif

        waveletLength_ = (int) (peak_cut_off_/avMZSpacing_);
        generateGammaValues_();

        std::vector<DPeakArray<1, PeakType > >* pwts = NULL;
        std::vector<double>* wt_thresholds = NULL;

        UnsignedInt nr_scans = traits_->getData().size();

        for (UnsignedInt i=0; i<nr_scans; ++i)
        {
            CoordinateType current_rt = traits_->getData()[i].getRetentionTime();

            DPeakArray<1, PeakType> current_scan = traits_->getData()[i].getContainer();

						#ifdef DEBUG_FEATUREFINDER

            String filename = "scan_" + String( traits_->getData()[i].getRetentionTime() );
            std::ofstream outfile(filename.c_str());
            for (UnsignedInt k=0; k<current_scan.size();++k)
            {
                outfile << current_scan[k].getPos() << " " << current_scan[k].getIntensity() << std::endl;
            }
            outfile.close();
						#endif

            std::cout << "Spectrum " << i << " (" << current_rt << ") of " << nr_scans << std::endl;

            // store peak data, once for each charge state
            pwts = new std::vector<DPeakArray<1, PeakType > > (max_charge_, traits_->getData()[i].getContainer() );
            wt_thresholds = new std::vector<double> (charges_.size(), 0);

            // compute cwt
// 						std::cout << "FastMultiCorrelate" << std::endl;
            fastMultiCorrelate(current_scan, pwts, wt_thresholds);
            // compute scores of charge states
// 						std::cout << "identifyCharge" << std::endl;
            identifyCharge(*pwts, wt_thresholds, i, current_rt);
// 						std::cout << "Deleting pointer" << std::endl;
            delete (pwts);
            delete (wt_thresholds);
        } // end of for (each scan)

// 				std::cout << "done. " << std::endl;
        // filter detected isotopic pattern
        filterHashByRTVotes();

				std::cout << "Searching for indizes. " << std::endl;
        CoordinateType min_mass = traits_->getData().getMin().Y();

        for (SweepLineHash::const_iterator citer = hash_.begin();
                citer != hash_.end();
                ++citer)
        {
            std::cout << "m/z range: ";
            std::cout << (min_mass + (citer->first-1)*avMZSpacing_) << " ";
            std::cout << (min_mass+ (citer->first)*avMZSpacing_) << " " << std::endl;

            for (std::list<double>::const_iterator iter_cl2 = citer->
                    second.first.begin();
                    iter_cl2 != citer->second.first.end();
                    ++iter_cl2)
            {
                std::cout << "'rt: " <<  *iter_cl2 << " |  ";

                for (std::list<double>::const_iterator iter_l = citer->
                        second.second.begin();
                        iter_l != citer->second.second.end();
                        ++iter_l)
                {
                    std::cout << *iter_l << " ";
                }
                std::cout << " | " << std::endl;
            }	// end of times

            std::cout << "----------------------------------------------------------------" << std::endl;

        }

				hash_iter_ = hash_.begin();		
				
				is_initialized_ = true;
				
    } // end of if (!is_initialized_)

		if (hash_iter_ == hash_.end() || hash_.size() == 0 ) throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, 1);
	 			
		
		// compute mass we are searching for
		CoordinateType min_mass = traits_->getData().getMin().Y();
		CoordinateType mass_to_find = min_mass + (hash_iter_->first-1)*avMZSpacing_;
		IndexSet region_;		
		
		ScanIndexMSExperiment<MapType> scan_index = traits_->getScanIndex();
	
		// check all scans that support this isotopic pattern
		for (std::list<double>::const_iterator iter_cl2 = hash_iter_->second.first.begin(); 
		  		iter_cl2 != hash_iter_->second.first.end(); 
		   		++iter_cl2)
				{
				
				CoordinateType rt_to_find = *iter_cl2;
				
// 				std::cout << "Searching for rt: " << rt_to_find << std::endl;
				unsigned int current_scan = scan_index.getRank(rt_to_find);
				
				if (current_scan >= (scan_index.size()-1) )
				{
// 					std::cout << "Wrong scan number:" << current_scan  << std::endl;
					break;
				}
				
// 				std::cout << "Searching for " << current_scan << std::endl;
// 				std::cout << "Scan index " << scan_index_.size() << std::endl;
							
				PeakIterator scan_begin = scan_index[current_scan];
				PeakIterator scan_end   = scan_index[current_scan+1];				
					
				PeakIterator insert_iter = std::lower_bound(scan_begin,scan_end,mass_to_find,MZless());	
				int peak_index = insert_iter.getPeakNumber();
				
	
// 				std::cout << "Adding peak at mass " << traits_->getPeakMz(peak_index) << std::endl;
				
				CoordinateType miso_mass = traits_->getPeakMz(peak_index);
				CoordinateType miso_rt       =	traits_->getPeakRt(peak_index);			
				
				// walk a bit to the left from the monoisotopic peak
				for (int p=0; p<=10; ++p)
				{
					if ( (peak_index - p) > 0 
					     && traits_->getPeakFlag( (peak_index - p) ) == FeaFiTraits::UNUSED
						 && traits_->getPeakRt( (peak_index - p) ) == miso_rt
						 && traits_->getPeakMz( peak_index) - traits_->getPeakMz( (peak_index - p) ) < 2  )
					{
					 	region_.add( (peak_index - p) );
					}
				}
		
				CoordinateType mass_distance = 0;							
				int nr_peaks = traits_->getNumberOfPeaks();
				
				//we added the first peak (hopefully the monoisotopic one), now
				// we want to walk for about 10 Thompson (or Dalton or whatever)
				// into positive  (increasing) m/z direction
				while (mass_distance < 7
				           && peak_index < (nr_peaks-2)
						   && traits_->getPeakRt(peak_index) == miso_rt)
				{
					++peak_index;
					if (traits_->getPeakFlag(peak_index) == FeaFiTraits::UNUSED) region_.add(peak_index);
// 					std::cout << "Adding peak " << peak_index << std::endl;
						mass_distance = (traits_->getPeakMz(peak_index+1) - miso_mass);
// 					std::cout << "Current mass distance : " << mass_distance << std::endl;
				} 
// 				std::cout << "This scan is done." << std::endl;
								
	}	// for (std::list...)
	
	std::cout << "Done. Size of region: " << region_.size() << std::endl;	
	
	region_.sort();
	 ++hash_iter_;
	
    return region_;
}

void IsotopeWaveletSeeder::computeSpacings_()
{
    CoordinateType MZspacing_sum = 0;
    UnsignedInt rt_counts                 = 0;
		min_spacing_                            = INT_MAX;

    PeakIterator piter = traits_->getData().peakBegin();
    double last_mz   = piter->getPos();
    double last_rt     = piter.getRt();
    ++piter;

    for (;	piter != traits_->getData().peakEnd(); ++piter)
    {

        if (piter.getRt() != last_rt) // check if a new scan has begun
        {
            ++rt_counts;
            last_mz = piter->getPos();
						last_rt   = piter.getRt();
            continue;
        }

        double current_spacing = (piter->getPos() - last_mz);

        if (fabs(current_spacing) < min_spacing_)
            min_spacing_ = current_spacing;

        MZspacing_sum += current_spacing;
        last_mz = piter->getPos();
    }
    // compute average spacing of points in m/z
    UnsignedInt nr_of_peaks = traits_->getData().getSize();
    avMZSpacing_ = (MZspacing_sum / (double)(nr_of_peaks - rt_counts - 1)); //-1, since there are n data points and hence n-1 spacings

    return;
}

void IsotopeWaveletSeeder::generateGammaValues_()
{
    std::cout << "Precomputing the Gamma function ...";
    preComputedGamma_.clear();

    double query = 0;
    UnsignedInt counter=0;
    while (query <= (max_charge_*peak_cut_off_ +1) )
    {
        preComputedGamma_[counter] = tgamma (query);
        query += min_spacing_;
        ++counter;
    }

    std::cout << " done." << std::endl;
}

void IsotopeWaveletSeeder::fastMultiCorrelate(const DPeakArray<1, PeakType >& signal,
        std::vector<DPeakArray<1, PeakType > >* pwts,
        std::vector<double>* wt_thresholds)
{
    std::vector<DPeakArray<1, PeakType > >* res = pwts;
    unsigned int signal_size = signal.size();

    WaveletCollection phis (charges_.size(), std::vector<double> (waveletLength_)); 		//all necessary wavelets (by rows)

    double cumSpacing=0, cSpacing=0, realMass=0, lambda=0, w_sum=0, w_s_sum=0, max_w_monoi_intens=0.25,
                                  align_offset, tmp_pos, tmp_pos1; //helping variables

    std::list<UnsignedInt>::const_iterator charge_iter;
    UnsignedInt k=0; //helping variables
    double max=0;

    for (UnsignedInt i=0; i<signal_size; ++i)
    {
        //Now, let's sample the wavelets
        for (charge_iter=charges_.begin(), k=0; charge_iter!=charges_.end(); ++charge_iter, ++k)
        {
            cumSpacing=0;
            w_sum=0;
            w_s_sum=0;
            realMass = signal[i].getPos() * (*charge_iter);
//         		std::cout << "realMass: " << realMass << std::endl;
            lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
            //integration = phiRawInt (lambda, 1.0/(double)(*charge_iter));
            //std::cout << "Exact: " << integration << "\t" << lambda << "\t" << 1.0/(double)(*charge_iter) << std::endl;

            max_w_monoi_intens=0.25/(*charge_iter);

            //Align the maximum monoisotopic peak of the wavelet with some signal point
            UnsignedInt j=0;
            double last=0;
            while (cumSpacing < max_w_monoi_intens)
            {
                cSpacing = signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos();
//            			std::cout << "cSpacing " << cSpacing << std::endl;
// 								std::cout << cumSpacing << " " << max_w_monoi_intens << std::endl;
                last=cumSpacing;
                if (cSpacing < 0)
                    cumSpacing += avMZSpacing_;
                else //The "normal" case
                    cumSpacing += signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos();
                ++j;
            }

            align_offset = max_w_monoi_intens-last;

            cumSpacing=align_offset;
            for (UnsignedInt j=0; j<waveletLength_; ++j)
            {
                tmp_pos = signal[(i+j)%signal_size].getPos();
                tmp_pos1 = signal[(i+j+1)%signal_size].getPos();

                realMass = tmp_pos1 * (*charge_iter);
                lambda = getLambda (realMass); //Lambda determines the distribution (the shape) of the wavelet
                //std::cout << "phi before: " << phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter))
                //<< "\tphi after" << phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter))- integration << std::endl;
                //std::cout << "cumSpacing: " << cumSpacing << std::endl;
                phis[k][j] = phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter));
                w_sum += phis[k][j];
                w_s_sum += phis[k][j] * phis[k][j];
                //std::cout << "t: " << cumSpacing << "\t" << lambda << "\t" << 1.0/(double)(*charge_iter)
                //			<< "\t" << integration << "\t" << phis[k][j] << std::endl;
                cSpacing = tmp_pos1 - tmp_pos;
                //cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal_size-1).
                //Since this case is only of theoretical interest (we do not expect any important signals at the very end of the
                //spectrum), we simlply use the average spacing in that case.
                if (cSpacing < 0)
                    cumSpacing += avMZSpacing_;
                else //The "normal" case
                    cumSpacing += tmp_pos1 - tmp_pos;
            }
            max=-INT_MAX;
            for (UnsignedInt j=0; j<waveletLength_; ++j)
            {
                phis[k][j] -= (w_sum/(double)waveletLength_);
                if (phis[k][j] > max)
                    max = phis[k][j];
            }
            for (UnsignedInt j=0; j<waveletLength_; ++j)
                phis[k][j] /= max;

            (*wt_thresholds)[k] = w_s_sum;
        }

        std::vector<double> sums (charges_.size());
        k=0;
        for (UnsignedInt j=i; j<signal_size && k<phis[0].size(); ++j, ++k)
        {	//Since all wavelet functions have the same length, we can simply use phis[0].size()

            for (UnsignedInt m=0; m<charges_.size(); ++m)
            {
                sums[m] += signal[j].getIntensity()*phis[m][k];
            }
        }

        for (UnsignedInt l=0; l<i && k<phis[0].size(); ++l, ++k)
        {	//Since all wavelet functions have the same length, we can simply use phis[0].size()

            for (UnsignedInt m=0; m<charges_.size(); ++m)
                sums[m] += signal[l].getIntensity()*phis[m][k];
        }

        //Store the current convolution result
        for (UnsignedInt m=0; m<charges_.size(); ++m)
        {
//   					std::cout << "Wavelet transformed data: " << sums[m] << std::endl;
            (*res)[m][i].setIntensity(sums[m]);
        }

    }

} // end of fastMultiCorrelate(...)


void IsotopeWaveletSeeder::identifyCharge (const std::vector<DPeakArray<1, PeakType > >& candidates,
        std::vector<double>* wt_thresholds,
        const UnsignedInt scan,
        const CoordinateType current_rt)
{
    std::vector<double> int_mins (candidates[0].size(),INT_MIN), zeros (candidates[0].size(),0);
    WaveletCollection scoresC (candidates.size(), zeros);

    std::vector<unsigned int> start_indices, end_indices;
    //In order to determine the start and end indices, we first need to know the width of the region one should consider
    //to estimate the mean and the sd of the pattern candidate.
    //That region is defined by the position of the heighest amplitude +/- waveletLength_.

    typedef MSSpectrum<DRawDataPoint<2> >::ContainerType containerType;
    containerType::iterator iter;
    unsigned int start_index, end_index, c_index, i_iter; //Helping variables
    double seed_mz, c_check_point, c_val, c_av_intens;
    std::vector<bool> processed (candidates[0].size(), false);
    std::pair<int, int> c_between;
    int start, end, goto_left;


    for (unsigned int c=0; c<candidates.size(); ++c)
    {
    		std::cout << "Checking charge state " << (c+1) << std::endl;
			  processed = std::vector<bool> (candidates[0].size(), false); 				//Reset
        containerType c_candidate(candidates[c].size());
        
				//Ugly, but do not how to do this in a better (and easy) way
        for (unsigned int i=0; i<candidates[c].size(); ++i)
        {
            c_candidate[i].setPosition(DPosition<2>( candidates[c][i].getPos(),i));
            c_candidate[i].setIntensity( candidates[c][i].getIntensity() );
        }

        sort (c_candidate.begin(), c_candidate.end(), 	ReverseComparator< DRawDataPoint<2>::IntensityLess>() );
        c_av_intens = getAbsMean (candidates[c], 0, candidates[c].size());

        std::cout << "Average intensity: " << c_av_intens << std::endl;
				std::cout << "Threshold for wt: " << ((*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) << std::endl;

        for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter)
        {
     	   		if (iter->getIntensity() <= (*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) break;
        }

// 				std::cout << "writing output. " << std::endl;
				#ifdef DEBUG_FEATUREFINDER
        String filename = "cwt_" + String(current_rt) + "_charge_" + String(c+1);
        std::ofstream outfile(filename.c_str());
				containerType::iterator write_iter;
        for (write_iter=c_candidate.begin(); write_iter != c_candidate.end(); ++write_iter)
        {
            outfile << write_iter->getPos() << " " << write_iter->getIntensity() << std::endl;
        }
        outfile.close();
				#endif
				
        c_candidate.erase (iter, c_candidate.end());

        i_iter=0;
        for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter, ++i_iter)
        {
            // Retrieve index
            c_index = (int) (iter->getPosition().Y());

            if (processed[c_index]) continue;               

            start_index = c_index-waveletLength_-1;
            end_index = c_index+waveletLength_+1;
            seed_mz=iter->getPosition().X();

            //Catch impossible cases
            if (end_index >= candidates[c].size() || start_index > end_index)  continue;
               
            //Mark as processed
            for (unsigned int z=start_index; z<=end_index; ++z) processed[z] = true;                

            start=(-2*(peak_cut_off_ - 1))+1, end=(2*(peak_cut_off_ - 1))-1;
            goto_left = c_index - waveletLength_ - 1;

            for (int v=start; v<=end; ++v)
            {
                c_check_point = seed_mz+v*0.5/((double)c+1);
                c_between = getNearBys (scan, c_check_point, goto_left);
								
                if (c_between.first < 0 || c_between.second < 0) break;
                
                c_val = getInterpolatedValue (candidates[c][c_between.first].getPos(),
                                              						c_check_point,
                                              						candidates[c][c_between.second].getPos(),
                                              						candidates[c][c_between.first].getIntensity(),
                                              						candidates[c][c_between.second].getIntensity());
                
								if (fabs(c_val) < c_av_intens)
								{	
// 									std::cout << "(fabs(c_val) < c_av_intens)" << std::endl;
									continue;               
								}
                // What the hell is happening here?
                if (abs(v)%2 == 1) //i.e. whole
                {
                    scoresC[c][c_index] -= c_val;
                }
                else //i.e. peak
                {
                    scoresC[c][c_index] += c_val;
                }
            }

            if (scoresC[c][c_index] <= intensity_factor_*iter->getIntensity())
            {
                scoresC[c][c_index] = 0;
            }
        }
    }

		std::cout << "Scoring regions. " << std::endl;
    //Now, since we computed all scores, we can hash all mz positions

    unsigned int numOfCharges = candidates.size(), numOfMZPositions = candidates[0].size();
    //This is a vector telling us the next mz position in charge i we have to hash
    std::vector<unsigned int> positions (numOfCharges, 0);
// 		std::cout << "numOfCharges " << numOfCharges << std::endl;
		
    unsigned int c_hash_key, count_finished_charges=0;
    double allZero;
    DoubleList c_pair;
    std::list<double> c_list, c_fill_list;
    std::list<double>::iterator iter_cl, iter_cl_hash;
    std::pair<SweepLineHash::iterator, SweepLineHash::iterator> iter_hash;
    while (1) // ;-)
    {
        //Termination criterion
        //Test for every charge ...
        for (unsigned int c=0; c<numOfCharges; ++c)
        {
            //... if we hashed already all possible mz coordinates
            if (positions[c] >= numOfMZPositions)
            {
                //if so ... goto FINISHED_HASHING
                if (++count_finished_charges >= numOfCharges)
                    goto FINISHED_HASHING;
                positions[c] = -1;
            }
        }
        //End of Termination criterion


        //The hashing
        for (unsigned int c=0; c<numOfCharges; ++c)
        {
            if (positions[c] >= numOfMZPositions) continue;               

            //positions[c] also tells us the next candidate to hash
            c_list.push_back (scoresC[c][positions[c]++]);
        }

        for (unsigned int c=0; c<numOfCharges-1; ++c)
				{
            if (positions[c+1] != positions[c])
            {
                std::cout << "Quadro Zack!" << std::endl;
            }
				}

// 				std::cout << "scan  " << scan << " positions[0]-1 " << (positions[0]-1) <<  std::endl; 
// 				std::cout << "getPos: " << traits_->getData()[scan].getContainer()[positions[0]-1].getPos() << std::endl;	
// 				std::cout << "getMin().Y() " << traits_->getData().getMin().Y() << std::endl;
        c_hash_key = (unsigned int) ((traits_->getData()[scan].getContainer()[positions[0]-1].getPos() - traits_->getData().getMin().Y()) / avMZSpacing_);
// 				std::cout << "c_hash_key " << c_hash_key << std::endl;
				
        allZero=true;
        for (iter_cl=c_list.begin(); iter_cl!=c_list.end(); ++iter_cl)
				{
            if (*iter_cl != 0) allZero=false;                
				}
// 				std::cout << "all zero done " << c_hash_key << std::endl;
        if (!c_list.empty() && !allZero)
        {
            iter_hash=hash_.equal_range(c_hash_key);
// 						std::cout << "equal range " << c_hash_key << std::endl;
            while (iter_hash.first != iter_hash.second)
            {
								// m/z already in hash table 
//                 std::cout << "M/Z already in hash table ... " << std::endl;

                if (scan != 0)
								{
                    if (find(iter_hash.first->second.first.begin(), iter_hash.first->second.first.end(),
                             traits_->getData()[scan-1].getRetentionTime()) == iter_hash.first->second.first.end())
                    {
                        //i.e. there is no neighbouring entry before this retention time
                        //i.e. we can treat this case as if no entry is present in the hash
                        ++iter_hash.first;
                        continue;
                    }
								}	
                c_fill_list = iter_hash.first->second.first;
                c_fill_list.push_back(current_rt);
// 								std::cout << "Adding rt " << current_rt << std::endl;
                c_fill_list.unique(); // It might be the case, that we have several votes for the same RT and MZ 
								                           // by different charges => unique the list.
                
								for (iter_cl = c_list.begin(), iter_cl_hash = iter_hash.first->second.second.begin();
                        iter_cl != c_list.end();
                        ++iter_cl, ++iter_cl_hash)
                {
                    *iter_cl += *iter_cl_hash;
                }

                hash_.erase (iter_hash.first);
                c_pair = DoubleList (c_fill_list, c_list);
                goto FINISH;
            }

						// new entry
// 						std::cout << "New entry . " << std::endl;
            //Store RT and the corresponding score
            c_fill_list.clear();
            c_fill_list.push_back(current_rt);
            c_pair = DoubleList (c_fill_list, c_list);

FINISH:
// 						std::cout << "Inserting entry at " << c_hash_key << std::endl;
            hash_.insert (SweepLineHash::value_type(c_hash_key, c_pair));
        }

        c_list.clear();
    }

FINISHED_HASHING:

// 		std::cout << "Scoring done. " << std::endl;

    return;
}

double IsotopeWaveletSeeder::getAbsMean (const DPeakArray<1, PeakType >& signal,
        const unsigned int startIndex,
        const unsigned int endIndex) const
{
    double res=0;
    for (unsigned int i=startIndex; i<endIndex; ++i)
        res += fabs(signal[i].getIntensity());

    return (res/(double)(endIndex-startIndex+1));
}

void IsotopeWaveletSeeder::filterHashByRTVotes ()
{
    std::cout << "Hash size before filtering: " << hash_.size() << std::endl << std::endl;

    SweepLineHash::iterator iter_f, iter_b;
    std::vector<SweepLineHash::iterator> toDelete;
    std::list<double>::iterator iter_l;

    for (SweepLineHash::iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
    {
        if (iter->second.first.size() <= rt_votes_cutoff_)
            toDelete.push_back(iter);
    }

    for (unsigned int i=0; i<toDelete.size(); ++i)
        hash_.erase (toDelete[i]);

    std::cout << "Hash size after filtering: " << hash_.size() << std::endl << std::endl;

}

} // end of namespace OpenMS
