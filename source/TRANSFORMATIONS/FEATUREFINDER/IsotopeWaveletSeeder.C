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
// $Maintainer: Ole Schulz-Trieglaff , Rene Hussong$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>
#include <OpenMS/SYSTEM/StopWatch.h>

using namespace std;

namespace OpenMS
{
	IsotopeWaveletSeeder::IsotopeWaveletSeeder()
		: BaseSeeder(), 
			is_initialized_(false),
	    peak_cut_off_(5),
	    waveletLength_(0),
	    avMZSpacing_(0),
	    min_spacing_(0),
	    rt_votes_cutoff_(3),
	    intensity_factor_(0), 
			avg_intensity_factor_(0)
	{
    setName(getProductName());

		// minimal number of scans for an isotopic pattern
		defaults_.setValue("rtvotes_cutoff",5);

		// max and min charge states examined
		defaults_.setValue("max_charge",4);
		defaults_.setValue("min_charge",1);
		
		// intensity threshold in wt
		defaults_.setValue("intensity_factor",1.5);
		defaults_.setValue("avg_intensity_factor",3.0);
		
		// determines width of feature box
		defaults_.setValue("mass_tolerance_right",2.5);
		defaults_.setValue("mass_tolerance_left",6.0);
		
		// number of scans used for alignment
    defaults_.setValue("scans_to_sumup",5);
		defaults_.setValue("tolerance_scansum",0.1);
		defaults_.setValue("min_samplingrate",0.05);
				
    defaultsToParam_();
	}
	
	IsotopeWaveletSeeder::~IsotopeWaveletSeeder()
	{
	}

  IsotopeWaveletSeeder::IsotopeWaveletSeeder(const IsotopeWaveletSeeder& rhs)
    : BaseSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  IsotopeWaveletSeeder& IsotopeWaveletSeeder::operator= (const IsotopeWaveletSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }

	void IsotopeWaveletSeeder::updateMembers_()
	{
		rt_votes_cutoff_ = param_.getValue("rtvotes_cutoff");
		intensity_factor_ = param_.getValue("intensity_factor");
		avg_intensity_factor_ = param_.getValue("avg_intensity_factor");
		mass_tolerance_right_ = param_.getValue("mass_tolerance_right");
		mass_tolerance_left_ = param_.getValue("mass_tolerance_left");
		tolerance_scansum_ = param_.getValue("tolerance_scansum");
		
    // store charge states
    for (UnsignedInt i=(UnsignedInt)param_.getValue("min_charge"); i<=(UnsignedInt)param_.getValue("max_charge"); ++i)
    {
    	charges_.push_back(i);           
		}
	}
	
	FeaFiModule::IndexSet IsotopeWaveletSeeder::nextSeed() throw (NoSuccessor)
	{
    if (!is_initialized_)
    {
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
			StopWatch watch;
			watch.start();
			for (UnsignedInt i=0; i<nr_scans; ++i)
			{
				CoordinateType current_rt = traits_->getData()[i].getRetentionTime();
	
				//A copy of the scan has to be made as the intensities are modified later
				MapType::SpectrumType current_scan = traits_->getData()[i];
	
#ifdef DEBUG_FEATUREFINDER
				String filename =  String("scan_") + traits_->getData()[i].getRetentionTime();
				std::ofstream outfile(filename.c_str());
				for (UnsignedInt k=0; k<current_scan.size();++k)
				{
					outfile << current_scan[k].getPos() << " " << current_scan[k].getIntensity() << std::endl;
				}
				outfile.close();
#endif
							
				// align and sum
				sumUp_(current_scan,i);
										
				std::cout << "Spectrum " << i << " (" << current_rt << ") of " << nr_scans << std::endl;
	
				// store peak data, once for each charge state
				pwts = new std::vector<DPeakArray<1, PeakType > > (charges_.size(), traits_->getData()[i].getContainer() );
				wt_thresholds = new std::vector<double> (charges_.size(), 0);
	
				// compute wavelet transform
				fastMultiCorrelate(current_scan, pwts, wt_thresholds);
				// compute scores of charge states
				identifyCharge(*pwts, wt_thresholds, i);
				delete (pwts);
				delete (wt_thresholds);
			} // end of for (each scan)
	
			// filter detected isotopic pattern
			filterHashByRTVotes();
	
			// the following lines print the masses of all detected peak clusters
			// it takes some time if many clusters were found, so we skip it.
#ifdef DEBUG_FEATUREFINDER
//		 CoordinateType min_mass = traits_->getData().getMin().Y();
// 
//		 for (SweepLineHash::const_iterator citer = hash_.begin();
//				 citer != hash_.end();
//				 ++citer)
//		 {
// 				
//			 std::cout << "m/z range: ";
//			 std::cout << (min_mass + (citer->first-1)*avMZSpacing_) << " ";
//			 std::cout << (min_mass+ (citer->first)*avMZSpacing_) << " " << std::endl;
// 												
//			 for (std::list<UnsignedInt>::const_iterator iter_cl2 = citer->second.first.begin();
//				   iter_cl2 != citer->second.first.end();
//				   ++iter_cl2)
//			 {
//				 std::cout << "rt: " <<  traits_->getData()[*iter_cl2].getRetentionTime() << " |  ";
// 
//				 for (std::list<double>::const_iterator iter_l = citer->
//						 second.second.begin();
//						 iter_l != citer->second.second.end();
//						 ++iter_l)
//				 {
//					 std::cout << *iter_l << " ";
//				 }
//				 std::cout << " | " << std::endl;
//			 }	// end of times
// 
//			 std::cout << "----------------------------------------------------------------" << std::endl;
// 
//		 }
#endif
				watch.stop();
				std::cout << "Time spent for cwt: " << watch.getClockTime() << " [s] " << std::endl;
				hash_iter_ = hash_.begin();		
				
				is_initialized_ = true;
					
			} // end of if (!is_initialized_)
	
			if (hash_iter_ == hash_.end() || hash_.size() == 0 )
			{
				throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(1,1));
		 	}	
			
			// compute mass we are searching for
			CoordinateType min_mass = traits_->getData().getMin().Y();
			CoordinateType mass_to_find = min_mass + (hash_iter_->first-1)*avMZSpacing_;
			IndexSet region;		
				
			// check all scans that support this isotopic pattern
			for (std::list<UnsignedInt>::const_iterator iter_cl2 = hash_iter_->second.first.begin(); 
			  		iter_cl2 != hash_iter_->second.first.end(); 
			   		++iter_cl2)
			{
				UnsignedInt current_scan = *iter_cl2;
				
				if (current_scan >= (traits_->getData().size() ) )
				{
					// this shouldn't happen ;-)
					break;
				}			
				
				SpectrumType::ConstIterator insert_iter = std::lower_bound(traits_->getData()[current_scan].begin(),traits_->getData()[current_scan].end(),mass_to_find,PeakType::NthPositionLess<0>());	
				
				CoordinateType miso_mass = insert_iter->getPos();
			
				// The isotope wavelet operates on mass bins and not actual masses in the spectrum
				// We therefore need to check a couple of surrounding peaks in order to find the monoisotopic one
				// walk to the left
				for (UnsignedInt p=0; p<=10; ++p)
				{
					if ( miso_mass - (insert_iter - p)->getPos() < mass_tolerance_right_  )
					{
				 		region.insert( make_pair(current_scan,insert_iter - traits_->getData()[current_scan].begin()) );
					}
					//abort when the left border is reached
					if (insert_iter - p == traits_->getData()[current_scan].begin())
					{
						break;
					}
				}
	
			CoordinateType mass_distance = 0;
			// walk to the right
			while (mass_distance < mass_tolerance_left_ && insert_iter != (traits_->getData()[current_scan].begin()-1) )
			{
				++insert_iter;
				region.insert( make_pair(current_scan,insert_iter - traits_->getData()[current_scan].begin()) );
				mass_distance = ((insert_iter + 1)->getPos() - miso_mass);
			} 
													
		}		// for (std::list...)
		
		std::cout << "Done. Size of region: " << region.size() << std::endl;	
		
		++hash_iter_;
		
		return region;
	}
	
	void IsotopeWaveletSeeder::computeSpacings_()
	{
		CoordinateType MZspacing_sum = 0;
		min_spacing_ = INT_MAX;
		
		double current_spacing;
		for (MapType::ConstIterator it = traits_->getData().begin(); it != traits_->getData().end(); ++it)
		{
			if (it->size() == 0) continue;
			
			for (SpectrumType::ConstIterator it2 = it->begin() + 1; it2 != it->end(); ++it2)
			{
				current_spacing = it2->getPos() - (it2-1)->getPos();
				if (current_spacing < min_spacing_)
				{
					min_spacing_ = current_spacing;
				}
				MZspacing_sum += current_spacing;
			}
		}	
		
		// compute average spacing of points in m/z
		avMZSpacing_ = MZspacing_sum / (double)(traits_->getData().getSize() - traits_->getData().size());
	
		CoordinateType min_sampling = param_.getValue("min_samplingrate");
		if (avMZSpacing_ < min_sampling)
		{
			avMZSpacing_ = min_sampling; 
		}
	}
	
	void IsotopeWaveletSeeder::generateGammaValues_()
	{
		std::cout << "Precomputing the Gamma function ...";
		preComputedGamma_.clear();
	
		double query = 0;
		UnsignedInt counter=0;
		UnsignedInt max_charge = param_.getValue("max_charge");
			
		while (query <= (max_charge*peak_cut_off_ +1) )
		{
			preComputedGamma_[counter] = tgamma (query);
			query += min_spacing_;
			++counter;
		}
		std::cout << " done." << std::endl;
	}
	
	void IsotopeWaveletSeeder::fastMultiCorrelate(const SpectrumType& signal, std::vector<DPeakArray<1, PeakType > >* pwts, std::vector<double>* wt_thresholds)
	{
		std::vector<DPeakArray<1, PeakType > >* res = pwts;
		unsigned int signal_size = signal.size();
	
		WaveletCollection phis (charges_.size(), std::vector<double> (waveletLength_)); 		//all necessary wavelets (by rows)
	
		//helping variables
		double cumSpacing=0, cSpacing=0, realMass=0, lambda=0, w_sum=0, w_s_sum=0, max_w_monoi_intens=0.25, align_offset, tmp_pos, tmp_pos1;
	
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
	
				lambda = getLambda (realMass); 	//Lambda determines the distribution (the shape) of the wavelet
				 max_w_monoi_intens=0.25/(*charge_iter);
	
				//Align the maximum monoisotopic peak of the wavelet with some signal point
				UnsignedInt j=0;
				double last=0;
	
				while (cumSpacing < max_w_monoi_intens)
				{
					cSpacing = signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos();
					last=cumSpacing;
					
					if (cSpacing <= 0)
					{
						cumSpacing += avMZSpacing_;
					}
					else //The "normal" case
					{
						cumSpacing += signal[(i+j+1)%signal_size].getPos() - signal[(i+j)%signal_size].getPos();
					}					
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
				   
					phis[k][j] = phiRaw (cumSpacing, lambda, 1.0/(double)(*charge_iter));
					w_sum += phis[k][j];
					w_s_sum += phis[k][j] * phis[k][j];
					cSpacing = tmp_pos1 - tmp_pos;
					
					//cSpacing might get negative, as soon as the wavelet approaches the end of the signal (if i=signal_size-1).
					//Since this case is only of theoretical interest (we do not expect any important signals at the very end of the
					//spectrum), we simlply use the average spacing in that case.
					if (cSpacing < 0)
					{
						cumSpacing += avMZSpacing_;
					}
					else //The "normal" case
					{
						cumSpacing += tmp_pos1 - tmp_pos;
					}
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
	
					//Since all wavelet functions have the same length, we can simply use phis[0].size()
			for (UnsignedInt j=i; j<signal_size && k<phis[0].size(); ++j, ++k)
			{	
				for (UnsignedInt m=0; m<charges_.size(); ++m)
				{
					sums[m] += signal[j].getIntensity()*phis[m][k];
				}
			}
	
			for (UnsignedInt l=0; l<i && k<phis[0].size(); ++l, ++k)
			{
	
				for (UnsignedInt m=0; m<charges_.size(); ++m)
				{
					sums[m] += signal[l].getIntensity()*phis[m][k];
				}
			}
	
			//Store the convolution result
			for (UnsignedInt m=0; m<charges_.size(); ++m)
			{
				(*res)[m][i].setIntensity(sums[m]);
			}
	
		}
	
	} // end of fastMultiCorrelate(...)
	
	
	void IsotopeWaveletSeeder::identifyCharge (const std::vector<DPeakArray<1, PeakType > >& candidates, std::vector<double>* wt_thresholds, UnsignedInt scan)
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
	
		ChargeVector::const_iterator charge_iter = charges_.begin();
	
		for (unsigned int c=0; c<candidates.size(); ++c)
		{		
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
					
			for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter)
			{
		 		if (iter->getIntensity() <= (*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) break;
			}
	
#ifdef DEBUG_FEATUREFINDER
			std::cout << "Checking charge state " << *charge_iter << std::endl;
			std::cout << "Average intensity: " << c_av_intens << std::endl;
			std::cout << "Threshold for wt: " << ((*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) << std::endl;
	
			// write debug output
			CoordinateType current_rt = traits_->getData()[scan].getRetentionTime();
			String filename = String("cwt_") + current_rt + "_charge_" + (c+1);
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
					
					if (fabs(c_val) < c_av_intens) continue;  
								
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
						
			++charge_iter;
		} // end of for (unsigned int c=0; c<candidates.size(); ++c)
			
			
		//Now, since we computed all scores, we can hash all mz positions
		UnsignedInt numOfCharges = candidates.size(), numOfMZPositions = candidates[0].size();
		//This is a vector telling us the next mz position in charge i we have to hash
		std::vector<UnsignedInt> positions (numOfCharges, 0);
			
		UnsignedInt c_hash_key, count_finished_charges=0;
		bool allZero;
		DoubleList c_pair;
	
		// c_list is now charge_scores
		// c_fill_list is scans
		std::list<double> charge_scores;
		std::list<UnsignedInt> scans;
			
		std::list<double>::iterator iter_cl, iter_cl_hash;
		std::pair<SweepLineHash::iterator, SweepLineHash::iterator> iter_hash;
			
		double score_sum = 0.0;
			
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
				charge_scores.push_back (scoresC[c][positions[c]++]);
			}
	
			for (unsigned int c=0; c<numOfCharges-1; ++c)
			{
				if (positions[c+1] != positions[c])
				{
					std::cout << "Quadro Zack!" << std::endl;
				}
			}
		   
			// generate hash key 
			c_hash_key = (UnsignedInt) ((traits_->getData()[scan].getContainer()[positions[0]-1].getPos() - traits_->getData().getMin().Y()) / avMZSpacing_);
	
			allZero=true;
			for (iter_cl=charge_scores.begin(); iter_cl!=charge_scores.end(); ++iter_cl)
			{
				if (*iter_cl != 0) allZero=false;				
			}
	
			if (!charge_scores.empty() && !allZero)
			{
				iter_hash=hash_.equal_range(c_hash_key);
	
				while (iter_hash.first != iter_hash.second)
				{
					// m/z already in hash table 
					if (scan != 0)
					{
						if (	find(iter_hash.first->second.first.begin(), iter_hash.first->second.first.end(), (scan-1) ) == iter_hash.first->second.first.end() )
						{
							//i.e. there is no neighbouring entry before this retention time
							//i.e. we can treat this case as if no entry is present in the hash
							++iter_hash.first;
							continue;
						}
					}	
					scans = iter_hash.first->second.first;
					scans.push_back(scan);
									
					// It might be the case, that we have several votes for the same RT and MZ 
					// by different charges => unique the list.
					scans.unique(); 
					
					for (iter_cl = charge_scores.begin(), iter_cl_hash = iter_hash.first->second.second.begin();
							iter_cl != charge_scores.end();
							++iter_cl, ++iter_cl_hash)
					{
						*iter_cl += *iter_cl_hash;
					}
	
					hash_.erase (iter_hash.first);
									
					// normalize scores
					score_sum = 0.0;
					for (std::list<double>::iterator it = charge_scores.begin();
						  it != charge_scores.end();
								++it)
					{
						score_sum += *it;
					}
					
					for (std::list<double>::iterator it = charge_scores.begin();
						  it != charge_scores.end();
								++it)
					{
						*it /= score_sum;
					}
																	
					c_pair = DoubleList (scans, charge_scores);
					goto FINISH;
				}
	
				// new hash entry
				scans.clear();
				scans.push_back(scan);
				
				// normalize scores
				score_sum = 0.0;
				for (std::list<double>::iterator it = charge_scores.begin();
						  it != charge_scores.end();
								++it)
				{
					score_sum += *it;
				}
						
				for (std::list<double>::iterator it = charge_scores.begin();
					  it != charge_scores.end();
							++it)
				{
					*it /= score_sum;
				}
													
				c_pair = DoubleList (scans, charge_scores);
	
				FINISH:
				hash_.insert (SweepLineHash::value_type(c_hash_key, c_pair));
			}
	
			charge_scores.clear();
		}
	
		FINISHED_HASHING:
	
		// done
		return;
	}

	double IsotopeWaveletSeeder::getAbsMean (const DPeakArray<1, PeakType >& signal, UnsignedInt startIndex, UnsignedInt endIndex) const
	{
	  double res=0;
	  for (unsigned int i=startIndex; i<endIndex; ++i)
	  {
	  	res += fabs(signal[i].getIntensity());
		}
	  return (res/(double)(endIndex-startIndex+1));
	}

	void IsotopeWaveletSeeder::filterHashByRTVotes()
	{
	  std::cout << "Hash size before filtering: " << hash_.size() << std::endl << std::endl;
	
	  SweepLineHash::iterator iter_f, iter_b;
	  std::vector<SweepLineHash::iterator> toDelete;
	  std::list<double>::iterator iter_l;
	
	  for (SweepLineHash::iterator iter = hash_.begin(); iter != hash_.end(); ++iter)
	  {
	    if (iter->second.first.size() <= rt_votes_cutoff_)
	    {
	        toDelete.push_back(iter);
	  	}
	  }
	
	  for (unsigned int i=0; i<toDelete.size(); ++i)
	  {
	  	hash_.erase (toDelete[i]);
		}
		
		std::cout << "Hash size after filtering: " << hash_.size() << std::endl << std::endl;
	}

	void IsotopeWaveletSeeder::sumUp_(SpectrumType& scan, UnsignedInt current_scan_index)
	{
		//Sum up those following scans that exist
		for ( UnsignedInt i=current_scan_index + 1
					; i < current_scan_index + 1 + (UnsignedInt)(param_.getValue("scans_to_sumup")) && i < traits_->getData().size()
				 	; ++i
				)
		{
			AlignAndSum_(scan,traits_->getData()[i]);
		}
	}
	
	void IsotopeWaveletSeeder::AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour)
	{
		if (scan.size() == 0 || neighbour.size() == 0) return;
	
		double mass_tolerance = 0.1;
	
		UnsignedInt index_newscan = 0;
		for (UnsignedInt k=0; k<neighbour.size(); ++k)
		{
			PeakType p			   = neighbour[k];
			CoordinateType mass = p.getPos();
	
			while (scan[index_newscan].getPos() < mass && index_newscan < scan.size())
				++index_newscan;
	
			// This seems to happen more frequently than expected -> quit the loop
			if (index_newscan >= scan.size() ) break;
	
			if (index_newscan > 0)
			{
				double left_diff   = fabs(scan[index_newscan-1].getPos() - mass);
				double right_diff = fabs(scan[index_newscan].getPos() - mass);
				// 					cout << "Checking neighbours: " << left_diff << " " << right_diff << endl;
	
				// check which neighbour is closer
				if (left_diff < right_diff && (left_diff < mass_tolerance) )
				{
					// 						cout << "Left. Old intensity: " << scan[ (index_newscan-1) ].getIntensity() << endl;
					scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() + p.getIntensity() );
					// 						cout << "Left. New intensity: " << scan[ (index_newscan-1) ].getIntensity() << endl;
				}
				else if (right_diff < mass_tolerance)
				{
					// 						cout << "Right. Old intensity: " << scan[ (index_newscan) ].getIntensity() << endl;
					scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() + p.getIntensity() );
					// 						cout << "Right. New intensity: " << scan[ (index_newscan) ].getIntensity() << endl;
				}
			}
			else // no left neighbour available
			{
				double right_diff = fabs(scan[index_newscan].getPos() - mass);
				if (right_diff < mass_tolerance)
				{
					scan[index_newscan].getIntensity() += p.getIntensity();
				}
			}
		} // end for (all peaks in neighbouring scan)
	}

} // end of namespace OpenMS
