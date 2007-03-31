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
// $Maintainer: Ole Schulz-Trieglaff, Rene Hussong$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletSeeder.h>

using namespace std;

namespace OpenMS
{
	IsotopeWaveletSeeder::IsotopeWaveletSeeder()
		: BaseSweepSeeder(), 
			wavelet_initialized_(false),
	    peak_cut_off_(5),
	    waveletLength_(0),
	    avMZSpacing_(0),
	    min_spacing_(0),
	    intensity_factor_(0), 
			avg_intensity_factor_(0)
	{
    setName(getProductName());

		// max and min charge states examined
		defaults_.setValue("max_charge",4);
		defaults_.setValue("min_charge",1);
		
		// intensity threshold in cwt
		defaults_.setValue("intensity_factor",1.5);
		defaults_.setValue("avg_intensity_factor",3.0);
				
    defaultsToParam_();
	}
	
	IsotopeWaveletSeeder::~IsotopeWaveletSeeder()
	{
	}

  IsotopeWaveletSeeder::IsotopeWaveletSeeder(const IsotopeWaveletSeeder& rhs)
    : BaseSweepSeeder(rhs),
    	wavelet_initialized_(false)
  {
    updateMembers_();
  }
  
  IsotopeWaveletSeeder& IsotopeWaveletSeeder::operator= (const IsotopeWaveletSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSweepSeeder::operator=(rhs);
    wavelet_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }

	void IsotopeWaveletSeeder::updateMembers_()
	{
		// update member of base class first
		BaseSweepSeeder::updateMembers_();

		intensity_factor_         = param_.getValue("intensity_factor");
		avg_intensity_factor_   = param_.getValue("avg_intensity_factor");
		
		// delete old charge states
		charges_.clear();
		
    // store charge states
    for (UInt i=(UInt)param_.getValue("min_charge"); i<=(UInt)param_.getValue("max_charge"); ++i)
    {
    	charges_.push_back(i);           
		}
	}
	
	IsotopeWaveletSeeder::ScoredMZVector IsotopeWaveletSeeder::detectIsotopicPattern_(SpectrumType& scan ) 
	{
			if (!wavelet_initialized_)
			{
      	// compute spacings
      	computeSpacings_();
			
				#ifdef DEBUG_FEATUREFINDER
      	std::cout << "Average m/z spacing: " << avMZSpacing_ << std::endl;
      	std::cout << "Minimal m/z spacing: " << min_spacing_ << std::endl;
				#endif
	
	    	waveletLength_ = (Int) (peak_cut_off_/avMZSpacing_);
				generateGammaValues_();
	
				wavelet_initialized_ = true;
			}
			
			std::vector<DPeakArray<1, PeakType > >* pwts = NULL;
			std::vector<double>* wt_thresholds = NULL;
	
			// store peak data, once for each charge state
			pwts = new std::vector<DPeakArray<1, PeakType > > (charges_.size(), scan.getContainer() );
			wt_thresholds = new std::vector<double> (charges_.size(), 0);
	
			// compute wavelet transform
			fastMultiCorrelate_(scan, pwts, wt_thresholds);
			// compute scores of charge states
			ScoredMZVector scmzvec = identifyCharge_(*pwts, wt_thresholds, scan);
			
			delete (pwts);
			delete (wt_thresholds);
	
			return scmzvec;
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
				current_spacing = it2->getMZ() - (it2-1)->getMZ();
				if (current_spacing < min_spacing_)
				{
					min_spacing_ = current_spacing;
				}
				MZspacing_sum += current_spacing;
			}
		}	
		
		// compute average spacing of points in m/z
		avMZSpacing_ = MZspacing_sum / (double)(traits_->getData().getSize() - traits_->getData().size());
	}
	
	void IsotopeWaveletSeeder::generateGammaValues_()
	{
		#ifdef DEBUG_FEATUREFINDER
		std::cout << "Precomputing the Gamma function ...";
		preComputedGamma_.clear();
		#endif
		
		double query = 0;
		UInt counter=0;
		UInt max_charge = param_.getValue("max_charge");
			
		while (query <= (max_charge*peak_cut_off_ +1) )
		{
			preComputedGamma_[counter] = tgamma (query);
			query += min_spacing_;
			++counter;
		}
		
		#ifdef DEBUG_FEATUREFINDER
		std::cout << " done." << std::endl;
		#endif
	}
	
	void IsotopeWaveletSeeder::fastMultiCorrelate_(const SpectrumType& signal, std::vector<DPeakArray<1, PeakType > >* pwts, std::vector<double>* wt_thresholds)
	{
		std::vector<DPeakArray<1, PeakType > >* res = pwts;
		unsigned int signal_size = signal.size();
	
		WaveletCollection phis (charges_.size(), std::vector<double> (waveletLength_)); 		//all necessary wavelets (by rows)
	
		//helping variables
		double cumSpacing=0, cSpacing=0, realMass=0, lambda=0, w_sum=0, w_s_sum=0, max_w_monoi_intens=0.25, align_offset, tmp_pos, tmp_pos1;
	
		std::list<UInt>::const_iterator charge_iter;
		UInt k=0; //helping variables
		double max=0;

		for (UInt i=0; i<signal_size; ++i)
		{
	
			//Now, let's sample the wavelets
			for (charge_iter=charges_.begin(), k=0; charge_iter!=charges_.end(); ++charge_iter, ++k)
			{
				cumSpacing=0;
				w_sum=0;
				w_s_sum=0;
				realMass = signal[i].getMZ() * (*charge_iter);
	
				lambda = getLambda_(realMass); 	//Lambda determines the distribution (the shape) of the wavelet
				max_w_monoi_intens=0.25/(*charge_iter);
	
				//Align the maximum monoisotopic peak of the wavelet with some signal point
				UInt j=0;
				double last=0;
	
				while (cumSpacing < max_w_monoi_intens)
				{
					cSpacing = signal[(i+j+1)%signal_size].getMZ() - signal[(i+j)%signal_size].getMZ();
					last=cumSpacing;
					
					if (cSpacing <= 0)
					{
						cumSpacing += avMZSpacing_;
					}
					else //The "normal" case
					{
						cumSpacing += signal[(i+j+1)%signal_size].getMZ() - signal[(i+j)%signal_size].getMZ();
					}					
					++j;
				}
	
				align_offset = max_w_monoi_intens-last;
	
				cumSpacing=align_offset;
	
				for (UInt j=0; j<waveletLength_; ++j)
				{
					tmp_pos = signal[(i+j)%signal_size].getMZ();
					tmp_pos1 = signal[(i+j+1)%signal_size].getMZ();
	
					realMass = tmp_pos1 * (*charge_iter);
					lambda = getLambda_(realMass); //Lambda determines the distribution (the shape) of the wavelet
				   
					phis[k][j] = phiRaw_(cumSpacing, lambda, 1.0/(double)(*charge_iter));
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
				for (UInt j=0; j<waveletLength_; ++j)
				{
					phis[k][j] -= (w_sum/(double)waveletLength_);
					if (phis[k][j] > max)
						max = phis[k][j];
				}
				for (UInt j=0; j<waveletLength_; ++j)
					phis[k][j] /= max;
	
				(*wt_thresholds)[k] = w_s_sum;
			}
	
			std::vector<double> sums (charges_.size());
			k=0;
	
			//Since all wavelet functions have the same length, we can simply use phis[0].size()
			for (UInt j=i; j<signal_size && k<phis[0].size(); ++j, ++k)
			{	
				for (UInt m=0; m<charges_.size(); ++m)
				{
					sums[m] += signal[j].getIntensity()*phis[m][k];
				}
			}
	
			for (UInt l=0; l<i && k<phis[0].size(); ++l, ++k)
			{
	
				for (UInt m=0; m<charges_.size(); ++m)
				{
					sums[m] += signal[l].getIntensity()*phis[m][k];
				}
			}
	
			//Store the convolution result
			for (UInt m=0; m<charges_.size(); ++m)
			{
				(*res)[m][i].setIntensity(sums[m]);
			}
	
		}
	
	} // end of fastMultiCorrelate(...)
	
	
	IsotopeWaveletSeeder::ScoredMZVector IsotopeWaveletSeeder::identifyCharge_(const std::vector<DPeakArray<1, PeakType > >& candidates, 
	                                                                                                                     std::vector<double>* wt_thresholds, 
																																								                                       SpectrumType& scan)
	{
		std::vector<IntensityType> int_mins (candidates[0].size(),INT_MIN), zeros (candidates[0].size(),0);
		WaveletCollection scoresC (candidates.size(), zeros);
		
		ScoredMZVector scmzvec;	// container of scored mz positions
			
		std::vector<UInt> start_indices, end_indices;
		
		//In order to determine the start and end indices, we first need to know the width of the region one should consider
		//to estimate the mean and the sd of the pattern candidate.
		//That region is defined by the position of the highest amplitude +/- waveletLength_.
			
		TempContainerType::iterator iter;
		UInt start_index, end_index, c_index, i_iter; //Helping variables
		CoordinateType seed_mz, c_check_point, c_val;
		IntensityType c_av_intens;
		std::vector<bool> processed (candidates[0].size(), false);
		std::pair<Int, Int> c_between;
		Int start, end, goto_left;
	
		ChargeVector::const_iterator charge_iter = charges_.begin();
	
		for (UInt c=0; c<candidates.size(); ++c)
		{		
			processed = std::vector<bool> (candidates[0].size(), false); 				//Reset
			TempContainerType c_candidate(candidates[c].size());
			
			//Ugly, but do not how to do this in a better (and easy) way
			for (UInt i=0; i<candidates[c].size(); ++i)
			{
				c_candidate[i].setPosition(DPosition<2>( candidates[c][i].getPos(),i));
				c_candidate[i].setIntensity( candidates[c][i].getIntensity() );
			}
				
			sort (c_candidate.begin(), c_candidate.end(), 	ReverseComparator< RawDataPoint2D::IntensityLess>() );
			c_av_intens = getAbsMean_(candidates[c], 0, candidates[c].size());
					
			for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter)
			{
		 		if (iter->getIntensity() <= (*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) break;
			}
	
			#ifdef DEBUG_FEATUREFINDER
			std::cout << "Checking charge state " << *charge_iter << std::endl;
			std::cout << "Average intensity: " << c_av_intens << std::endl;
			std::cout << "Threshold for wt: " << ((*wt_thresholds)[c]*avg_intensity_factor_*c_av_intens) << std::endl;
	
			// write debug output
			CoordinateType current_rt = scan.getRetentionTime();
			String filename = String("isowavcwt_") + current_rt + "_charge_" + (c+1);
			std::ofstream outfile(filename.c_str());
			TempContainerType::iterator write_iter;
			for (write_iter=c_candidate.begin(); write_iter != c_candidate.end(); ++write_iter)
			{
				outfile << write_iter->getPosition().getX() << " " << write_iter->getIntensity() << std::endl;
			}
			outfile.close();
			#endif
					
			c_candidate.erase (iter, c_candidate.end());
			
			i_iter=0;
			for (iter=c_candidate.begin(); iter != c_candidate.end(); ++iter, ++i_iter)
			{
				// Retrieve index
				c_index = (Int) (iter->getPosition().getY());	
	
				if (processed[c_index])
				{
				 continue;			   
				}
				
				start_index = c_index-waveletLength_-1;
				end_index = c_index+waveletLength_+1;
				seed_mz=iter->getPosition().getX();	// m/z is still X coordinate
	
				//Catch impossible cases
				if (end_index >= candidates[c].size() || start_index > end_index) 
				{
				 	continue;
				}  
				//Mark as processed
				for (UInt z=start_index; z<=end_index; ++z) processed[z] = true;				
	
				start=(-2*(peak_cut_off_ - 1))+1, end=(2*(peak_cut_off_ - 1))-1;
				goto_left = c_index - waveletLength_ - 1;
					
				for (Int v=start; v<=end; ++v)
				{
					c_check_point = seed_mz+v*0.5/((double)c+1);
					c_between = getNearBys_(scan, c_check_point, goto_left);
					
					if (c_between.first < 0 || c_between.second < 0) 
					{
						break;
					}
					
					c_val = getInterpolatedValue_(candidates[c][c_between.first].getPos(),
													  						c_check_point,
													  						candidates[c][c_between.second].getPos(),
													  						candidates[c][c_between.first].getIntensity(),
													  						candidates[c][c_between.second].getIntensity());
																				
					if (fabs(c_val) < c_av_intens)
					{
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
						
			++charge_iter;
		} // end of for (unsigned int c=0; c<candidates.size(); ++c)
			
			
		//Now, since we computed all scores, we can hash all mz positions
		UInt num_charges = candidates.size(); 
		Int num_mzpos = candidates[0].size();
		//This vector tells us the next mz position in charge i we have to hash
		std::vector<Int> positions (num_charges, 0);
			
 		UInt count_finished_charges=0;
		bool allZero;
		DoubleList c_pair;
	
		// c_list is now charge_scores
		// c_fill_list is scans
		std::vector<ProbabilityType> charge_scores;
		std::list<UInt> scans;
			
		std::vector<double>::iterator iter_cl, iter_cl_hash;

		while (1) // ;-)
		{
			//Termination criterion
			//Test for every charge ...
			for (UInt c=0; c<num_charges; ++c)
			{
				//... if we hashed already all possible mz coordinates
				if (positions[c] >= num_mzpos)
				{
					if (++count_finished_charges >= num_charges) return scmzvec;
			
					positions[c] = -1;
				}
			}
			//End of Termination criterion
		
			//The hashing
			for (UInt c=0; c<num_charges; ++c)
			{
				if (positions[c] >= num_mzpos) continue;			   
	
				//positions[c] also tells us the next candidate to hash
				charge_scores.push_back (scoresC[c][positions[c]++]);
			}
	
			for (UInt c=0; c<num_charges-1; ++c)
			{
				OPENMS_PRECONDITION(positions[c+1] == positions[c], "positions[c+1] != positions[c]");
			}
		
			allZero=true;
			for (iter_cl=charge_scores.begin(); iter_cl!=charge_scores.end(); ++iter_cl)
			{
				if (*iter_cl != 0) allZero=false;				
			}
	
			if (!charge_scores.empty() && !allZero)
			{
			
				ScoredChargeType sc_charge;	
				for (UInt i = 0; i < charge_scores.size(); ++i)
				{				
					if (charge_scores[i] >= sc_charge.second)
					{
						sc_charge.first       = (i+1);
						sc_charge.second	= charge_scores[i];
					}
				}
				scmzvec.push_back( make_pair( (positions[0]-1),sc_charge) );	
				
			}
	
			charge_scores.clear();
		}

		// done
		return scmzvec;
	}

	double IsotopeWaveletSeeder::getAbsMean_(const DPeakArray<1, PeakType >& signal, UInt startIndex, UInt endIndex) const
	{
	  double res=0;
	  for (unsigned int i=startIndex; i<endIndex; ++i)
	  {
	  	res += fabs(signal[i].getIntensity());
		}
	  return (res/(double)(endIndex-startIndex+1));
	}

} // end of namespace OpenMS
