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
	    signal_avg_factor_(0), 
			cwt_avg_factor_(0),
			null_var_(),
			n_null_()
	{
    setName(getProductName());

		// max and min charge states examined
		defaults_.setValue("max_charge",4);
		defaults_.setValue("min_charge",1);
		
		// intensity threshold in cwt
		defaults_.setValue("signal_avg_factor",3.0);
		defaults_.setValue("cwt_avg_factor",3.0);
				
    defaultsToParam_();
	}
	
	IsotopeWaveletSeeder::~IsotopeWaveletSeeder()
	{
	}

  IsotopeWaveletSeeder::IsotopeWaveletSeeder(const IsotopeWaveletSeeder& rhs)
    : BaseSweepSeeder(rhs),
    	wavelet_initialized_(false),
			peak_cut_off_(rhs.peak_cut_off_),
			waveletLength_(rhs.peak_cut_off_),
			avMZSpacing_(rhs.avMZSpacing_),
			min_spacing_(rhs.min_spacing_),
			signal_avg_factor_(rhs.signal_avg_factor_),
			cwt_avg_factor_(rhs.cwt_avg_factor_),
			null_var_(rhs.null_var_),
			n_null_(rhs.n_null_)
  {
    updateMembers_();
  }
  
  IsotopeWaveletSeeder& IsotopeWaveletSeeder::operator= (const IsotopeWaveletSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSweepSeeder::operator=(rhs);
		
    wavelet_initialized_ = false;
    peak_cut_off_         = rhs.peak_cut_off_;
		waveletLength_       = rhs.waveletLength_;
		avMZSpacing_        = rhs.avMZSpacing_;
		min_spacing_         = rhs.min_spacing_;
		signal_avg_factor_  = rhs.signal_avg_factor_;
		cwt_avg_factor_     = rhs.cwt_avg_factor_;
		null_var_               = rhs.null_var_;
		n_null_                 = rhs.n_null_;
		
    updateMembers_();
    
    return *this;
  }

	void IsotopeWaveletSeeder::updateMembers_()
	{
		// update member of base class first
		BaseSweepSeeder::updateMembers_();

		signal_avg_factor_  = param_.getValue("signal_avg_factor");
		cwt_avg_factor_     = param_.getValue("cwt_avg_factor");
		
		// delete old charge states
		charges_.clear();
		
    // store charge states
    for (UInt i=(UInt)param_.getValue("min_charge"); i<=(UInt)param_.getValue("max_charge"); ++i)
    {
    	charges_.push_back(i);           
		}
		
		null_var_.resize(charges_.size(),0.0);
		n_null_.resize(charges_.size(),0);
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
			std::vector<DPeakArray<PeakType > >* pwts = NULL;
			// store peak data, once for each charge state
			pwts = new std::vector<DPeakArray<PeakType > > (charges_.size(), scan.getContainer() );
			
			// compute wavelet transform
			fastMultiCorrelate_(scan, pwts);
			// compute scores of charge states
			ScoredMZVector scmzvec = identifyCharge_(*pwts,scan);

			delete (pwts);
			
			return scmzvec;
	}
	
	void IsotopeWaveletSeeder::computeSpacings_()
	{
		CoordinateType MZspacing_sum = 0;
		min_spacing_ = numeric_limits<CoordinateType>::max();
		
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
		#endif
		
		preComputedGamma_.clear();
		
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
			
	void IsotopeWaveletSeeder::computeNullVariance_(const DPeakArray<PeakType >& cwt, const UInt charge_index )
	{
			
		CoordinateType first_mass = (cwt.size() > 0) ? cwt[ (cwt.size()-1) ].getMZ() : 0.0 ;
		CoordinateType mass_diff  = 0.0;
		UInt j = (cwt.size() > 0) ? (cwt.size() - 1) : 0;
		
		IntensityType mean = 0.0;
		IntensityType S       = 0.0;
		UInt n                     = 0;
						
		for (; mass_diff < 6.0 && j > 0; --j)
		{
			++n;
  		IntensityType delta = cwt[j].getIntensity() - mean;
  		mean = mean + delta/n;
  		S += delta*( cwt[j].getIntensity()  - mean);
			
			mass_diff    = cwt[j].getMZ() - first_mass;
			
		}
		n_null_[charge_index] = n;		
		null_var_[charge_index] = S / (n-1);
	}
		
	void IsotopeWaveletSeeder::fastMultiCorrelate_(const SpectrumType& signal, std::vector<DPeakArray<PeakType > >* pwts)
	{
		std::vector<DPeakArray<PeakType > >* res = pwts;
		UInt signal_size = signal.size();
	
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
				max = numeric_limits<double>::min();
				for (UInt j=0; j<waveletLength_; ++j)
				{
					phis[k][j] -= (w_sum/(double)waveletLength_);
					if (phis[k][j] > max)
						max = phis[k][j];
				}
				
				for (UInt j=0; j<waveletLength_; ++j)
					phis[k][j] /= max;
				
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
	
	
	IsotopeWaveletSeeder::ScoredMZVector IsotopeWaveletSeeder::identifyCharge_(std::vector<DPeakArray<PeakType > >& candidates,
																																								                                       SpectrumType& scan)
	{
	
		ScoredMZVector scmzvec;	 // scored positions
				
		std::vector<IntensityType> cwt_thresholds(candidates.size(),0.0);		// threshold for cwt intensities (one for each charge state)
		std::vector<UInt> last_pattern(candidates.size(),0);													// the last index where a pattern of was found (one for each charge)

		for (UInt c = 0; c < candidates.size(); ++c)
		{		
			IntensityType avg_cwt  = 0;
			
			computeNullVariance_(candidates[c],c);
						
			// compute average intensity in cwt			
			for (UInt i =0; i< candidates[c].size(); ++i)
			{
				avg_cwt   += candidates[c][i].getIntensity();			
			}						
			avg_cwt 	/= candidates[c].size();
			cwt_thresholds.at(c) =  avg_cwt * cwt_avg_factor_;		
			
			#ifdef DEBUG_FEATUREFINDER
			//write debug output
			CoordinateType current_rt = scan.getRT();
			String filename = String("isowavcwt_") + current_rt + "_charge_" + (c+1);
			std::ofstream outfile(filename.c_str());
			for (DPeakArray<PeakType >::const_iterator write_iter=candidates[c].begin(); 
						write_iter != candidates[c].end(); 
						++write_iter)
			{
				outfile << write_iter->getMZ() << " " << write_iter->getIntensity() << std::endl;
			}
			outfile.close();
			#endif
			
		}
	
		IntensityType avg_scan          = 0;		// average intensity in signal		
		IntensityType scan_threshold = 0;		// threshold for signal intensity
		
		for (UInt z=0; z<scan.size();++z)
		{
			avg_scan += scan[z].getIntensity();
		}
		avg_scan /= scan.size();
		scan_threshold = avg_scan * signal_avg_factor_;
		
		//std::cout << "Average intensity in scan: " << avg_scan << std::endl;
		//std::cout << "Intensity threshold for signal: " << scan_threshold << std::endl;
				
		for (UInt i = 1; i < candidates[0].size(); ++i) 			// cwt's for all charge states have the same length....
		{																						
								
				if (scan[i].getIntensity() < scan_threshold)	continue;	// ignore low intensity signals
		
				// vector of p-values
				std::vector<ProbabilityType> charge_scores( candidates.size(),numeric_limits<ProbabilityType>::max() );

				for (UInt c=0; c<candidates.size(); ++c)		// for all charge states
				{											
					// test if :
					// intensity in cwt is higher than threshold 
					// we are at a local max
					// the last pattern is already behind us
					if ( candidates[c][i].getIntensity() > cwt_thresholds[c] && 
							 i > last_pattern[c] && 
							(i - last_pattern[c]) > peak_cut_off_/ (c+1)  &&                                   
							(scan[i-1].getIntensity() - scan[i].getIntensity() < 0.0) && 
						  (scan[i+1].getIntensity() - scan[i].getIntensity() < 0.0) )						  	  
					{						
							UInt max = findNextMax_(candidates[c],i);
							ProbabilityType pvalue = testLocalVariance_(candidates[c],max,c);	
							charge_scores.at(c) = pvalue;
							last_pattern.at(c)     = max; 																	// store index of last pattern
					} // end if (local max...)
					
				}  // end for all (charge states)
				
				// determine highest scoring charge state
				ProbabilityType best_score = numeric_limits<ProbabilityType>::max();
				UInt best_charge                = 0;
				
				for (UInt z=0;z<charge_scores.size();++z)
				{
					if (charge_scores.at(z) < best_score)	
					{
						best_score  = charge_scores.at(z);
						best_charge = (z+1);
					}				
				}
				
				if (best_score == numeric_limits<ProbabilityType>::max()) continue;
					
				if (best_score == 0.0) 
				{
					best_score += 0.00000001; // add pseudo count for zero p-values
				}
				
				ScoredChargeType sc_charge;
				sc_charge.first 			= best_charge;		// charge
				sc_charge.second = best_score;			// score
					
				scmzvec.push_back( make_pair( last_pattern.at(best_charge-1) ,sc_charge) );	// store scored m/z position
					
			}	// end for (UInt c = 0; c < candidates[0].size(); ++c)

		// done
		return scmzvec;
	}

	double IsotopeWaveletSeeder::getAbsMean_(const DPeakArray<PeakType >& signal, UInt startIndex, UInt endIndex) const
	{
	  double res=0;
	  for (UInt i=startIndex; i<endIndex; ++i)
	  {
	  	res += fabs(signal[i].getIntensity());
		}
	  return (res/(double)(endIndex-startIndex+1));
	}
	
	UInt IsotopeWaveletSeeder::findNextMax_(const DPeakArray<PeakType >& cwt, const UInt index)
	{
	
		UInt max_index = index;
		IntensityType max_intensity = cwt[index].getIntensity();
		
		CoordinateType first_mass =  cwt[index].getMZ();
		CoordinateType mass_diff  = 0.0;
		
		// check to the left
		UInt i = index;
		while (mass_diff < 3.0 && i >= 1)
		{
			if (cwt[i].getIntensity() > 	max_intensity)
			{
				max_intensity = cwt[i].getIntensity();
				max_index     = i;
			}
			mass_diff = (first_mass - cwt[i].getMZ());
			--i;						
		}
		
		// check to the right
		i = index;
		while (mass_diff < 3.0 && i<cwt.size())
		{
			if (cwt[i].getIntensity() > 	max_intensity)
			{
				max_intensity = cwt[i].getIntensity();
				max_index     = i;
			}
			mass_diff = (cwt[i].getMZ() - first_mass);
			++i;						
		}		
		
		return max_index;
	}

	IsotopeWaveletSeeder::ProbabilityType IsotopeWaveletSeeder::testLocalVariance_(const DPeakArray<PeakType >& cwt, const UInt& start, const UInt charge_index)
	{			
		IntensityType cwt_sum    = 0.0;
		IntensityType cwt_sqsum = 0.0;
	
		CoordinateType first_mass =  cwt[start].getMZ();
		CoordinateType mass_diff  = 0.0;
		UInt i = start;
		for (; mass_diff < 5.0 && i < cwt.size() ;++i)
		{
			cwt_sum    += cwt[i].getIntensity();
			cwt_sqsum += (cwt[i].getIntensity() * cwt[i].getIntensity());		
			mass_diff    = cwt[i].getMZ() - first_mass;
		} 			
		UInt N_local = i-start;		
		
		first_mass =  cwt[start].getMZ();
		mass_diff  = 0.0;
		i = start;
		for (; mass_diff < 5.0 && i >= 1 ;--i)
		{
			cwt_sum    += cwt[i].getIntensity();
			cwt_sqsum += (cwt[i].getIntensity() * cwt[i].getIntensity());		
			mass_diff    = cwt[i].getMZ() - first_mass;
		} 			
		N_local += start-i;		
		
		IntensityType local_var = ( N_local * cwt_sqsum - ( cwt_sum * cwt_sum) ) / ( N_local * (N_local-1) );
		
		IntensityType f_stat  = local_var/null_var_[charge_index];
		ProbabilityType pval = (1 - gsl_cdf_fdist_P(f_stat, (N_local-1) , (n_null_[ charge_index ]-1) )) ;
		
		return pval; 			
	}

} // end of namespace OpenMS
