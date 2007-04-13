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
			: BaseSweepSeeder(), 
				is_initialized_(false)
	{
		setName(getProductName());
	
		// lower and upper bounds for distances between isotopic peaks (defaults)
		// we check for charges states up to three
		// charge 1
		defaults_.setValue("charge1_ub",1.3f);
		defaults_.setValue("charge1_lb",0.70f);
		// charge 2
		defaults_.setValue("charge2_ub",0.70f);
		defaults_.setValue("charge2_lb",0.40f);
		// charge 3
		defaults_.setValue("charge3_ub",0.40f);
		defaults_.setValue("charge3_lb",0.1f);
			
		// params for the cwt
		defaults_.setValue("cwt_scale",0.1f);
		defaults_.setValue("noise_level_signal",1000);
		defaults_.setValue("noise_level_cwt",6000);
		
		// minimum number of maxima in cwt
		defaults_.setValue("min_peaks_per_scan",2);
		
		defaultsToParam_();
	}
	
	MarrWaveletSeeder::~MarrWaveletSeeder()
	{
	}

  MarrWaveletSeeder::MarrWaveletSeeder(const MarrWaveletSeeder& rhs)
    : BaseSweepSeeder(rhs),
    	is_initialized_(false)
  {
    updateMembers_();
  }
  
  MarrWaveletSeeder& MarrWaveletSeeder::operator= (const MarrWaveletSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSweepSeeder::operator=(rhs);
    is_initialized_ = false;
    
    updateMembers_();
    
    return *this;
  }
	
	void MarrWaveletSeeder::updateMembers_()
	{
		// update member of base class first
		BaseSweepSeeder::updateMembers_();
	
		// retrieve values for accepted peaks distances
		charge1_ub_	= param_.getValue("charge1_ub");
		charge1_lb_	 = param_.getValue("charge1_lb");
	
		charge2_ub_	= param_.getValue("charge2_ub");
		charge2_lb_	 = param_.getValue("charge2_lb");
	
		charge3_ub_	= param_.getValue("charge3_ub");
		charge3_lb_	 = param_.getValue("charge3_lb");
				
		// thresholds for signal and cwt
		noise_level_signal_ 							 = param_.getValue("noise_level_signal");
		noise_level_cwt_                   = param_.getValue("noise_level_cwt");
	
		// scale of Marr wavelet
		cwt_scale_                          = param_.getValue("cwt_scale");
		
		// minimum number of maxima per scan
		min_peaks_ = param_.getValue("min_peaks_per_scan");
	}
	
	MarrWaveletSeeder::ScoredMZVector MarrWaveletSeeder::detectIsotopicPattern_(SpectrumType& current_scan )
	{	
			UInt current_charge = 0;							// charge state of the current isotopic cluster
			
			double spacing_cwt = 0.0001; 		// spacing between sampling points of the wavelet
			double resolution_cwt = 1.0;			 // compute convolution for each data point
					
			ScoredMZVector scored_positions;
			
			cwt_.init(cwt_scale_, spacing_cwt);
			cwt_.transform(current_scan.begin(), current_scan.end(),resolution_cwt);
	
			#ifdef DEBUG_FEATUREFINDER
			String fname = String("cwt_") + current_scan.getRT();
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
				return scored_positions;
			}
	
			// search for maximal positions in the cwt and extract potential peaks
			vector<Int> local_maxima;
			getMaxPositions_(current_scan.begin(), current_scan.end(), cwt_, local_maxima 
																	#ifdef DEBUG_FEATUREFINDER
											 						, current_scan.getRT() 
																	#endif
											 						);
	
			UInt nr_maxima = local_maxima.size();
			cout << "# local maxima in cwt : " << nr_maxima << endl;
						
			// compute score and charge estimate for each position
			for (UInt i=0; i < ( nr_maxima - 1 ) && nr_maxima > 0; ++i)
			{
						
				// distance in m/z to next local max.
				CoordinateType curr_mz = current_scan[ local_maxima[i] ].getMZ();
				CoordinateType dist2nextpeak = (current_scan[ local_maxima[i+1] ].getMZ() - curr_mz);
	
				// test for different charge states
				current_charge = distanceToCharge_(dist2nextpeak);			
							
				if (current_charge > 0) 	// 0 => no pattern
				{
					ScoredChargeType sc_charge;
					sc_charge.first = current_charge;
			
					// count number of local maxima supporting this score
					UInt count = 1; 
					for (UInt c = (i+1); c < nr_maxima; ++c)
					{	
						CoordinateType this_mz = 	current_scan[ local_maxima[ c ] ].getMZ();
						CoordinateType prev_mz = current_scan[ local_maxima[ (c-1) ] ].getMZ();
									
						UInt next_charge = 	distanceToCharge_( (this_mz - prev_mz) );
											
						if (next_charge != current_charge) 	break;
					
						++count;	
					}
					
					if (count >= min_peaks_)  										
					{
						
						#ifdef DEBUG_FEATUREFINDER
						cout << "Pattern at " << current_scan[ local_maxima[i] ].getMZ() << endl;
						cout	<< "There are " << count << " supporting local maxima. " << endl;
						#endif
						
						sc_charge.second = testLocalVariance_(local_maxima,i);
						scored_positions.push_back( make_pair(local_maxima[i],sc_charge) );	
					}
				}
				
			}
			
			return scored_positions;
			
	} // end of void detectIsotopicPattern_(SpectrumType& )
	
	void MarrWaveletSeeder::getMaxPositions_( const SpectrumType::const_iterator&  first, 
																						const SpectrumType::const_iterator&  last, 
																						const ContinuousWaveletTransform& wt, 
																						vector<Int>& localmax
																						#ifdef DEBUG_FEATUREFINDER 
																							, CoordinateType current_rt
																						#endif
																						)
	{
		if (wt.getSize() == 0) return;
	
		Int zeros_left_index  = wt.getLeftPaddingIndex();
		Int zeros_right_index = wt.getRightPaddingIndex();
			
		// Points to most intensive data point in the signal
		SpectrumType::const_iterator it_max_pos;
			
		IntensityType max_value = 0.0;	// intensity of the signal peak
		UInt max_index   = 0;		// index of the signal peak (counted from the very first data point)

		Int start = zeros_left_index + 2;
		Int end  = zeros_right_index - 1;
		
		#ifdef DEBUG_FEATUREFINDER 
		// remove debug file (below we append entries only,
		// so if the file exists, we will get some quite confusing results)
		String fname = String("cwt_localmax_") + current_rt;
		if (File::exists(fname) )
		{
			File::remove(fname);
		}		
		#endif
	
		Int i=0, j=0;
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
				Int radius = 3;  // search radius for peaks in data
				Int start_intervall = (( i - (Int)radius) < 0 ) ? 0 : ( i - (Int)radius);
				Int end_intervall  = (( i + (Int)radius) >= distance(first,last)) ? 0 : (i + (Int)radius);
												
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
					localmax.push_back(max_index);
				}
			}
		}
	}
	
	MarrWaveletSeeder::ProbabilityType MarrWaveletSeeder::testLocalVariance_(const vector<Int>& local_maxima , const UInt max_index)
	{	
		IntensityType cwt_sum    = 0.0;
		IntensityType cwt_sqsum = 0.0;
		
		UInt N = 30; // number of samples: 10 to the left, 20 to the right 
					
		// compute local variance and test against variance in end or beginning of spectrum
		
		Int start = ( local_maxima[ max_index ] < 10 ) ? 0  : (local_maxima[ max_index ] - 10);
		Int end  = (local_maxima[ max_index ] + 20) < (Int) cwt_.getSize() ? (local_maxima[ max_index ] + 20) : cwt_.getSize();

		for (Int j = start ; j < end ; ++j)	// walk to the left in cwt
		{
			cwt_sum    += cwt_[j];
			cwt_sqsum += (cwt_[j] * cwt_[j]);							
		} 
					
		
		IntensityType local_var = ( N * cwt_sqsum - ( cwt_sum * cwt_sum) ) / ( N * (N-1) );
		
		// we need to check for several cases:
		// a) are we at the very beginning of the scan => compute variance at the end
		// b) end of scan => other way around
		// c) scan to short => do what ?
		
		if ( cwt_.getSize() < 60 )
		{
			// cwt too small, return very small p-value
			return 0.0001;	
		}
	
		IntensityType null_var = 0;
		cwt_sum  =  cwt_sqsum  = 0.0;

		// check if we are at the beginning or end of the cwt
		if (start < 30 )
		{
			// take end	
			for (Int j = (cwt_.getSize() - 30); j < cwt_.getSize(); ++j)	// walk to the left in cwt
			{
				cwt_sum    += cwt_[j];
				cwt_sqsum += (cwt_[j] * cwt_[j]);							
			} 
		}
		else
		{
			// take begin
			for (UInt j = 0; j < 30; ++j)	// walk to the left in cwt
			{
				cwt_sum    += cwt_[j];
				cwt_sqsum += (cwt_[j] * cwt_[j]);							
			} 
		}
		
		null_var = ( N * cwt_sqsum - ( cwt_sum * cwt_sum) ) / ( N * (N-1) );
		
		IntensityType f_stat = local_var / null_var;		
		
// 		cout << "local_var " << local_var << endl;
// 		cout << "null_var " << null_var << endl;
// 		cout << "Value of f_stat " << f_stat << endl;
// 		cout << "p-value is " << (1 - gsl_cdf_fdist_P(f_stat,29,29)) << endl;

		return f_stat; 
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
		return 0;
	}

} // end of namespace OpenMS
