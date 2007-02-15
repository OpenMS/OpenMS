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
// $Maintainer: Alexandra Zerck $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>

using namespace std;
using namespace OpenMS::TwoDOptimizationFunctions;

namespace OpenMS
{
  namespace TwoDOptimizationFunctions
  {
   
    
			
    typedef DRawDataPoint<1> RawDataPointType;
    unsigned int total_nr_peaks;
    // a vector holding iterators of the raw data of the regions involved in the optimization
    std::vector<std::pair<int,int> > signal2D; 
    //iterator to the beginning of the raw data
    MSExperiment<RawDataPointType>::ConstIterator raw_data_first;
    // the picked peaks
    MSExperiment<DPickedPeak<1> >::Iterator picked_peaks_iter;
    // a vector storing the information which peaks are matching in the different scans involved in the
    // 2d optimization
    std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > > matching_peaks;
    std::multimap<double,IsotopeCluster>::iterator iso_map_iter;

    // Evaluation of the target function for nonlinear optimization.
    int residual2D(const gsl_vector* x, void* params , gsl_vector* f)
    {
      // According to the gsl conventions, x contains the parameters to be optimized.
      // In our case these are the weighted average mz-position and left and right width for all
      // matching peaks in all considered scans. Additionally the intensities of all peaks are stored.
      // Params might contain any additional parameters. We handle these using class members
      // instead.
      // The vector f is supposed to contain the result when we return from this function.
      // Note: GSL wants the values for each data point i as one component of the results vector
			double computed_signal, current_position, experimental_signal;
			double p_height, p_position, p_width;
			int count =0;
			int counter_posf=0;
			unsigned int num_scans = TwoDOptimizationFunctions::signal2D.size()/2;
			std::set<std::pair<UnsignedInt,UnsignedInt> >::iterator peak_iter = iso_map_iter->second.peaks_.begin();
			gsl_vector_set_zero(f);
   
			//iterate over all scans
			for (size_t current_scan = 0; current_scan < num_scans; ++current_scan)
				{
					unsigned int curr_scan_idx = current_scan + iso_map_iter->second.peaks_.begin()->first;
					//iterate over all points of the signal
					for (int current_point = 0;
							 current_point +  TwoDOptimizationFunctions::signal2D[2*current_scan].second
								 <= TwoDOptimizationFunctions::signal2D[2*current_scan+1].second;
							 ++current_point)
						{
							computed_signal   = 0.;
							current_position  = ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)->getContainer()
																	 .begin() + TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getPos();
							experimental_signal = ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)->getContainer().begin()
																		 + TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getIntensity();
	
#ifdef DEBUG_2D
							std::cout << "experimental signal rt "<<(raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
								->getRetentionTime()
												<< "\tmz " << ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
																			 ->getContainer().begin()+ TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getPos()
												<< "\tint " << ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
																				->getContainer().begin()+ TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)
								->getIntensity()<<std::endl;
#endif

							size_t current_peak = 0;
							peak_iter = iso_map_iter->second.peaks_.begin();
							while(peak_iter != iso_map_iter->second.peaks_.end() && peak_iter->first != curr_scan_idx) ++peak_iter;
							//iterate over all peaks of the current scan
							while(peak_iter != iso_map_iter->second.peaks_.end() && peak_iter->first == curr_scan_idx)
								{
									int peak_idx = distance(iso_map_iter->second.peaks_.begin(),peak_iter);
									MSSpectrum<DPickedPeak<1> >::Iterator p_peak_iter =
										(picked_peaks_iter + peak_iter->first)->begin() + peak_iter->second;
									double mz_in_hash = p_peak_iter->getPos() * 5;
									std::map<int,std::vector<MSSpectrum<DPickedPeak<1> >::Iterator> >::iterator  m_spec_iter =
										matching_peaks.begin();
									int map_idx=0;
									while(m_spec_iter->first != (int)(mz_in_hash+0.5) )
										{
											++map_idx;
											++m_spec_iter;
										}
									std::vector<MSSpectrum<DPickedPeak<1> >::Iterator>::iterator m_peak_iter =
										m_spec_iter->second.begin();

									while(*m_peak_iter != p_peak_iter && m_peak_iter !=m_spec_iter->second.end())
										{
											++m_peak_iter;
										}
									// if the current peak is in the reference scan take all parameters from the vector x
#ifdef DEBUG_2D
									std::cout << "ref_scan : "<< current_peak << "\t "
														<<gsl_vector_get(x,3*current_peak) << "\t" <<gsl_vector_get(x,3*current_peak+1)
														<< "\t"<<gsl_vector_get(x,3*current_peak+2)<<std::endl;
#endif
									//Store the current parameters for this peak
									p_position    = gsl_vector_get(x,total_nr_peaks+3*map_idx);
									p_height 	    = gsl_vector_get(x,peak_idx);
									p_width  	    = (current_position <= p_position) ?
										gsl_vector_get(x,total_nr_peaks+3*map_idx+1) :
										gsl_vector_get(x,total_nr_peaks+3*map_idx+2);
									++count;
		  


									//is it a Lorentz or a Sech - Peak?
									if (((picked_peaks_iter + peak_iter->first)->begin()+peak_iter->second)->getPeakShape()
											== PeakShapeType::LORENTZ_PEAK)
										{
#ifdef DEBUG_2D
											std::cout << "p_height "<< p_height << "\tp_position "<< p_position << "\tcurrent_position "
																<< current_position << std::endl;
#endif
											computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
										}
									// if it's a Sech - Peak
									else 
										{
#ifdef DEBUG_2D
											std::cout << "p_height "<< p_height << "\tp_position "<< p_position << "\tcurrent_position "
																<< current_position << std::endl;
#endif
											computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
										}
									++current_peak;
									++peak_iter;
		  
								}// end while
#ifdef DEBUG_2D
							std::cout << "computed vs experimental signal: "<< computed_signal << "\t"
												<< experimental_signal<<std::endl;
#endif
							gsl_vector_set(f, counter_posf,fabs(computed_signal - experimental_signal));
							++counter_posf;
						}
	  
				}
      
			// penalties : especially negative heights have to be penalised
			double penalty = 0.;
			OptimizationFunctions::PenaltyFactorsInt* penalties = (OptimizationFunctions::PenaltyFactorsInt *)params;
      
      
			//iterate over all peaks again to compute the penalties
			// first look at all positions and width parameters
			unsigned int peak=0,current_peak=0;
			std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > >::iterator map_iter=matching_peaks.begin();
			for (;map_iter != matching_peaks.end(); ++map_iter)
				{
					std::vector<MSSpectrum<DPickedPeak<1> >::Iterator >::iterator vec_iter
						= map_iter->second.begin();
					double old_position = 0,old_width_l=0,old_width_r=0;
					double weight =0;
					double old_height,p_height;
					for(;vec_iter != map_iter->second.end();++vec_iter)
						{
							old_height = (*vec_iter)->getIntensity();
							weight += old_height;
							old_position += (*vec_iter)->getPos() * old_height;
							old_width_l += (*vec_iter)->getLeftWidthParameter() * old_height;
							old_width_r += (*vec_iter)->getRightWidthParameter() * old_height;

							p_height     = gsl_vector_get(x, peak);
							++peak;
	      
							if(p_height < 1)
								{
									penalty += 1000000*penalties->height*pow(fabs(p_height - old_height),2);
								}
	      
						}
					old_position /= weight;
					old_width_l /= weight;
					old_width_r /= weight;
 	  
					double p_position   = gsl_vector_get(x, total_nr_peaks+3*current_peak);
					double p_width_l    = gsl_vector_get(x, total_nr_peaks+3*current_peak +1);
					double p_width_r    = gsl_vector_get(x, total_nr_peaks+3*current_peak +2);
					if(p_width_l < 0 )
						{
							penalty += 1e7*penalties->lWidth*pow(fabs(p_width_l - old_width_l),2);
						}
					else if(p_width_l < 1) penalty += 1000*penalties->lWidth*pow(fabs(p_width_l - old_width_l),2);
					if(p_width_r < 0 )
						{
							penalty += 1e7*penalties->rWidth*pow(fabs(p_width_r - old_width_r),2);
						}
					else if(p_width_r < 1) penalty += 1000*penalties->rWidth*pow(fabs(p_width_r - old_width_r),2);
					if(p_position < 0) 
						{
							penalty +=100*penalties->pos * pow(p_position - old_position, 2);
						}
					++current_peak;
				}
       
			gsl_vector_set(f, f->size -1, penalty);
			return GSL_SUCCESS; 
    }


    /** Compute the Jacobian of the residual, where each row of the matrix corresponds to a 
     *  point in the data.
     */
    int jacobian2D(const gsl_vector* x, void* params, gsl_matrix* J)
    {
			// For the conventions on x and params c.f. the commentary in residual()
			//
			// The matrix J is supposed to contain the result when we return from this function.
			// Note: GSL expects the Jacobian as follows:
			// 					- each row corresponds to one data point
			// 					- each column corresponds to one parameter

			double computed_signal, current_position, experimental_signal;
			double p_height, p_position, p_width;
			double diff, denom_inv,ddl_left,ddl_right,ddx0,sinh_term;
			int count =0;
			int counter_posf=0;
			unsigned int num_scans = TwoDOptimizationFunctions::signal2D.size()/2;
			std::vector<double> ov_weight(matching_peaks.size(),0);
			std::set<std::pair<UnsignedInt,UnsignedInt> >::iterator peak_iter = iso_map_iter->second.peaks_.begin();
   
			//iterate over all scans
			for (size_t current_scan = 0; current_scan < num_scans; ++current_scan)
				{
					unsigned int curr_scan_idx = current_scan + iso_map_iter->second.peaks_.begin()->first;
					// iterate over all points of the signal
					for (int current_point = 0;
							 current_point +  TwoDOptimizationFunctions::signal2D[2*current_scan].second
								 <= TwoDOptimizationFunctions::signal2D[2*current_scan+1].second;
							 ++current_point)
						{
							computed_signal   = 0.;
							current_position  = ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
																	 ->getContainer().begin() + TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getPos();
							experimental_signal = ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)->getContainer().begin()
																		 + TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getIntensity();
	
#ifdef DEBUG_2D
							std::cout << "experimental signal rt "<<(raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
								->getRetentionTime()
												<< "\tmz " << ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
																			 ->getContainer().begin()+ TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)->getPos()
												<< "\tint " << ((raw_data_first + TwoDOptimizationFunctions::signal2D[2*current_scan].first)
																				->getContainer().begin()+ TwoDOptimizationFunctions::signal2D[2*current_scan].second+current_point)
								->getIntensity()<<std::endl;
#endif

							size_t current_peak = 0;
							peak_iter = iso_map_iter->second.peaks_.begin();
							while(peak_iter != iso_map_iter->second.peaks_.end() && peak_iter->first != curr_scan_idx) ++peak_iter;
							//iterate over all peaks of the current scan
							while(peak_iter != iso_map_iter->second.peaks_.end() && peak_iter->first == curr_scan_idx)
								{
									int peak_idx = distance(iso_map_iter->second.peaks_.begin(),peak_iter);
									MSSpectrum<DPickedPeak<1> >::Iterator p_peak_iter =
										(picked_peaks_iter + peak_iter->first)->begin() + peak_iter->second;
									double mz_in_hash = p_peak_iter->getPos() * 5;
									std::map<int,std::vector<MSSpectrum<DPickedPeak<1> >::Iterator> >::iterator  m_spec_iter =
										matching_peaks.begin();
									int map_idx=0;
									while(m_spec_iter->first != (int)(mz_in_hash+0.5) )
										{
											++map_idx;
											++m_spec_iter;
										}
									std::vector<MSSpectrum<DPickedPeak<1> >::Iterator>::iterator m_peak_iter =
										m_spec_iter->second.begin();
		  
									while(*m_peak_iter != p_peak_iter && m_peak_iter !=m_spec_iter->second.end())
										{
											++m_peak_iter;
										}
									// if the current peak is in the reference scan take all parameters from the vector x
#ifdef DEBUG_2D
									std::cout << "ref_scan : "<< current_peak << "\t "
														<<gsl_vector_get(x,3*current_peak) << "\t" <<gsl_vector_get(x,3*current_peak+1)
														<< "\t"<<gsl_vector_get(x,3*current_peak+2)<<std::endl;
#endif
									// Store the current parameters for this peak
									p_position    = gsl_vector_get(x,total_nr_peaks+3*map_idx);
									p_height 	    = gsl_vector_get(x,peak_idx);
									p_width  	    = (current_position <= p_position) ?
										gsl_vector_get(x,total_nr_peaks+3*map_idx+1) :
										gsl_vector_get(x,total_nr_peaks+3*map_idx+2);
									++count;
									double weight = ((picked_peaks_iter + peak_iter->first)->
																	 begin()+peak_iter->second)->getIntensity();
									ov_weight[map_idx] += weight;
									double ddx0_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx);
									double ddl_left_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx+1);
									double ddl_right_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx+2);
									//is it a Lorentz or a Sech - Peak?
									if (((picked_peaks_iter + peak_iter->first)->begin()+peak_iter->second)->getPeakShape()
											== PeakShapeType::LORENTZ_PEAK)
										{
											diff      = current_position - p_position;
											// partial derivative with respect to the height,...
											denom_inv = 1./(1. + pow(p_width * diff, 2));
											// left width,...
											ddl_left  = (current_position <= p_position)
												? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2)
												: 0;
											// right width ...
											ddl_right = (current_position  > p_position)
												? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2)
												: 0;
		      
											// and position
											ddx0 = 2*p_height*pow(p_width,2)*diff*pow(denom_inv, 2);
		      
		      
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx    , ddx0* weight  + ddx0_old);
											gsl_matrix_set(J, counter_posf, peak_idx                     , denom_inv);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +1 , ddl_left * weight + ddl_left_old);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +2 , ddl_right * weight + ddl_right_old);
										}
									// if it's a Sech - Peak
									else
										{
											diff      = current_position - p_position;
											denom_inv = 1./cosh(p_width*diff);
											// The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
											// and can assume that all derivatives vanish
											sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width*diff);
		      
		      
											ddl_left  = (current_position <= p_position)
												? -2 * p_height * sinh_term*diff * pow(denom_inv,3) : 0;
											ddl_right = (current_position  > p_position)
												? -2 * p_height * sinh_term*diff * pow(denom_inv,3) : 0;
		      
											ddx0      = 2*p_height*p_width*sinh_term * pow(denom_inv,3);
		      
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx     , ddx0 * weight + ddx0_old);
											gsl_matrix_set(J, counter_posf, peak_idx                      , pow(denom_inv, 2));
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 1, ddl_left * weight + ddl_left_old);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 2, ddl_right * weight + ddl_right_old);

										}
									++current_peak;
									++peak_iter;
		  
								}// end while
	     
							++counter_posf;
						}
	  
				}

			for(unsigned int cluster=0; cluster < matching_peaks.size();++cluster)
				{
					for(unsigned int j=0; j < J->size1-1;++j)
						{
							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster  )/ov_weight[cluster]);
							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster + 1,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster+1)/ov_weight[cluster]);

							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster + 2,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster+2)/ov_weight[cluster]);
						}
				}
			// penalties : especially negative heights have to be penalised
			OptimizationFunctions::PenaltyFactorsInt* penalties = (OptimizationFunctions::PenaltyFactorsInt *)params;
      
      
			//iterate over all peaks again to compute the penalties
			// first look at all positions and width parameters
			unsigned int peak=0,current_peak=0;
			std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > >::iterator map_iter=matching_peaks.begin();
       
			for (;map_iter != matching_peaks.end(); ++map_iter)
				{
					std::vector<MSSpectrum<DPickedPeak<1> >::Iterator >::iterator vec_iter
						= map_iter->second.begin();
					double old_position = 0,old_width_l=0,old_width_r=0;
					double weight =0;
					double old_height,p_height;
					double penalty_h=0, penalty_l=0, penalty_r=0,penalty_p=0; 
					for(;vec_iter != map_iter->second.end();++vec_iter)
						{
							old_height = (*vec_iter)->getIntensity();
							weight += old_height;
							old_position += (*vec_iter)->getPos() * old_height;
							old_width_l += (*vec_iter)->getLeftWidthParameter() * old_height;
							old_width_r += (*vec_iter)->getRightWidthParameter() * old_height;

							p_height     = gsl_vector_get(x, peak);
	      
	      
							double penalty_height = 2.*penalties->height*fabs(p_height-old_height);
							if(p_height < 1)
								{
									penalty_h += 1000000*penalty_height;
								}
							gsl_matrix_set(J,counter_posf,peak,penalty_h);
							++peak;
						}
					old_position /= weight;
					old_width_l /= weight;
					old_width_r /= weight;


					double p_position   = gsl_vector_get(x, total_nr_peaks+3*current_peak);
					double p_width_l    = gsl_vector_get(x, total_nr_peaks+3*current_peak +1);
					double p_width_r    = gsl_vector_get(x, total_nr_peaks+3*current_peak +2);
					double penalty_lwidth = 2.*penalties->lWidth*fabs(p_width_l - old_width_l); 
					double penalty_rwidth = 2.*penalties->rWidth*fabs(p_width_r - old_width_r);
					double penalty_pos    = 2.*penalties->pos*fabs(p_position-old_position);
#ifdef DEBUG_2D
					std::cout << "penalty_lwidth " << penalty_lwidth << "penalty_rwidth " << penalty_rwidth
										<< "penalty_pos " << penalty_pos << std::endl;
#endif
					if(p_width_l < 0 )
						{
							penalty_l += 1e7*penalty_lwidth;
						}
					else if(p_width_l < 1) penalty_l += 2000*penalties->lWidth*(fabs(p_width_l - old_width_l));
					if(p_width_r < 0 )
						{
							penalty_r += 1e7*penalty_rwidth;
						}
					else if(p_width_r < 1) penalty_r += 2000*penalties->rWidth*(fabs(p_width_r - old_width_r));
					if(p_position < 0)
						{
							penalty_p +=200*penalty_pos;
						}


       	  gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak+1, penalty_l);
					gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak+2, penalty_r);
       	  gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak  , penalty_p);
       	

	  
					++current_peak;
				}
			return GSL_SUCCESS;       
    }

    //Driver function for the evaluation of function and jacobian.
    int evaluate2D(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      residual2D(x, params, f);
      jacobian2D(x, params, J);

      return GSL_SUCCESS;
    }

  }// namespace


	void TwoDOptimization::findMatchingPeaks_(std::multimap<double, IsotopeCluster>::iterator& it,
																						MSExperiment< DPickedPeak<1> >& ms_exp)
	{
		IndexSet::const_iterator iter = it->second.peaks_.begin();
		for(; iter != it->second.peaks_.end(); ++iter)
			{
				
				double mz = (ms_exp[iter->first][iter->second]).getPos();
				mz *= 5;
				matching_peaks_[(int)(mz+0.5)].push_back(ms_exp[iter->first].begin()+iter->second);
			}

		std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > >::iterator it2 = matching_peaks_.begin();
#ifdef DEBUG_2D
		for(;it2 != matching_peaks_.end();++it2)
			{
				std::cout << it2->first << " has "<<it2->second.size()<<" elements:"<<std::endl;
				for(unsigned int i=0;i<it2->second.size();++i) std::cout << it2->second[i]->getPos()<<"\t";
				std::cout<<std::endl;
			}
#endif
	}

	template < >
	void TwoDOptimization::optimizeRegions_<MSExperiment<DRawDataPoint<1> >::const_iterator, DPickedPeak<1> >
	(MSExperiment<DRawDataPoint<1> >::const_iterator& first,
	 MSExperiment<DRawDataPoint<1> >::const_iterator& last,
	 MSExperiment<DPickedPeak<1> >& ms_exp)
	{
      
		//#ifdef DEBUG_2D
		int counter =0;
		//#endif
		// go through the clusters
		for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
				 it != iso_map_.end();
				 ++it)
			{
#ifdef DEBUG_2D
				std::cout << "element: " << counter<< std::endl;
				std::cout << "mz: "<< it->first <<"\tcharge: " << it->second.charge_ << std::endl<<"rts: ";
				for(unsigned int i=0;i<it->second.scans_.size();++i) std::cout << it->second.scans_[i] << "\n";
				std::cout<<std::endl<<"peaks: ";
				for(unsigned int i=0;i<it->second.peaks_.size();++i) std::cout << it->second.peaks_[i].first
																																			 << "\t" <<  it->second.peaks_[i].second<<std::endl;
				std::cout << std::endl << std::endl;
	
#endif

				// prepare for optimization:
				// determine the matching peaks
				matching_peaks_.clear();
				findMatchingPeaks_(it,ms_exp);
				matching_peaks = matching_peaks_;
				// and the endpoints of each isotope pattern in the cluster
				getRegionEndpoints_(ms_exp,first,last,counter,400);
	  
				// peaks have to be stored globally
				iso_map_iter = it; 
	 
				picked_peaks_iter = ms_exp.begin();
				raw_data_first =  first;
	 
				int nr_diff_peaks = matching_peaks_.size();
				total_nr_peaks = it->second.peaks_.size();
	  
				int nr_parameters = nr_diff_peaks*3 + total_nr_peaks;
	  
				gsl_vector *start_value=gsl_vector_alloc(nr_parameters);
				gsl_vector_set_zero(start_value);

				// initialize parameters for optimization
				std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > >::iterator m_peaks_it
					= matching_peaks.begin();
				double av_mz=0,av_lw=0,av_rw=0,avr_height=0,height;
				int peak_counter = 0;
				int diff_peak_counter =0;
				// go through the matching peaks
				for(;m_peaks_it != matching_peaks.end();++m_peaks_it)
					{
						av_mz=0,av_lw=0,av_rw=0,avr_height=0;
						std::vector<MSSpectrum<DPickedPeak<1> >::Iterator>::iterator iter_iter = (m_peaks_it)->second.begin();
						for(;iter_iter != m_peaks_it->second.end();++iter_iter)
							{
								height = (*iter_iter)->getIntensity();
								avr_height += height;
								av_mz += (*iter_iter)->getPos() * height;
								av_lw += (*iter_iter)->getLeftWidthParameter() * height;
								av_rw += (*iter_iter)->getRightWidthParameter() * height;
								gsl_vector_set(start_value,peak_counter,height);
								++peak_counter;
							}
						gsl_vector_set(start_value,total_nr_peaks+3*diff_peak_counter, av_mz/avr_height);
						gsl_vector_set(start_value,total_nr_peaks+3*diff_peak_counter+1, av_lw/avr_height);
						gsl_vector_set(start_value,total_nr_peaks+3*diff_peak_counter+2, av_rw/avr_height);
						++diff_peak_counter;
					}

#ifdef DEBUG_2D
				std::cout << "----------------------------\n\nstart_value: "<<std::endl;
				for(unsigned int k=0;k<start_value->size;++k)
					{
						std::cout << gsl_vector_get(start_value,k)<<std::endl;
					}
#endif
	  
				int num_positions = 0;
				for(unsigned int i=0; i<TwoDOptimizationFunctions::signal2D.size(); i+=2)
					{
						num_positions += (TwoDOptimizationFunctions::signal2D[i+1].second - TwoDOptimizationFunctions::signal2D[i].second +1);
#ifdef DEBUG_2D	      
						std::cout << TwoDOptimizationFunctions::signal2D[i+1].second << " - "<< TwoDOptimizationFunctions::signal2D[i].second <<" +1 "<< std::endl;
#endif
					}
#ifdef DEBUG_2D	  
				std::cout << "num_positions : "<< num_positions << std::endl;
#endif	  
				// The gsl algorithms require us to provide function pointers for the evaluation of
				// the target function.
				gsl_multifit_function_fdf fit_function;
				fit_function.f   = TwoDOptimizationFunctions::residual2D;
				fit_function.df  = TwoDOptimizationFunctions::jacobian2D;
				fit_function.fdf = TwoDOptimizationFunctions::evaluate2D;
				// dunno why, but gsl crashes when n is smaller than p!!!!!!!! ????
				fit_function.n   = std::max(num_positions+1,(int)(nr_parameters));
				fit_function.p   = nr_parameters;
				fit_function.params = &penalties_;
#ifdef DEBUG_2D
				std::cout << "fit_function.n "<<fit_function.n
									<< "\tfit_function.p "<<fit_function.p<<std::endl;
#endif
				const gsl_multifit_fdfsolver_type *type = gsl_multifit_fdfsolver_lmsder;
	  
				gsl_multifit_fdfsolver *fit=gsl_multifit_fdfsolver_alloc(type,
																																 std::max(num_positions+1,(int)(nr_parameters)),
																																 nr_parameters);
	  
				gsl_multifit_fdfsolver_set(fit, &fit_function, start_value);

	 
	  
				// initial norm
#ifdef DEBUG_2D
				std::cout << "Before optimization: ||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
#endif
				// Iteration
				int iteration = 0;
				int status;
      
				do 
					{
						iteration++;
						status = gsl_multifit_fdfsolver_iterate(fit);
#ifdef DEBUG_2D
						std::cout << "Iteration " << iteration << "; Status " << gsl_strerror(status) << "; " << std::endl;
						std::cout << "||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
						std::cout << "Number of parms: " << nr_parameters << std::endl;
						std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
#endif
						if (isnan(gsl_blas_dnrm2(fit->dx)))
							break;
	      
						status = gsl_multifit_test_delta(fit->dx, fit->x, eps_abs_, eps_rel_);
						if (status != GSL_CONTINUE)
							break;
	      
					}
				while (status == GSL_CONTINUE && iteration < max_iteration_);
	
#ifdef DEBUG_2D		
				std::cout << "Finished! No. of iterations" << iteration << std::endl;
				std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
				double chi = gsl_blas_dnrm2(fit->f);
				std::cout << "After optimization: || f || = " << gsl_blas_dnrm2(fit->f) << std::endl;
				std::cout << "chisq/dof = " << pow(chi, 2.0) / (num_positions - nr_parameters);	


				std::cout << "----------------------------------------------\n\nnachher"<<std::endl;
				for(unsigned int k=0;k<fit->x->size;++k)
					{
						std::cout << gsl_vector_get(fit->x,k)<<std::endl;
					}
#endif
				int peak_idx =0;
				std::map<int, std::vector<MSSpectrum<DPickedPeak<1> >::Iterator > >::iterator it=matching_peaks.begin();
				for(; it != matching_peaks.end(); ++it)
					{
						int i = distance(matching_peaks.begin(),it);
						for(unsigned int j=0;j<it->second.size();++j)
							{

#ifdef DEBUG_2D
								std::cout << "pos: "<<it->second[j]->getPos()<<"\nint: "<<it->second[j]->getIntensity()
													<<"\nlw: "<<it->second[j]->getLeftWidthParameter()
													<<"\nrw: "<<it->second[j]->getRightWidthParameter() << "\n";
		  
#endif
		  
								it->second[j]->setPos(gsl_vector_get(fit->x,total_nr_peaks+3*i));

								it->second[j]->setIntensity(gsl_vector_get(fit->x,peak_idx));
		  
								it->second[j]->setLeftWidthParameter(gsl_vector_get(fit->x,total_nr_peaks+3*i+1));
		  
								it->second[j]->setRightWidthParameter(gsl_vector_get(fit->x,total_nr_peaks+3*i+2));

#ifdef DEBUG_2D
								std::cout << "pos: "<<it->second[j]->getPos()<<"\nint: "<<it->second[j]->getIntensity()
													<<"\nlw: "<<it->second[j]->getLeftWidthParameter()
													<<"\nrw: "<<it->second[j]->getRightWidthParameter() << "\n";
		  
#endif
		  
								++peak_idx;

		  
							}
					}

				gsl_multifit_fdfsolver_free(fit);
				gsl_vector_free(start_value);
				++counter; 
			}// end for
      
	}


	template < >
	void TwoDOptimization::optimizeRegionsScanwise_<MSExperiment<DRawDataPoint<1> >::const_iterator, DPickedPeak<1> >
	(MSExperiment<DRawDataPoint<1> >::const_iterator& first,
	 MSExperiment<DRawDataPoint<1> >::const_iterator& last,
	 MSExperiment<DPickedPeak<1> >& ms_exp)
	{

		int counter =0;
		picked_peaks_iter = ms_exp.begin();
		raw_data_first =  first;


		struct OpenMS::OptimizationFunctions::PenaltyFactors penalties;
	        
	      
		DataValue dv = param_.getValue("2D_optimization:penalties:position");
		if (dv.isEmpty() || dv.toString() == "")
			penalties.pos = 0.;
		else
			penalties.pos = (float)dv;
	      
		dv = param_.getValue("2D_optimization:penalties:left_width");
		if (dv.isEmpty() || dv.toString() == "")
			penalties.lWidth = 1.;
		else
			penalties.lWidth = (float)dv;
	      
		dv = param_.getValue("2D_optimization:penalties:right_width");
		if (dv.isEmpty() || dv.toString() == "")
			penalties.rWidth = 1.;
		else
			penalties.rWidth = (float)dv;
	      
		unsigned int max_iteration;
		dv = param_.getValue("2D_optimization:iterations");
		if (dv.isEmpty() || dv.toString() == "")
			max_iteration = 15;
		else
			max_iteration = (unsigned int)dv;

		double eps_abs;
		dv = param_.getValue("2D_optimization:delta_abs_error");
		if (dv.isEmpty() || dv.toString() == "")
			eps_abs = 1e-04f;
		else
			eps_abs = (double)dv;
	      
		double eps_rel;
		dv = param_.getValue("2D_optimization:delta_rel_error");
		if (dv.isEmpty() || dv.toString() == "")
			eps_rel = 1e-04f;
		else
			eps_rel = (double)dv;

		std::vector<PeakShape> peak_shapes;
       

		// go through the clusters
		for (std::multimap<double, IsotopeCluster>::iterator it = iso_map_.begin();
				 it != iso_map_.end();
				 ++it)
			{
				iso_map_iter = it; 
#ifdef DEBUG_2D
				std::cout << "element: " << counter<< std::endl;
				std::cout << "mz: "<< it->first <<"\tcharge: " << it->second.charge_ << std::endl<<"rts: ";
				for(unsigned int i=0;i<it->second.scans_.size();++i) std::cout << it->second.scans_[i] << "\n";
				std::cout<<std::endl<<"peaks: ";
				for(unsigned int i=0;i<it->second.peaks_.size();++i) std::cout << it->second.peaks_[i].first
																																			 << "\t" <<  it->second.peaks_[i].second<<std::endl;
				std::cout << std::endl << std::endl;
	
#endif
	  
	  
				// prepare for optimization:
				// determine the matching peaks
				// and the endpoints of each isotope pattern in the cluster
				getRegionEndpoints_(ms_exp,first,last,counter,400);
	  
				unsigned int idx = 0;
				for(unsigned int i=0; i < TwoDOptimizationFunctions::signal2D.size()/2; ++i)
					{
						OpenMS::OptimizationFunctions::positions_.clear();
						OpenMS::OptimizationFunctions::signal_.clear();
	      
						MSSpectrum<DRawDataPoint<1> >::const_iterator ms_it =
							(raw_data_first + TwoDOptimizationFunctions::signal2D[2*i].first)->begin()+TwoDOptimizationFunctions::signal2D[2*i].second;
						int size = distance(ms_it,(raw_data_first + TwoDOptimizationFunctions::signal2D[2*i].first)->begin()+TwoDOptimizationFunctions::signal2D[2*i+1].second);
						OpenMS::OptimizationFunctions::positions_.reserve(size);
						OpenMS::OptimizationFunctions::signal_.reserve(size);

						while(ms_it != (raw_data_first + TwoDOptimizationFunctions::signal2D[2*i].first)->begin()+TwoDOptimizationFunctions::signal2D[2*i+1].second)
							{
								OpenMS::OptimizationFunctions::positions_.push_back(ms_it->getPos());
								OpenMS::OptimizationFunctions::signal_.push_back(ms_it->getIntensity());
								++ms_it;
							}

						
						Idx pair;
						pair.first =  iso_map_iter->second.peaks_.begin()->first + idx;
						
						IndexSet::const_iterator set_iter = lower_bound(iso_map_iter->second.peaks_.begin(),
																														iso_map_iter->second.peaks_.end(),
																														pair,IndexLess());
						
						
						// find the last entry with this rt-value
						++pair.first;
						IndexSet::const_iterator set_iter2 = lower_bound(iso_map_iter->second.peaks_.begin(),
																														 iso_map_iter->second.peaks_.end(),
																														 pair,IndexLess());
						
						while(set_iter != set_iter2)
							{
								DPickedPeak<1> peak = *(ms_exp[set_iter->first].begin()+set_iter->second);
								PeakShape shape(peak.getIntensity(),
																peak.getPos(),
																peak.getLeftWidthParameter(),
																peak.getRightWidthParameter(),
																peak.getArea(),
																peak.getPeakShape());
								peak_shapes.push_back(shape);
								++set_iter;
								
							}
						std::cout	<< "rt "<<(raw_data_first + TwoDOptimizationFunctions::signal2D[2*i].first)->getRetentionTime() <<"\n";
						OptimizePick opt(penalties,max_iteration,eps_abs,eps_rel);
						std::cout	<< "vorher\n";
						for(unsigned int p =0; p < peak_shapes.size();++p)
							{
								std::cout	<< peak_shapes[p].mz_position <<"\t" << peak_shapes[p].height
													<<"\t" << peak_shapes[p].left_width <<"\t" << peak_shapes[p].right_width 	<<std::endl;
							}
						opt.optimize(peak_shapes);
						std::cout	<< "nachher\n";
						for(unsigned int p =0; p < peak_shapes.size();++p)
							{
								std::cout	<< peak_shapes[p].mz_position <<"\t" << peak_shapes[p].height
													<<"\t" << peak_shapes[p].left_width <<"\t" << peak_shapes[p].right_width 	<<std::endl;
							}
	     
						std::sort(peak_shapes.begin(),peak_shapes.end(),PeakShape::PositionLess());
						pair.first =  iso_map_iter->second.peaks_.begin()->first + idx;
						
						set_iter = lower_bound(iso_map_iter->second.peaks_.begin(),
																	 iso_map_iter->second.peaks_.end(),
																	 pair,IndexLess());
						unsigned int i=0;
						while(i < peak_shapes.size())
							{
								(ms_exp[set_iter->first][set_iter->second]) . setPos(peak_shapes[i].mz_position);

								(ms_exp[set_iter->first][set_iter->second]) . setIntensity(peak_shapes[i].height);

								(ms_exp[set_iter->first][set_iter->second]) . setLeftWidthParameter(peak_shapes[i].left_width);

								(ms_exp[set_iter->first][set_iter->second]) . setRightWidthParameter(peak_shapes[i].right_width);
							
								++set_iter;
								++i;
							}
						++idx;
						peak_shapes.clear();
					}

				++counter;
			}



      
	}

	/// Finds the neighbour of the peak denoted by @p current_mz in the previous scan
  std::vector<double>::iterator TwoDOptimization::searchInScan_(std::vector<double>::iterator scan_begin,
																																std::vector<double>::iterator scan_end ,
																																double current_mz)
  {

    // perform binary search to find the neighbour in rt dimension
    // 	lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
    std::vector<double>::iterator insert_iter = lower_bound(scan_begin,scan_end,current_mz);

    // the peak found by lower_bound does not have to be the closest one, therefore we have
    // to check both neighbours
    if ( insert_iter == scan_end ) // we are at the and have only one choice
      {
				return --insert_iter;
      }
    else
      {
				// if the found peak is at the beginning of the spectrum,
				// there is not much we can do.
				if ( insert_iter == scan_begin )
					{
						return insert_iter;
					}
				else // see if the next smaller one fits better
					{
						double delta_mz = fabs(*insert_iter - current_mz);
						--insert_iter;

						if ( fabs(*insert_iter - current_mz) < delta_mz )
							{
								return insert_iter; // peak to the left is closer (in m/z dimension)
							}
						else
							{
								return ++insert_iter;    // peak to the right is closer
							}
					}
      }
		
  } // end of searchInScan_

      
    
}
