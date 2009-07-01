// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/TwoDOptimization.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{

    // Evaluation of the target function for nonlinear optimization.
    Int TwoDOptimization::residual2D_(const gsl_vector* x, void* params , gsl_vector* f)
    {
      // According to the gsl conventions, x contains the parameters to be optimized.
      // In our case these are the weighted average mz-position and left and right width for all
      // matching peaks in all considered scans. Additionally the intensities of all peaks are stored.
      // Params might contain any additional parameters. We handle these using class members
      // instead.
      // The vector f is supposed to contain the result when we return from this function.
      // Note: GSL wants the values for each data point i as one component of the results vector
      DoubleReal computed_signal, current_position, experimental_signal,step,last_position;
      DoubleReal p_height, p_position, p_width;
      Int count =0;
      Int counter_posf=0;
			std::vector<std::pair<SignedSize,SignedSize> >& signal2D = static_cast<TwoDOptimization::Data*> (params) ->signal2D; 
			std::multimap<DoubleReal,IsotopeCluster>::iterator iso_map_iter=static_cast<TwoDOptimization::Data*> (params) ->iso_map_iter;
			Size total_nr_peaks=static_cast<TwoDOptimization::Data*> (params) ->total_nr_peaks;
			std::map<Int, std::vector<PeakIndex> >& matching_peaks=static_cast<TwoDOptimization::Data*> (params) ->matching_peaks;
			MSExperiment<> &picked_peaks = static_cast<TwoDOptimization::Data*> (params) ->picked_peaks;
			MSExperiment<Peak1D>::ConstIterator raw_data_first = static_cast<TwoDOptimization::Data*> (params) ->raw_data_first;
			OptimizationFunctions::PenaltyFactorsIntensity& penalties=static_cast<TwoDOptimization::Data*> (params) ->penalties;
// 			std::vector<DoubleReal> &positions=static_cast<TwoDOptimization::Data*> (params) ->positions;
// 			std::vector<DoubleReal> &signal=static_cast<TwoDOptimization::Data*> (params) ->signal;

      Size num_scans = signal2D.size()/2;
      IsotopeCluster::ChargedIndexSet::iterator peak_iter = iso_map_iter->second.peaks.begin();
      gsl_vector_set_zero(f);

      //iterate over all scans
      for (Size current_scan = 0; current_scan < num_scans; ++current_scan)
				{
					Size curr_scan_idx = current_scan + iso_map_iter->second.peaks.begin()->first;
					current_position = ((raw_data_first 
															 + signal2D[2*current_scan].first)->begin() 
															+ signal2D[2*current_scan].second)->getMZ();
					//iterate over all points of the signal
					for (Int current_point = 1;
							 current_point +  signal2D[2*current_scan].second
								 <= signal2D[2*current_scan+1].second;
							 ++current_point)
						{
							last_position = current_position;
							
							computed_signal   = 0.;
							current_position  = ((raw_data_first 
																		+ signal2D[2*current_scan].first)->begin() 
																	 + signal2D[2*current_scan].second+current_point)->getMZ();
							experimental_signal = ((raw_data_first 
																			+ signal2D[2*current_scan].first)->begin()
																		 + signal2D[2*current_scan].second+current_point)->getIntensity();
							step = current_position - last_position;
#ifdef DEBUG_2D
          std::cout << "experimental signal rt "<<(raw_data_first 
                + signal2D[2*current_scan].first)->getRT()
                << "\tmz " << ((raw_data_first 
                + signal2D[2*current_scan].first)->begin()
															 + signal2D[2*current_scan].second+current_point)->getMZ()
                << "\tint " << ((raw_data_first 
                + signal2D[2*current_scan].first)->begin()
                + signal2D[2*current_scan].second+current_point)->getIntensity()<<std::endl;
#endif

							size_t current_peak = 0;
							peak_iter = iso_map_iter->second.peaks.begin();
							while(peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first != curr_scan_idx) ++peak_iter;
							//iterate over all peaks of the current scan
							while(peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first == curr_scan_idx)
								{
									Int peak_idx = distance(iso_map_iter->second.peaks.begin(),peak_iter);
									DoubleReal mz_in_hash = ((picked_peaks[peak_iter->first]).begin() + peak_iter->second)->getMZ() * 10;
									std::map<Int,std::vector<PeakIndex> >::iterator  m_spec_iter = matching_peaks.begin();
									Int map_idx=0;
									while(m_spec_iter->first != (Int)(mz_in_hash+0.5) )
										{
											++map_idx;
											++m_spec_iter;
										}
									// if the current peak is in the reference scan take all parameters from the vector x
#ifdef DEBUG_2D
									std::cout << "ref_scan : "<< current_peak << "\t "
														<<gsl_vector_get(x,3*current_peak) << "\t" <<gsl_vector_get(x,3*current_peak+1)
														<< "\t"<<gsl_vector_get(x,3*current_peak+2)<<std::endl;
#endif
									//Store the current parameters for this peak
									p_position    = gsl_vector_get(x,total_nr_peaks+3*map_idx);
									p_height      = gsl_vector_get(x,peak_idx);
									p_width       = (current_position <= p_position) ?
										gsl_vector_get(x,total_nr_peaks+3*map_idx+1) :
										gsl_vector_get(x,total_nr_peaks+3*map_idx+2);
									++count;



									//is it a Lorentz or a Sech - Peak?
									if ((PeakShape::Type)(Int)Math::round((picked_peaks[peak_iter->first]).getFloatDataArrays()[5][peak_iter->second]) == PeakShape::LORENTZ_PEAK)
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
							gsl_vector_set(f, counter_posf,step*(computed_signal - experimental_signal));
							++counter_posf;
						}

				}

      // penalties : especially negative heights have to be penalised
      DoubleReal penalty = 0.;
      


      //iterate over all peaks again to compute the penalties
      // first look at all positions and width parameters
      UInt peak=0,current_peak=0;
      std::map<Int, std::vector<PeakIndex> >::iterator map_iter=matching_peaks.begin();
      for (;map_iter != matching_peaks.end(); ++map_iter)
				{
					std::vector<PeakIndex >::iterator vec_iter = map_iter->second.begin();
					DoubleReal old_position = 0,old_width_l=0,old_width_r=0;
					DoubleReal weight =0;
					DoubleReal old_height,p_height;
					for(;vec_iter != map_iter->second.end();++vec_iter)
						{
						
							old_height = (vec_iter)->getPeak(picked_peaks).getIntensity();
							weight += old_height;
							old_position += (vec_iter)->getPeak(picked_peaks).getMZ() * old_height;
							old_width_l += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[3][vec_iter->peak]* old_height;
							old_width_r += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[4][vec_iter->peak]* old_height;
	
							p_height     = gsl_vector_get(x, peak);
							++peak;

							if(p_height < 1)
								{
									penalty += 1000000*penalties.height*pow(fabs(p_height - old_height),2);
								}

						}
					old_position /= weight;
					old_width_l /= weight;
					old_width_r /= weight;

					DoubleReal p_position   = gsl_vector_get(x, total_nr_peaks+3*current_peak);
					DoubleReal p_width_l    = gsl_vector_get(x, total_nr_peaks+3*current_peak +1);
					DoubleReal p_width_r    = gsl_vector_get(x, total_nr_peaks+3*current_peak +2);
					if(p_width_l < 0 )
						{
							penalty += 1e7*penalties.lWidth*pow(fabs(p_width_l - old_width_l),2);
						}
					else if(p_width_l < 1) penalty += 1000*penalties.lWidth*pow(fabs(p_width_l - old_width_l),2);
					if(p_width_r < 0 )
						{
							penalty += 1e7*penalties.rWidth*pow(fabs(p_width_r - old_width_r),2);
						}
					else if(p_width_r < 1) penalty += 1000*penalties.rWidth*pow(fabs(p_width_r - old_width_r),2);
					if(p_position < 0)
						{
							penalty +=100*penalties.pos * pow(p_position - old_position, 2);
						}
					if(fabs(old_width_r-p_width_r) > 1)
						{
							penalty +=1000*penalties.rWidth * pow(old_width_r-p_width_r, 2);
						}
					if(fabs(old_width_l-p_width_l) > 1)
						{
							penalty +=1000*penalties.lWidth * pow(old_width_l-p_width_l, 2);
						}
					if(fabs(old_position-p_position) > 0.2)
						{
							penalty +=1000*penalties.pos * pow(p_position - old_position, 2);
						}
				
					++current_peak;
				}

      gsl_vector_set(f, f->size -1, penalty);
      return GSL_SUCCESS;
    }


    /** Compute the Jacobian of the residual, where each row of the matrix corresponds to a
     *  point in the data.
     */
    Int TwoDOptimization::jacobian2D_(const gsl_vector* x, void* params, gsl_matrix* J)
    {
      // For the conventions on x and params c.f. the commentary in residual()
      //
      // The matrix J is supposed to contain the result when we return from this function.
      // Note: GSL expects the Jacobian as follows:
      //          - each row corresponds to one data point
      //          - each column corresponds to one parameter

      DoubleReal computed_signal, current_position, experimental_signal,last_position,step;
      DoubleReal p_height, p_position, p_width;
      DoubleReal diff, denom_inv,ddl_left,ddl_right,ddx0,sinh_term;
      Int count =0;
      Int counter_posf=0;

			std::vector<std::pair<SignedSize,SignedSize> >& signal2D = static_cast<TwoDOptimization::Data*> (params) ->signal2D; 
			std::multimap<DoubleReal,IsotopeCluster>::iterator iso_map_iter=static_cast<TwoDOptimization::Data*> (params) ->iso_map_iter;
			Size total_nr_peaks=static_cast<TwoDOptimization::Data*> (params) ->total_nr_peaks;
			std::map<Int, std::vector<PeakIndex> >& matching_peaks=static_cast<TwoDOptimization::Data*> (params) ->matching_peaks;
			std::vector<DoubleReal> ov_weight(matching_peaks.size(),0);
			MSExperiment<> &picked_peaks = static_cast<TwoDOptimization::Data*> (params) ->picked_peaks;
			MSExperiment<Peak1D>::ConstIterator raw_data_first = static_cast<TwoDOptimization::Data*> (params) ->raw_data_first;
			OptimizationFunctions::PenaltyFactorsIntensity& penalties=static_cast<TwoDOptimization::Data*> (params) ->penalties;
// 			std::vector<DoubleReal> &positions=static_cast<TwoDOptimization::Data*> (params) ->positions;
// 			std::vector<DoubleReal> &signal=static_cast<TwoDOptimization::Data*> (params) ->signal;
      IsotopeCluster::ChargedIndexSet::iterator peak_iter = iso_map_iter->second.peaks.begin();
			Size num_scans = signal2D.size()/2;
      //iterate over all scans
      for (Size current_scan = 0; current_scan < num_scans; ++current_scan)
				{
					Size curr_scan_idx = current_scan + iso_map_iter->second.peaks.begin()->first;
					current_position = ((raw_data_first 
															 + signal2D[2*current_scan].first)->begin() 
															+ signal2D[2*current_scan].second)->getMZ();
					// iterate over all points of the signal
					for (Int current_point = 1;
							 current_point +  signal2D[2*current_scan].second
								 <= signal2D[2*current_scan+1].second;
							 ++current_point)
						{
							last_position = current_position;
							computed_signal   = 0.;
							current_position  = ((raw_data_first 
																		+ signal2D[2*current_scan].first)->begin() 
																	 + signal2D[2*current_scan].second+current_point)->getMZ();
							experimental_signal = ((raw_data_first + signal2D[2*current_scan].first)->begin()
																		 + signal2D[2*current_scan].second+current_point)->getIntensity();

							step = current_position - last_position;
							
#ifdef DEBUG_2D
          std::cout << "experimental signal rt "<<(raw_data_first 
                        + signal2D[2*current_scan].first)->getRT()
                    << "\tmz " << ((raw_data_first 
                        + signal2D[2*current_scan].first)->begin()
																	 + signal2D[2*current_scan].second+current_point)->getMZ()
                    << "\tint " << ((raw_data_first 
                        + signal2D[2*current_scan].first)->begin()
                        + signal2D[2*current_scan].second+current_point)->getIntensity()<<std::endl;

#endif

							size_t current_peak = 0;
							peak_iter = iso_map_iter->second.peaks.begin();
							while(peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first != curr_scan_idx) ++peak_iter;
							//iterate over all peaks of the current scan
							while(peak_iter != iso_map_iter->second.peaks.end() && peak_iter->first == curr_scan_idx)
								{
									Int peak_idx = distance(iso_map_iter->second.peaks.begin(),peak_iter);
									DoubleReal mz_in_hash = ((picked_peaks[peak_iter->first]).begin() + peak_iter->second)->getMZ() * 10;
									std::map<Int,std::vector<PeakIndex> >::iterator  m_spec_iter =	matching_peaks.begin();
									Int map_idx=0;
									while(m_spec_iter->first != (Int)(mz_in_hash+0.5) )
										{
											++map_idx;
											++m_spec_iter;
										}
									// if the current peak is in the reference scan take all parameters from the vector x
#ifdef DEBUG_2D
									std::cout << "ref_scan : "<< current_peak << "\t "
														<<gsl_vector_get(x,3*current_peak) << "\t" <<gsl_vector_get(x,3*current_peak+1)
														<< "\t"<<gsl_vector_get(x,3*current_peak+2)<<std::endl;
#endif
									// Store the current parameters for this peak
									p_position    = gsl_vector_get(x,total_nr_peaks+3*map_idx);
									p_height      = gsl_vector_get(x,peak_idx);
									p_width       = (current_position <= p_position) ?
										gsl_vector_get(x,total_nr_peaks+3*map_idx+1) :
										gsl_vector_get(x,total_nr_peaks+3*map_idx+2);
									++count;
									DoubleReal weight = step*((picked_peaks[peak_iter->first]).begin()+peak_iter->second)->getIntensity();
									ov_weight[map_idx] += weight;
									DoubleReal ddx0_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx);
									DoubleReal ddl_left_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx+1);
									DoubleReal ddl_right_old = gsl_matrix_get(J, counter_posf, total_nr_peaks +3*map_idx+2);
									//is it a Lorentz or a Sech - Peak?
									
									if ((PeakShape::Type)(Int)Math::round((picked_peaks[peak_iter->first]).getFloatDataArrays()[5][peak_iter->second]) == PeakShape::LORENTZ_PEAK)
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
											gsl_matrix_set(J, counter_posf, peak_idx                     , step*denom_inv);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +1 , ddl_left * weight + ddl_left_old);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +2 , ddl_right * weight + ddl_right_old);
							
											
// 											diff      = current_position - p_position;

// 											DoubleReal denom1 = 1. + pow(p_width*diff,2);
// 											DoubleReal help = p_height / denom1 - experimental_signal;
// 											// partial derivative with respect to the height,...
// 											denom_inv = (2.*help)/denom1;
// 											// left width,...
// 											ddl_left  = (current_position <= p_position)
// 												? (-4 * help * p_height*p_width*pow(diff,2) ) / pow(denom1, 2)
// 												: 0;
// 											// right width ...
// 											ddl_right = (current_position  > p_position)
// 												? (-4 * help * p_height*p_width*pow(diff,2) ) / pow(denom1, 2)
// 												: 0;

// 											// and position
// 											ddx0 = (4 * help * p_height*pow(p_width,2)*diff);

// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx    , ddx0* weight  + ddx0_old);
// 											gsl_matrix_set(J, counter_posf, peak_idx                     , denom_inv);
// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +1 , ddl_left * weight + ddl_left_old);
// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx +2 , ddl_right * weight + ddl_right_old);
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
											gsl_matrix_set(J, counter_posf, peak_idx                      , step*pow(denom_inv, 2));
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 1, ddl_left * weight + ddl_left_old);
											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 2, ddl_right * weight + ddl_right_old);

											
// 											diff      = current_position - p_position;
// 											DoubleReal cosh_term = cosh(p_width*diff);
// 											DoubleReal enum_term = p_height / pow(cosh_term,2) - experimental_signal;
// 											denom_inv = 1./cosh(p_width*diff);
											
// 											// The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
// 											// and can assume that all derivatives vanish
// 											sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width*diff);


// 											ddl_left  = (current_position <= p_position)
// 												? -4 * enum_term * p_height * sinh_term*diff * pow(denom_inv,3) : 0;
// 											ddl_right = (current_position  > p_position)
// 												? -4 * enum_term * p_height * sinh_term*diff * pow(denom_inv,3) : 0;

// 											ddx0      = 4*enum_term*p_height*p_width*sinh_term * pow(denom_inv,3);

// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx     , ddx0 * weight + ddx0_old);
// 											gsl_matrix_set(J, counter_posf, peak_idx                      , 2*enum_term*pow(denom_inv, 2));
// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 1, ddl_left * weight + ddl_left_old);
// 											gsl_matrix_set(J, counter_posf, total_nr_peaks +3*map_idx  + 2, ddl_right * weight + ddl_right_old);

										}
									++current_peak;
									++peak_iter;

								}// end while

							++counter_posf;
						}

				}

      for(Size cluster=0; cluster < matching_peaks.size();++cluster)
				{
					for(UInt j=0; j < J->size1-1;++j)
						{
							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster  )/ov_weight[cluster]);
							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster + 1,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster+1)/ov_weight[cluster]);

							gsl_matrix_set(J, j, total_nr_peaks + 3*cluster + 2,
														 gsl_matrix_get(J, j, total_nr_peaks + 3*cluster+2)/ov_weight[cluster]);
						}
				}

      //iterate over all peaks again to compute the penalties
      // first look at all positions and width parameters
      UInt peak=0,current_peak=0;
      std::map<Int, std::vector<PeakIndex> >::iterator map_iter=matching_peaks.begin();

      for (;map_iter != matching_peaks.end(); ++map_iter)
				{
					std::vector<PeakIndex>::iterator vec_iter	= map_iter->second.begin();
					DoubleReal old_position = 0,old_width_l=0,old_width_r=0;
					DoubleReal weight =0;
					DoubleReal old_height,p_height;
					DoubleReal penalty_h=0, penalty_l=0, penalty_r=0,penalty_p=0;
					for(;vec_iter != map_iter->second.end();++vec_iter)
						{
							old_height = (vec_iter)->getPeak(picked_peaks).getIntensity();
							weight += old_height;
							old_position += (vec_iter)->getPeak(picked_peaks).getMZ() * old_height;
							old_width_l += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[3][vec_iter->peak]* old_height;
							old_width_r += picked_peaks[vec_iter->spectrum].getFloatDataArrays()[4][vec_iter->peak]* old_height;
	
							p_height     = gsl_vector_get(x, peak);


							DoubleReal penalty_height = 2.*penalties.height*fabs(p_height-old_height);
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

					//					std::cout << old_position << "vs. ";
					DoubleReal p_position   = gsl_vector_get(x, total_nr_peaks+3*current_peak);
					DoubleReal p_width_l    = gsl_vector_get(x, total_nr_peaks+3*current_peak +1);
					DoubleReal p_width_r    = gsl_vector_get(x, total_nr_peaks+3*current_peak +2);
					DoubleReal penalty_lwidth = 2.*penalties.lWidth*fabs(p_width_l - old_width_l);
					DoubleReal penalty_rwidth = 2.*penalties.rWidth*fabs(p_width_r - old_width_r);
					DoubleReal penalty_pos    = 2.*penalties.pos*fabs(p_position-old_position);
					//std::cout << p_position<<std::endl;
#ifdef DEBUG_2D
					std::cout << "penalty_lwidth " << penalty_lwidth << "penalty_rwidth " << penalty_rwidth
										<< "penalty_pos " << penalty_pos << std::endl;
#endif
					if(p_width_l < 0 )
						{
							penalty_l += 1e7*penalty_lwidth;
						}
					else if(p_width_l < 1) penalty_l += 2000*penalties.lWidth*(fabs(p_width_l - old_width_l));
					if(p_width_r < 0 )
						{
							penalty_r += 1e7*penalty_rwidth;
						}
					else if(p_width_r < 1) penalty_r += 2000*penalties.rWidth*(fabs(p_width_r - old_width_r));
					if(p_position < 0)
						{
							penalty_p +=200*penalty_pos;
						}

					if(fabs(old_position-p_position) > 0.2)
						{
							
							penalty_p +=2000*penalties.pos * fabs(p_position - old_position);
						
						}
					if(fabs(old_width_r-p_width_r) > 1)
						{
							penalty_r +=1000*penalty_rwidth;
						}
					if(fabs(old_width_l-p_width_l) > 1)
						{
							penalty_l +=1000*penalty_lwidth;
						}
					
					gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak+1, penalty_l);
					gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak+2, penalty_r);
					gsl_matrix_set(J, counter_posf, total_nr_peaks+3*current_peak  , penalty_p);



					++current_peak;
				}
      return GSL_SUCCESS;
    }

    //Driver function for the evaluation of function and jacobian.
    Int TwoDOptimization::evaluate2D_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      residual2D_(x, params, f);
      jacobian2D_(x, params, J);

      return GSL_SUCCESS;
    }


	TwoDOptimization::TwoDOptimization()
		:DefaultParamHandler("TwoDOptimization")
	{
		// 2D optimization parameters
		defaults_.setValue("penalties:position",0.0,"If the position changes more than 0.2Da during the fitting it can be penalized");
		defaults_.setValue("penalties:height",1.0,"penalty term for the fitting of the intensity:"\
											 "If it gets negative during the fitting it can be penalized.");
		defaults_.setValue("penalties:left_width",0.0,"penalty term for the fitting of the left width:"\
											 "If the left width gets too broad or negative during the fitting it can be penalized.");
		defaults_.setValue("penalties:right_width",0.0,"penalty term for the fitting of the right width:"\
											 "If the right width gets too broad or negative during the fitting it can be penalized.");
		defaults_.setValue("2d:tolerance_mz",2.2,"mz tolerance for cluster construction", StringList::create("advanced"));
		defaults_.setValue("2d:max_peak_distance",1.2,"maximal peak distance in mz in a cluster", StringList::create("advanced"));
		defaults_.setValue("delta_abs_error",1e-05f,"if the absolute error gets smaller than this value the fitting is stopped.", StringList::create("advanced"));
		defaults_.setValue("delta_rel_error",1e-05f,"if the relative error gets smaller than this value the fitting is stopped.", StringList::create("advanced"));
		defaults_.setValue("iterations",10,"maximal number of iterations for the fitting step");


		defaultsToParam_();
		updateMembers_();
	}
	
	TwoDOptimization::TwoDOptimization(const TwoDOptimization& opt)
		:DefaultParamHandler(opt)
	{
		updateMembers_();
	}
	
	
	TwoDOptimization& TwoDOptimization::operator=(const TwoDOptimization& opt)
	{
		if(&opt == this) return *this; 
		DefaultParamHandler::operator=(opt);
		updateMembers_();

		return *this;
	}

	
  void TwoDOptimization::findMatchingPeaks_(std::multimap<DoubleReal, IsotopeCluster>::iterator& it, MSExperiment<>& ms_exp)
  {
  	IsotopeCluster::ChargedIndexSet::const_iterator iter = it->second.peaks.begin();
    for(; iter != it->second.peaks.end(); ++iter)
			{

				DoubleReal mz = (ms_exp[iter->first][iter->second]).getMZ();
				mz *= 10;
				matching_peaks_[(Int)(mz+0.5)].push_back(PeakIndex(iter->first,iter->second));
			}

		
#ifdef DEBUG_2D
    std::map<Int, PeakIndex >::iterator it2 = matching_peaks_.begin();
    for(;it2 != matching_peaks_.end();++it2)
			{
				std::cout << it2->first << " has "<<it2->second.size()<<" elements:"<<std::endl;
				for(Size i=0;i<it2->second.size();++i) std::cout << it2->second[i]->getPeak(ms_exp).getMZ()<<"\t";
				std::cout<<std::endl;
			}
#endif

  }

  // Finds the neighbour of the peak denoted by @p current_mz in the previous scan
  std::vector<DoubleReal>::iterator TwoDOptimization::searchInScan_(std::vector<DoubleReal>::iterator scan_begin,
																																std::vector<DoubleReal>::iterator scan_end ,
																																DoubleReal current_mz)
  {

    // perform binary search to find the neighbour in rt dimension
    //  lower_bound finds the peak with m/z current_mz or the next larger peak if this peak does not exist.
    std::vector<DoubleReal>::iterator insert_iter = lower_bound(scan_begin,scan_end,current_mz);

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
						DoubleReal delta_mz = fabs(*insert_iter - current_mz);
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

	void TwoDOptimization::updateMembers_()
	{
		penalties_.height = (DoubleReal)param_.getValue("penalties:height");
		penalties_.pos = (DoubleReal)param_.getValue("penalties:position");
		penalties_.lWidth = (DoubleReal)param_.getValue("penalties:left_width");
		penalties_.rWidth = (DoubleReal)param_.getValue("penalties:right_width");
		max_peak_distance_ = (DoubleReal)param_.getValue("2d:max_peak_distance");
		tolerance_mz_ = (DoubleReal)param_.getValue("2d:tolerance_mz");
		eps_abs_= (DoubleReal)param_.getValue("delta_abs_error");
    eps_rel_= (DoubleReal)param_.getValue("delta_rel_error");
		max_iteration_= (UInt)param_.getValue("iterations");;
		
	}

}
