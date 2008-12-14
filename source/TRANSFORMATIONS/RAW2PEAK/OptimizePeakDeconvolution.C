// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <boost/math/special_functions/acosh.hpp>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePeakDeconvolution.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{

	
  const double OptimizePeakDeconvolution::dist_ = 1.003;  
  namespace OptimizationFunctions
  {
    int charge;

   
    
    std::vector<PeakShape> peaks_DC_;
    std::vector<double> positions_DC_;
    std::vector<double> signal_DC_;

    
    

		// Evaluation of the target function for nonlinear optimization.
    int residualDC(const gsl_vector* x, void* params , gsl_vector* f)
    {
      // According to the gsl conventions, x contains the parameters to be optimized.
      // The first two entries are the left and right width, respectively.They are equal
      // for all peaks. Then the height and position of all peaks are stored.
      //
      // Params might contain any additional parameters. We handle these using class members
      // instead.
      // The vector f is supposed to contain the result when we return from this function.
      // Note: GSL wants the values for each data point i as one component of the results vector


      double leftwidth = gsl_vector_get(x,0);
      double rightwidth = gsl_vector_get(x,1);
      //double posP1 = gsl_vector_get(x,2);
      
      // iterate over all points of the signal
      for (size_t current_point = 0; current_point < positions_DC_.size(); current_point++)
				{
					double computed_signal     = 0.;
					double current_position    = positions_DC_[current_point];
					double experimental_signal = signal_DC_[current_point];	
    
					//iterate over all peaks
					for (size_t current_peak = 0; current_peak < peaks_DC_.size(); current_peak++)
						{
							//Store the current parameters for this peak
							double p_height 	   = gsl_vector_get(x,2+ 2*current_peak);
							double p_position    = gsl_vector_get(x,2+ 2*current_peak +1);
							double p_width  	   = (current_position <= p_position) ? leftwidth : rightwidth;

							//is it a Lorentz or a Sech - Peak?
							if (peaks_DC_[current_peak].type == PeakShape::LORENTZ_PEAK)
								{
									computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
								}
							else // It's a Sech - Peak 
								{
									computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
								}
						}
					gsl_vector_set(f, current_point, computed_signal - experimental_signal);
				}

      // penalties : especially negative heights have to be penalised
      double penalty = 0.;
      PenaltyFactorsIntensity* penalties = (PenaltyFactorsIntensity *)params;
      double penalty_pos    = penalties->pos;
      double penalty_lwidth = penalties->lWidth;
      double penalty_rwidth = penalties->rWidth;
      double penalty_intensity = penalties->height;


      //iterate over all peaks again to compute the penalties
      for (size_t current_peak = 0; current_peak < peaks_DC_.size(); current_peak++)
      	{
 					double p_position = gsl_vector_get(x, 2+2*current_peak+1);
					if(current_peak < peaks_DC_.size()-1)
						{
	      
							double next_p_position  = gsl_vector_get(x, 2+2*current_peak+3);
							// if distance between peaks does not match the peptide mass rule
							if( fabs(fabs(p_position - next_p_position) - 1.003/OptimizationFunctions::charge) > 0.05)
								{
									// penalize it
									penalty +=  penalty_pos * 10000
										* pow(fabs(fabs(p_position - next_p_position) - 1.003/OptimizationFunctions::charge),2);
								}
						}
					double old_position   = peaks_DC_[current_peak].mz_position;
					double old_width_l 	= peaks_DC_[current_peak].left_width;
      	  double old_width_r 	= peaks_DC_[current_peak].right_width;
					double old_height    = peaks_DC_[current_peak].height;
	  
					double p_width_l    = gsl_vector_get(x, 0);
      	  double p_width_r    = gsl_vector_get(x, 1);
					double p_height     = gsl_vector_get(x,2+2*current_peak);

					if(p_height <  1)
						{
							penalty += 100000*penalty_intensity*pow(fabs(p_height - old_height),2);

						}
					if(p_width_l < 0 )
						{
							penalty += penalty_lwidth * peaks_DC_.size()*10000*pow(fabs(p_width_l - old_width_l),2);
						}
					else if (p_width_l < 1.5 ) penalty += 10000*pow(fabs(p_width_l - old_width_l),2);
					if(p_width_r < 0 )
						{
							penalty += penalty_rwidth *peaks_DC_.size()*10000*pow(fabs(p_width_r - old_width_r),2);
						}
					else if (p_width_r < 1.5 ) penalty += 10000*pow(fabs(p_width_r - old_width_r),2);
					if(fabs(old_position - p_position) > 0.5)
						{
							penalty += 10000*penalty_pos*pow(fabs(old_position - p_position),2);
						}


	  
				}
      gsl_vector_set(f, f->size -1, penalty);
      return GSL_SUCCESS;
    }

    /** Compute the Jacobian of the residual, where each row of the matrix corresponds to a 
     *  point in the data.
     */
    int jacobianDC(const gsl_vector* x, void* params, gsl_matrix* J)
    {
      // For the conventions on x and params c.f. the commentary in residual()
      //  
      // The matrix J is supposed to contain the result when we return from this function.
      // Note: GSL expects the Jacobian as follows:
      // 					- each row corresponds to one data point
      // 					- each column corresponds to one parameter


      double leftwidth = gsl_vector_get(x,0);
      double rightwidth = gsl_vector_get(x,1);
 

      gsl_matrix_set_zero(J);
    
      
      // iterate over all points of the signal
      for (size_t current_point = 0; current_point < positions_DC_.size(); current_point++)
				{
					double current_position    = positions_DC_[current_point];

					// iterate over all peaks
					for (size_t current_peak = 0; current_peak < peaks_DC_.size(); current_peak++)
						{

	      
							//Store the current parameters for this peak
							double p_height 	   = gsl_vector_get(x, 2+2*current_peak);
							double p_position    = gsl_vector_get(x, 2+2*current_peak+1);
							double p_width  	    = (current_position <= p_position) ? leftwidth : rightwidth;

							//is it a Lorentz or a Sech - Peak?
							if (peaks_DC_[current_peak].type == PeakShape::LORENTZ_PEAK)
								{
									double diff      = current_position - p_position;
									double denom_inv = 1./(1. + pow(p_width * diff, 2));
		  
									double ddl_left  = (current_position <= p_position) 
										? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) 
										: 0;
		  
									double ddl_right = (current_position  > p_position) 
										? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2) 
										: 0;

									// left and right width are the same for all peaks,
									// the sums of the derivations over all peaks are stored in the first two columns
									gsl_matrix_set(J, current_point, 0, gsl_matrix_get(J,current_point,0) + ddl_left);
									gsl_matrix_set(J, current_point, 1, gsl_matrix_get(J,current_point,1) + ddl_right); 
		  
									double ddx0	   = 2*p_height*pow(p_width,2)*diff*pow(denom_inv, 2);

									// partial derivation with respect to intensity
									gsl_matrix_set(J, current_point, 2 + 2*current_peak,     denom_inv);

									// partial derivation with respect to the mz-position
									gsl_matrix_set(J, current_point, 2 + 2*current_peak+1 , ddx0);
								}
							else // It's a Sech - Peak 
								{
									double diff      = current_position - p_position;
									double denom_inv = 1./cosh(p_width*diff);

									// The remaining computations are not stable if denom_inv == 0. In that case, we are far away from the peak
									// and can assume that all derivatives vanish
									double sinh_term = (fabs(denom_inv) < 1e-6) ? 0.0 : sinh(p_width*diff);


									double ddl_left  = (current_position <= p_position) 
										? -2 * p_height * sinh_term*diff * pow(denom_inv,3) : 0;
									double ddl_right = (current_position  > p_position) 
										? -2 * p_height * sinh_term*diff * pow(denom_inv,3) : 0;
		  
									gsl_matrix_set(J, current_point, 0, gsl_matrix_get(J,current_point,0) + ddl_left);
									gsl_matrix_set(J, current_point, 1, gsl_matrix_get(J,current_point,1) + ddl_right); 
		 
									double ddx0      = 2*p_height*p_width*sinh_term * pow(denom_inv,3);

									gsl_matrix_set(J, current_point, 2 + 2*current_peak, pow(denom_inv, 2));
									gsl_matrix_set(J, current_point, 2 + 2*current_peak+1,ddx0 );
								}
						}
				}

    
      /** Now iterate over all peaks again to compute the
       *  penalties.
       */
      PenaltyFactorsIntensity* penalties = (PenaltyFactorsIntensity *)params;
      for (size_t current_peak = 0; current_peak < peaks_DC_.size(); current_peak++)
      	{

	  
					double penalty_p = 0;
		 		 	double p_position = gsl_vector_get(x, 2+2*current_peak+1);
					if(current_peak < peaks_DC_.size()-1)
						{
	      
							double next_p_position  = gsl_vector_get(x, 2+2*current_peak+3);
							// if distance between peaks does not match the peptide mass rule
							if( fabs(fabs(p_position - next_p_position) - 1.003/OptimizationFunctions::charge) > 0.05)
								{
									// penalize it
									penalty_p += penalties->pos * 20000
										* fabs(fabs(p_position - next_p_position) - 1.003/OptimizationFunctions::charge);
		  
								}
						}
					//  std::cout << "penalty_p "<<penalty_p<<std::endl;
      	  double p_width_left = gsl_vector_get(x,0); 
      	  double p_width_right = gsl_vector_get(x,1); 
					double p_height   = gsl_vector_get(x,2+2*current_peak);

					double old_position    = peaks_DC_[current_peak].mz_position;
      	  double old_width_left  = peaks_DC_[current_peak].left_width;
      	  double old_width_right = peaks_DC_[current_peak].right_width;
					double old_height      = peaks_DC_[current_peak].height;

					double penalty_h =0., penalty_l=0., penalty_r=0.;
					if(p_height < 1)
						{
							penalty_h += 100000*2*penalties->height *(fabs(p_height) - fabs(old_height));
						}

					if(p_width_left < 0 )
						{
							penalty_l += peaks_DC_.size()*2*penalties->lWidth*10000*(fabs(p_width_left - old_width_left));
						}
					else if (p_width_left < 1.5 ) penalty_l += 2*penalties->lWidth*10000*pow(fabs(p_width_left - old_width_left),2);
					if(p_width_right < 0 )
						{
							penalty_r += peaks_DC_.size()*2*penalties->rWidth*10000*(fabs(p_width_right - old_width_right));
						}
					else if (p_width_right < 1.5 ) penalty_r += 2*penalties->rWidth*10000*pow(fabs(p_width_right - old_width_right),2);
					if(fabs(old_position - p_position) > 0.5)
						{
							penalty_p += 10000*penalties->pos*2*fabs(old_position - p_position);
						}
	
	  
	  
					gsl_matrix_set(J, positions_DC_.size(), 2+2*current_peak, 100*penalty_h);
      	  gsl_matrix_set(J, positions_DC_.size(), 0, 100*penalty_l);
					gsl_matrix_set(J, positions_DC_.size(), 1, 100*penalty_r);
      	  gsl_matrix_set(J, positions_DC_.size(), 2+2*current_peak+1, 100*penalty_p);
      	}
	
      return GSL_SUCCESS;
    }

    // Driver function for the evaluation of function and jacobian.
    int evaluateDC(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      residualDC(x, params, f);
      jacobianDC(x, params, J);

      return GSL_SUCCESS;
    }
    

  }// namespace OptimizationFunctions

	
 
	OptimizePeakDeconvolution::OptimizePeakDeconvolution( )
		: DefaultParamHandler("OptimizePeakDeconvolution"),charge_(1)
	{
		
		defaults_.setValue("max_iteration",10,"maximal number of iterations for the fitting step");
		defaults_.setValue("eps_abs",1e-04,"if the absolute error gets smaller than this value the fitting is stopped", StringList::create("advanced"));
		defaults_.setValue("eps_rel",1e-04,"if the relative error gets smaller than this value the fitting is stopped", StringList::create("advanced"));

		defaults_.setValue("penalties:left_width",0.0,"penalty term for the fitting of the left width:"\
											 "If the left width gets too broad or negative during the fitting it can be penalized.");
		defaults_.setValue("penalties:right_width",0.0,"penalty term for the fitting of the right width:"\
											 "If the right width gets too broad or negative during the fitting it can be penalized.");
		defaults_.setValue("penalties:height",0.0,"penalty term for the fitting of the intensity:"\
											 "If it gets negative during the fitting it can be penalized.");
		defaults_.setValue("penalties:position",0.0,"penalty term for the fitting of the peak position:"\
											 "If the position changes more than 0.5Da during the fitting it can be penalized as well as "\
											 "discrepancies of the peptide mass rule.");

		defaults_.setValue("fwhm_threshold",1.0,"If a peaks is broader than fwhm_threshold, it is assumed that it contains another peaks and an additional peak is added.");

		defaultsToParam_();
	}


	void OptimizePeakDeconvolution::updateMembers_()
	{
		penalties_.rWidth = (float)param_.getValue("penalties:right_width");
		penalties_.lWidth = (float)param_.getValue("penalties:left_width");
		penalties_.height = (float)param_.getValue("penalties:height");
		penalties_.pos    = (float)param_.getValue("penalties:position");
	
	}

	bool OptimizePeakDeconvolution::optimize(std::vector<PeakShape>& peaks, int failure)
	{
      
		if (peaks.size() == 0)	return true;


#ifdef DEBUG_DECONV
		std::cout<<"peaksanzahl:"<<peaks.size();
		std::cout<<"\tpeaks[0].mz_position:"<<peaks[0].mz_position<<std::endl;
      
		for(unsigned int j=0;j<peaks.size();++j)
			{
				std::cout<<"\tpeaks[j].mz_position:"<<peaks[j].mz_position;
				std::cout<<"\tpeaks[j].height:"<<peaks[j].height<<std::endl;
				std::cout<<"\tpeaks[j].left_width:"<<peaks[j].left_width;
				std::cout<<"\tpeaks[j].right_width:"<<peaks[j].right_width<<std::endl<<std::endl;
			}
      
		for(unsigned int j=0;j<OptimizationFunctions::positions_DC_.size();++j)
			{
				std::cout<<"positions_DC_["<<j<<"]="<<OptimizationFunctions::positions_DC_[j]<<std::endl;
			}
      
#endif
      
		// the input peaks are stored in a temporary vector
		std::vector<PeakShape> temp_shapes = peaks;
      
      
		int global_peak_number=0;

		double min = DBL_MAX;
		int best_charge,num_peaks;
		int best_num_peaks;
		gsl_vector *best_result=gsl_vector_alloc(2+2*OptimizationFunctions::peaks_DC_.size());;


		// try three different charge states : charge-1, charge, charge +1
		// take the best solution
		int l = (charge_-1 > 1) ? charge_-1 : charge_;
		int start_l = l;
#ifdef DEBUG_DECONV
		std::cout<<"charge "<<l<<" max_charge" << charge_+1
						 <<"\tpeaks.size() "<<peaks.size()<<std::endl;
#endif
		best_charge = l;
		best_num_peaks = peaks.size();
		num_peaks = peaks.size();
		for(;l<charge_+2;++l)
			{


				num_peaks = getNumberOfPeaks_(l, temp_shapes);
#ifdef DEBUG_DECONV
				std::cout<<"charge "<<l<<" #peaks "<<num_peaks<<"\tpeaks.size() "
								 <<OptimizationFunctions::peaks_DC_.size()<<std::endl;
#endif
				gsl_vector *start_value;
				// the vector storing the start values for the parameters has to be filled
				// differently depending on the usage of the peptide mass rule
				start_value=gsl_vector_alloc(2+2*OptimizationFunctions::peaks_DC_.size());
				for (size_t i = 0; i < OptimizationFunctions::peaks_DC_.size(); i++)
					{
						gsl_vector_set(start_value, 2+2*i, OptimizationFunctions::peaks_DC_[i].height);
						gsl_vector_set(start_value, 3+2*i, OptimizationFunctions::peaks_DC_[i].mz_position);
					}
	    
	  
				// Initialize the parameters for the optimization 

				// all peaks shall have the same width
				double wl = OptimizationFunctions::peaks_DC_[0].left_width;
				double wr = OptimizationFunctions::peaks_DC_[0].right_width;
				if (isnan(wl))
					{
						for(unsigned int i=0;i<OptimizationFunctions::peaks_DC_.size();++i)
							{
								OptimizationFunctions::peaks_DC_[i].left_width = 1;
							}
						wl = 1.;
					}
				if (isnan(wr))
					{
						for(unsigned int i=0;i<OptimizationFunctions::peaks_DC_.size();++i)
							{
								OptimizationFunctions::peaks_DC_[i].right_width = 1;
							}
						wr = 1.;
					}
	  
				gsl_vector_set(start_value, 0, wl);
				gsl_vector_set(start_value, 1, wr);
	 	  
 		
				// The gsl algorithms require us to provide function pointers for the evaluation of
				// the target function.
				gsl_multifit_function_fdf fit_function;

				fit_function.f      = OptimizationFunctions::residualDC;
				fit_function.df     = OptimizationFunctions::jacobianDC;
				fit_function.fdf    = OptimizationFunctions::evaluateDC;
	  
				fit_function.n      = std::max(OptimizationFunctions::positions_DC_.size()+1,
																			 2+2*OptimizationFunctions::peaks_DC_.size());
	  
				fit_function.p	  = 2+2*OptimizationFunctions::peaks_DC_.size();

				fit_function.params = &penalties_;
#ifdef DEBUG_DECONV
				std::cout<<"fit_function.p "<<fit_function.p<<"\t fit_function.n "<<fit_function.n<<std::endl;
				std::cout<<"peaks.size() "<<OptimizationFunctions::peaks_DC_.size()<<std::endl;
#endif
				const gsl_multifit_fdfsolver_type *type = gsl_multifit_fdfsolver_lmsder;
				gsl_multifit_fdfsolver *fit;
				fit= gsl_multifit_fdfsolver_alloc(type,
																					std::max(OptimizationFunctions::positions_DC_.size()+1,
																									 2+2*OptimizationFunctions::peaks_DC_.size()),
																					2+2*OptimizationFunctions::peaks_DC_.size());
	  

				OptimizationFunctions::charge = l;
	  

				gsl_multifit_fdfsolver_set(fit, &fit_function, start_value);

#ifdef DEBUG_DECONV
				// initial norm 
				std::cout << "Before optimization: ||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
#endif
				// Iteration
				int iteration = 0;
				int status;

				do 
					{
						iteration++;
						status = gsl_multifit_fdfsolver_iterate(fit);

#ifdef DEBUG_DECONV
						std::cout << "Iteration " << iteration << "; Status " << gsl_strerror(status) << "; " << std::endl;
						std::cout << "||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
						std::cout << "Number of parms: " << OptimizationFunctions::peaks_DC_.size() + 3 << std::endl;
						std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
#endif
						if (isnan(gsl_blas_dnrm2(fit->dx)))
							break;
				
						status = gsl_multifit_test_delta(fit->dx, fit->x,(float)param_.getValue("eps_abs"),
																						 (float)param_.getValue("eps_rel"));
						if (status != GSL_CONTINUE)
							break;
						if(!checkFWHM_(peaks,fit) && failure <1) return false;
					}
				while (status == GSL_CONTINUE && iteration < (int)param_.getValue("max_iteration"));

	
				double chi = gsl_blas_dnrm2(fit->f);
#ifdef DEBUG_DECONV
				std::cout << "Finished! Charge " << l << "\tIterations: "<<iteration<< std::endl;
				std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
	  
				std::cout << "chisq/dof = " << pow(chi, 2.0) / (OptimizationFunctions::positions_DC_.size()
																												- (3+OptimizationFunctions::peaks_DC_.size()));
				std::cout << "\nAfter optimization: ||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
#endif
				if((l == start_l) || (chi < min))
					{
						if(l!= start_l)  gsl_vector_free(best_result);
						best_result=gsl_vector_alloc(2+2*OptimizationFunctions::peaks_DC_.size());
						gsl_vector_memcpy(best_result,fit->x);
						min = chi;
						best_charge = l;
						best_num_peaks = OptimizationFunctions::peaks_DC_.size();
					}
				iteration = 0;

	   
	   
				gsl_multifit_fdfsolver_free(fit);
				gsl_vector_free(start_value);
	   
			}
		global_peak_number += best_num_peaks;
		// iterate over all peaks and store the optimized values in peaks
		if(best_num_peaks > 0)
			{
				peaks.resize(best_num_peaks);
				for (int current_peak = 0; current_peak < best_num_peaks; current_peak++)
					{
	      
						// Store the current parameters for this peak
	      
						peaks[current_peak].left_width  = gsl_vector_get(best_result, 0);
						peaks[current_peak].right_width = gsl_vector_get(best_result, 1);
	      
						peaks[current_peak].height      = gsl_vector_get(best_result,2+ 2*current_peak);
						peaks[current_peak].mz_position = gsl_vector_get(best_result, 2+ 2*current_peak+1);

	    
	  
						//	compute the area
						//  is it a Lorentz or a Sech - Peak? 
						if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
							{
								PeakShape p = peaks[current_peak];
								double x_left_endpoint=p.mz_position+1/p.left_width*sqrt(p.height/1-1);
								double x_right_endpoint=p.mz_position+1/p.right_width*sqrt(p.height/1-1);
#ifdef DEBUG_DECONV
								std::cout<<"x_left_endpoint "<<x_left_endpoint<<" x_right_endpoint "<<x_right_endpoint<<std::endl;
								std::cout<<"p.height"<<p.height<<std::endl;
#endif
								double area_left=-p.height/p.left_width*atan(p.left_width*(x_left_endpoint-p.mz_position));
								double area_right=-p.height/p.right_width*atan(p.right_width*(p.mz_position-x_right_endpoint));
								peaks[current_peak].area=area_left+area_right;
		  
							}
						else  //It's a Sech - Peak 
							{
								PeakShape p = peaks[current_peak];
								double x_left_endpoint=p.mz_position+1/p.left_width*boost::math::acosh(sqrt(p.height/0.001));
								double x_right_endpoint=p.mz_position+1/p.right_width*boost::math::acosh(sqrt(p.height/0.001));
#ifdef DEBUG_DECONV
								std::cout<<"x_left_endpoint "<<x_left_endpoint<<" x_right_endpoint "<<x_right_endpoint<<std::endl;
								std::cout<<"p.height"<<p.height<<std::endl;
#endif
								double area_left=-p.height/p.left_width*(sinh(p.left_width*(p.mz_position-x_left_endpoint))
																												 /cosh(p.left_width*(p.mz_position-x_left_endpoint)));
								double area_right=-p.height/p.right_width*(sinh(p.right_width*(p.mz_position-x_right_endpoint))
																													 /cosh(p.right_width*(p.mz_position-x_right_endpoint)));
								peaks[current_peak].area=area_left+area_right;
	      
							}
	      
					}
			}
		charge_ = best_charge;
		gsl_vector_free(best_result);
      
		return true;
	}


  bool OptimizePeakDeconvolution::checkFWHM_(std::vector<PeakShape>& peaks, gsl_multifit_fdfsolver *& fit)
  {
    double fwhm_threshold = (double)param_.getValue("fwhm_threshold");
   
    PeakShape p;
    for (unsigned int current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
				p.left_width  = gsl_vector_get(fit->x, 0);
				p.right_width = gsl_vector_get(fit->x, 1);
				p.type        = peaks[current_peak].type;

				if(p.getFWHM() > fwhm_threshold) return false;
      }
    
    return true;
  }

	int OptimizePeakDeconvolution::getNumberOfPeaks_(int charge, std::vector<PeakShape>& temp_shapes)
	{
		double dist = dist_/charge;

		OptimizationFunctions::peaks_DC_.clear();
      
		unsigned int shape=0;
#ifdef DEBUG_DECONV
		std::cout<<"temp_shapes[0].mz_position "<<temp_shapes[0].mz_position
						 <<"\t dist "<<dist<<"\tp_index "<<shape<<std::endl;
#endif
		// while the peak's position is smaller than the last considered position
		// take the peak for optimization
		while( (temp_shapes[0].mz_position + shape*dist <
						OptimizationFunctions::positions_DC_[OptimizationFunctions::positions_DC_.size()-1]) &&
					 (shape < temp_shapes.size() ) )
			{
				OptimizationFunctions::peaks_DC_.push_back(temp_shapes[shape]);
#ifdef DEBUG_DECONV
				std::cout<<"temp_shapes[0].mz_position + p_index*dist = "<<temp_shapes[0].mz_position + shape*dist<<std::endl;
#endif
				++shape;
			}
      
		return shape;
      
	}

    
	
  
}
