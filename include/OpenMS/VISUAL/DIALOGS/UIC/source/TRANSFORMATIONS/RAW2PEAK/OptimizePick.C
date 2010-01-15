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
// $Maintainer: Eva Lange $
// $Authors: $
// --------------------------------------------------------------------------

#include <cmath>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/OptimizePick.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>


namespace OpenMS
{
  namespace OptimizationFunctions
  {



    // Print the computed signal
			void printSignal(const gsl_vector* x,void* param, float resolution)
    {

			std::vector<DoubleReal>& positions = static_cast<OptimizePick::Data*> (param) ->positions;
  		std::vector<PeakShape>& peaks = static_cast<OptimizePick::Data*> (param) ->peaks;
      std::cout << "Printing Signal" << std::endl;
      if (resolution == 1.)
      {
        // iterate over all points of the signal
        for (size_t current_point = 0; current_point < positions.size(); current_point++)
        {
          double computed_signal     = 0.;
          double current_position    = positions[current_point];

          // iterate over all peaks
          for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
          {
            // Store the current parameters for this peak
            double p_height 		 = gsl_vector_get(x, 4*current_peak);
            double p_position    = gsl_vector_get(x, 4*current_peak + 3);
            double p_width  		 = (current_position <= p_position) ? gsl_vector_get(x, 4*current_peak + 1)
                                 : gsl_vector_get(x, 4*current_peak + 2);

            // is it a Lorentz or a Sech - Peak?
            if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
            {
              computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
            }
            else // It's a Sech - Peak
            {
              computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
            }
          }
          std::cerr << positions[current_point] << " " << computed_signal << std::endl;
        }
      }
      else
      {
        // Compute step width
        float sw = (positions[1] - positions[0])/resolution;
        for (int i=0; i<positions.size() * resolution; i++)
        {
          double computed_signal     = 0.;
          double current_position    = positions[0] + i*sw;

          // iterate over all peaks
          for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
          {
            // Store the current parameters for this peak
            double p_height 		 = gsl_vector_get(x, 4*current_peak);
            double p_position    = gsl_vector_get(x, 4*current_peak + 3);
            double p_width  		 = (current_position <= p_position) ? gsl_vector_get(x, 4*current_peak + 1)
                                 : gsl_vector_get(x, 4*current_peak + 2);

            // is it a Lorentz or a Sech - Peak?
            if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
            {
              computed_signal += p_height / (1. + pow(p_width * (current_position - p_position), 2));
            }
            else // It's a Sech - Peak
            {
              computed_signal += p_height / pow(cosh(p_width * (current_position - p_position)), 2);
            }
          }

          std::cerr.precision(writtenDigits<DoubleReal>());

          std::cerr << positions[0] + i*sw << " " << computed_signal << std::endl;
        }
      }
    }


    // Evaluation of the target function for nonlinear optimization.
    int residual(const gsl_vector* x, void* params , gsl_vector* f)
    {
      // According to the gsl conventions, x contains the parameters to be optimized.
      // In our case, this means that we store for each peak four consecutive values:
      //  - its height
      //  - its left width
      //  - its right width
      //  - its position
      //
      // Params might contain any additional parameters. We handle these using class members
      // instead.
      // The vector f is supposed to contain the result when we return from this function.
      // Note: GSL wants the values for each data point i as one component of the results vector
				std::vector<DoubleReal>& signal = static_cast<OptimizePick::Data*> (params) ->signal;
				std::vector<DoubleReal>& positions = static_cast<OptimizePick::Data*> (params) ->positions;
				std::vector<PeakShape>& peaks = static_cast<OptimizePick::Data*> (params) ->peaks;
        OptimizationFunctions::PenaltyFactors& penalties=static_cast<OptimizePick::Data*> (params) ->penalties;
      // iterate over all points of the signal
      for (size_t current_point = 0; current_point < positions.size(); current_point++)
      {
        double computed_signal     = 0.;
        double current_position    = positions[current_point];
        double experimental_signal = signal[current_point];

        // iterate over all peaks
        for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
        {
          // Store the current parameters for this peak
          double p_height 		 = gsl_vector_get(x, 4*current_peak);
          double p_position    = gsl_vector_get(x, 4*current_peak + 3);
          double p_width  		 = (current_position <= p_position) ? gsl_vector_get(x, 4*current_peak + 1)
                               : gsl_vector_get(x, 4*current_peak + 2);

          // is it a Lorentz or a Sech - Peak?
          if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
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

      double penalty = 0.;
//      struct PenaltyFactors* penalties = (struct PenaltyFactors *)params;
      double penalty_pos    = penalties.pos;
      double penalty_lwidth = penalties.lWidth;
      double penalty_rwidth = penalties.rWidth;

      // iterate over all peaks again to compute the penalties
      for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
        double old_position = peaks[current_peak].mz_position;
        double old_width_l 	= peaks[current_peak].left_width;
        double old_width_r 	= peaks[current_peak].right_width;
        double p_position   = gsl_vector_get(x, 4*current_peak + 3);
        double p_width_l 		= gsl_vector_get(x, 4*current_peak + 1);
        double p_width_r    = gsl_vector_get(x, 4*current_peak + 2);

        //penalty += pow(p_position - old_position, 2) + pow(p_width_l - old_width_l, 2) + pow(p_width_r - old_width_r, 2);
        penalty += 		penalty_pos    * pow(p_position - old_position, 2)
                     + penalty_lwidth * pow(p_width_l - old_width_l, 2)
                     + penalty_rwidth * pow(p_width_r - old_width_r, 2);
      }

      gsl_vector_set(f, positions.size(), 100*penalty);

      return GSL_SUCCESS;
    }

    /** Compute the Jacobian of the residual, where each row of the matrix corresponds to a
     *  point in the data.
     */
    int jacobian(const gsl_vector* x, void* params, gsl_matrix* J)
    {
      // For the conventions on x and params c.f. the commentary in residual()
      //
      // The matrix J is supposed to contain the result when we return from this function.
      // Note: GSL expects the Jacobian as follows:
      // 	- each row corresponds to one data point
      // 	- each column corresponds to one parameter
			//	std::vector<DoubleReal>& signal = static_cast<OptimizePick::Data*> (params) ->signal;
				std::vector<DoubleReal>& positions = static_cast<OptimizePick::Data*> (params) ->positions;
				std::vector<PeakShape>& peaks = static_cast<OptimizePick::Data*> (params) ->peaks;
        OptimizationFunctions::PenaltyFactors& penalties=static_cast<OptimizePick::Data*> (params) ->penalties;
      // iterate over all points of the signal
      for (size_t current_point = 0; current_point < positions.size(); current_point++)
      {
        double current_position    = positions[current_point];

        // iterate over all peaks
        for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
        {
          // Store the current parameters for this peak
          double p_height 		 = gsl_vector_get(x, 4*current_peak);
          double p_position    = gsl_vector_get(x, 4*current_peak + 3);
          double p_width  		 = (current_position <= p_position) ? gsl_vector_get(x, 4*current_peak + 1)
                               : gsl_vector_get(x, 4*current_peak + 2);

          // is it a Lorentz or a Sech - Peak?
          if (peaks[current_peak].type == PeakShape::LORENTZ_PEAK)
          {
            double diff      = current_position - p_position;
            double denom_inv = 1./(1. + pow(p_width * diff, 2));

            double ddl_left  = (current_position <= p_position)
                               ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2)
                               : 0;

            double ddl_right = (current_position  > p_position)
                               ? -2 * p_height * pow(diff, 2) * p_width * pow(denom_inv, 2)
                               : 0;

            double ddx0			 = -2*p_height*pow(p_width,2)*diff*pow(denom_inv, 2);

            gsl_matrix_set(J, current_point, 4*current_peak,     denom_inv);
            gsl_matrix_set(J, current_point, 4*current_peak + 1, ddl_left);
            gsl_matrix_set(J, current_point, 4*current_peak + 2, ddl_right);
            gsl_matrix_set(J, current_point, 4*current_peak + 3, ddx0     );
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
            double ddx0      = 2*p_height*p_width*sinh_term * pow(denom_inv,3);

            gsl_matrix_set(J, current_point, 4*current_peak, pow(denom_inv, 2));
            gsl_matrix_set(J, current_point, 4*current_peak + 1, ddl_left         );
            gsl_matrix_set(J, current_point, 4*current_peak + 2, ddl_right        );
            gsl_matrix_set(J, current_point, 4*current_peak + 3, ddx0             );
          }
        }
      }

      // Now iterate over all peaks again to compute the penalties.
//   struct PenaltyFactors* penalties = (struct PenaltyFactors *)params;
      for (size_t current_peak = 0; current_peak < peaks.size(); current_peak++)
      {
        double p_width_left = gsl_vector_get(x, 4*current_peak + 1);
        double p_width_right = gsl_vector_get(x, 4*current_peak + 2);
        double p_position   = gsl_vector_get(x, 4*current_peak + 3);

        double old_width_left  = peaks[current_peak].left_width;
        double old_width_right = peaks[current_peak].right_width;
        double old_position    = peaks[current_peak].mz_position;


        double penalty_l = 2.*penalties.lWidth*(p_width_left - old_width_left);
        double penalty_r = 2.*penalties.rWidth*(p_width_right - old_width_right);
        double penalty_p=0;
        if (fabs(p_position - old_position) < 0.2)
        {
          penalty_p = 2.*penalties.pos*(p_position - old_position);
        }

        gsl_matrix_set(J, positions.size(), 4*current_peak,   0.);
        gsl_matrix_set(J, positions.size(), 4*current_peak+1, 100*penalty_l);
        gsl_matrix_set(J, positions.size(), 4*current_peak+2, 100*penalty_r);
        gsl_matrix_set(J, positions.size(), 4*current_peak+3, 100*penalty_p);
      }

      return GSL_SUCCESS;
    }

    // Driver function for the evaluation of function and jacobian.
    int evaluate(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      residual(x, params, f);
      jacobian(x, params, J);

      return GSL_SUCCESS;
    }
  }


  OptimizePick::OptimizePick(const struct OptimizationFunctions::PenaltyFactors& penalties,
                             const int max_iteration,
                             const double eps_abs,
                             const double eps_rel )
  {

    penalties_=penalties;

    max_iteration_=max_iteration;
    eps_abs_=eps_abs;
    eps_rel_=eps_rel;

#ifdef DEBUG_PEAK_PICKING
    std::cout << "max iteration " << max_iteration_
    << "\n eps abs " << eps_abs
    << "\n eps rel " << eps_rel
    << "\n penalty factor pos " << penalties.pos
    << "\n penalty factor left width " << penalties.lWidth
    << "\n penalty factor right width " << penalties.rWidth
    << std::endl;
#endif

  }

  OptimizePick::~OptimizePick()
  {}

void OptimizePick::optimize(std::vector<PeakShape>& peaks,Data& data)
  {
    if (peaks.size() == 0)
      return;

    size_t global_peak_number=0;
    data.peaks.assign(peaks.begin(),peaks.end());

    gsl_vector *start_value=gsl_vector_alloc(4*data.peaks.size());
    // We have to initialize the parameters for the optimization
    for (size_t i = 0; i < data.peaks.size(); i++)
    {
      PeakShape current_peak = data.peaks[i];
      double h  = current_peak.height;
      double wl = current_peak.left_width;
      double wr = current_peak.right_width;
      double p  = current_peak.mz_position;
      if (boost::math::isnan(wl))
      {
        data.peaks[i].left_width = 1;
        wl = 1.;
      }
      if (boost::math::isnan(wr))
      {
        data.peaks[i].right_width = 1;
        wr = 1.;
      }

      gsl_vector_set(start_value, 4*i  , h);
      gsl_vector_set(start_value, 4*i+1, wl);
      gsl_vector_set(start_value, 4*i+2, wr);
      gsl_vector_set(start_value, 4*i+3, p);
    }


    // The gsl algorithms require us to provide function pointers for the evaluation of
    // the target function.
    gsl_multifit_function_fdf fit_function;
    fit_function.f    	= OptimizationFunctions::residual;
    fit_function.df  	  = OptimizationFunctions::jacobian;
    fit_function.fdf		= OptimizationFunctions::evaluate;
    fit_function.n			= std::max(data.positions.size()+1,4*data.peaks.size());
    fit_function.p			= 4*data.peaks.size();
//    fit_function.params = &penalties_;
		data.penalties = penalties_;
    fit_function.params = &data;

    const gsl_multifit_fdfsolver_type *type = gsl_multifit_fdfsolver_lmsder;

    gsl_multifit_fdfsolver *fit=gsl_multifit_fdfsolver_alloc(type, std::max(data.positions.size()+1, 4*data.peaks.size()), 4*data.peaks.size());

    gsl_multifit_fdfsolver_set(fit, &fit_function, start_value);

    // initial norm
    // std::cout << "Before optimization: ||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;

    // Iteration
    unsigned int iteration = 0;
    int status;

    do
    {
      iteration++;
      status = gsl_multifit_fdfsolver_iterate(fit);
#ifdef DEBUG_PEAK_PICKING
      std::cout << "Iteration " << iteration << "; Status " << gsl_strerror(status) << "; " << std::endl;
      std::cout << "||f|| = " << gsl_blas_dnrm2(fit->f) << std::endl;
      std::cout << "Number of parms: " << data.peaks.size() * 4 << std::endl;
      std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
#endif
      if (boost::math::isnan(gsl_blas_dnrm2(fit->dx)))
        break;

      // We use the gsl function gsl_multifit_test_delta to decide if we can finish the iteration.
      // We only finish if all new parameters deviates only by a small amount from the parameters of the last iteration
      status = gsl_multifit_test_delta(fit->dx, fit->x, eps_abs_, eps_rel_);
      if (status != GSL_CONTINUE)
        break;

    }
    while (status == GSL_CONTINUE && iteration < max_iteration_);

#ifdef DEBUG_PEAK_PICKING
    std::cout << "Finished!" << std::endl;
    std::cout << "Delta: " << gsl_blas_dnrm2(fit->dx) << std::endl;
    double chi = gsl_blas_dnrm2(fit->f);
    std::cout << "chisq/dof = " << pow(chi, 2.0) / (data.positions.size() - 4*data.peaks.size());
#endif

    //		OptimizationFunctions::printSignal(fit->x, 5.,param);

    // iterate over all peaks and store the optimized values in peaks
   for (size_t current_peak = 0; current_peak <data.peaks.size(); current_peak++)
    {
      // Store the current parameters for this peak
      peaks[global_peak_number+current_peak].height 		 = gsl_vector_get(fit->x, 4*current_peak);
      peaks[global_peak_number+current_peak].mz_position = gsl_vector_get(fit->x, 4*current_peak + 3);
      peaks[global_peak_number+current_peak].left_width  = gsl_vector_get(fit->x, 4*current_peak + 1);
      peaks[global_peak_number+current_peak].right_width = gsl_vector_get(fit->x, 4*current_peak + 2);

      //	compute the area
      // is it a Lorentz or a Sech - Peak?
      if (peaks[global_peak_number+current_peak].type == PeakShape::LORENTZ_PEAK)
      {
        PeakShape p = peaks[global_peak_number+current_peak];
        double x_left_endpoint=p.mz_position+1/p.left_width*sqrt(p.height/1-1);
        double x_rigth_endpoint=p.mz_position+1/p.right_width*sqrt(p.height/1-1);
        double area_left=-p.height/p.left_width*atan(p.left_width*(x_left_endpoint-p.mz_position));
        double area_right=-p.height/p.right_width*atan(p.right_width*(p.mz_position-x_rigth_endpoint));
        peaks[global_peak_number+current_peak].area=area_left+area_right;
      }
      else  //It's a Sech - Peak
      {
        PeakShape p = peaks[global_peak_number+current_peak];
				double x_left_endpoint=p.mz_position+1/p.left_width* boost::math::acosh(sqrt(p.height/0.001));
        double x_rigth_endpoint=p.mz_position+1/p.right_width* boost::math::acosh(sqrt(p.height/0.001));
        double area_left=-p.height/p.left_width*(sinh(p.left_width*(p.mz_position-x_left_endpoint))/cosh(p.left_width*(p.mz_position-x_left_endpoint)));
        double area_right=-p.height/p.right_width*(sinh(p.right_width*(p.mz_position-x_rigth_endpoint))/cosh(p.right_width*(p.mz_position-x_rigth_endpoint)));
        peaks[global_peak_number+current_peak].area=area_left+area_right;
      }
    }
    global_peak_number += data.peaks.size();

    gsl_multifit_fdfsolver_free(fit);
    gsl_vector_free(start_value);
  }


  // double OptimizePick::correlate_(const PeakShape& peak,
//                                  double left_endpoint,
// 																	double right_endpoint,Data& data)
//   {
//     double SSxx = 0., SSyy = 0., SSxy = 0.;

//     // compute the averages
//     double data_average=0., fit_average=0.;
//     double data_sqr=0., fit_sqr=0.;
//     double cross=0.;

//     int number_of_points = 0;

//     int first=0;
//     int last=data.positions.size()-1;

//     // search for the left endpoint position
//     while (data.positions[first] < left_endpoint)
//     {
//       ++first;
//     }

//     // search for the right endpoint position
//     while (data.positions[last] > right_endpoint)
//     {
//       --last;
//     }


//     // for separate overlapping peak correlate until the max position...
//     for (int i=first; i <= last; i++)
//     {
//       double data_val = data.signal[i];
//       double peak_val = peak(data.positions[i]);

//       data_average += data_val;
//       fit_average  += peak_val;

//       data_sqr += data_val * data_val;
//       fit_sqr  += peak_val * peak_val;

//       cross += data_val * peak_val;

//       number_of_points++;
//     }

//     if (number_of_points == 0)
//       return 0.;

//     data_average /= number_of_points;
//     fit_average  /= number_of_points;

//     SSxx = data_sqr - number_of_points * (data_average * data_average);
//     SSyy = fit_sqr - number_of_points * (fit_average * fit_average);
//     SSxy = cross - number_of_points * (data_average * fit_average);

//     return (SSxy * SSxy) / (SSxx * SSyy);
//   }

}

