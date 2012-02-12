// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>

#include <algorithm>
#include <cstdlib>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>


namespace OpenMS
{
  LowessSmoothing::LowessSmoothing()
    : DefaultParamHandler("LowessSmoothing")
  {
    defaults_.setValue("window_size", 10, "The number of peaks to be included for local fitting in one window.");
    defaultsToParam_();
  }

  LowessSmoothing::~LowessSmoothing()
  {
  }

  void LowessSmoothing::smoothData(const DoubleVector& input_x, const DoubleVector& input_y, DoubleVector& smoothed_output)
  {
    if (input_x.size() != input_y.size()) {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Sizes of x and y values not equal! Aborting... ", String(input_x.size()));
    }

    // unable to smooth over 2 or less data points (we need at least 3)
    if (input_x.size()<=2)
    {
      smoothed_output = input_y;
      return;
    }

    Size input_size = input_y.size();
    // alpha_ = 1 / ( input_size / window_size_ );

    // const Size q = floor( input_size * alpha );
    const Size q = (window_size_ < input_size) ? window_size_ : input_size-1;

    DoubleVector distances(input_size, 0.0);
    DoubleVector sortedDistances(input_size, 0.0);

    // DoubleVector smooth_yvals_Vec(input_size);

    gsl_matrix *X = gsl_matrix_alloc(input_size, 3);
    gsl_matrix *cov = gsl_matrix_alloc(3, 3);

    gsl_vector *weights = gsl_vector_alloc(input_size);
    gsl_vector *c = gsl_vector_alloc(3);
    gsl_vector *x = gsl_vector_alloc(3);
    gsl_vector *yvals_ = gsl_vector_alloc(input_size);

    DoubleReal y, yErr, chisq;

    // Setup the model matrix X for a quadratic fit.
    // yvals_ = new double[input_size];
    // Size idx = 0;

    for(Size p_idx = 0; p_idx < input_y.size(); ++p_idx)
    {
      DoubleReal rt = input_x[p_idx];

      gsl_matrix_set(X, p_idx, 0, 1.0);
      gsl_matrix_set(X, p_idx, 1, rt);
      gsl_matrix_set(X, p_idx, 2, rt * rt);

      gsl_vector_set(yvals_, p_idx, input_y[p_idx]);

      // ++idx;
    }



    //for(DoubleVector::const_iterator outer_peak_it = input_y.begin(); outer_peak_it != input_y.end(); ++outer_peak_it )
    for(Size outer_idx = 0; outer_idx < input_y.size(); ++outer_idx)
    {
      // Compute distances.
      // Size inner_idx = 0;
      for(Size inner_idx = 0; inner_idx < input_y.size(); ++inner_idx)
      {
        distances[inner_idx] = std::fabs(input_x[outer_idx] - input_x[inner_idx]);
        sortedDistances[inner_idx] = distances[inner_idx];
        // ++inner_idx;
      }

      // Sort distances in order from smallest to largest.
      // std::sort(sortedDistances, sortedDistances + input_size);
      std::sort(sortedDistances.begin(), sortedDistances.end());


      // Compute weights.
      for (Size inner_idx = 0; inner_idx < input_size; ++inner_idx)
      {
        gsl_vector_set(weights, inner_idx, tricube_(distances[inner_idx], sortedDistances[q]));
      }


      gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(input_size, 3);
      gsl_multifit_wlinear( X, weights, yvals_, c, cov, &chisq, work);
      gsl_multifit_linear_free(work);


      DoubleReal rt = input_x[outer_idx];
      gsl_vector_set(x, 0, 1.0);
      gsl_vector_set(x, 1, rt);
      gsl_vector_set(x, 2, rt * rt);

      gsl_multifit_linear_est (x, c, cov, &y, &yErr);

      smoothed_output.push_back(y);
    }

    gsl_matrix_free(X);
    gsl_vector_free(weights);
    gsl_vector_free(c);
    gsl_matrix_free(cov);

    return ;
  }

  DoubleReal LowessSmoothing::tricube_( DoubleReal u, DoubleReal t )
  {
    // In our case, u represents a distance and hence should be strictly positive.
    if (u < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Value of u must be strictly positive! Aborting...", String(u));
    }

    // 0 <= u < t; u is regarded as 0.0 if fabs(u) falls below epsilon
    if ( (fabs(u) < std::numeric_limits<double>::epsilon() || (0.0 < u)) && (u < t) )
    {
      // (1 - (u/t)^3)^3
      return pow( ( 1.0 - pow(u/t, 3.0)), 3.0 );
    }
    // u >= t
    else
    {
      return 0.0;
    }
  }

  void LowessSmoothing::updateMembers_()
  {
    window_size_ = (Size)param_.getValue("window_size");
  }

} //namespace OpenMS
