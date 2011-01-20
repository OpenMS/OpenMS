// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/SIMULATION/EGHFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <boost/math/special_functions/fpclassify.hpp>

#define DEBUG_EGHFITTER

namespace OpenMS
{
  EGHFitter1D::EGHFitter1D()
  : LevMarqFitter1D()
  {
    setName(getProductName());
    defaults_.setValue("statistics:variance",1.0,"Variance of the model.", StringList::create("advanced"));
    defaultsToParam_();
  }

  EGHFitter1D::EGHFitter1D(const EGHFitter1D& source)
  : LevMarqFitter1D(source)
  {
    setParameters( source.getParameters() );
    updateMembers_();
  }

  EGHFitter1D::~EGHFitter1D()
  {
  }

  EGHFitter1D& EGHFitter1D::operator = (const EGHFitter1D& source)
  {
    if (&source == this) return *this;

    LevMarqFitter1D::operator = (source);
    setParameters( source.getParameters() );
    updateMembers_();

    return *this;
  }

  Int EGHFitter1D::residual_(const gsl_vector* x, void* params, gsl_vector* f)
  {
    Size n = static_cast<EGHFitter1D::Data*> (params) ->n;
    RawDataArrayType set = static_cast<EGHFitter1D::Data*> (params) ->set;

    CoordinateType H  = gsl_vector_get( x, 0 );
    CoordinateType tR = gsl_vector_get( x, 1 );
    CoordinateType sigma_square = gsl_vector_get( x, 2 );
    CoordinateType tau = gsl_vector_get( x, 3 );

    CoordinateType t_diff,t_diff2,denominator = 0.0;

    CoordinateType fegh = 0.0;

    // iterate over all points of the signal
    for ( Size i = 0; i < n; i++ )
    {
      DoubleReal t = set[i].getPos();

      t_diff = t - tR;
      t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

      denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

      if(denominator > 0.0)
      {
        fegh = H * exp(- t_diff2 / denominator);
      }
      else
      {
        fegh = 0.0;
      }

      gsl_vector_set( f, i, ( fegh - set[i].getIntensity() ) );
    }

    return GSL_SUCCESS;
  }

  Int EGHFitter1D::jacobian_(const gsl_vector* x, void* params, gsl_matrix* J)
  {
    Size n =  static_cast<EGHFitter1D::Data*> (params) ->n;
    RawDataArrayType set = static_cast<EGHFitter1D::Data*> (params) ->set;

    CoordinateType H  = gsl_vector_get( x, 0 );
    CoordinateType tR = gsl_vector_get( x, 1 );
    CoordinateType sigma_square = gsl_vector_get( x, 2 );
    CoordinateType tau = gsl_vector_get( x, 3 );

    CoordinateType derivative_H, derivative_tR, derivative_sigma_square, derivative_tau = 0.0;
    CoordinateType t_diff,t_diff2,exp1,denominator = 0.0;


    // iterate over all points of the signal
    for ( Size i = 0; i < n; i++ )
    {
      CoordinateType t = set[i].getPos();

      t_diff = t - tR;
      t_diff2 = t_diff * t_diff; // -> (t - t_R)^2

      denominator = 2 * sigma_square + tau * t_diff; // -> 2\sigma_{g}^{2} + \tau \left(t - t_R\right)

      if( denominator > 0)
      {
        exp1 = exp(- t_diff2 / denominator);

        // \partial H f_{egh}(t) = \exp\left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right)
        derivative_H = exp1;

        // \partial t_R f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{\left( 4 \sigma_{g}^{2} + \tau \left(t-t_R \right) \right) \left(t-t_R \right)}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
        derivative_tR = H * exp1 * ( ((4 * sigma_square + tau * t_diff) * t_diff) / (denominator * denominator));

        // \partial \sigma_{g}^{2} f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ 2 \left(t - t_R\right)^2}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
        derivative_sigma_square = H * exp1 * ((2 * t_diff2) / (denominator * denominator));

        // \partial \tau f_{egh}(t) &=& H \exp \left( \frac{-\left(t-t_R \right)^2}{2\sigma_{g}^{2} + \tau \left(t - t_R\right)} \right) \left( \frac{ \left(t - t_R\right)^3}{\left( 2\sigma_{g}^{2} + \tau \left(t - t_R\right) \right)^2} \right)
        derivative_tau = H * exp1 * ((t_diff * t_diff2) / (denominator * denominator));
      }
      else
      {
        derivative_H = 0.0;
        derivative_tR = 0.0;
        derivative_sigma_square = 0.0;
        derivative_tau = 0.0;
      }

      // set the jacobian matrix
      gsl_matrix_set( J, i, 0, derivative_H );
      gsl_matrix_set( J, i, 1, derivative_tR );
      gsl_matrix_set( J, i, 2, derivative_sigma_square );
      gsl_matrix_set( J, i, 3, derivative_tau );
    }

    return GSL_SUCCESS;
  }

  Int EGHFitter1D::evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
  {
    EGHFitter1D::residual_( x, params, f );
    EGHFitter1D::jacobian_( x, params, J );

    return GSL_SUCCESS;
  }

  void EGHFitter1D::printState_(Int iter, gsl_multifit_fdfsolver * s)
  {
    printf ( "iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
        gsl_vector_get( s->x, 0 ),
        gsl_vector_get( s->x, 1 ),
        gsl_vector_get( s->x, 2 ),
        gsl_vector_get( s->x, 3 ),
        gsl_blas_dnrm2( s->f ) );
  }

  EGHFitter1D::QualityType EGHFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
  {
    // Calculate bounding box
    min_ = max_ = set[0].getPos();
    for ( Size pos=1; pos < set.size(); ++pos)
    {
      CoordinateType tmp = set[pos].getPos();
      if ( min_ > tmp ) min_ = tmp;
      if ( max_ < tmp ) max_ = tmp;
    }

    // Enlarge the bounding box by a few multiples of the standard deviation
    {
      stdev1_ = sqrt ( statistics_.variance() ) * tolerance_stdev_box_;
      min_ -= stdev1_;
      max_ += stdev1_;
    }

    // Set advanced parameters for residual_  und jacobian_ method
    EGHFitter1D::Data d;
    d.n = set.size();
    d.set = set;

    // Compute start parameters
    setInitialParameters_(set);

    // Optimize parameter with Levenberg-Marquardt algorithm (GLS)
    CoordinateType x_init[ 4 ] = { height_, retention_, sigma_square_, tau_ };
    optimize_(set, 4, x_init, &(residual_), &(jacobian_), &(evaluate_), &d);

    // Set optimized parameters
    height_ = x_init[0];
    retention_ = x_init[1];
    sigma_square_ = x_init[2];
    tau_ = x_init[3];

#ifdef DEBUG_EGHFITTER
    LOG_DEBUG << "Fitter returned \n";
    LOG_DEBUG << "height:       " << height_ << "\n";
    LOG_DEBUG << "retention:    " << retention_ << "\n";
    LOG_DEBUG << "sigma_square: " << sigma_square_ << "\n";
    LOG_DEBUG << "tau:          " << tau_ << std::endl;

    if ( getGslStatus_() != "success" )
    {
      std::cout << "status: " << getGslStatus_() << std::endl;
    }
#endif

    // build model
    model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("EGHModel"));
    model->setInterpolationStep( interpolation_step_ );

    Param tmp;
    tmp.setValue( "statistics:variance", statistics_.variance() );
    tmp.setValue( "statistics:mean", statistics_.mean() );

    tmp.setValue( "bounding_box:compute", "false" ); // disable auto computation of bounding box
    tmp.setValue( "bounding_box:min", min_ );
    tmp.setValue( "bounding_box:max", max_ );

    tmp.setValue( "egh:height", height_ );
    tmp.setValue( "egh:retention", retention_);

    tmp.setValue("egh:guess_parameter", "false"); // disable guessing of parameters from A/B
    tmp.setValue("egh:tau", tau_);
    tmp.setValue("egh:sigma_square", sigma_square_);

    model->setParameters( tmp );


    // calculate pearson correlation
    std::vector<Real> real_data;
    real_data.reserve(set.size());
    std::vector<Real> model_data;
    model_data.reserve(set.size());

    for (Size i=0; i < set.size(); ++i)
    {
      real_data.push_back(set[i].getIntensity());
      model_data.push_back( model->getIntensity( DPosition<1>(set[i].getPosition()) ) );
    }

    QualityType correlation = Math::pearsonCorrelationCoefficient(real_data.begin(), real_data.end(), model_data.begin(), model_data.end());
    if (boost::math::isnan(correlation)) correlation = -1.0;

    return correlation;
  }

  void EGHFitter1D::setInitialParameters_(const RawDataArrayType& set)
  {
    // sum over all intensities
    CoordinateType sum = 0.0;
    for (Size i=0; i<set.size(); ++i) sum += set[i].getIntensity();

    // calculate the median
    //Size median = 0;
    //Real count = 0.0;
    Size apex_rt = 0;
    CoordinateType apex = 0.0;
    for ( Size i = 0; i < set.size(); ++i )
    {
      //count += set[i].getIntensity();
      //if ( count <= sum / 2 ) median = i;

      if (set[i].getIntensity() > apex)
      {
        apex = set[i].getIntensity();
        apex_rt = i;
      }

    }

    // calculate the height of the peak
    height_ = set[apex_rt].getIntensity();

    // calculate retention time
    retention_ = set[apex_rt].getPos();


    // guess A / B for alpha = 0.5 -> left/right half max distance

    Size i = apex_rt;
    while(i > 0)
    {
      if(set[i].getIntensity() / height_ < 0.5) break;
      else --i;
    }
    CoordinateType A = retention_ - set[i+1].getPos();

    i = apex_rt;
    while(i < set.size())
    {
      if(set[i].getIntensity() / height_ < 0.5) break;
      else ++i;
    }
    CoordinateType B = set[i-1].getPos() - retention_;

    // compute estimates for tau / sigma_square based on A/B
    CoordinateType log_alpha = log(0.5);

    tau_ = ( -1 / log_alpha ) * (B - A);
    sigma_square_ = ( -1 / (2 * log_alpha) ) * (B * A);

#ifdef DEBUG_EGHFITTER
    LOG_DEBUG << "Initial parameters\n";
    LOG_DEBUG << "height:       " << height_ << "\n";
    LOG_DEBUG << "retention:    " << retention_ << "\n";
    LOG_DEBUG << "A:            " << A << "\n";
    LOG_DEBUG << "B:            " << B << "\n";
    LOG_DEBUG << "sigma_square: " << sigma_square_ << "\n";
    LOG_DEBUG << "tau:          " << tau_ << std::endl;
#endif
  }

  void EGHFitter1D::updateMembers_()
  {
    LevMarqFitter1D::updateMembers_();
    statistics_.setVariance(param_.getValue("statistics:variance"));
  }

}
