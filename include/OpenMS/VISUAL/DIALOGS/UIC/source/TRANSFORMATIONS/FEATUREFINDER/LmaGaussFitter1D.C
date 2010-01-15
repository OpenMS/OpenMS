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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LmaGaussFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{
    LmaGaussFitter1D::LmaGaussFitter1D()
    : LevMarqFitter1D()
    {
      setName(getProductName());
      defaults_.setValue("statistics:variance",1.0,"Variance of the model.", StringList::create("advanced"));
      defaultsToParam_();
    }

    LmaGaussFitter1D::LmaGaussFitter1D(const LmaGaussFitter1D& source)
    : LevMarqFitter1D(source)
    {
        setParameters( source.getParameters() );
        updateMembers_();
    }

    LmaGaussFitter1D::~LmaGaussFitter1D()
    {
    }

    LmaGaussFitter1D& LmaGaussFitter1D::operator = (const LmaGaussFitter1D& source)
    {
        if (&source == this) return *this;

        LevMarqFitter1D::operator = (source);
        setParameters( source.getParameters() );
        updateMembers_();

        return *this;
    }

    Int LmaGaussFitter1D::residual_(const gsl_vector* x, void* params, gsl_vector* f)
    {
        Size n = static_cast<LmaGaussFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LmaGaussFitter1D::Data*> (params) ->set;

        CoordinateType normal_s = gsl_vector_get( x, 0 );
        CoordinateType normal_m = gsl_vector_get( x, 1 );
        CoordinateType normal_scale = gsl_vector_get( x, 2 );

        CoordinateType Yi = 0.0;

        for ( Size i = 0; i < n; i++ )
        {
          CoordinateType t = set[i].getPos();

          Yi = ( 1 / ( sqrt( 2 * Constants::PI ) * normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

          gsl_vector_set( f, i, ( Yi - set[i].getIntensity() ) );
        }

      return GSL_SUCCESS;
    }

    Int LmaGaussFitter1D::jacobian_(const gsl_vector* x, void* params, gsl_matrix* J)
    {
        Size n = static_cast<LmaGaussFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LmaGaussFitter1D::Data*> (params) ->set;

        CoordinateType normal_s = gsl_vector_get( x, 0 );
        CoordinateType normal_m = gsl_vector_get( x, 1 );
        CoordinateType normal_scale = gsl_vector_get( x, 2 );

        CoordinateType derivative_normal_s, derivative_normal_m, derivative_normal_scale = 0.0;

        for ( Size i = 0; i < n; i++ )
        {
          CoordinateType t = set[i].getPos();

          // f'(normal_s)
          derivative_normal_s = -( ( 1 / sqrt( 2 * Constants::PI ) ) / ( normal_s * normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale + ( ( 1 / sqrt( 2 * Constants::PI ) ) / ( normal_s * normal_s * normal_s * normal_s ) ) * ( ( t - normal_m ) * ( t - normal_m ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

	  			// f'(normal_m)
          derivative_normal_m = ( ( 1 / sqrt( 2 * Constants::PI ) ) / ( normal_s * normal_s * normal_s ) ) * ( t - normal_m ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) ) * normal_scale;

          // f'(normal_scale)
          derivative_normal_scale = ( ( 1 / sqrt( 2 * Constants::PI ) ) / ( normal_s ) ) * exp( -( ( t - normal_m ) * ( t - normal_m ) ) / ( 2 * normal_s * normal_s ) );

          // set the jacobian matrix of the normal distribution
          gsl_matrix_set( J, i, 0, derivative_normal_s );
          gsl_matrix_set( J, i, 1, derivative_normal_m );
          gsl_matrix_set( J, i, 2, derivative_normal_scale );
        }

        return GSL_SUCCESS;
    }

    Int LmaGaussFitter1D::evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      LmaGaussFitter1D::residual_( x, params, f );
      LmaGaussFitter1D::jacobian_( x, params, J );

      return GSL_SUCCESS;
    }

    void LmaGaussFitter1D::printState_(Int iter, gsl_multifit_fdfsolver * s)
    {
      printf ( "in loop iter: %4u x = % 15.8f % 15.8f % 15.8f |f(x)| = %g\n", iter,
               gsl_vector_get( s->x, 0 ),
               gsl_vector_get( s->x, 1 ),
               gsl_vector_get( s->x, 2 ),
               gsl_blas_dnrm2( s->f ) );
    }

    LmaGaussFitter1D::QualityType LmaGaussFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
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
        LmaGaussFitter1D::Data d;
        d.n= set.size();
        d.set = set;

         // Compute start parameter
        setInitialParameters_(set);

        // Optimize parameter with Levenberg-Marquardt algorithm (GLS)
        CoordinateType x_init[ 3 ] = { standard_deviation_, expected_value_, scale_factor_ };
        if ( symmetric_ == false )
        {
          optimize_(set, 3, x_init, &(residual_), &(jacobian_), &(evaluate_), &d);
        }

        // Set optimized parameter
        standard_deviation_ = x_init[0];
        expected_value_ = x_init[1];
        scale_factor_ = x_init[2];

#ifdef DEBUG_FEATUREFINDER
        if ( getGslStatus_() != "success" )
        {
          std::cout << "status: " << getGslStatus_() << std::endl;
        }
#endif
        // build model
        model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("LmaGaussModel"));
	      model->setInterpolationStep( interpolation_step_ );

        Param tmp;
        tmp.setValue( "bounding_box:min", min_ );
        tmp.setValue( "bounding_box:max", max_ );
        tmp.setValue( "statistics:variance", statistics_.variance() );
        tmp.setValue( "statistics:mean", statistics_.mean() );
        tmp.setValue( "lma:scale_factor", scale_factor_ );
        tmp.setValue( "lma:standard_deviation", standard_deviation_ );
        tmp.setValue( "lma:expected_value", expected_value_ );
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

    void LmaGaussFitter1D::setInitialParameters_(const RawDataArrayType& set)
    {
      // sum over all intensities
      CoordinateType sum = 0.0;
      for (Size i=0; i<set.size(); ++i) sum += set[i].getIntensity();

      // calculate the median
      SignedSize median = 0;
      Real count = 0.0;
      for ( Size i = 0; i < set.size(); ++i )
      {
        count += set[i].getIntensity();
        if ( count <= sum * 0.5 ) median = i;
      }

      CoordinateType sumS = 0.0;
      for ( Size i = 0; i < set.size(); ++i )
      {
        sumS += pow( ( set[i].getPos() - set[median].getPos() ), 2 );
      }

      // calculate the stardard deviation
      standard_deviation_ = sqrt( sumS / ( set.size() - 1 ) );

      // set expeceted value
      expected_value_ = set[median].getPos();

      // scaling factor of the peak
      scale_factor_ = set[median].getIntensity();

      // default is an asymmetric peak
      symmetric_ = false;
    }

    void LmaGaussFitter1D::updateMembers_()
    {
        LevMarqFitter1D::updateMembers_();
        statistics_.setVariance(param_.getValue("statistics:variance"));
    }

}
