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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <boost/math/special_functions/fpclassify.hpp>

namespace OpenMS
{
    EmgFitter1D::EmgFitter1D()
    : LevMarqFitter1D()
    {
      setName(getProductName());
      defaults_.setValue("statistics:variance",1.0,"Variance of the model.", StringList::create("advanced"));
      defaultsToParam_();
    }

    EmgFitter1D::EmgFitter1D(const EmgFitter1D& source)
    : LevMarqFitter1D(source)
    {
        setParameters( source.getParameters() );
        updateMembers_();
    }

    EmgFitter1D::~EmgFitter1D()
    {
    }

    EmgFitter1D& EmgFitter1D::operator = (const EmgFitter1D& source)
    {
        if (&source == this) return *this;

        LevMarqFitter1D::operator = (source);
        setParameters( source.getParameters() );
        updateMembers_();

        return *this;
    }

    Int EmgFitter1D::residual_(const gsl_vector* x, void* params, gsl_vector* f)
    {
        Size n = static_cast<EmgFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<EmgFitter1D::Data*> (params) ->set;

        CoordinateType h = gsl_vector_get( x, 0 );
        CoordinateType w = gsl_vector_get( x, 1 );
        CoordinateType s = gsl_vector_get( x, 2 );
        CoordinateType z = gsl_vector_get( x, 3 );

        CoordinateType Yi = 0.0;

        // iterate over all points of the signal
        for ( Size i = 0; i < n; i++ )
        {
            DoubleReal t = set[i].getPos();

            // Simplified EMG
            Yi = ( h * w / s ) * sqrt( 2.0 * Constants::PI ) * exp( ( pow( w, 2 ) / ( 2 * pow( s, 2 ) ) ) - ( ( t - z ) / s ) ) / ( 1 + exp( ( -2.4055 / sqrt( 2.0 ) ) * ( ( ( t - z ) / w ) - w / s ) ) );

            gsl_vector_set( f, i, ( Yi - set[i].getIntensity() ) );
        }

      return GSL_SUCCESS;
    }

    Int EmgFitter1D::jacobian_(const gsl_vector* x, void* params, gsl_matrix* J)
    {
        Size n =  static_cast<EmgFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<EmgFitter1D::Data*> (params) ->set;

        CoordinateType h = gsl_vector_get( x, 0 );
        CoordinateType w = gsl_vector_get( x, 1 );
        CoordinateType s = gsl_vector_get( x, 2 );
        CoordinateType z = gsl_vector_get( x, 3 );

        const CoordinateType emg_const = 2.4055;
        const CoordinateType sqrt_2pi = sqrt( 2 * Constants::PI );
        const CoordinateType sqrt_2 = sqrt( 2.0 );

        CoordinateType exp1, exp2, exp3 = 0.0;
        CoordinateType derivative_height, derivative_width, derivative_symmetry, derivative_retention = 0.0;

        // iterate over all points of the signal
        for ( Size i = 0; i < n; i++ )
        {
            CoordinateType t = set[i].getPos();

            exp1 = exp( ( ( w * w ) / ( 2 * s * s ) ) - ( ( t - z ) / s ) );
            exp2 = ( 1 + exp( ( -emg_const / sqrt_2 ) * ( ( ( t - z ) / w ) - w / s ) ) );
            exp3 = exp( ( -emg_const / sqrt_2 ) * ( ( ( t - z ) / w ) - w / s ) );

            // f'(h)
            derivative_height = w / s * sqrt_2pi * exp1 / exp2;

            // f'(h)
            derivative_width = h / s * sqrt_2pi * exp1 / exp2 + ( h * w * w ) / ( s * s * s ) * sqrt_2pi * exp1 / exp2 + ( emg_const * h * w ) / s * sqrt_2pi * exp1 * ( -( t - z ) / ( w * w ) - 1 / s ) * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

            // f'(s)
            derivative_symmetry = - h * w / ( s * s ) * sqrt_2pi * exp1 / exp2 + h * w / s * sqrt_2pi * ( -( w * w ) / ( s * s * s ) + ( t - z ) / ( s * s ) ) * exp1 / exp2 + ( emg_const * h * w * w ) / ( s * s * s ) * sqrt_2pi * exp1 * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

            // f'(z)
            derivative_retention = h * w / ( s * s ) * sqrt_2pi * exp1 / exp2 - ( emg_const * h ) / s * sqrt_2pi * exp1 * exp3 / ( ( exp2 * exp2 ) * sqrt_2 );

            // set the jacobian matrix
            gsl_matrix_set( J, i, 0, derivative_height );
            gsl_matrix_set( J, i, 1, derivative_width );
            gsl_matrix_set( J, i, 2, derivative_symmetry );
            gsl_matrix_set( J, i, 3, derivative_retention );
        }

        return GSL_SUCCESS;
    }

    Int EmgFitter1D::evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
    {
      EmgFitter1D::residual_( x, params, f );
      EmgFitter1D::jacobian_( x, params, J );

      return GSL_SUCCESS;
    }

    void EmgFitter1D::printState_(Int iter, gsl_multifit_fdfsolver * s)
    {
      printf ( "iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
               gsl_vector_get( s->x, 0 ),
               gsl_vector_get( s->x, 1 ),
               gsl_vector_get( s->x, 2 ),
               gsl_vector_get( s->x, 3 ),
               gsl_blas_dnrm2( s->f ) );
    }

    EmgFitter1D::QualityType EmgFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
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
        EmgFitter1D::Data d;
        d.n = set.size();
        d.set = set;

        // Compute start parameters
        setInitialParameters_(set);

        // Optimize parameter with Levenberg-Marquardt algorithm (GLS)
        CoordinateType x_init[ 4 ] = { height_, width_, symmetry_, retention_ };
        if ( symmetric_ == false )
        {
          optimize_(set, 4, x_init, &(residual_), &(jacobian_), &(evaluate_), &d);
        }

        // Set optimized parameters
        height_ = x_init[0];
        width_ = x_init[1];
        symmetry_ = x_init[2];
        retention_ = x_init[3];

#ifdef DEBUG_FEATUREFINDER
        if ( getGslStatus_() != "success" )
        {
          std::cout << "status: " << getGslStatus_() << std::endl;
        }
#endif

	      // build model
      	model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("EmgModel"));
	      model->setInterpolationStep( interpolation_step_ );

        Param tmp;
        tmp.setValue( "bounding_box:min", min_ );
        tmp.setValue( "bounding_box:max", max_ );
        tmp.setValue( "statistics:variance", statistics_.variance() );
        tmp.setValue( "statistics:mean", statistics_.mean() );
        tmp.setValue( "emg:height", height_ );
        tmp.setValue( "emg:width", width_ );
        tmp.setValue( "emg:symmetry", symmetry_ );
        tmp.setValue( "emg:retention", retention_ );
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

    void EmgFitter1D::setInitialParameters_(const RawDataArrayType& set)
    {
      // sum over all intensities
      CoordinateType sum = 0.0;
      for (Size i=0; i<set.size(); ++i) sum += set[i].getIntensity();

      // calculate the median
      Size median = 0;
      Real count = 0.0;
      for ( Size i = 0; i < set.size(); ++i )
      {
        count += set[i].getIntensity();
        if ( count <= sum / 2 ) median = i;
      }

      // calculate the height of the peak
      height_ = set[median].getIntensity();

      // calculate retention time
      retention_ = set[median].getPos();

      // default is an asymmetric peak
      symmetric_ = false;

      // calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
      symmetry_ = fabs( set[set.size()-1].getPos() - set[median].getPos() ) / fabs( set[median].getPos() - set[0].getPos() );

      // check the symmetry
      if ( boost::math::isinf( symmetry_ ) || boost::math::isnan( symmetry_ ) )
      {
        symmetric_ = true;
        symmetry_ = 10;
      }

      // optimize the symmetry
      // The computations can lead to an overflow error at very low values of symmetry (s~0).
      // For s~5 the parameter can be aproximized by the Levenberg-Marquardt argorithms.
      // (the other parameters are much greater than one)
      if ( symmetry_ < 1 ) symmetry_ += 5;

      // calculate the width of the peak
      // rt-values with intensity zero are not allowed for calculation of the width
      // normally: width_ = fabs( set[set.size() - 1].getPos() - set[0].getPos() );
      // but its better for the emg function to proceed from narrow peaks
      width_ = symmetry_;
    }

    void EmgFitter1D::updateMembers_()
    {
        LevMarqFitter1D::updateMembers_();
        statistics_.setVariance(param_.getValue("statistics:variance"));
    }

}
