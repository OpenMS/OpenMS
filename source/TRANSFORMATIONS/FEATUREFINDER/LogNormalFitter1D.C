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
// $Maintainer: Marcel Grunert $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LogNormalFitter1D.h>
    
namespace OpenMS
{
      LogNormalFitter1D::LogNormalFitter1D()
      : LevMarqFitter1D()
      {
        setName(getProductName());
        defaults_.setValue("statistics:variance",1.0,"Variance of the model", true);
        defaultsToParam_();
      }
    
      LogNormalFitter1D::LogNormalFitter1D(const LogNormalFitter1D& source)
      : LevMarqFitter1D(source)
      {
        setParameters( source.getParameters() );
        updateMembers_();
      }
    
      LogNormalFitter1D::~LogNormalFitter1D()
      {
      }
    
      LogNormalFitter1D& LogNormalFitter1D::operator = (const LogNormalFitter1D& source)
      {
        if (&source == this) return *this;
    
        LevMarqFitter1D::operator = (source);
        setParameters( source.getParameters() );
        updateMembers_();
    
        return *this;
      }
                    
      Int LogNormalFitter1D::residual_(const gsl_vector* x, void* params, gsl_vector* f)
      {
        UInt n = static_cast<LogNormalFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LogNormalFitter1D::Data*> (params) ->set;
            
        CoordinateType h = gsl_vector_get( x, 0 );
        CoordinateType w = gsl_vector_get( x, 1 );
        CoordinateType s = gsl_vector_get( x, 2 );
        CoordinateType z = gsl_vector_get( x, 3 );
        CoordinateType r = 2; 
  
        CoordinateType Yi = 0.0;
  
        for ( UInt i = 0; i < n; i++ )
        {
          CoordinateType t = set[i].getPos();
  
          Yi = h * exp( -log( r ) / ( log( s ) * log( s ) ) * pow( log( ( t - z ) * ( s * s - 1 ) / ( w * s ) + 1 ), 2 ) );
  
          gsl_vector_set( f, i, ( Yi - set[i].getIntensity()) );
        }
    
        return GSL_SUCCESS;
      }
                  
      Int LogNormalFitter1D::jacobian_(const gsl_vector* x, void* params, gsl_matrix* J)
      {
        UInt n = static_cast<LogNormalFitter1D::Data*> (params) ->n;
        RawDataArrayType set = static_cast<LogNormalFitter1D::Data*> (params) ->set;
            
        CoordinateType h = gsl_vector_get( x, 0 );
        CoordinateType w = gsl_vector_get( x, 1 );
        CoordinateType s = gsl_vector_get( x, 2 );
        CoordinateType z = gsl_vector_get( x, 3 );
        CoordinateType r = 2;
  
        CoordinateType derivative_height, derivative_width, derivative_symmetry, derivative_retention, derivative_r = 0.0;
  
        // iterate over all points of the signal
        for ( UInt i = 0; i < n; i++ )
        {
          CoordinateType t = set[i].getPos();
  
          CoordinateType exp1 = exp( -log( r ) / ( log( s ) * log( s ) ) * pow( log( ( t - z ) * ( s * s - 1 ) / ( w * s ) + 1 ), 2 ) );
          CoordinateType term1 = ( ( ( t - z ) * ( s * s - 1 ) ) / ( w * s ) ) + 1;
          CoordinateType log_s = log( s );
          CoordinateType log_term1 = log( term1 );
          CoordinateType log_r = log( r );
  
          derivative_height = exp1;
  
          derivative_width = 2 * h * log_r / ( log_s * log_s ) * log_term1 * ( t - z ) * ( s * s - 1 ) / ( w * w ) / s / term1 * exp1;
  
          derivative_symmetry = h * ( 2 * log_r / ( log_s * log_s * log_s ) * ( log_term1 * log_term1 ) / s - 2 * log_r / ( log_s * log_s ) * log_term1 * ( 2 * ( t - z ) / w - ( t - z ) * ( s * s - 1 ) / ( w * s * s ) ) / term1 ) * exp1;
  
          derivative_retention = 2 * h * log_r / ( log_s * log_s ) * log_term1 * ( s * s - 1 ) / ( w * s ) / term1 * exp1;
  
          derivative_r = -h / r / ( log_s * log_s ) * ( log_term1 * log_term1 ) * exp1;
  
            // set the jacobian matrix
          gsl_matrix_set( J, i, 0, derivative_height );
          gsl_matrix_set( J, i, 1, derivative_width );
          gsl_matrix_set( J, i, 2, derivative_symmetry );
          gsl_matrix_set( J, i, 3, derivative_retention );
         }
        
          return GSL_SUCCESS;
       }
        
      Int LogNormalFitter1D::evaluate_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
      {
        LogNormalFitter1D::residual_( x, params, f );
        LogNormalFitter1D::jacobian_( x, params, J );
    
        return GSL_SUCCESS;
      }
          
      void LogNormalFitter1D::printState_(Int iter, gsl_multifit_fdfsolver * s)
      {
        printf ( "iter: %4u x = % 15.8f % 15.8f  % 15.8f  % 15.8f |f(x)| = %g\n", iter,
                gsl_vector_get( s->x, 0 ),
                gsl_vector_get( s->x, 1 ),
                gsl_vector_get( s->x, 2 ),
                gsl_vector_get( s->x, 3 ),
                gsl_blas_dnrm2( s->f ) );
      }
      
    
      LogNormalFitter1D::QualityType LogNormalFitter1D::fit1d(const RawDataArrayType& set, InterpolationModel*& model)
      {
        // Calculate bounding box
        min_ = max_ = set[0].getPos();
        for ( UInt pos=1; pos < set.size(); ++pos)
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
        LogNormalFitter1D::Data d;
        d.n = set.size();;
        d.set = set;
               
        // Compute start parameter
        setInitialParameters_(set);
  
        // Optimize parameter with Levenberg-Marquardt algorithm (GLS)                
        CoordinateType x_init[ 4 ] = { height_, width_, symmetry_, retention_ };
        if ( symmetric_ == false )
        {
          optimize_(set, 4, x_init, &(residual_), &(jacobian_), &(evaluate_), &d);
        }
          
        // Set optimized parameter
        height_ = x_init[0];
        width_ = x_init[1];
        symmetry_ = x_init[2];   
        retention_ = x_init[3];
        r_ = r_;
          
  #ifdef DEBUG_FEATUREFINDER                
        if ( getGslStatus_() != "success" )
        {
          std::cout << "status: " << getGslStatus_() << std::endl;
        } 
  #endif 
                  
       	// build model
        model = static_cast<InterpolationModel*> (Factory<BaseModel<1> >::create("LogNormalModel"));
				model->setInterpolationStep( interpolation_step_ );
  
        Param tmp;
        tmp.setValue( "bounding_box:min", min_ );
        tmp.setValue( "bounding_box:max", max_ );
        tmp.setValue( "statistics:variance", statistics_.variance() );
        tmp.setValue( "statistics:mean", statistics_.mean() );
        tmp.setValue( "lognormal:height", height_ );
        tmp.setValue( "lognormal:width", width_ );
        tmp.setValue( "lognormal:symmetry", symmetry_ );
        tmp.setValue( "lognormal:retention", retention_ );
        tmp.setValue( "lognormal:r", r_ );
        model->setParameters( tmp );
        
        return 1.0;
      }
        
      void LogNormalFitter1D::setInitialParameters_(const RawDataArrayType& set)
      {
        // sum over all intensities
        CoordinateType sum = 0.0;
        for (UInt i=0; i<set.size(); ++i) sum += set[i].getIntensity();
      
        // calculate the median
        Int median = 0;
        Real count = 0.0;
        for ( UInt i = 0; i < set.size(); ++i )
        {
          count += set[i].getIntensity();
          if ( count <= sum * 0.5 ) median = i;
        }
      
        // calculate the height of the peak
        height_ = set[median].getIntensity();
  
        // calculate the width of the peak
        // rt-values with intensity zero are not allowed for calculation of the width
        width_ = fabs( set[set.size() - 1].getPos() - set[0].getPos() );
      
        // calculate retention time
        retention_ = set[median].getPos();
      
        // default is an asymmetric peak
        symmetric_ = false;
      
        // calculate the symmetry (fronted peak: s<1 , tailed peak: s>1)
        symmetry_ = fabs( set[set.size()-1].getPos() - set[median].getPos() ) / fabs( set[median].getPos() - set[0].getPos() );
        
        // check the symmetry
        if ( isinf( symmetry_ ) || isnan( symmetry_ ) )
        {
          symmetric_ = true;
          symmetry_ = 10;
        }
      
        // optimize the symmetry
        // The computations can lead to an overflow error at very low values of symmetry (s~0).
        if ( symmetry_ <= 0.8 ) symmetry_ = 0.8;
        if ( symmetry_ == 1 ) symmetry_ = 1.1;
        if ( symmetry_ >= 1.5 ) symmetry_ = 1.4;
  
        // it is better to proceed from narrow peaks
        width_ *= 0.5;
        
        /* set the parameter r of the log normal function;
        r is the ratio between h and the height at which w and s are computed;
        r = 2, see "Mathematical functions for representation of chromatographic peaks", V.B. Di Marco(2001) */
        r_ = 2;
        
      }
         
      void LogNormalFitter1D::updateMembers_()
      {
        LevMarqFitter1D::updateMembers_();
        statistics_.setVariance(param_.getValue("statistics:variance"));
      }
    
}
