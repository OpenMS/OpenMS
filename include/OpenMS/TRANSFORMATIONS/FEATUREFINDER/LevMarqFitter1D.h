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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace OpenMS
{
    namespace Internal
    {
        /// Helper struct ... the structure contains two components, the size of an area and a raw data container
        struct Data
        {
            typedef RawDataPoint1D RawDataPointType;
            typedef DPeakArray<RawDataPointType > RawDataArrayType;
            
            size_t n;
            RawDataArrayType set;
        };
    }

/** 
  	@brief Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization
	@ingroup FeatureFinder
*/
    class LevMarqFitter1D
    : public Fitter1D
    {

      public:
        
        /// Default constructor
        LevMarqFitter1D()
        : Fitter1D()
        {
          this->defaults_.setValue( "max_iteration", 500, "Maximum number of iterations fitting.", true );
          this->defaults_.setValue( "deltaAbsError", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithm.", true );
          this->defaults_.setValue( "deltaRelError", 0.0001, "Relative error used by the Levenberg-Marquardt algorithm.", true );
        }
  
        /// copy constructor
        LevMarqFitter1D(const LevMarqFitter1D& source)
        : Fitter1D(source),
            max_iteration_(source.max_iteration_),
            abs_error_(source.abs_error_),
            rel_error_(source.rel_error_)
        {
        }
                          
        /// destructor
        virtual ~LevMarqFitter1D()
        {
        }
  
        /// assignment operator
        virtual LevMarqFitter1D& operator = (const LevMarqFitter1D& source)
        {
            if (&source ==this) return *this;
      
            Fitter1D::operator = (source);
            max_iteration_ = source.max_iteration_;
            abs_error_ = source.abs_error_;
            rel_error_ = source.rel_error_;
      
            return *this;
        }
           
    protected:
	
        /// GSL status		
        int gsl_status_;
        /// Parameter indicates symmetric peaks
        bool symmetric_;
        /// Maximum number of iterations
        CoordinateType max_iteration_;
        /** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
        /// Absolute error
        CoordinateType abs_error_;
        /// Relative error
        CoordinateType rel_error_;

        /** Diplay the intermediate state of the solution. The solver state contains 
            the vector s->x which is the current position, and the vector s->f with 
            corresponding function values */
        virtual void printState_(size_t iter, gsl_multifit_fdfsolver * s) = 0;  
     
        /// Return GSL status as string
        const String getGslStatus_()
        {
          return gsl_strerror( gsl_status_ );
        }
           
        /// Optimize start parameter
        void optimize_(const RawDataArrayType& set, int num_params, CoordinateType x_init[],
                       int (* residual)(const gsl_vector * x, void * params, gsl_vector * f),
                       int (* jacobian)(const gsl_vector * x, void * params, gsl_matrix * J),
                       int (* evaluate)(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J) )
        {
          const gsl_multifit_fdfsolver_type * T;
          gsl_multifit_fdfsolver *s;
    
          int status;
          size_t iter = 0;
          const size_t n = set.size();
    
          // number of parameter to be optimize
          unsigned int p = num_params;
          
          // gsl always expects N>=p or default gsl error handler invoked, 
          // cause Jacobian be rectangular M x N with M>=N
          if ( n < p ) throw UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, gsl always expects N>=p" );
    
          gsl_matrix *covar = gsl_matrix_alloc( p, p );
          gsl_multifit_function_fdf f;
    
          gsl_vector_view x;
          x = gsl_vector_view_array( x_init, p );
        
          const gsl_rng_type * type;
          gsl_rng * r;
          gsl_rng_env_setup();
          type = gsl_rng_default;
          r = gsl_rng_alloc ( type );
    
          struct Internal::Data d = { n, set };
        
          f.f = (residual);
          f.df = (jacobian);
          f.fdf = (evaluate);
          f.n = set.size();
          f.p = p;
          f.params = &d;
    
          T = gsl_multifit_fdfsolver_lmsder;
          s = gsl_multifit_fdfsolver_alloc( T, n, p );
          gsl_multifit_fdfsolver_set( s, &f, &x.vector );

#ifdef DEBUG_FEATUREFINDER
          printState_(iter, s);
#endif
  
          // this is the loop for fitting
          do
          {
            iter++;
            status = gsl_multifit_fdfsolver_iterate ( s );
          
#ifdef DEBUG_FEATUREFINDER
            printState_(iter, s);
#endif
            
            /* check if solver is stuck */
            if ( status ) break;
            status = gsl_multifit_test_delta( s->dx, s->x, abs_error_, rel_error_ );
          }
          while ( status == GSL_CONTINUE && iter < max_iteration_ );
          
          // This function uses Jacobian matrix J to compute the covariance matrix of the best-fit parameters, covar. The parameter epsrel (0.0) is used to remove linear-dependent columns when J is rank deficient.
          gsl_multifit_covar( s->J, 0.0, covar );
  
#ifdef DEBUG_FEATUREFINDER
          gsl_matrix_fprintf( stdout, covar, "covar %g" );
#endif
  
#define FIT(i) gsl_vector_get(s->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
  
          // Set GSl status
          gsl_status_ = status;
        
#ifdef DEBUG_FEATUREFINDER
          { 
            // chi-squared value
            double chi = gsl_blas_dnrm2(s->f);
            double dof = n - p;
            double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
                                
            printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
                              
            for (size_t i=0; i<p; ++i)
            {
              printf( i".Parameter = %.5f +/- %.5f\n", FIT( i ), c*ERR( i ) );
            }
                      
          }
#endif
        
          // set optimized parameter  
          for (size_t i=0; i<p; ++i)
          {
            x_init[i] = FIT( i );
          }
          
          gsl_multifit_fdfsolver_free (s);
          gsl_matrix_free (covar);
          gsl_rng_free (r);
        }       

        void updateMembers_()
        {
            Fitter1D::updateMembers_();
            max_iteration_ = this->param_.getValue("max_iteration");
            abs_error_ = this->param_.getValue("deltaAbsError");
            rel_error_ = this->param_.getValue("deltaRelError");
        }
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H
