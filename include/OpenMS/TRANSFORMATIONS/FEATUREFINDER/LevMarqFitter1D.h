// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

#include <gsl/gsl_rng.h> // gsl random number generators
#include <gsl/gsl_randist.h> // gsl random number distributions
#include <gsl/gsl_vector.h> // gsl vector and matrix definitions
#include <gsl/gsl_multifit_nlin.h> // gsl multidimensional fitting
#include <gsl/gsl_blas.h> // gsl linear algebra stuff

namespace OpenMS
{
  
    /** 
      @brief Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization.
		*/
    class OPENMS_DLLAPI LevMarqFitter1D
    : public Fitter1D
    {

      public:
      
      	typedef std::vector < double > ContainerType;
        
        /// Default constructor
        LevMarqFitter1D()
        : Fitter1D()
        {
          this->defaults_.setValue( "max_iteration", 500, "Maximum number of iterations using by Levenberg-Marquardt algorithm.", StringList::create("advanced") );
          this->defaults_.setValue( "deltaAbsError", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithm.", StringList::create("advanced") );
          this->defaults_.setValue( "deltaRelError", 0.0001, "Relative error used by the Levenberg-Marquardt algorithm.", StringList::create("advanced") );
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
        Int gsl_status_;
        /// Parameter indicates symmetric peaks
        bool symmetric_;
        /// Maximum number of iterations
        Int max_iteration_;
        /** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
        /// Absolute error
        CoordinateType abs_error_;
        /// Relative error
        CoordinateType rel_error_;
       
        /** Display the intermediate state of the solution. The solver state contains 
            the vector s->x which is the current position, and the vector s->f with 
            corresponding function values */
        virtual void printState_(Int iter, gsl_multifit_fdfsolver * s) = 0;  
     
        /// Return GSL status as string
        const String getGslStatus_()
        {
          return gsl_strerror( gsl_status_ );
        }
           

				/**
					@brief Optimize start parameter
					
					@exception Exception::UnableToFit is thrown if fitting cannot be performed
				*/
        void optimize_(const RawDataArrayType& set, Int num_params, CoordinateType x_init[],
                       Int (* residual)(const gsl_vector * x, void * params, gsl_vector * f),
                       Int (* jacobian)(const gsl_vector * x, void * params, gsl_matrix * J),
                       Int (* evaluate)(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J),
											 void* advanced_params
											)
        {
          
          const gsl_multifit_fdfsolver_type *T;
          gsl_multifit_fdfsolver *s;
    
          Int status;
          Int iter = 0;
          const UInt n = (UInt)set.size();
    
          // number of parameters to be optimized
          UInt p = num_params;
          
          // gsl always expects N>=p or default gsl error handler invoked, 
          // cause Jacobian be rectangular M x N with M>=N
          if ( n < p ) throw Exception::UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, gsl always expects N>=p" );
    
		 		 // allocate space for a covariance matrix of size p by p
          gsl_matrix *covar = gsl_matrix_alloc( p, p );
          gsl_multifit_function_fdf f;
    
          gsl_vector_view x = gsl_vector_view_array( x_init, p );
        
          gsl_rng_env_setup();
          
          // set up the function to be fit
          f.f = (residual); // the function of residuals
          f.df = (jacobian); // the gradient of this function
          f.fdf = (evaluate); // combined function and gradient
          f.n = set.size(); // number of points in the data set
          f.p = p; // number of parameters in the fit function
          f.params = advanced_params; // // structure with the data and error bars
          
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
            
            // perform a single iteration of the fitting routine
            status = gsl_multifit_fdfsolver_iterate ( s );
          
#ifdef DEBUG_FEATUREFINDER
						// customized routine to print out current parameters
						printState_(iter, s);
#endif
           
            /* check if solver is stuck */
            if ( status ) break;
            
            // test for convergence with an absolute and relative error
            status = gsl_multifit_test_delta( s->dx, s->x, abs_error_, rel_error_ );
          }
          while ( status == GSL_CONTINUE && iter < max_iteration_ );
        
          // This function uses Jacobian matrix J to compute the covariance matrix of the best-fit parameters, covar. 
          // The parameter epsrel (0.0) is used to remove linear-dependent columns when J is rank deficient.
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
            DoubleReal chi = gsl_blas_dnrm2(s->f);
            DoubleReal dof = n - p;
            DoubleReal c = GSL_MAX_DBL(1, chi / sqrt(dof)); 
                                
            printf("chisq/dof = %g\n",  pow(chi, 2.0) / dof);
                              
            for (Size i=0; i<p; ++i)
            {
              std::cout << i;
			  			printf(".Parameter = %.5f +/- %.5f\n", FIT( i ), c*ERR( i ) );
            }           
          }
#endif
        
          // set optimized parameters
          for (Size i = 0; i < p; ++i)
          {
            x_init[i] = FIT( i );
          }
          
          gsl_multifit_fdfsolver_free (s);
          gsl_matrix_free (covar);
       
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
