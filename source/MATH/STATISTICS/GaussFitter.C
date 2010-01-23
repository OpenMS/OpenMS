// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/MATH/STATISTICS/GaussFitter.h>
#include <sstream>

using namespace std;

#define GAUSS_FITTER_VERBOSE
#undef  GAUSS_FITTER_VERBOSE

#include <iostream>

namespace OpenMS
{
	namespace Math
	{
		GaussFitter::GaussFitter()
		{
			init_param_.A = 0.06;
			init_param_.x0 = 3.0;
			init_param_.sigma = 0.5;
		}
	
		GaussFitter::~GaussFitter()
		{
		}
	
		void GaussFitter::setInitialParameters(const GaussFitResult& param)
		{
			init_param_.A = param.A;
			init_param_.x0 = param.x0;
			init_param_.sigma = param.sigma;
		}
	
		const String& GaussFitter::getGnuplotFormula() const
		{
			return gnuplot_formula_;
		}
	
		int GaussFitter::gaussFitterf_(const gsl_vector* x, void* params, gsl_vector* f)
		{
			vector<DPosition<2> >* data = static_cast<vector<DPosition<2> >*>(params);
			
			double A = gsl_vector_get (x, 0);
			double x0 = gsl_vector_get (x, 1);
			double sig = gsl_vector_get (x, 2);
		
			UInt i = 0;
			for (vector<DPosition<2> >::iterator it = data->begin(); it != data->end(); ++it)
			{
				gsl_vector_set(f, i++,  A * exp(-1.0 * pow(it->getX() - x0, 2) / (2 * pow(sig, 2))) - it->getY());
			}
		
			return GSL_SUCCESS;
		}
	
		int GaussFitter::gaussFitterdf_(const gsl_vector* x, void* params, gsl_matrix* J)
		{
			vector<DPosition<2> >* data = static_cast<vector<DPosition<2> >*>(params);
	
			double A = gsl_vector_get (x, 0);
			double x0 = gsl_vector_get (x, 1);
			double sig = gsl_vector_get (x, 2);
		
			UInt i = 0;
			for (vector<DPosition<2> >::iterator it = data->begin(); it != data->end(); ++it)
			{
				gsl_matrix_set (J, i, 0, exp(-1.0 * pow(it->getX() - x0, 2) / (2 * pow(sig, 2))));
				gsl_matrix_set (J, i, 1, (A * exp(-1.0 * pow(it->getX() - x0, 2) / (2 * pow(sig, 2))) * (-1 * (-2 * it->getX() + 2.0 * x0) / (2 * pow(sig, 2)))));
				gsl_matrix_set (J, i++, 2, (A * exp(-1.0 * pow(it->getX() - x0, 2) / (2 * pow(sig, 2))) * (pow(it->getX() - x0, 2) / (4 * pow(sig, 3)))));				
			}
		  return GSL_SUCCESS;
		}
	
		int GaussFitter::gaussFitterfdf_(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)
		{
		  gaussFitterf_(x, params, f);
		  gaussFitterdf_(x, params, J);
		  return GSL_SUCCESS;
		}
	
#ifdef GAUSS_FITTER_VERBOSE
		void GaussFitter::printState_(size_t iter, gsl_multifit_fdfsolver * s)
		{
		  printf ("iter: %3u x = % 15.8f % 15.8f % 15.8f "
		          "|f(x)| = %g\n",
		          iter,
		          gsl_vector_get(s->x, 0), 
		          gsl_vector_get(s->x, 1),
		          gsl_vector_get(s->x, 2), 
		          /*gsl_blas_dnrm2(s->f)*/0);
		}
#endif
	
		GaussFitter::GaussFitResult GaussFitter::fit(vector<DPosition<2> >& input)
		{
		  const gsl_multifit_fdfsolver_type *T = NULL;
		  gsl_multifit_fdfsolver *s = NULL;
		
		  int status(0);
		  size_t iter = 0;
		
		  const size_t p = 3;
		
		  //gsl_matrix *covar = gsl_matrix_alloc (p, p);
		  gsl_multifit_function_fdf f;
		  double x_init[3] = { init_param_.A, init_param_.x0, init_param_.sigma };
		  gsl_vector_view x = gsl_vector_view_array (x_init, p);
		  const gsl_rng_type * type = NULL;
		  gsl_rng * r = NULL;
		
		  gsl_rng_env_setup();
		
		  type = gsl_rng_default;
		  r = gsl_rng_alloc (type);
		
		  f.f = &gaussFitterf_;
		  f.df = &gaussFitterdf_;
		  f.fdf = &gaussFitterfdf_;
		  f.n = input.size();
		  f.p = p;
		  f.params = &input;
		
		  T = gsl_multifit_fdfsolver_lmsder;
		  s = gsl_multifit_fdfsolver_alloc (T, input.size(), p);
		  gsl_multifit_fdfsolver_set (s, &f, &x.vector);
	
			#ifdef GAUSS_FITTER_VERBOSE
		  printState_(iter, s);
			#endif
		
		  do
		  {
		    iter++;
		    status = gsl_multifit_fdfsolver_iterate (s);
		
				#ifdef GAUSS_FITTER_VERBOSE
		    printf ("status = %s\n", gsl_strerror (status));
		    printState_(iter, s);
	
				cerr << "f(x)=" << gsl_vector_get(s->x, 0) << " * exp(-(x - " << gsl_vector_get(s->x, 1) << ") ** 2 / 2 / (" << gsl_vector_get(s->x, 2) << ") ** 2)";
				#endif
	
		    if (status)
				{
			    break;
				}
	
		    status = gsl_multifit_test_delta (s->dx, s->x,
		                                        1e-4, 1e-4);
		  }
		  while (status == GSL_CONTINUE && iter < 500);
	
			if (status!=GSL_SUCCESS)
			{
				gsl_rng_free(r);
				gsl_multifit_fdfsolver_free(s);
	
				throw Exception::UnableToFit(__FILE__,__LINE__,__PRETTY_FUNCTION__,"UnableToFit-GaussFitter","Could not fit the gaussian to the data");
			}
		  
			// write the result in a GaussFitResult struct
			GaussFitResult result;
			result.A = gsl_vector_get(s->x, 0);
			result.x0 = gsl_vector_get(s->x, 1);
			result.sigma = gsl_vector_get(s->x, 2);
	
			// build a formula with the fitted parameters for gnuplot
			stringstream formula;
			formula << "f(x)=" << result.A << " * exp(-(x - " << result.x0 << ") ** 2 / 2 / (" << result.sigma << ") ** 2)";
			gnuplot_formula_ = formula.str();
			#ifdef GAUSS_FITTER_VERBOSE
			cout << gnuplot_formula_ << endl;
			#endif
			
			gsl_rng_free(r);
			gsl_multifit_fdfsolver_free (s);
	
			return result;
		}
		
	} //namespace Math
} // namespace OpenMS

