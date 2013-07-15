// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//

#include <sstream>

#include <cmath>

#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>

using namespace std;

#define GUMBEL_DISTRIBUTION_FITTER_VERBOSE
#undef  GUMBEL_DISTRIBUTION_FITTER_VERBOSE

namespace OpenMS
{
  namespace Math
  {
    GumbelDistributionFitter::GumbelDistributionFitter()
    {
      init_param_.a = 0.25;
      init_param_.b = 0.1;
    }

    GumbelDistributionFitter::~GumbelDistributionFitter()
    {
    }

    void GumbelDistributionFitter::setInitialParameters(const GumbelDistributionFitResult & param)
    {
      init_param_.a = param.a;
      init_param_.b = param.b;
    }

    const String & GumbelDistributionFitter::getGnuplotFormula() const
    {
      return gnuplot_formula_;
    }

    int GumbelDistributionFitter::gumbelDistributionFitterf_(const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f)
    {
      vector<DPosition<2> > * data = static_cast<vector<DPosition<2> > *>(params);

      double a = deprecated_gsl_vector_get(x, 0);
      double b = deprecated_gsl_vector_get(x, 1);

      UInt i = 0;
      for (vector<DPosition<2> >::iterator it = data->begin(); it != data->end(); ++it)
      {
        double the_x = it->getX();
        double z = exp((a - the_x) / b);
        deprecated_gsl_vector_set(f, i++, (z * exp(-1 * z)) / b - it->getY());
      }

      return deprecated_gsl_SUCCESS;
    }

    // compute Jacobian matrix for the different parameters
    int GumbelDistributionFitter::gumbelDistributionFitterdf_(const deprecated_gsl_vector * x, void * params, deprecated_gsl_matrix * J)
    {
      vector<DPosition<2> > * data = static_cast<vector<DPosition<2> > *>(params);

      double a = deprecated_gsl_vector_get(x, 0);
      double b = deprecated_gsl_vector_get(x, 1);
      UInt i(0);
      for (vector<DPosition<2> >::iterator it = data->begin(); it != data->end(); ++it, ++i)
      {
        double the_x = it->getX();
        double z = exp((a - the_x) / b);
        double f = z * exp(-1 * z);
        double part_dev_a = (f - pow(z, 2) * exp(-1 * z)) / pow(b, 2);
        deprecated_gsl_matrix_set(J, i, 0, part_dev_a);
        double dev_z =  ((the_x - a) / pow(b, 2));
        double cum = f * dev_z;
        double part_dev_b = ((cum - z * cum) * b - f) / pow(b, 2);
        deprecated_gsl_matrix_set(J, i, 1, part_dev_b);




      }
      return deprecated_gsl_SUCCESS;
    }

    int GumbelDistributionFitter::gumbelDistributionFitterfdf_(const deprecated_gsl_vector * x, void * params, deprecated_gsl_vector * f, deprecated_gsl_matrix * J)
    {
      gumbelDistributionFitterf_(x, params, f);
      gumbelDistributionFitterdf_(x, params, J);
      return deprecated_gsl_SUCCESS;
    }

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
    void GumbelDistributionFitter::printState_(size_t iter, deprecated_gsl_multifit_fdfsolver * s)
    {
      printf("iter: %3u x = % 15.8f % 15.8f "
             "|f(x)| = %g\n",
             (unsigned int)iter,
             deprecated_gsl_vector_get(s->x, 0),
             deprecated_gsl_vector_get(s->x, 1),
             deprecated_gsl_blas_dnrm2(s->f));
    }

#endif

    GumbelDistributionFitter::GumbelDistributionFitResult GumbelDistributionFitter::fit(vector<DPosition<2> > & input)
    {
      const deprecated_gsl_multifit_fdfsolver_type * T = NULL;
      deprecated_gsl_multifit_fdfsolver * s = NULL;

      int status = 0;
      size_t iter = 0;

      const size_t p = 2;

      double x_init[2] = { init_param_.a, init_param_.b };
      deprecated_gsl_vector_view_ptr x = deprecated_gsl_vector_view_array(x_init, p);
      const deprecated_gsl_rng_type * type = NULL;
      deprecated_gsl_rng * r = NULL;

      deprecated_gsl_rng_env_setup();

      type = deprecated_wrapper_get_gsl_rng_default();
      r = deprecated_gsl_rng_alloc(type);
      deprecated_gsl_multifit_function_fdf_ptr f
      		  = deprecated_wrapper_gsl_multifit_fdfsolver_lmsder_new (
      				&gumbelDistributionFitterf_,
      				&gumbelDistributionFitterdf_,
      				&gumbelDistributionFitterfdf_,
      			    input.size(),
      				p,
      				&input );

      T = deprecated_wrapper_get_multifit_fdfsolver_lmsder();
      s = deprecated_gsl_multifit_fdfsolver_alloc(T, input.size(), p);
      deprecated_gsl_multifit_fdfsolver_set(s, f.get(),
    		  deprecated_wrapper_gsl_vector_view_get_vector( x) );

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
      printState_(iter, s);
#endif

      do
      {
        ++iter;
        status = deprecated_gsl_multifit_fdfsolver_iterate(s);

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
        printf("status = %s\n", deprecated_gsl_strerror(status));
        printState_(iter, s);
#endif
        if (status)
        {
          break;
        }

        status = deprecated_gsl_multifit_test_delta(
        		deprecated_wrapper_gsl_multifit_fdfsolver_get_dx(s),
        		deprecated_wrapper_gsl_multifit_fdfsolver_get_x(s), 1e-4, 1e-4);
#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
        printf("Status = '%s'\n", deprecated_gsl_strerror(status));
#endif
      }
      while (status == deprecated_gsl_CONTINUE && iter < 1000);

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
      printf("Final status = '%s'\n", deprecated_gsl_strerror(status));
#endif

      if (status != deprecated_gsl_SUCCESS)
      {
        deprecated_gsl_rng_free(r);
        deprecated_gsl_multifit_fdfsolver_free(s);

        throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-GumbelDistributionFitter", "Could not fit the gumbel distribution to the data");
      }

      // write the result in a GumbelDistributionFitResult struct
      GumbelDistributionFitResult result;
      result.a = deprecated_gsl_vector_get(
    		  deprecated_wrapper_gsl_multifit_fdfsolver_get_x(s), 0);
      result.b = deprecated_gsl_vector_get(
    		  deprecated_wrapper_gsl_multifit_fdfsolver_get_x(s), 1);

      // build a formula with the fitted parameters for gnuplot
      stringstream formula;
      formula << "f(x)=" << "(1/" << result.b << ") * " << "exp(( " << result.a << "- x)/" << result.b << ") * exp(-exp((" << result.a << " - x)/" << result.b << "))";
      gnuplot_formula_ = formula.str();

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
      cout << gnuplot_formula_ << endl;
#endif

      deprecated_gsl_rng_free(r);
      deprecated_gsl_multifit_fdfsolver_free(s);

      return result;
    }

  }   //namespace Math
} // namespace OpenMS
