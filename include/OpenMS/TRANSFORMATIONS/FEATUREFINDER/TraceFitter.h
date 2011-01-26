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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_TRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_TRACEFITTER_H

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

namespace OpenMS
{

  /**
   * @brief Abstract fitter for RT profile fitting
   *
   * This class provides the basic interface and some functionality to fit multiple mass traces to
   * a given RT shape model using the Levenberg-Marquardt algorithm.
   *
   * @htmlinclude OpenMS_TraceFitter.parameters
   *
   * @todo docu needs update
   *
   */
  template <class PeakType>
  class TraceFitter
    : public DefaultParamHandler
  {

  public:
    /// default constructor.
    TraceFitter()
      : DefaultParamHandler("TraceFitter")
    {
      this->defaults_.setValue( "max_iteration", 500, "Maximum number of iterations using by Levenberg-Marquardt algorithm.", StringList::create("advanced") );
      this->defaults_.setValue( "epsilon_abs", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithm.", StringList::create("advanced") );
      this->defaults_.setValue( "epsilon_rel", 0.0001, "Relative error used by the Levenberg-Marquardt algorithm.", StringList::create("advanced") );
    }

    /// copy constructor
    TraceFitter(const TraceFitter& source)
      : DefaultParamHandler(source),
        epsilon_abs_(source.epsilon_abs_),
        epsilon_rel_(source.epsilon_rel_),
        max_iterations_(source.max_iterations_)
    {
    }


    /// assignment operator
    virtual TraceFitter& operator = (const TraceFitter& source)
    {
      DefaultParamHandler::operator=(source);
      max_iterations_ = source.max_iterations_;
      epsilon_abs_  = source.epsilon_abs_;
      epsilon_rel_ = source.epsilon_rel_;
      updateMembers_();

      return *this;
    }

    /// destructor
    virtual ~TraceFitter()
    {
    }

    /**
     * Main method of the TraceFitter which triggers the actual fitting.
     */
    virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType> & /* traces */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the lower bound of the fitted RT model
     */
    virtual DoubleReal getLowerRTBound() const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the upper bound of the fitted RT model
     */
    virtual DoubleReal getUpperRTBound() const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the height of the fitted model
     */
    virtual DoubleReal getHeight() const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the center position of the fitted model
     */
    virtual DoubleReal getCenter() const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the mass trace width at half max (FWHM)
     */
    virtual DoubleReal getFWHM() const
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns the theoretical value of the fitted model at position k in the passed Mass Trace
     */
    virtual DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> & /* trace */, Size /* k */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Checks if the fitted model fills out at least 'min_rt_span' of the RT span
     *
     * @param rt_bounds RT boundaries of the fitted model
     * @param min_rt_span Minimum RT span in relation to extended area that has to remain after model fitting
     */
    virtual bool checkMinimalRTSpan(const std::pair<DoubleReal,DoubleReal> & /* rt_bounds */, const DoubleReal /* min_rt_span */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Checks if the fitted model is not to big
     *
     * @param max_rt_span Maximum RT span in relation to extended area that the model is allowed to have
     */
    virtual bool checkMaximalRTSpan(const DoubleReal /* max_rt_span */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * ???
     * @todo docu needs update
     */
    virtual DoubleReal getFeatureIntensityContribution()
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Returns a textual representation of the fitted model function, that can be plotted using Gnuplot
     *
     * @parm trace The MassTrace that should be plotted
     * @param function_name The name of the function (e.g. f(x) -> function_name = f)
     * @param baseline The intensity of the baseline
     * @param rt_shift A shift value, that allows to plot all RT profiles side by side, even if they would overlap in reality.
     *                 This should be 0 for the first mass trace and increase by a fixed value for each mass trace.
     */
    virtual String getGnuplotFormula(FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType> const & /* trace */, const char /* function_name */, const DoubleReal /* baseline */, const DoubleReal /* rt_shift */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

  protected:

    /**
     *
     */
    virtual void printState_(SignedSize /* iter */, gsl_multifit_fdfsolver * /* s */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    virtual void updateMembers_()
    {
      max_iterations_ = this->param_.getValue("max_iteration");
      epsilon_abs_ = this->param_.getValue("epsilon_abs");
      epsilon_rel_ = this->param_.getValue("epsilon_rel");
    }

    /**
     *
     */
    virtual void getOptimizedParameters_(gsl_multifit_fdfsolver * /* s */)
    {
      throw Exception::NotImplemented(__FILE__,__LINE__,__PRETTY_FUNCTION__);
    }

    /**
     * Optimize the given parameters using the Levenberg-Marquardt algorithm.
     */
    void optimize_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType> & traces, const Size num_params, double x_init[],
        Int (* residual)(const gsl_vector * x, void * params, gsl_vector * f),
        Int (* jacobian)(const gsl_vector * x, void * params, gsl_matrix * J),
        Int (* evaluate)(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J))
    {
      const gsl_multifit_fdfsolver_type *T;
      gsl_multifit_fdfsolver *s;

      const size_t data_count = traces.getPeakCount();

      // gsl always expects N>=p or default gsl error handler invoked,
      // cause Jacobian be rectangular M x N with M>=N
      if ( data_count < num_params ) throw Exception::UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, gsl always expects N>=p" );

      gsl_multifit_function_fdf func;
      gsl_vector_view x = gsl_vector_view_array(x_init, num_params);
      const gsl_rng_type * type;
      gsl_rng * r;
      gsl_rng_env_setup();
      type = gsl_rng_default;
      r = gsl_rng_alloc(type);
      func.f = (residual);
      func.df = (jacobian);
      func.fdf = (evaluate);
      func.n = data_count;
      func.p = num_params;
      func.params = &traces;
      T = gsl_multifit_fdfsolver_lmsder;
      s = gsl_multifit_fdfsolver_alloc(T, data_count, num_params);
      gsl_multifit_fdfsolver_set(s, &func, &x.vector);
      SignedSize iter = 0;
      Int gsl_status_;
			do
      {
        iter++;
        gsl_status_ = gsl_multifit_fdfsolver_iterate(s);
        printState_(iter, s);
        if (gsl_status_) break;
        gsl_status_ = gsl_multifit_test_delta(s->dx, s->x, epsilon_abs_, epsilon_rel_);
      }
      while (gsl_status_ == GSL_CONTINUE && iter < max_iterations_);

      // get the parameters out of the fdfsolver
      getOptimizedParameters_(s);

      gsl_multifit_fdfsolver_free(s);
    }

    /** Test for the convergence of the sequence by comparing the last iteration step dx with the absolute error epsabs and relative error epsrel to the current position x */
    /// Absolute error
    DoubleReal epsilon_abs_;
    /// Relative error
    DoubleReal epsilon_rel_;
    /// Maximum number of iterations
    SignedSize max_iterations_;
    /// Current status of the gsl fitting
    //Int gsl_status_;
  };

}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_RTFITTING_H
