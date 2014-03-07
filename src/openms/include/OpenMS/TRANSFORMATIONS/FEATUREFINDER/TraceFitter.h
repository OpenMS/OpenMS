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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_TRACEFITTER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_TRACEFITTER_H

#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

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
   * @todo docu needs update
   *
   */
  template <class PeakType>
  class TraceFitter :
    public DefaultParamHandler
  {

public:
    /// default constructor.
    TraceFitter() :
      DefaultParamHandler("TraceFitter")
    {
      defaults_.setValue("max_iteration", 500, "Maximum number of iterations used by the Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
      defaults_.setValue("epsilon_abs", 0.0001, "Absolute error used by the Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
      defaults_.setValue("epsilon_rel", 0.0001, "Relative error used by the Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
      defaults_.setValue("weighted", "false", "Weight mass traces according to their theoretical intensities.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("weighted", ListUtils::create<String>("true,false"));
      defaultsToParam_();
    }

    /// copy constructor
    TraceFitter(const TraceFitter& source) :
      DefaultParamHandler(source),
      epsilon_abs_(source.epsilon_abs_),
      epsilon_rel_(source.epsilon_rel_),
      max_iterations_(source.max_iterations_),
      weighted_(source.weighted_)
    {
      updateMembers_();
    }

    /// assignment operator
    virtual TraceFitter& operator=(const TraceFitter& source)
    {
      DefaultParamHandler::operator=(source);
      max_iterations_ = source.max_iterations_;
      epsilon_abs_  = source.epsilon_abs_;
      epsilon_rel_ = source.epsilon_rel_;
      weighted_ = source.weighted_;
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
    virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>& traces) = 0;

    /**
     * Returns the lower bound of the fitted RT model
     */
    virtual DoubleReal getLowerRTBound() const = 0;

    /**
     * Returns the upper bound of the fitted RT model
     */
    virtual DoubleReal getUpperRTBound() const = 0;

    /**
     * Returns the height of the fitted model
     */
    virtual DoubleReal getHeight() const = 0;

    /**
     * Returns the center position of the fitted model
     */
    virtual DoubleReal getCenter() const = 0;

    /**
     * Returns the mass trace width at half max (FWHM)
     */
    virtual DoubleReal getFWHM() const = 0;

    /**
     * Evaluate the fitted model at a time point
     */
    virtual DoubleReal getValue(DoubleReal rt) const = 0;

    /**
     * Returns the theoretical value of the fitted model at position k in the passed mass trace
     *
     * @param trace the mass trace for which the value should be computed
     * @param k  use the position of the k-th peak to compute the value
     */
    DoubleReal computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, Size k)
    {
      double rt = trace.peaks[k].first;

      return trace.theoretical_int * getValue(rt);
    };

    /**
     * Checks if the fitted model fills out at least 'min_rt_span' of the RT span
     *
     * @param rt_bounds RT boundaries of the fitted model
     * @param min_rt_span Minimum RT span in relation to extended area that has to remain after model fitting
     */
    virtual bool checkMinimalRTSpan(const std::pair<DoubleReal, DoubleReal>& rt_bounds, const DoubleReal min_rt_span) = 0;

    /**
     * Checks if the fitted model is not to big
     *
     * @param max_rt_span Maximum RT span in relation to extended area that the model is allowed to have
     */
    virtual bool checkMaximalRTSpan(const DoubleReal max_rt_span) = 0;

    /**
     * Returns the peak area of the fitted model
     */
    virtual DoubleReal getArea() = 0;

    /**
     * Returns a textual representation of the fitted model function, that can be plotted using Gnuplot
     *
     * @param trace The mass trace that should be plotted
     * @param function_name The name of the function (e.g. f(x) -> function_name = f)
     * @param baseline The intensity of the baseline
     * @param rt_shift A shift value, that allows to plot all RT profiles side by side, even if they would overlap in reality.
     *                 This should be 0 for the first mass trace and increase by a fixed value for each mass trace.
     */
    virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<PeakType>& trace, const char function_name, const DoubleReal baseline, const DoubleReal rt_shift) = 0;

protected:

    /// Structure for passing data to GSL functions
    struct ModelData
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType>* traces_ptr;
      bool weighted;
    };

    /**
     * Prints the state of the current iteration (e.g., values of the parameters)
     *
     * @param iter Number of current iteration.
     * @param s The solver that also contains all the parameters.
     */
    virtual void printState_(SignedSize iter, gsl_multifit_fdfsolver* s) = 0;

    virtual void updateMembers_()
    {
      max_iterations_ = this->param_.getValue("max_iteration");
      epsilon_abs_ = this->param_.getValue("epsilon_abs");
      epsilon_rel_ = this->param_.getValue("epsilon_rel");
      weighted_ = this->param_.getValue("weighted") == "true";
    }

    /**
     * Updates all member variables to the fitted values stored in the solver.
     *
     * @param s The solver containing the fitted parameter values.
     */
    virtual void getOptimizedParameters_(gsl_multifit_fdfsolver* s) = 0;

    /**
     * Optimize the given parameters using the Levenberg-Marquardt algorithm.
     */
    void optimize_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<PeakType> & traces, const Size num_params, double x_init[],
                   Int (* residual)(const gsl_vector* x, void* params, gsl_vector* f),
                   Int (* jacobian)(const gsl_vector* x, void* params, gsl_matrix* J),
                   Int (* evaluate)(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix * J))
    {
      const gsl_multifit_fdfsolver_type* T;
      gsl_multifit_fdfsolver* s;

      const size_t data_count = traces.getPeakCount();

      // gsl always expects N>=p or default gsl error handler invoked,
      // cause Jacobian be rectangular M x N with M>=N
      if (data_count < num_params) throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, gsl always expects N>=p");

      gsl_multifit_function_fdf func;
      gsl_vector_view x = gsl_vector_view_array(x_init, num_params);
      gsl_rng_env_setup();
      func.f = (residual);
      func.df = (jacobian);
      func.fdf = (evaluate);
      func.n = data_count;
      func.p = num_params;
      ModelData params = {&traces, weighted_};
      func.params = &params;
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
    /// Whether to weight mass traces by theoretical intensity during the optimization
    bool weighted_;

  };

}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_RTFITTING_H
