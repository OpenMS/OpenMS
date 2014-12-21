// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <unsupported/Eigen/NonLinearOptimization>

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
  class TraceFitter :
    public DefaultParamHandler
  {

public:
    /** Generic functor for LM-Optimization */
    //TODO: This is copy and paste from LevMarqFitter1d.h. Make a generic wrapper for LM optimization
    class GenericFunctor
    {
    public:
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GenericFunctor(int dimensions, int num_data_points)
      : m_inputs(dimensions), m_values(num_data_points) {}

      virtual ~GenericFunctor() {}

      virtual int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) = 0;
      // compute Jacobian matrix for the different parameters
      virtual int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J) = 0;

    protected:
      const int m_inputs, m_values;
    };

    /// default constructor
    TraceFitter() :
      DefaultParamHandler("TraceFitter")
    {
      defaults_.setValue("max_iteration", 500, "Maximum number of iterations used by the Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
      defaults_.setValue("weighted", "false", "Weight mass traces according to their theoretical intensities.", ListUtils::create<String>("advanced"));
      defaults_.setValidStrings("weighted", ListUtils::create<String>("true,false"));
      defaultsToParam_();
    }

    /// copy constructor
    TraceFitter(const TraceFitter& source) :
      DefaultParamHandler(source),
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
    virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D>& traces) = 0;

    /**
     * Returns the lower bound of the fitted RT model
     */
    virtual double getLowerRTBound() const = 0;

    /**
     * Returns the upper bound of the fitted RT model
     */
    virtual double getUpperRTBound() const = 0;

    /**
     * Returns the height of the fitted model
     */
    virtual double getHeight() const = 0;

    /**
     * Returns the center position of the fitted model
     */
    virtual double getCenter() const = 0;

    /**
     * Returns the mass trace width at half max (FWHM)
     */
    virtual double getFWHM() const = 0;

    /**
     * Evaluate the fitted model at a time point
     */
    virtual double getValue(double rt) const = 0;

    /**
     * Returns the theoretical value of the fitted model at position k in the passed mass trace
     *
     * @param trace the mass trace for which the value should be computed
     * @param k  use the position of the k-th peak to compute the value
     */
    double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D>& trace, Size k)
    {
      double rt = trace.peaks[k].first;

      return trace.theoretical_int * getValue(rt);
    }

    /**
     * Checks if the fitted model fills out at least 'min_rt_span' of the RT span
     *
     * @param rt_bounds RT boundaries of the fitted model
     * @param min_rt_span Minimum RT span in relation to extended area that has to remain after model fitting
     */
    virtual bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span) = 0;

    /**
     * Checks if the fitted model is not to big
     *
     * @param max_rt_span Maximum RT span in relation to extended area that the model is allowed to have
     */
    virtual bool checkMaximalRTSpan(const double max_rt_span) = 0;

    /**
     * Returns the peak area of the fitted model
     */
    virtual double getArea() = 0;

    /**
     * Returns a textual representation of the fitted model function, that can be plotted using Gnuplot
     *
     * @param trace The mass trace that should be plotted
     * @param function_name The name of the function (e.g. f(x) -> function_name = f)
     * @param baseline The intensity of the baseline
     * @param rt_shift A shift value, that allows to plot all RT profiles side by side, even if they would overlap in reality.
     *                 This should be 0 for the first mass trace and increase by a fixed value for each mass trace.
     */
    virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace<Peak1D>& trace, const char function_name, const double baseline, const double rt_shift) = 0;

protected:
    struct ModelData
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces<Peak1D>* traces_ptr;
      bool weighted;
    };

    virtual void updateMembers_()
    {
      max_iterations_ = this->param_.getValue("max_iteration");
      weighted_ = this->param_.getValue("weighted") == "true";
    }

    /**
     * Updates all member variables to the fitted values stored in the solver.
     *
     * @param s The solver containing the fitted parameter values.
     */
    virtual void getOptimizedParameters_(const Eigen::VectorXd&) = 0;
    /**
     * Optimize the given parameters using the Levenberg-Marquardt algorithm.
     */
    void optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor)
    {
      //TODO: this function is copy&paste from LevMarqFitter1d.h. Make a generic wrapper for
      //LM optimization
      int data_count = functor.values();
      int num_params = functor.inputs();

      // LM always expects N>=p, cause Jacobian be rectangular M x N with M>=N
      if (data_count < num_params) throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Skipping feature, we always expects N>=p");


      Eigen::LevenbergMarquardt<GenericFunctor> lmSolver (functor);
      lmSolver.parameters.maxfev = max_iterations_;
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
          throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-FinalSet", "Could not fit the gaussian to the data: Error " + String(status));
      }

      getOptimizedParameters_(x_init);
    }

    /// Maximum number of iterations
    SignedSize max_iterations_;
    /// Whether to weight mass traces by theoretical intensity during the optimization
    bool weighted_;

  };

}

#endif // #ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMPICKED_RTFITTING_H
