// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <vector>

namespace OpenMS
{

  /**
    @brief Abstract base class for fitting an RT-shape model to one or more mass traces.

    Provides the shared parameter handling and a Levenberg-Marquardt
    driver (@ref optimize_) for subclasses that supply a concrete RT
    profile (Gaussian, EGH, ...). Subclasses implement @ref fit, the
    geometric query methods (@ref getLowerRTBound,
    @ref getUpperRTBound, @ref getHeight, @ref getCenter,
    @ref getFWHM, @ref getValue, @ref getArea), the validation hooks
    (@ref checkMinimalRTSpan, @ref checkMaximalRTSpan), the plot
    helper @ref getGnuplotFormula, and @ref getOptimizedParameters_.

    @par Parameters
    - @c "max_iteration" (singular): maximum number of LM iterations
      (default @c 500).
    - @c "weighted" (@c "true"/@c "false"): whether mass traces are
      weighted by their theoretical intensity during fitting (default
      @c "false"). The flag is exposed as @ref weighted_ for use by
      subclasses.

    @ingroup FeatureFinder
  */
  class OPENMS_DLLAPI TraceFitter :
    public DefaultParamHandler
  {

public:
    /**
      @brief Generic LM functor with a raw-pointer interface, used to keep Eigen out of the public API.

      Subclasses implement the @c operator() to fill the residual
      vector and @ref df to fill the Jacobian; the cpp wraps each side
      in an @c Eigen::Map before handing it to the LM solver.
    */
    //TODO: This is copy and paste from LevMarqFitter1d.h. Make a generic wrapper for LM optimization
    class GenericFunctor
    {
public:
      /// Number of free parameters (input dimensionality).
      int inputs() const;
      /// Number of residuals (data points).
      int values() const;

      /**
        @brief Construct with the given problem dimensions.

        @param[in] dimensions      Number of free parameters.
        @param[in] num_data_points Number of residuals.
      */
      GenericFunctor(int dimensions, int num_data_points);

      virtual ~GenericFunctor();

      /**
        @brief Compute residuals at the current parameter vector.

        @param[in]  x    Parameter vector of size @ref inputs.
        @param[out] fvec Residual vector of size @ref values.
        @return @c 0 on success.
      */
      virtual int operator()(const double* x, double* fvec) = 0;

      /**
        @brief Compute the Jacobian at the current parameter vector.

        @param[in]  x Parameter vector of size @ref inputs.
        @param[out] J Jacobian matrix, @ref values rows by @ref inputs
                      columns, stored column-major.
        @return @c 0 on success.
      */
      virtual int df(const double* x, double* J) = 0;

protected:
      const int m_inputs, m_values;
    };

    /// Default constructor; registers the @c "max_iteration" and @c "weighted" defaults on @ref DefaultParamHandler.
    TraceFitter();

    /// Copy constructor.
    TraceFitter(const TraceFitter& source);

    /// Assignment operator.
    TraceFitter& operator=(const TraceFitter& source);

    /// Destructor.
    ~TraceFitter() override;

    /**
      @brief Run the LM fit on @p traces.

      Subclasses fill the model parameters and (typically) call
      @ref optimize_ on a problem-specific functor.

      @param[in,out] traces Mass traces to fit; subclasses may set
                            additional bookkeeping on @p traces.
    */
    virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces) = 0;

    /**
      @brief Lower RT bound of the fitted RT model.

      @return Lower RT bound (retention-time units).
    */
    virtual double getLowerRTBound() const = 0;

    /**
      @brief Upper RT bound of the fitted RT model.

      @return Upper RT bound (retention-time units).
    */
    virtual double getUpperRTBound() const = 0;

    /**
      @brief Height of the fitted RT model at its centre.

      @return Model intensity at the centre RT.
    */
    virtual double getHeight() const = 0;

    /**
      @brief Centre of the fitted RT model.

      @return Centre retention time of the fitted model.
    */
    virtual double getCenter() const = 0;

    /**
      @brief Full width at half maximum of the fitted RT model.

      @return FWHM in retention-time units.
    */
    virtual double getFWHM() const = 0;

    /**
      @brief Evaluate the fitted model at retention time @p rt.

      @param[in] rt Retention time at which to evaluate the model.
      @return Model intensity at @p rt.
    */
    virtual double getValue(double rt) const = 0;

    /**
      @brief Evaluate the fitted model at the RT of the @p k-th peak of @p trace, scaled by the trace's theoretical intensity.

      @param[in] trace Mass trace; only the RT of its @p k-th peak is read.
      @param[in] k     Index of the peak in @p trace.peaks.
      @return @c trace.theoretical_int * getValue(trace.peaks[k].first).
    */
    double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k) const;

    /**
      @brief Whether the fitted RT model covers at least @p min_rt_span of the extended search area.

      @param[in] rt_bounds   RT boundaries of the extended search area
                             (lower, upper).
      @param[in] min_rt_span Minimum RT span the fitted model must cover
                             relative to @p rt_bounds, as a fraction.
      @return @c true when the model spans at least @p min_rt_span of
              the search area.
    */
    virtual bool checkMinimalRTSpan(const std::pair<double, double>& rt_bounds, const double min_rt_span) = 0;

    /**
      @brief Whether the fitted RT model stays within @p max_rt_span of the extended search area.

      @param[in] max_rt_span Maximum allowed RT span of the fitted
                             model relative to the extended search
                             area, as a fraction.
      @return @c true when the model does not exceed @p max_rt_span of
              the search area.
    */
    virtual bool checkMaximalRTSpan(const double max_rt_span) = 0;

    /**
      @brief Area under the fitted RT model.

      @return Integrated intensity of the fitted model.
    */
    virtual double getArea() = 0;

    /**
      @brief Gnuplot expression of the fitted model for @p trace.

      @param[in] trace         Mass trace whose theoretical intensity scales the expression.
      @param[in] function_name Name of the function in the resulting expression (e.g. @c 'f' for @c "f(x)").
      @param[in] baseline      Intensity of the plotted baseline.
      @param[in] rt_shift      RT offset added to the expression so several traces can be plotted side by side
                               (typically @c 0 for the first trace and incremented for subsequent traces).
      @return Gnuplot expression as a @c String.
    */
    virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift) = 0;

protected:
    /// Helper bundle passed to functors so they can read both the traces and the weighting flag.
    struct ModelData
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces* traces_ptr; ///< Mass traces being fitted.
      bool weighted;                                                     ///< Mirror of @ref weighted_ at fit time.
    };

    void updateMembers_() override;

    /**
      @brief Subclass hook: copy the optimised parameter vector back into the subclass's model fields.

      Called by @ref optimize_ after the LM solver has converged.

      @param[in] s Optimised parameter vector (size matches the functor's @c inputs).
    */
    virtual void getOptimizedParameters_(const std::vector<double>& s) = 0;

    /**
      @brief Run Levenberg-Marquardt on @p functor starting from @p x_init.

      On return, @p x_init holds the optimised parameter vector and
      @ref getOptimizedParameters_ has been called. The LM iteration
      limit is taken from the @c "max_iteration" parameter.

      @param[in,out] x_init  Initial parameter vector; replaced with
                             the optimised vector on success.
      @param[in,out] functor LM functor providing residuals and the
                             Jacobian.

      @throws Exception::UnableToFit when the number of residuals
                                     reported by @p functor is smaller
                                     than its number of free
                                     parameters, or when the LM solver
                                     reports a non-recoverable state
                                     (@c NotStarted / @c Running /
                                     @c ImproperInputParameters).
    */
    void optimize_(std::vector<double>& x_init, GenericFunctor& functor);

    /// Maximum number of LM iterations; refreshed from @c "max_iteration" in @ref updateMembers_.
    SignedSize max_iterations_;
    /// Whether mass traces are weighted by theoretical intensity during fitting; refreshed from @c "weighted" in @ref updateMembers_.
    bool weighted_;

  };

}

