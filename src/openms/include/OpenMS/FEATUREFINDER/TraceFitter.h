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

// forward decl
namespace Eigen
{
    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    class Matrix;
    using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
    using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;
}

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
  class OPENMS_DLLAPI TraceFitter :
    public DefaultParamHandler
  {

public:
    /** Generic functor for LM-Optimization */
    //TODO: This is copy and paste from LevMarqFitter1d.h. Make a generic wrapper for LM optimization
    class GenericFunctor
    {
public:
      int inputs() const;
      int values() const;

      GenericFunctor(int dimensions, int num_data_points);

      virtual ~GenericFunctor();

      virtual int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) = 0;
      // compute Jacobian matrix for the different parameters
      virtual int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) = 0;

protected:
      const int m_inputs, m_values;
    };

    /// default constructor
    TraceFitter();

    /// copy constructor
    TraceFitter(const TraceFitter& source);

    /// assignment operator
    TraceFitter& operator=(const TraceFitter& source);

    /// destructor
    ~TraceFitter() override;

    /**
     * Main method of the TraceFitter which triggers the actual fitting.
     */
    virtual void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces) = 0;

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
    double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k) const;

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
    virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift) = 0;

protected:
    struct ModelData
    {
      FeatureFinderAlgorithmPickedHelperStructs::MassTraces* traces_ptr;
      bool weighted;
    };

    void updateMembers_() override;

    /**
     * Updates all member variables to the fitted values stored in the solver.
     *
     * @param s The solver containing the fitted parameter values.
     */
    virtual void getOptimizedParameters_(const Eigen::VectorXd& s) = 0;
    /**
     * Optimize the given parameters using the Levenberg-Marquardt algorithm.
     */
    void optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor);

    /// Maximum number of iterations
    SignedSize max_iterations_;
    /// Whether to weight mass traces by theoretical intensity during the optimization
    bool weighted_;

  };

}

