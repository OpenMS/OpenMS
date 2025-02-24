// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/TraceFitter.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Core>

namespace OpenMS
{

  int TraceFitter::GenericFunctor::inputs() const
  {
    return m_inputs;
  }

  int TraceFitter::GenericFunctor::values() const
  {
    return m_values;
  }

  TraceFitter::GenericFunctor::GenericFunctor(int dimensions, int num_data_points) :
    m_inputs(dimensions), m_values(num_data_points)
  {
  }

  TraceFitter::GenericFunctor::~GenericFunctor() = default;

  TraceFitter::TraceFitter() :
    DefaultParamHandler("TraceFitter")
  {
    defaults_.setValue("max_iteration", 500, "Maximum number of iterations used by the Levenberg-Marquardt algorithm.", {"advanced"});
    defaults_.setValue("weighted", "false", "Weight mass traces according to their theoretical intensities.", {"advanced"});
    defaults_.setValidStrings("weighted", {"true","false"});
    defaultsToParam_();
  }

  TraceFitter::TraceFitter(const TraceFitter& source) :
    DefaultParamHandler(source),
    max_iterations_(source.max_iterations_),
    weighted_(source.weighted_)
  {
    updateMembers_();
  }

  TraceFitter& TraceFitter::operator=(const TraceFitter& source)
  {
    DefaultParamHandler::operator=(source);
    max_iterations_ = source.max_iterations_;
    weighted_ = source.weighted_;
    updateMembers_();

    return *this;
  }

  TraceFitter::~TraceFitter() = default;

  double TraceFitter::computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k) const
  {
    double rt = trace.peaks[k].first;

    return trace.theoretical_int * getValue(rt);
  }

  void TraceFitter::updateMembers_()
  {
    max_iterations_ = this->param_.getValue("max_iteration");
    weighted_ = this->param_.getValue("weighted") == "true";
  }

  void TraceFitter::optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor)
  {
    //TODO: this function is copy&paste from LevMarqFitter1d.h. Make a generic wrapper for
    //LM optimization
    int data_count = functor.values();
    int num_params = functor.inputs();

    // LM always expects N>=p, cause Jacobian be rectangular M x N with M>=N
    if (data_count < num_params)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Skipping feature, we always expects N>=p");
    }

    Eigen::LevenbergMarquardt<GenericFunctor> lmSolver(functor);
    lmSolver.parameters.maxfev = max_iterations_;
    Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

    //the states are poorly documented. after checking the source, we believe that
    //all states except NotStarted, Running and ImproperInputParameters are good
    //termination states.
    if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Could not fit the gaussian to the data: Error " + String(status));
    }

    getOptimizedParameters_(x_init);
  }

} // namespace OpenMS
