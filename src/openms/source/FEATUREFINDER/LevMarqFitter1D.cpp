// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <fstream>

namespace OpenMS
{

    void LevMarqFitter1D::optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor) const
    {
      //TODO: this function is copy&paste from TraceFitter.h. Make a generic wrapper for
      //LM optimization
      int data_count = functor.values();
      int num_params = functor.inputs();

      // LM always expects N>=p, cause Jacobian be rectangular M x N with M>=N
      if (data_count < num_params)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Skipping feature, we always expects N>=p");
      }
      Eigen::LevenbergMarquardt<GenericFunctor> lmSolver (functor);
      lmSolver.parameters.maxfev = max_iteration_;
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Could not fit the gaussian to the data: Error " + String(status));
      }
    }

    void LevMarqFitter1D::updateMembers_()
    {
      Fitter1D::updateMembers_();
      max_iteration_ = this->param_.getValue("max_iteration");
    }


} // namespace OpenMS

