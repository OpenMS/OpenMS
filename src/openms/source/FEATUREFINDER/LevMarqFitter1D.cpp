// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/LevMarqFitter1D.h>

#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Core>
#include <fstream>

namespace OpenMS
{
    /// Internal adapter to wrap GenericFunctor for Eigen's LM solver (const version)
    class GenericFunctorConstEigenAdapter
    {
    public:
      GenericFunctorConstEigenAdapter(const LevMarqFitter1D::GenericFunctor& functor)
        : functor_(functor)
      {
      }

      int inputs() const { return functor_.inputs(); }
      int values() const { return functor_.values(); }

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
      {
        return functor_(x.data(), fvec.data());
      }

      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const
      {
        return functor_.df(x.data(), J.data());
      }

    private:
      const LevMarqFitter1D::GenericFunctor& functor_;
    };

    void LevMarqFitter1D::optimize_(std::vector<double>& x_init, GenericFunctor& functor) const
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

      // Create Eigen vector and copy data from std::vector
      Eigen::VectorXd x_eigen = Eigen::Map<Eigen::VectorXd>(x_init.data(), x_init.size());

      // Create adapter to wrap our functor for Eigen's LM solver
      GenericFunctorConstEigenAdapter adapter(functor);

      Eigen::LevenbergMarquardt<GenericFunctorConstEigenAdapter> lmSolver(adapter);
      lmSolver.parameters.maxfev = max_iteration_;
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_eigen);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Could not fit the gaussian to the data: Error " + String(status));
      }
      // Copy results back to x_init
      std::copy(x_eigen.data(), x_eigen.data() + x_eigen.size(), x_init.begin());
    }

    void LevMarqFitter1D::updateMembers_()
    {
      Fitter1D::updateMembers_();
      max_iteration_ = this->param_.getValue("max_iteration");
    }


} // namespace OpenMS

