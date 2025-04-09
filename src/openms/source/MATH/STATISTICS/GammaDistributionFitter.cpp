// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <sstream>
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <unsupported/Eigen/NonLinearOptimization>

#include <OpenMS/MATH/STATISTICS/GammaDistributionFitter.h>

// #define GAMMA_DISTRIBUTION_FITTER_VERBOSE
// #undef  GAMMA_DISTRIBUTION_FITTER_VERBOSE

namespace OpenMS::Math
{

    GammaDistributionFitter::GammaDistributionFitter() :
      init_param_(1.0, 5.0)
    {
    }

    GammaDistributionFitter::~GammaDistributionFitter()
    = default;

    void GammaDistributionFitter::setInitialParameters(const GammaDistributionFitResult& param)
    {
      init_param_ = param;
    }

    struct GammaFunctor
    {
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GammaFunctor(unsigned dimensions, const std::vector<DPosition<2> >* data) :
        m_inputs(dimensions), 
        m_values(static_cast<int>(data->size())), 
        m_data(data) 
      {
      }

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
      {

        double b = x(0);
        double p = x(1);

        UInt i = 0;

        // gamma distribution is only defined for positive parameter values
        if (b > 0.0 && p > 0.0)
        {
          for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
          {
            double the_x = it->getX();
            fvec(i) =  std::pow(b, p) / std::tgamma(p) * std::pow(the_x, p - 1) * std::exp(-b * the_x) - it->getY();
          }
        }
        else
        {
          for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
          {
            fvec(i) = -it->getY();
          }
        }


        return 0;
      }

      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const
      {

        double b = x(0);
        double p = x(1);

        UInt i = 0;
        // gamma distribution is only defined for positive parameter values
        if (b > 0.0 && p > 0.0)
        {
          for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
          {
            double the_x = it->getX();

            // partial deviation regarding b
            double part_dev_b = std::pow(the_x, p - 1) * std::exp(-the_x * b) / std::tgamma(p) * (p * std::pow(b, p - 1) - the_x * std::pow(b, p));
            J(i, 0) = part_dev_b;

            // partial deviation regarding p
            double factor = std::exp(-b * the_x) * std::pow(the_x, p - 1) * std::pow(b, p) / std::pow(std::tgamma(p), 2);
            double argument = (std::log(b) + std::log(the_x)) * std::tgamma(p) - std::tgamma(p) * boost::math::digamma(p);
            double part_dev_p = factor * argument;
            J(i, 1) = part_dev_p;
          }
        }
        else
        {
          for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
          {
            J(i, 0) = 0.0;
            J(i, 1) = 0.0;
          }
        }
        return 0;
      }

      const int m_inputs, m_values;
      const std::vector<DPosition<2> >* m_data;
    };

    GammaDistributionFitter::GammaDistributionFitResult GammaDistributionFitter::fit(const std::vector<DPosition<2> >& input) const
    {
      Eigen::VectorXd x_init(2);
      x_init << init_param_.b, init_param_.p;
      GammaFunctor functor(2, &input);
      Eigen::LevenbergMarquardt<GammaFunctor> lmSolver(functor);
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GammaDistributionFitter", "Could not fit the gamma distribution to the data");
      }

#ifdef GAMMA_DISTRIBUTION_FITTER_VERBOSE
      std::stringstream formula;
      formula << "f(x)=" << "(" << x_init(0) << " ** " << x_init(1) << ") / gamma(" << x_init(1) << ") * x ** (" << x_init(1) << " - 1) * exp(- " << x_init(0) << " * x)";
      std::cout << formula.str() << std::endl;
#endif

      return GammaDistributionFitResult(x_init(0), x_init(1));
    }

} // namespace OpenMS //namespace Math
