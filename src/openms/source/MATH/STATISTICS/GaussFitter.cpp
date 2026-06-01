// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <boost/math/distributions/normal.hpp>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace std;

// #define GAUSS_FITTER_VERBOSE
// #undef  GAUSS_FITTER_VERBOSE

namespace OpenMS::Math
{
    GaussFitter::GaussFitter()
    : init_param_(0.06, 3.0, 0.5)
    {
    }

    GaussFitter::~GaussFitter()
    = default;

    void GaussFitter::setInitialParameters(const GaussFitResult & param)
    {
      init_param_ = param;
    }

    struct GaussFunctor
    {
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GaussFunctor(int dimensions, const std::vector<DPosition<2> >* data)
      : m_inputs(dimensions), 
        m_values(static_cast<int>(data->size())),
        m_data(data)
      {}

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
      {
        const double A = x(0);
        const double x0 = x(1);
        const double sig = x(2);
        const double sig2 = 2 * sig * sig;

        UInt i = 0;
        for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
        {
          fvec(i) = A * std::exp(- (it->getX() - x0) * (it->getX() - x0) / sig2) - it->getY();
        }

        return 0;
      }
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J) const
      {
        const double A = x(0);
        const double x0 = x(1);
        const double sig = x(2);
        const double sig2 = 2 * sig * sig;
        const double sig3 = 2 * sig2 * sig;

        UInt i = 0;
        for (std::vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
        {
          const double xd = (it->getX() - x0);
          const double xd2 = xd*xd;
          double j0 = std::exp(-1.0 * xd2 / sig2);
          J(i,0) = j0;
          J(i,1) = (A * j0 * (-(-2 * it->getX() + 2.0 * x0) / sig2));
          J(i,2) = (A * j0 * (xd2 / sig3));
        }
        return 0;
      }

      const int m_inputs, m_values;
      const std::vector<DPosition<2> >* m_data;
    };

    GaussFitter::GaussFitResult GaussFitter::fit(vector<DPosition<2> > & input) const
    {
      Eigen::VectorXd x_init (3);
      x_init(0) = init_param_.A;
      x_init(1) = init_param_.x0;
      x_init(2) = init_param_.sigma;
      GaussFunctor functor (3, &input);
      Eigen::LevenbergMarquardt<GaussFunctor> lmSolver (functor);
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      // the states are poorly documented. after checking the source and
      // http://www.ultimatepp.org/reference%24Eigen_demo%24en-us.html we believe that
      // all states except TooManyFunctionEvaluation and ImproperInputParameters are good
      // termination states.
      if (status == Eigen::LevenbergMarquardtSpace::ImproperInputParameters ||
          status == Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation)
      {
          throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GaussFitter", "Could not fit the Gaussian to the data: Error " + String(status));
      }
      
      x_init(2) = fabs(x_init(2)); // sigma can be negative, but |sigma| would actually be the correct solution

#ifdef GAUSS_FITTER_VERBOSE
      std::stringstream formula;
      formula << "f(x)=" << result.A << " * exp(-(x - " << result.x0 << ") ** 2 / 2 / (" << result.sigma << ") ** 2)";
      std::cout << formular.str() << std::endl;
#endif
      
      return GaussFitResult (x_init(0), x_init(1), x_init(2));
    }

    // static
    std::vector<double> GaussFitter::eval(const std::vector<double>& evaluation_points, const GaussFitter::GaussFitResult& model)
    {
      std::vector<double> out;
      out.reserve(evaluation_points.size());
      boost::math::normal_distribution<> ndf(model.x0, model.sigma);
      double int0 = model.A / boost::math::pdf(ndf, model.x0); // intensity normalization factor of the max @ x0 (simply multiplying the CDF with A is wrong!)
      for (Size i = 0; i < evaluation_points.size(); ++i)
      {
        out.push_back(boost::math::pdf(ndf, evaluation_points[i]) * int0 );
      }
      return out;
    }

    double GaussFitter::GaussFitResult::eval(const double x) const
    {
      boost::math::normal_distribution<> ndf(x0, sigma);
      double int0 = A / boost::math::pdf(ndf, x0); // intensity normalization factor of the max @ x0 (simply multiplying the CDF with A is wrong!)
      return (boost::math::pdf(ndf, x) * int0 );
    }

    double GaussFitter::GaussFitResult::log_eval_no_normalize(const double x) const
    {
      //TODO we could cache log sigma but then we would need to make the members private and update log sigma whenever
      // sigma is reset
      //TODO for likelihood maximization also the halflogtwopi constant could be removed
      return -log(sigma) - halflogtwopi - 0.5 * pow((x - x0) / sigma, 2.0);
    }
  
} // namespace OpenMS  //namespace Math
