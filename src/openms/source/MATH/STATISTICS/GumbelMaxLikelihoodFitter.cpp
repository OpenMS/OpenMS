// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
//

#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <unsupported/Eigen/NonLinearOptimization>

using namespace std;

namespace OpenMS::Math
{
    namespace // anonymous namespace to prevent name clashes with GumbleDistributionFitter
    {
      // Generic functor
      template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
      struct Functor
      {
        typedef _Scalar Scalar;
        enum {
          InputsAtCompileTime = NX,
          ValuesAtCompileTime = NY
        };
        typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
        typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
        typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

        int m_inputs, m_values;

        Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
        Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

        int inputs() const { return m_inputs; }
        int values() const { return m_values; }

      };

      struct GumbelDistributionFunctor : Functor<double>
      {

        GumbelDistributionFunctor(const std::vector<double>& data, const std::vector<double>& weights):
            Functor<double>(2,2),
            m_data(data), m_weights(weights)
        {
        }

        int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
        {
          fvec(0) = 0.0;
          double sigma = fabs(x(1));
          double logsigma = log(sigma);
          auto wit = m_weights.cbegin();
          for (auto it = m_data.cbegin(); it != m_data.cend(); ++it, ++wit)
          {
            double diff = (*it - x(0)) / sigma;
            fvec(0) += *wit * (-logsigma - diff - exp(-diff));
          }
          double foo = -fvec(0);
          fvec(0) = foo;
          fvec(1) = 0.0;
          return 0;
        }
        const std::vector<double>& m_data;
        const std::vector<double>& m_weights;
      };
    }

    GumbelMaxLikelihoodFitter::GumbelDistributionFitResult GumbelMaxLikelihoodFitter::fitWeighted(const std::vector<double> & x, const std::vector<double> & w)
    {
      Eigen::VectorXd x_init (2);
      x_init(0) = init_param_.a;
      x_init(1) = init_param_.b;
      GumbelDistributionFunctor functor (x, w);
      Eigen::NumericalDiff<GumbelDistributionFunctor> numDiff(functor);
      Eigen::LevenbergMarquardt<Eigen::NumericalDiff<GumbelDistributionFunctor>,double> lm(numDiff);
      Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::Status::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GumbelMaxLikelihoodFitter", "Could not fit the gumbel distribution to the data");
      }

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
      // build a formula with the fitted parameters for gnuplot
      stringstream formula;
      formula << "f(x)=" << "(1/" << x_init(1) << ") * " << "exp(( " << x_init(0) << "- x)/" << x_init(1) << ") * exp(-exp((" << x_init(0) << " - x)/" << x_init(1) << "))";
      cout << formula.str() << endl;
#endif
      init_param_.a = x_init(0);
      init_param_.b = fabs(x_init(1));

      return {x_init(0), fabs(x_init(1))};
    }

    double GumbelMaxLikelihoodFitter::GumbelDistributionFitResult::log_eval_no_normalize(const double x) const
    {
      // -log b is a constant again
      double diff = (x - a)/b;
      return -log(b) - diff - exp(- diff);
    }

    GumbelMaxLikelihoodFitter::GumbelMaxLikelihoodFitter():
    init_param_({0.25, 0.1})
    {}

    GumbelMaxLikelihoodFitter::GumbelMaxLikelihoodFitter(GumbelMaxLikelihoodFitter::GumbelDistributionFitResult init):
        init_param_(init)
    {}

    GumbelMaxLikelihoodFitter::~GumbelMaxLikelihoodFitter() = default;

    void GumbelMaxLikelihoodFitter::setInitialParameters(const GumbelDistributionFitResult & param)
    {
      init_param_ = param;
    }

} // namespace OpenMS   //namespace Math
