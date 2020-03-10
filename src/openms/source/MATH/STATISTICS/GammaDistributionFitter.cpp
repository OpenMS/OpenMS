// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

namespace OpenMS
{
  namespace Math
  {
    GammaDistributionFitter::GammaDistributionFitter() :
      init_param_(1.0, 5.0)
    {
    }

    GammaDistributionFitter::~GammaDistributionFitter()
    {
    }

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

      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec)
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
            fvec(i) =  std::pow(b, p) / boost::math::tgamma(p) * std::pow(the_x, p - 1) * std::exp(-b * the_x) - it->getY();
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
      int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J)
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
            double part_dev_b = std::pow(the_x, p - 1) * std::exp(-the_x * b) / boost::math::tgamma(p) * (p * std::pow(b, p - 1) - the_x * std::pow(b, p));
            J(i, 0) = part_dev_b;

            // partial deviation regarding p
            double factor = std::exp(-b * the_x) * std::pow(the_x, p - 1) * std::pow(b, p) / std::pow(boost::math::tgamma(p), 2);
            double argument = (std::log(b) + std::log(the_x)) * boost::math::tgamma(p) - boost::math::tgamma(p) * boost::math::digamma(p);
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

    GammaDistributionFitter::GammaDistributionFitResult GammaDistributionFitter::fit(const std::vector<DPosition<2> >& input)
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

  } //namespace Math
} // namespace OpenMS
