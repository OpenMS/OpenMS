// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <boost/math/distributions/normal.hpp>
#include <unsupported/Eigen/NonLinearOptimization>

#include <sstream>
#include <iostream>

using namespace std;

// #define GAUSS_FITTER_VERBOSE
// #undef  GAUSS_FITTER_VERBOSE

namespace OpenMS
{
  namespace Math
  {
    GaussFitter::GaussFitter()
    : init_param_(0.06, 3.0, 0.5)
    {
    }

    GaussFitter::~GaussFitter()
    {
    }

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
        m_values(data->size()), 
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
          throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-GaussFitter", "Could not fit the Gaussian to the data: Error " + String(status));
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

  }   //namespace Math
} // namespace OpenMS
