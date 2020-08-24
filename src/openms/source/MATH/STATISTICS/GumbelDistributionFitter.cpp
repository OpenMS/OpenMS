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
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
//


#include <unsupported/Eigen/NonLinearOptimization>

#include <OpenMS/MATH/STATISTICS/GumbelDistributionFitter.h>

using namespace std;

// #define GUMBEL_DISTRIBUTION_FITTER_VERBOSE
// #undef  GUMBEL_DISTRIBUTION_FITTER_VERBOSE

namespace OpenMS
{
  namespace Math
  {
    double GumbelDistributionFitter::GumbelDistributionFitResult::log_eval_no_normalize(const double x) const
    {
      // -log b is a constant again
      double diff = (x - a)/b;
      return -log(b) - diff - exp(- diff);
    }

    GumbelDistributionFitter::GumbelDistributionFitter()
    {
      init_param_ = GumbelDistributionFitResult(0.25, 0.1);
    }

    GumbelDistributionFitter::~GumbelDistributionFitter() = default;

    void GumbelDistributionFitter::setInitialParameters(const GumbelDistributionFitResult & param)
    {
      init_param_ = param;
    }

    struct GumbelDistributionFunctor
    {
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GumbelDistributionFunctor(unsigned dimensions, const std::vector<DPosition<2> >* data)
      : m_inputs(dimensions), 
        m_values(static_cast<int>(data->size())), 
        m_data(data) 
      {
      }

      int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec)
      {
        double a = x(0); //location
        double b = x(1); //scale

        UInt i = 0;
        for (vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
        {
          double the_x = it->getX();
          double z = exp((a - the_x) / b);
          fvec(i) = (z * exp(-1 * z)) / b - it->getY();
        }
        return 0;
      }
      // compute Jacobian matrix for the different parameters
      int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J)
      {
        double a = x(0);
        double b = x(1);
        UInt i = 0;
        for (vector<DPosition<2> >::const_iterator it = m_data->begin(); it != m_data->end(); ++it, ++i)
        {
          double the_x = it->getX();
          double z = exp((a - the_x) / b);
          double f = z * exp(-1 * z);
          double part_dev_a = (f - pow(z, 2) * exp(-1 * z)) / pow(b, 2);
          J(i,0) = part_dev_a;
          double dev_z =  ((the_x - a) / pow(b, 2));
          double cum = f * dev_z;
          double part_dev_b = ((cum - z * cum) * b - f) / pow(b, 2);
          J(i,1) = part_dev_b;
        }
        return 0;
      }

      const int m_inputs, m_values;
      const std::vector<DPosition<2> >* m_data;
    };

    GumbelDistributionFitter::GumbelDistributionFitResult GumbelDistributionFitter::fit(vector<DPosition<2> > & input)
    {

      Eigen::VectorXd x_init (2);
      x_init(0) = init_param_.a;
      x_init(1) = init_param_.b;
      GumbelDistributionFunctor functor (2, &input);
      Eigen::LevenbergMarquardt<GumbelDistributionFunctor> lmSolver (functor);
      Eigen::LevenbergMarquardtSpace::Status status = lmSolver.minimize(x_init);

      //the states are poorly documented. after checking the source, we believe that
      //all states except NotStarted, Running and ImproperInputParameters are good
      //termination states.
      if (status <= Eigen::LevenbergMarquardtSpace::ImproperInputParameters)
      {
        throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-GumbelDistributionFitter", "Could not fit the gumbel distribution to the data");
      }

#ifdef GUMBEL_DISTRIBUTION_FITTER_VERBOSE
      // build a formula with the fitted parameters for gnuplot
      stringstream formula;
      formula << "f(x)=" << "(1/" << x_init(1) << ") * " << "exp(( " << x_init(0) << "- x)/" << x_init(1) << ") * exp(-exp((" << x_init(0) << " - x)/" << x_init(1) << "))";
      cout << formula.str() << endl;
#endif

      return {x_init(0), x_init(1)};
    }

  }   //namespace Math
} // namespace OpenMS
