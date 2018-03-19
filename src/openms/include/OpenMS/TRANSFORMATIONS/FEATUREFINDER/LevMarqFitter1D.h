// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

#include <unsupported/Eigen/NonLinearOptimization>

#include <algorithm>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{

  /**
    @brief Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization.
      */
  class OPENMS_DLLAPI LevMarqFitter1D :
    public Fitter1D
  {

public:

    typedef std::vector<double> ContainerType;

    /** Generic functor for LM-Optimization */
    //TODO: This is copy and paste from TraceFitter.h. Make a generic wrapper for LM optimization
    class GenericFunctor
    {
    public:
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GenericFunctor(int dimensions, int num_data_points)
      : m_inputs(dimensions), m_values(num_data_points) {}

      virtual ~GenericFunctor() {}

      virtual int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) = 0;

      // compute Jacobian matrix for the different parameters
      virtual int df(const Eigen::VectorXd &x, Eigen::MatrixXd &J) = 0;

    protected:
      const int m_inputs, m_values;
    };

    /// Default constructor
    LevMarqFitter1D() :
      Fitter1D()
    {
      this->defaults_.setValue("max_iteration", 500, "Maximum number of iterations using by Levenberg-Marquardt algorithm.", ListUtils::create<String>("advanced"));
    }

    /// copy constructor
    LevMarqFitter1D(const LevMarqFitter1D & source) :
      Fitter1D(source),
      max_iteration_(source.max_iteration_)
    {
    }

    /// destructor
    ~LevMarqFitter1D() override
    {
    }

    /// assignment operator
    virtual LevMarqFitter1D & operator=(const LevMarqFitter1D & source)
    {
      if (&source == this) return *this;

      Fitter1D::operator=(source);
      max_iteration_ = source.max_iteration_;

      return *this;
    }

protected:

    /// Parameter indicates symmetric peaks
    bool symmetric_;
    /// Maximum number of iterations
    Int max_iteration_;

    /**
        @brief Optimize start parameter

        @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    void optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor)
    {
      //TODO: this function is copy&paste from TraceFitter.h. Make a generic wrapper for
      //LM optimization
      int data_count = functor.values();
      int num_params = functor.inputs();

      // LM always expects N>=p, cause Jacobian be rectangular M x N with M>=N
      if (data_count < num_params) throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "UnableToFit-FinalSet", "Skipping feature, we always expects N>=p");

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

    void updateMembers_() override
    {
      Fitter1D::updateMembers_();
      max_iteration_ = this->param_.getValue("max_iteration");
    }

  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_LEVMARQFITTER1D_H
