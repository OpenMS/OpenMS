// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FEATUREFINDER/Fitter1D.h>
#include <algorithm>

// forward decl
namespace Eigen
{
  template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
  class Matrix;
  using MatrixXd = Matrix<double, -1, -1, 0, -1, -1>;
  using VectorXd = Matrix<double, -1, 1, 0, -1, 1>;
} // namespace Eigen

namespace OpenMS
{

  /**
    @brief Abstract class for 1D-model fitter using Levenberg-Marquardt algorithm for parameter optimization.
      */
  class OPENMS_DLLAPI LevMarqFitter1D : public Fitter1D
  {
  public:
    typedef std::vector<double> ContainerType;

    /** Generic functor for LM-Optimization */
    // TODO: This is copy and paste from TraceFitter.h. Make a generic wrapper for LM optimization
    class GenericFunctor
    {
    public:
      int inputs() const
      {
        return m_inputs;
      }
      int values() const
      {
        return m_values;
      }

      GenericFunctor(int dimensions, int num_data_points) : m_inputs(dimensions), m_values(num_data_points)
      {
      }

      virtual ~GenericFunctor()
      {
      }

      virtual int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const = 0;

      // compute Jacobian matrix for the different parameters
      virtual int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const = 0;

    protected:
      const int m_inputs, m_values;
    };

    /// Default constructor
    LevMarqFitter1D() : Fitter1D()
    {
      this->defaults_.setValue("max_iteration", 500, "Maximum number of iterations using by Levenberg-Marquardt algorithm.", {"advanced"});
    }

    /// copy constructor
    LevMarqFitter1D(const LevMarqFitter1D& source) : Fitter1D(source), max_iteration_(source.max_iteration_)
    {
    }

    /// destructor
    ~LevMarqFitter1D() override
    {
    }

    /// assignment operator
    LevMarqFitter1D& operator=(const LevMarqFitter1D& source)
    {
      if (&source == this)
        return *this;

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
    void optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor) const;

    void updateMembers_() override;
  };
} // namespace OpenMS
