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
}

namespace OpenMS
{
  class OPENMS_DLLAPI LevMarqFitter1D : public Fitter1D
  {
  public:
    typedef std::vector<double> ContainerType;

    class GenericFunctor
    {
    public:
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }

      GenericFunctor(int dimensions, int num_data_points) 
        : m_inputs(dimensions), m_values(num_data_points) {}

      virtual ~GenericFunctor() {}

      virtual int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const = 0;
      virtual int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const = 0;

    protected:
      const int m_inputs, m_values;
    };

    LevMarqFitter1D();
    LevMarqFitter1D(const LevMarqFitter1D& source);
    ~LevMarqFitter1D() override;
    LevMarqFitter1D& operator=(const LevMarqFitter1D& source);

  protected:
    bool symmetric_;
    Int max_iteration_;

    void optimize_(Eigen::VectorXd& x_init, GenericFunctor& functor) const;
    void updateMembers_() override;
  };
}
