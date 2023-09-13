// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Christian Ehrlich $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  /** Matrix utility functions.
   *
   */

  typedef boost::shared_ptr< const Eigen::MatrixXd > EigenMatrixXdPtr;
  typedef boost::shared_ptr< Eigen::MatrixXd > MutableEigenMatrixXdPtr;

  inline EigenMatrixXdPtr
  convertOpenMSMatrix2EigenMatrixXd( const Matrix<double>& m )
  {
    MutableEigenMatrixXdPtr em ( new Eigen::MatrixXd(m.rows(), m.cols()) );
    for (unsigned i=0; i<m.rows(); ++i)
    {
      for (unsigned j=0; j<m.cols(); ++j)
      {
        (*em)(i,j) = m(i,j);
      }
    }
    return em;
  }

  inline bool
  matrixIsIdentityMatrix(const Matrix<double>& channel_frequency)
  {
    for (Matrix<double>::SizeType i = 0; i < channel_frequency.rows(); ++i)
    {
      for (Matrix<double>::SizeType j = 0; j < channel_frequency.rows(); ++j)
      {
        // check if the entries are those of a identity matrix;
        // i==j -> m(i,j) == 1.0 && i!=j -> m(i,j) == 0.0
        if ((i == j && channel_frequency(i, j) != 1.0) || channel_frequency(i, j) != 0.0)
        {
          return false;
        }
      }
    }
    return true;
  }
}//namespace

