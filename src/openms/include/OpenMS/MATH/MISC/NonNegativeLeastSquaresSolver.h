// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/Utils/MatrixUtils.h>

namespace OpenMS
{
  /**
      @brief Wrapper for a non-negative least squares (NNLS) solver.

      It solves Ax=b, where x>0 in the least squares sense (i.e. minimum residual)
  */
  class OPENMS_DLLAPI NonNegativeLeastSquaresSolver
  {
public:

    enum RETURN_STATUS
    {
      SOLVED,
      ITERATION_EXCEEDED
    };

    /**
      @brief This is a wrapper for the external nnls library for the non-negative least square problem Ax=b, where x>0

      @param A Input matrix A of size m x n
      @param b Input vector (OpenMS::Matrix with one column) b of size m x 1
      @param x Output vector (OpenMS::Matrix with one column) with non-negative least square solution of size n x 1
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(const Matrix<double> & A, const Matrix<double> & b, Matrix<double> & x);

    /**
      @brief This is a wrapper for the external nnls library for the non-negative least square problem Ax=b, where x>0. Works without copies but inputs will be modified.

      @param A Input pointer to Eigen::MatrixXd A of size m x n (Note: due to an in-place algorithm, A will be modified!)
      @param b Input vector b of size m (Note: due to an in-place algorithm, b will be modified!)‚
      @param x Output vector with non-negative least square solution of size n. Contents will be overwritten!‚
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(MutableEigenMatrixXdPtr & A, std::vector<double> & b, std::vector<double> & x);
  };

} // namespace OpenMS

