// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <vector>

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
      @brief Solve the non-negative least square problem Ax=b, where x>0

      @param[in] A Input matrix A of size m x n
      @param[in] b Input vector (OpenMS::Matrix with one column) b of size m x 1
      @param[out] x Output vector (OpenMS::Matrix with one column) with non-negative least square solution of size n x 1
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(const Matrix<double>& A, const Matrix<double>& b, Matrix<double>& x);

    /**
      @brief Solve the non-negative least square problem Ax=b, where x>0. Works in-place.

      This version works directly on raw buffers and modifies inputs in-place for efficiency.

      @param A Input matrix A of size m x n (column-major storage). Will be modified!
      @param A_rows Number of rows in A
      @param A_cols Number of columns in A
      @param b Input vector b of size m. Will be modified!
      @param x Output vector with non-negative least square solution of size n. Contents will be overwritten!
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(double* A, int A_rows, int A_cols,
                     std::vector<double>& b, std::vector<double>& x);

    /**
      @brief Solve the non-negative least square problem Ax=b using Matrix and vectors.

      This is a convenience overload that works with OpenMS Matrix and std::vector.
      Note: The matrix A will be modified during computation!

      @param A Input matrix A of size m x n. Will be modified!
      @param b Input vector b of size m. Will be modified!
      @param x Output vector with non-negative least square solution of size n.
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(Matrix<double>& A, std::vector<double>& b, std::vector<double>& x);
  };

} // namespace OpenMS
