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

      This overload copies the input matrix and vector before solving, preserving the originals.

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

      @param[in,out] A Input matrix A of size m x n (column-major storage). Modified in-place by the solver.
      @param[in] A_rows Number of rows in A
      @param[in] A_cols Number of columns in A
      @param[in,out] b Input vector b of size m. Modified in-place by the solver.
      @param[out] x Output vector with non-negative least square solution of size n.
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(double* A, int A_rows, int A_cols,
                     std::vector<double>& b, std::vector<double>& x);

    /**
      @brief Solve the non-negative least square problem Ax=b using Matrix and vectors.

      This overload works with OpenMS Matrix and std::vector, modifying inputs in-place for efficiency.

      @param[in,out] A Input matrix A of size m x n. Modified in-place by the solver.
      @param[in,out] b Input vector b of size m. Modified in-place by the solver.
      @param[out] x Output vector with non-negative least square solution of size n.
      @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

      @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(Matrix<double>& A, std::vector<double>& b, std::vector<double>& x);
  };

} // namespace OpenMS
