// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/Matrix.h>

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

        @param A Input matrix A of size mxn
        @param b Input vector (matrix with one column) b of size mx1
        @param x Output vector with non-negative least square solution of size mx1
        @return status of solution (either NonNegativeLeastSquaresSolver::SOLVED, NonNegativeLeastSquaresSolver::ITERATION_EXCEEDED)

        @throws Exception::InvalidParameters if Matrix dimensions do not fit
    */
    static Int solve(const Matrix<double> & A, const Matrix<double> & b, Matrix<double> & x);
  };

} // namespace OpenMS

