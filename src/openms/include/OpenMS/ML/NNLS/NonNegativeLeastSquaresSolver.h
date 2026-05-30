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
    @brief Wrapper around the Lawson-Hanson non-negative least squares (NNLS)
           Fortran routine.

    Solves the constrained least-squares problem
    @c argmin_{x >= 0} ||Ax - b||_2 with all entries of @c x required
    to be non-negative.

    Three overloads are provided:
      - a non-destructive @c Matrix-based one that copies its inputs;
      - an in-place pointer-based one for callers that already own a
        column-major buffer;
      - an in-place @c Matrix / @c std::vector convenience wrapper
        around the pointer overload.

    @ingroup Math
  */
  class OPENMS_DLLAPI NonNegativeLeastSquaresSolver
  {
public:

    /// Return status of @ref solve.
    enum RETURN_STATUS
    {
      SOLVED,             ///< NNLS converged within the iteration limit; @p x holds the solution.
      ITERATION_EXCEEDED  ///< NNLS hit the iteration limit; @p x holds the best result so far.
    };

    /**
      @brief Solve @c argmin_{x >= 0} ||Ax - b||_2 without modifying the inputs.

      Copies @p A and @p b internally, then runs NNLS on the copies and
      writes the result into @p x (resized to @c A.cols() rows, 1 column).

      @param[in]  A Input matrix of size @c m x @c n.
      @param[in]  b Input vector as a column matrix of size @c m x 1.
      @param[out] x Receives the non-negative solution as a column matrix of size @c n x 1.
      @return @ref SOLVED on convergence, @ref ITERATION_EXCEEDED when
              the iteration limit was reached.

      @throws Exception::InvalidParameter when @c A.rows() does not
                                          match @c b.rows(), or when
                                          the underlying NNLS routine
                                          reports an invalid dimension
                                          (should not happen in
                                          practice).
    */
    static Int solve(const Matrix<double>& A, const Matrix<double>& b, Matrix<double>& x);

    /**
      @brief In-place pointer overload of @ref solve.

      Operates directly on the caller-owned buffer for @p A and the
      @p b vector, so neither needs to be copied. After the call the
      contents of @p A and @p b are clobbered by the NNLS workspace
      and should not be read again.

      @param[in,out] A      Pointer to a column-major @c double buffer of
                            length @c A_rows * @c A_cols. The caller
                            retains ownership, must keep the buffer
                            valid for the duration of the call, and
                            must not pass @c nullptr or a buffer that
                            is shorter than required. The contents are
                            unspecified on return.
      @param[in]     A_rows Number of rows of @p A (@c > @c 0).
      @param[in]     A_cols Number of columns of @p A (@c > @c 0).
      @param[in,out] b      Right-hand side of size @c A_rows; contents
                            are unspecified on return.
      @param[out]    x      Receives the non-negative solution of size
                            @c A_cols (resized internally).
      @return @ref SOLVED on convergence, @ref ITERATION_EXCEEDED when
              the iteration limit was reached.

      @throws Exception::InvalidParameter when @c A_rows does not match
                                          @c b.size(), or when the
                                          underlying NNLS routine
                                          reports an invalid dimension.

      @note Prefer the @c Matrix overload for additional bounds
            checking; a future revision may adopt @c std::span.
    */
    static Int solve(double* A, int A_rows, int A_cols,
                     std::vector<double>& b, std::vector<double>& x);

    /**
      @brief In-place @c Matrix / @c std::vector overload of @ref solve.

      Delegates to the pointer overload using @c A.data(), so @p A and
      @p b are clobbered by the NNLS workspace and should not be read
      after the call.

      @param[in,out] A Input matrix of size @c m x @c n; contents are
                       unspecified on return.
      @param[in,out] b Right-hand side of size @c m; contents are
                       unspecified on return.
      @param[out]    x Receives the non-negative solution of size @c n
                       (resized internally).
      @return @ref SOLVED on convergence, @ref ITERATION_EXCEEDED when
              the iteration limit was reached.

      @throws Exception::InvalidParameter when @c A.rows() does not
                                          match @c b.size(), or when
                                          the underlying NNLS routine
                                          reports an invalid dimension.
    */
    static Int solve(Matrix<double>& A, std::vector<double>& b, std::vector<double>& x);
  };

} // namespace OpenMS
