// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/ML/NNLS/NNLS.h>

namespace OpenMS
{
  Int NonNegativeLeastSquaresSolver::solve(const Matrix<double> & A, const Matrix<double> & b, Matrix<double> & x)
  {

    if (A.rows() != b.rows())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "NNSL::solve() #rows of A does not match #rows of b !");
    }

    // translate A to array a (column major order)
    double * a_vec = new double[A.rows() * A.cols()];
    size_t idx = 0;
    for (long int col = 0; col < A.cols(); ++col)
    {
      for (long int row = 0; row < A.rows(); ++row)
      {
        a_vec[idx] = A(row, col);
        idx++;
      }
    }

#ifdef NNLS_DEBUG
    //std::cout << "A:\n" << A << std::endl;
#endif

    // this needs to be int (not Int, Size or anything else), because the external nnls constructor expects it this way!
    int a_rows = (int)A.rows();
    int a_cols = (int)A.cols();

    // translate b
    double * b_vec = new double[a_rows];
    for (long int row = 0; row < b.rows(); ++row)
    {
      b_vec[row] = b(row, 0);
    }

#ifdef NNLS_DEBUG
    std::cout << "b:\n" << b << std::endl;
#endif

    // prepare solution array (directly copied from example)
    double * x_vec = new double[a_cols + 1];
    double rnorm;
    double * w = new double[a_cols + 1];
    double * zz = new double[a_rows + 1];
    int * indx = new int[a_cols + 1];
    int mode;

#ifdef NNLS_DEBUG
    std::cout << "solving ..." << std::endl;
#endif

    NNLS::nnls_(a_vec, &a_rows, &a_rows, &a_cols, b_vec, x_vec, &rnorm, w, zz, indx, &mode);

    // translate solution back to Matrix:
    x.getEigenMatrix().resize(a_cols, 1);
    x.getEigenMatrix().setZero();
    for (Int row = 0; row < a_cols; ++row)
    {
      x(row, 0) = x_vec[row];
    }

#ifdef NNLS_DEBUG
    std::cout << "done" << std::endl;
    std::cout << "solution x:\n" << x << std::endl;
#endif

    // clean up
    delete[] a_vec;
    delete[] b_vec;
    delete[] x_vec;
    delete[] w;
    delete[] zz;
    delete[] indx;

    if (mode == 1)
    {
      return SOLVED;
    }
    else if (mode == 2) // this should not happen (dimensions are bad)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "NonNegativeLeastSquaresSolver::solve() Bad dimension reported!");
    }
    else     /*if (mode==3)*/
    {
      return ITERATION_EXCEEDED;
    }
  }

} // namespace OpenMS
