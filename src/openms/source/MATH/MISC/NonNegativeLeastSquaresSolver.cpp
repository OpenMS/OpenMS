// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/MATH/MISC/NNLS/NNLS.h>

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
    for (size_t col = 0; col < A.cols(); ++col)
    {
      for (size_t row = 0; row < A.rows(); ++row)
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
    for (size_t row = 0; row < b.rows(); ++row)
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
    x.resize(a_cols, 1);
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
      return SOLVED;
    else if (mode == 2) // this should not happen (dimensions are bad)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "NonNegativeLeastSquaresSolver::solve() Bad dimension reported!");
    }
    else     /*if (mode==3)*/
      return ITERATION_EXCEEDED;

  }

} // namespace OpenMS
