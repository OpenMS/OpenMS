// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
  SavitzkyGolayFilter::SavitzkyGolayFilter() :
    ProgressLogger(),
    DefaultParamHandler("SavitzkyGolayFilter"),
    coeffs_()
  {
    defaults_.setValue("frame_length", 11, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("polynomial_order", 4, "Order or the polynomial that is fitted.");

    defaultsToParam_();
  }

  SavitzkyGolayFilter::~SavitzkyGolayFilter()
  {
  }

  void SavitzkyGolayFilter::updateMembers_()
  {
    frame_size_ = (UInt)param_.getValue("frame_length");
    order_ = (UInt)param_.getValue("polynomial_order");

    //recalculate coefficients
    if (!Math::isOdd(frame_size_))
    {
      frame_size_ += 1;
    }
    coeffs_.resize(frame_size_ * (frame_size_ / 2 + 1));

    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The degree of the polynomial has to be less than the frame length.", String(order_));
    }

    int nr, m = frame_size_ / 2;
    for (int nl = 0; nl <= m; ++nl)
    {
      nr = frame_size_ - 1 - nl;

      int i, j;
      double help;

      gsl_vector * sv = gsl_vector_alloc((int)order_ + 1);
      gsl_vector * work = gsl_vector_alloc((int)order_ + 1);
      gsl_matrix * A = gsl_matrix_calloc(frame_size_, (int)order_ + 1);
      gsl_matrix * V = gsl_matrix_calloc((int)order_ + 1, (int)order_ + 1);


      // compute a vandermonde matrix whose columns are powers of the vector [-nL,...,nR]
      for (i = -nl; i <= nr; i++)
      {
        for (j = 0; j <= (int)order_; j++)
        {
          gsl_matrix_set(A, i + nl, j, gsl_pow_int(i, j));
        }
      }

      // compute the singular-value decomposition of A
      if (gsl_linalg_SV_decomp(A, V, sv, work) != 1)
      {
        // compute B=V*inv(D)
        for (i = 0; i <= (int)order_; ++i)
        {
          gsl_vector_set(sv, i, (gsl_matrix_get(V, 0, i) / gsl_vector_get(sv, i)));
        }

        // compute B*transpose(U)*b, where b is the unit vector b=[1 0 ... 0]
        for (i = 0; i < (int)frame_size_; ++i)
        {
          help = 0;
          for (j = 0; j <= (int)order_; ++j)
          {
            help += gsl_vector_get(sv, j) * gsl_matrix_get(A, i, j);
          }
          coeffs_[(nl + 1) * frame_size_ - i - 1] = help;
        }
      }
      gsl_vector_free(sv);
      gsl_vector_free(work);
      gsl_matrix_free(A);
      gsl_matrix_free(V);
    }
  }

}
