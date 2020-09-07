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
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <Eigen/Core>
#include <Eigen/SVD>

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
    if (frame_size_ <= order_)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The degree of the polynomial has to be less than the frame length.", String(order_));
    }
    coeffs_.resize(frame_size_ * (frame_size_ / 2 + 1));

    for (int nl = 0; nl <= (int) (frame_size_ / 2); ++nl)
    {
      int nr = frame_size_ - 1 - nl;

      // compute a Vandermonde matrix whose columns are powers of the vector [-nL,...,nR]
      Eigen::MatrixXd A (frame_size_, order_ + 1);
      for (int i = -nl; i <= nr; i++)
      {
        for (int j = 0; j <= static_cast<int>(order_); j++)
        {
          A(i + nl, j) = std::pow((float)i, j); // pow(int, int) is not defined
        }
      }

      // compute the singular-value decomposition of A
      Eigen::JacobiSVD<Eigen::MatrixXd> svd (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

      Eigen::VectorXd B (order_ + 1);
      for (UInt i = 0; i <= order_; ++i)
      {
        B(i) = svd.matrixV()(0, i) / svd.singularValues()(i);
      }

      // compute B*transpose(U)*b, where b is the unit vector b=[1 0 ... 0]
      for (UInt i = 0; i < frame_size_; ++i)
      {
        coeffs_[(nl + 1) * frame_size_ - i - 1] = 0;
        for (UInt j = 0; j <= order_; ++j)
        {
          coeffs_[(nl + 1) * frame_size_ - i - 1] += B(j) * svd.matrixU()(i, j);
        }
      }
    }
  }
}
